"""Vorwärtsmodell und inverse Abel-Rekonstruktion im SID4-Phasenraster."""

from pathlib import Path

import numpy as np
from scipy.optimize import nnls

PRESSURES = (10, 20, 30, 40)
HEIGHTS_MM = (0.5, 1.0, 2.0, 3.0, 4.0)
MM_PER_PIXEL = 5 / 118
WAVELENGTH_UM = 0.635
HISTORICAL_CONVERSION_CM3 = 1.66541e22
RADIUS_PX = 70.0
RADIAL_BINS = 32
REGULARIZATION = 0.15
ROI = np.s_[5:123, 28:179]


def read_phase(path: Path, delimiter=None):
    frame = np.loadtxt(path, delimiter=delimiter, ndmin=2)
    if frame.shape != (128, 200) or not np.isfinite(frame).all():
        raise ValueError(f"{path.name}: erwartet werden 128 × 200 endliche Phasenwerte.")
    return frame


def differential_frames(root: Path, pressure: int):
    paths = sorted((root / "Abbildungen" / f"{pressure}bar").glob("substracted *.txt"))
    if len(paths) != 10:
        raise ValueError(f"{pressure} bar: zehn Differenzbilder erwartet, {len(paths)} gefunden.")
    return np.stack([read_phase(p, delimiter=",") for p in paths])


def signal_radians(frame):
    """Sensorphase in Wellen, historisches negatives Vorzeichen des Gasjets."""
    image = frame[ROI]
    x = np.arange(image.shape[1], dtype=float)
    background = np.r_[0:15, image.shape[1]-15:image.shape[1]]
    design = np.column_stack((x, np.ones_like(x)))
    coefficients = np.linalg.lstsq(design[background], image[:, background].T, rcond=None)[0]
    return (coefficients.T @ design.T - image) * (2 * np.pi)


def estimate_center(signal):
    """Eine gemeinsame Subpixelachse, aus der Spiegelung über alle Höhen."""
    x = np.arange(signal.shape[1], dtype=float)
    middle = (x.size - 1) / 2
    profile = signal.mean(axis=0)
    distances = np.arange(1, 61, dtype=float)
    candidates = middle + np.linspace(-3, 3, 49)
    losses = [
        np.mean((np.interp(c - distances, x, profile)
                 - np.interp(c + distances, x, profile)) ** 2)
        for c in candidates
    ]
    return float(candidates[np.argmin(losses)])


def shell_projection(offsets, edges):
    """Sehnenlängen in Pixeln für konzentrische, stückweise konstante Schalen."""
    edges = np.asarray(edges, dtype=float)
    offsets = np.asarray(offsets, dtype=float)
    if (edges.ndim != 1 or len(edges) < 2 or edges[0] != 0
            or not np.isfinite(edges).all() or (np.diff(edges) <= 0).any()
            or offsets.ndim != 1 or not np.isfinite(offsets).all()):
        raise ValueError("Das Radialgitter muss bei null beginnen und streng ansteigen.")
    radius = np.sqrt(np.maximum(edges[None, :] ** 2 - offsets[:, None] ** 2, 0))
    return 2 * np.diff(radius, axis=1)


def solve_shells(matrix, signal, strength=REGULARIZATION):
    """Nichtnegative Schalen mit Glättung der zweiten radialen Differenz."""
    matrix = np.asarray(matrix, dtype=float)
    signal = np.asarray(signal, dtype=float)
    if (matrix.ndim != 2 or signal.ndim != 2 or signal.shape[1] != matrix.shape[0]
            or not np.isfinite(matrix).all() or not np.isfinite(signal).all()
            or not np.isfinite(strength) or strength < 0):
        raise ValueError("Projektionsmatrix, Phasendaten oder Regularisierung sind ungültig.")
    bins = matrix.shape[1]
    penalty = strength * np.linalg.norm(matrix) / np.sqrt(bins)
    smoothness = np.diff(np.eye(bins), n=2, axis=0)
    design = np.vstack((matrix, penalty * smoothness))
    return np.stack([
        nnls(design, np.r_[row, np.zeros(smoothness.shape[0])])[0]
        for row in signal
    ])


def profiles_at_heights(field):
    output = []
    for height in HEIGHTS_MM:
        index = field.shape[0] - 1 - round(height / MM_PER_PIXEL)
        lo, hi = max(0, index - 2), min(len(field), index + 3)
        output.append(field[lo:hi].mean(axis=0))
    return np.asarray(output)


def diameter_at_half_peak(radius_mm, profile):
    """Doppelter äußerer Halbwertsradius relativ zur gemeinsamen Jetachse."""
    peak = np.max(profile)
    if peak <= 0 or profile[-1] >= peak / 2:
        return np.nan
    idx = np.flatnonzero(profile >= peak / 2)[-1]
    return float(2 * np.interp(
        peak / 2, profile[idx:idx+2][::-1], radius_mm[idx:idx+2][::-1]))


def relative_rmse(observed, predicted):
    denominator = np.linalg.norm(observed)
    return float(np.linalg.norm(observed - predicted) / denominator) if denominator else np.nan


def reconstruct(frames):
    signals = np.stack([signal_radians(frame) for frame in frames])
    signal = signals.mean(axis=0)
    center = estimate_center(signal)
    x = np.arange(signal.shape[1], dtype=float) - center
    usable = np.abs(x) < RADIUS_PX
    offsets = x[usable]
    y = signal[:, usable]
    edges = np.linspace(0, RADIUS_PX, RADIAL_BINS + 1)
    radii = (edges[:-1] + edges[1:]) / 2
    projection = shell_projection(offsets, edges)
    shells = solve_shells(projection, y)
    basis = np.column_stack([np.ones_like(radii)] + [
        1 - (-1)**n * np.cos(n * np.pi * radii / RADIUS_PX) for n in range(1, 11)
    ])
    projected_basis = projection @ basis
    cosine = np.linalg.lstsq(projected_basis, y.T, rcond=None)[0].T @ basis.T
    replicates = np.stack([
        profiles_at_heights(solve_shells(projection, a[:, usable])) for a in signals
    ])
    shifts = np.stack([
        profiles_at_heights(solve_shells(shell_projection(offsets - delta, edges), y))
        for delta in (-1, 1)
    ])
    regularization_variants = np.stack([
        profiles_at_heights(solve_shells(projection, y, strength))
        for strength in (REGULARIZATION / 2, 2 * REGULARIZATION)
    ])
    left = offsets <= 0
    right = offsets >= 0
    halves = np.stack([
        profiles_at_heights(solve_shells(projection[half], y[:, half]))
        for half in (left, right)
    ])
    scale = HISTORICAL_CONVERSION_CM3 * WAVELENGTH_UM / (1000 * MM_PER_PIXEL)
    return {
        "radii_mm": radii * MM_PER_PIXEL,
        "density": shells * scale,
        "cosine_density": cosine * scale,
        "replicate_profiles": replicates * scale,
        "shift_profiles": shifts * scale,
        "regularization_profiles": regularization_variants * scale,
        "half_profiles": halves * scale,
        "center_offset_px": center - (signal.shape[1] - 1) / 2,
        "shell_projection_error": relative_rmse(y, shells @ projection.T),
        "cosine_projection_error": relative_rmse(y, cosine @ projection.T),
        "left_right_profile_difference": relative_rmse(*halves),
    }


def gaussian_benchmark():
    edges = np.linspace(0, RADIUS_PX, RADIAL_BINS + 1)
    radii = (edges[:-1] + edges[1:]) / 2
    offsets = np.arange(-69, 70, dtype=float)
    sigma = 14.0
    profile = np.exp(-0.5 * (radii / sigma) ** 2)
    line_integral = np.sqrt(2 * np.pi) * sigma * np.exp(-0.5 * (offsets / sigma) ** 2)
    noise = np.random.default_rng(2026).normal(0, 0.01 * line_integral.max(), offsets.shape)
    matrix = shell_projection(offsets, edges)
    recovered = solve_shells(matrix, (line_integral + noise)[None, :])[0]
    return radii, profile, recovered
