"""Erzeuge wissenschaftliche Ergebnisgrafiken ohne Ausgabe der Eingangsmatrizen."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from gasjet import (
    HEIGHTS_MM, MM_PER_PIXEL, PRESSURES, diameter_at_half_peak,
    differential_frames, gaussian_benchmark, profiles_at_heights, read_phase,
    reconstruct, relative_rmse, signal_radians,
)


def save(fig, path):
    fig.savefig(path, dpi=135, facecolor="white")
    plt.close(fig)


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def benchmark(output):
    radius, truth, recovered = gaussian_benchmark()
    fig, ax = plt.subplots(figsize=(6.9, 3.4), layout="constrained")
    ax.plot(radius * MM_PER_PIXEL, truth, label="Analytische Dichte")
    ax.plot(radius * MM_PER_PIXEL, recovered, "--", label="Rekonstruktion")
    ax.set(xlabel="Radius (mm)", ylabel="Normierte Dichte", title="Synthetischer Gauß-Test · 1 % Phasenrauschen")
    ax.legend(frameon=False)
    ax.grid(alpha=.15)
    save(fig, output / "referenztest.png")
    return relative_rmse(truth, recovered)


def reference_comparison(root):
    folder = root / "Processed_data" / "10barGasDruck" / "Results"
    def phase(ms):
        return read_phase(folder / f"PHA SID4 13h42m03s{ms}ms.txt")
    before, gas, after = [phase(ms) for ms in (225, 265, 305)]
    left, right = [signal_radians(gas - ref) for ref in (before, after)]
    denominator = np.linalg.norm((left + right) / 2)
    change = float(np.linalg.norm(left - right) / denominator) if denominator else np.nan
    return {"pressure_bar": 10, "relative_rms_reference_change": change}


def analysis(root, output):
    results = {p: reconstruct(differential_frames(root, p)) for p in PRESSURES}
    fig_map, axes_map = plt.subplots(2, 2, figsize=(8.4, 6.3), layout="constrained", sharex=True, sharey=True)
    fig_profiles, axes_profiles = plt.subplots(2, 2, figsize=(8.4, 6.1), layout="constrained", sharex=True, sharey=True)
    fig_map.suptitle("Regularisierte Rekonstruktion · gemeinsamer Dichtemaßstab", fontweight="bold")
    fig_profiles.suptitle("Radialprofile bei 2 mm · Methoden und Streuung", fontweight="bold")
    density_max = max(r["density"].max() for r in results.values()) / 1e18
    quality, metrics = [], []
    height_index = HEIGHTS_MM.index(2.0)
    for (p, result), map_ax, profile_ax in zip(results.items(), axes_map.flat, axes_profiles.flat):
        r = result["radii_mm"]
        field = result["density"] / 1e18
        profiles = profiles_at_heights(result["density"])
        cosine = profiles_at_heights(result["cosine_density"])
        frame_profiles = result["replicate_profiles"]
        shifts = result["shift_profiles"]
        strengths = result["regularization_profiles"]
        vmax_r = r[-1] + (r[1] - r[0]) / 2
        im = map_ax.imshow(field[:, np.r_[len(r)-1:-1:-1, 0:len(r)]], aspect="auto",
                           extent=(-vmax_r, vmax_r, -MM_PER_PIXEL / 2, (field.shape[0] - .5) * MM_PER_PIXEL),
                           cmap="viridis", vmin=0, vmax=density_max)
        map_ax.set_title(f"{p} bar Vordruck", loc="left")
        samples = frame_profiles[:, height_index] / 1e18
        profile_ax.fill_between(r, samples.min(axis=0), samples.max(axis=0),
                                color="#006aa2", alpha=.17, label="Spannweite der zehn Differenzbilder")
        profile_ax.plot(r, profiles[height_index] / 1e18, color="#006aa2", label="Regularisierte Schalen")
        profile_ax.plot(r, cosine[height_index] / 1e18, "--", color="#ae5122", label="Kosinus-Basis")
        profile_ax.axhline(0, color="#555555", linewidth=.6)
        profile_ax.set_title(f"{p} bar Vordruck", loc="left")
        profile_ax.grid(alpha=.15)
        quality.append({"pressure_bar": p, "axis_shift_px": result["center_offset_px"],
                        "shell_projection_relative_rmse": result["shell_projection_error"],
                        "cosine_projection_relative_rmse": result["cosine_projection_error"],
                        "left_right_profile_relative_difference": result["left_right_profile_difference"]})
        for i, height in enumerate(HEIGHTS_MM):
            profile = profiles[i]
            peak = float(profile.max())
            metrics.append({
                "pressure_bar": p, "height_above_roi_bottom_mm": height,
                "center_shell_density_assumed_cm3": float(profile[0]),
                "peak_density_assumed_cm3": peak,
                "radius_of_peak_mm": float(r[np.argmax(profile)]),
                "half_peak_outer_diameter_mm": diameter_at_half_peak(r, profile),
                "pair_peak_min_assumed_cm3": float(frame_profiles[:, i].max(axis=1).min()),
                "pair_peak_max_assumed_cm3": float(frame_profiles[:, i].max(axis=1).max()),
                "axis_shift_peak_change_fraction": float(np.max(np.abs(shifts[:, i].max(axis=1) / peak - 1))) if peak else np.nan,
                "regularization_peak_change_fraction": float(np.max(np.abs(strengths[:, i].max(axis=1) / peak - 1))) if peak else np.nan,
            })
    fig_map.colorbar(im, ax=axes_map, label="Teilchendichte (10¹⁸ cm⁻³, historischer Faktor)")
    fig_map.supxlabel("Abstand von der Jetachse (mm)")
    fig_map.supylabel("Höhe im Auswertebereich (mm)")
    fig_profiles.supxlabel("Radius (mm)")
    fig_profiles.supylabel("Teilchendichte (10¹⁸ cm⁻³, historischer Faktor)")
    axes_profiles[0, 1].legend(frameon=False, fontsize=7)
    save(fig_map, output / "dichtekarten.png")
    save(fig_profiles, output / "profilvergleich.png")
    write_csv(output / "profile.csv", metrics)
    write_csv(output / "diagnostik.csv", quality)
    write_csv(output / "referenzwahl.csv", [reference_comparison(root)])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.data_root and args.output_dir.resolve().is_relative_to(args.data_root.resolve()):
        parser.error("Der Ausgabeordner muss außerhalb der Originaldaten liegen.")
    try:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        error = benchmark(args.output_dir)
        if args.data_root is not None:
            analysis(args.data_root, args.output_dir)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.exit(1, f"Auswertung fehlgeschlagen: {exc}\n")
    print(f"Synthetischer relativer Rekonstruktionsfehler: {error:.3%}")


if __name__ == "__main__":
    main()
