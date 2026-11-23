import numpy as np
import pytest

import gasjet


def test_projection_matches_cylinder_and_symmetry():
    x = np.array([-6, -3, 0, 3, 6])
    operator = gasjet.shell_projection(x, [0, 2, 5])
    np.testing.assert_allclose(operator @ np.ones(2), [0, 8, 10, 8, 0])
    assert np.all(operator >= 0)


def test_index_and_phase_units():
    x = np.arange(200)
    frame = np.broadcast_to(0.02 + .001 * x, (128, 200)).copy()
    frame[:, 90:110] -= .1
    result = gasjet.signal_radians(frame)
    np.testing.assert_allclose(result[:, :10], 0, atol=1e-12)
    np.testing.assert_allclose(result[:, 62:82], .2 * np.pi)


def test_nonnegative_gaussian_reconstruction():
    _, exact, reconstructed = gasjet.gaussian_benchmark()
    assert np.all(reconstructed >= 0)
    assert gasjet.relative_rmse(exact, reconstructed) < .08


def test_half_width_and_truncated_profile():
    r = np.arange(4, dtype=float)
    assert gasjet.diameter_at_half_peak(r, [1., .75, .25, 0.]) == pytest.approx(3)
    assert np.isnan(gasjet.diameter_at_half_peak(r, [1., .9, .8, .7]))


def test_rejects_invalid_phase_and_missing_frames(tmp_path):
    p = tmp_path / "bad.txt"
    p.write_text("1 2\n3 4", encoding="utf-8")
    with pytest.raises(ValueError, match="128"):
        gasjet.read_phase(p)
    with pytest.raises(ValueError, match="zehn Differenzbilder"):
        gasjet.differential_frames(tmp_path, 10)
