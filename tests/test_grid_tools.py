# =============================================================================
# Author: Luigi Sante Zampa
# Affiliation: NllocPy project (institutional affiliation not specified)
# Created: 2026-11-03
# Last modified: 2026-11-10
#
# Purpose:
#   Round-trip regression tests for native NonLinLoc headers, buffers and grid metadata.
# =============================================================================

"""Regression tests for NonLinLoc-compatible grid I/O."""
# -----------------------------------------------------------------------------
# Import external libraries
# -----------------------------------------------------------------------------
import numpy as np
import pytest

from NllocPy.grid_tools import NLLGrid, velocity_to_slow_len


@pytest.mark.parametrize("float_type", ["FLOAT", "DOUBLE"])
def test_model_grid_roundtrip_and_buffer_order(tmp_path, float_type):
    velocity = np.arange(1, 25, dtype=float).reshape(2, 3, 4) + 2
    slow_len = velocity_to_slow_len(velocity, 0.5)
    grid = NLLGrid()
    grid.array = slow_len
    grid.type = "SLOW_LEN"
    grid.float_type = float_type
    grid.x_orig, grid.y_orig, grid.z_orig = -10.0, 15.0, 2.0
    grid.dx = grid.dy = grid.dz = 0.5
    grid.proj_name = "SIMPLE"
    grid.orig_lat, grid.orig_lon, grid.map_rot = 46.0, 13.0, 15.0
    grid.basename = str(tmp_path / "layer.P.mod")
    grid.write_hdr_file()
    grid.write_buf_file()

    output = NLLGrid(str(tmp_path / "layer.P.mod.buf"))
    assert (output.nx, output.ny, output.nz) == (2, 3, 4)
    assert output.proj_name == "SIMPLE"
    assert (output.orig_lat, output.orig_lon, output.map_rot) == (46, 13, 15)
    np.testing.assert_allclose(output.array, slow_len, rtol=1e-6)
    dtype = np.float32 if float_type == "FLOAT" else np.float64
    np.testing.assert_array_equal(
        np.fromfile(tmp_path / "layer.P.mod.buf", dtype=dtype),
        slow_len.astype(dtype).ravel(order="C"))
    assert output.get_xyz(1, 2, 3) == (-9.5, 16.0, 3.5)


def test_legacy_header_no_precision_and_station(tmp_path):
    base = tmp_path / "layer.P.time.STA"
    (tmp_path / "layer.P.time.STA.hdr").write_text(
        "2 2 2 0 0 0 1 1 1 TIME\n"
        "STA 1.1 2.2 0.3\n"
        "TRANSFORM LAMBERT RefEllipsoid WGS-84 LatOrig 46 LongOrig 13 "
        "FirstStdParal 45 SecondStdParal 47 RotCW 2\n")
    np.arange(8, dtype=np.float32).tofile(tmp_path / "layer.P.time.STA.buf")
    grid = NLLGrid(str(base))
    assert grid.float_type == "FLOAT"
    assert grid.station == "STA"
    assert grid.proj_name == "LAMBERT"
    assert grid.first_std_paral == 45
    assert grid.second_std_paral == 47
    np.testing.assert_equal(grid.array, np.arange(8).reshape(2, 2, 2))


def test_invalid_velocity_and_truncated_buffer(tmp_path):
    with pytest.raises(ValueError, match="finite positive"):
        velocity_to_slow_len(np.array([[[-1.0]]]), 1)
    grid = NLLGrid()
    grid.array = np.ones((2, 2, 3))
    grid.type = "SLOW_LEN"
    grid.basename = str(tmp_path / "model")
    grid.write_hdr_file()
    grid.array[0, 0, 0] = np.nan
    with pytest.raises(ValueError, match="NaN/Inf"):
        grid.write_buf_file()
    np.ones(7, dtype=np.float32).tofile(tmp_path / "model.buf")
    with pytest.raises(ValueError, match="expected 12"):
        NLLGrid(str(tmp_path / "model"))


def test_array_crop_updates_geometry():
    grid = NLLGrid()
    grid.array = np.zeros((4, 5, 6))
    grid.array = grid.array[1:3, 2:4, 3:5]
    assert (grid.nx, grid.ny, grid.nz) == (2, 2, 2)


def test_packed_angle_read(tmp_path):
    p = tmp_path / "a.P.angle.STA"
    (tmp_path / "a.P.angle.STA.hdr").write_text(
        "1 1 2 0 0 0 1 1 1 ANGLE FLOAT\nTRANSFORM NONE\n")
    np.array([16 * 305 + 1, 1234, 0, 0], dtype=np.uint16).tofile(
        tmp_path / "a.P.angle.STA.buf")
    g = NLLGrid(str(p))
    assert g.dip[0, 0, 0] == 30.5
    assert g.azimuth[0, 0, 0] == 123.4
    assert g.quality[0, 0, 0] == 1
    assert np.isnan(g.dip[0, 0, 1])
    assert g.array is None
