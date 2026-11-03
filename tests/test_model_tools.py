"""Coordinate sorting, native binary layout and XY/XZ/YZ plotting regressions."""
import numpy as np
import pytest

from NllocPy.model_tools import (build_regular_model, create_model_file,
                                read_model_file, plot_model_grids)
from NllocPy.grid_io import NLLGrid


def _shuffled_model(seed=11, z_positive='up', vs=True):
    x = np.array([10., 11., 12., 13.])
    y = np.array([20., 21., 22.])
    depth = np.array([1., 2., 3., 4., 5.])
    X, Y, Z = np.meshgrid(x, y, depth, indexing='ij')
    vp = 4 + (X - 10) + 10 * (Y - 20) + 100 * (Z - 1)
    d = {'x': X.ravel(), 'y': Y.ravel(),
         'z': (-Z if z_positive == 'up' else Z).ravel(),
         'vp': vp.ravel()}
    if vs:
        d['vs'] = (vp / 2).ravel()
    ix = np.random.default_rng(seed).permutation(vp.size)
    return {key: val[ix] for key, val in d.items()}, vp, x, y, depth


@pytest.mark.parametrize('vertical', ['up', 'down'])
@pytest.mark.parametrize('seed', [7, 18])
def test_unordered_model_exact_p_s_and_binary(tmp_path, vertical, seed):
    model, expected, x, y, depth = _shuffled_model(seed, vertical)
    pf, sf, grd = create_model_file(model, path=tmp_path, step=1,
                                    z_positive=vertical, plot=False, proj_name='NONE')
    assert pf.endswith('.P.mod.buf')
    assert sf.endswith('.S.mod.buf')
    assert grd is not None
    vp = read_model_file(pf, return_type='array')[0]
    vs = read_model_file(sf, return_type='array')[0]
    np.testing.assert_allclose(vp, expected, rtol=1e-6)
    np.testing.assert_allclose(vs, expected / 2, rtol=1e-6)
    g = NLLGrid(pf)
    assert (g.x_orig, g.y_orig, g.z_orig) == (10, 20, 1)
    assert (g.nx, g.ny, g.nz) == (4, 3, 5)
    ix, iy, iz = (3, 1, 4)
    flat = ((ix * g.ny) + iy) * g.nz + iz
    raw = np.fromfile(pf, dtype=np.float32)
    assert raw[flat] == pytest.approx(1 / expected[ix, iy, iz], rel=1e-6)
    asdict = read_model_file(pf, return_type='dict')[0]
    assert asdict['values'][flat] == pytest.approx(expected[ix, iy, iz])
    assert (asdict['x'][flat], asdict['y'][flat], asdict['z'][flat]) == (13, 21, 5)


def test_no_interpolation_when_exact_even_if_requested():
    model, v, *_ = _shuffled_model()
    *_, vp, vs = build_regular_model(model, step=1, interpolate=True)
    np.testing.assert_array_equal(vp, v)
    np.testing.assert_array_equal(vs, v / 2)


def test_duplicate_missing_and_nonpositive_nodes():
    model, *_ = _shuffled_model()
    for k in model:
        model[k][1] = model[k][0]
    with pytest.raises(ValueError, match='Duplicate'):
        build_regular_model(model)
    model, *_ = _shuffled_model()
    model = {k: val[:-1] for k, val in model.items()}
    with pytest.raises(ValueError, match='incomplete'):
        build_regular_model(model)
    model, *_ = _shuffled_model()
    model['vp'][0] = np.nan
    with pytest.raises(ValueError, match='NaN'):
        build_regular_model(model)


def test_nonuniform_and_interpolation():
    model, *_ = _shuffled_model()
    model = {k: v[1:] for k, v in model.items()}
    out = build_regular_model(model, step=1, interpolate=True)
    assert out[-2].shape == (4, 3, 5)
    assert np.isfinite(out[-2]).all()


def test_limits_and_cropping(tmp_path):
    model, v, *_ = _shuffled_model()
    p, s, _ = create_model_file(model, path=tmp_path, step=1, proj_name='NONE')
    sub = read_model_file(p, return_type='array', lim=(11, 12, 21, 22, 2, 4))[0]
    np.testing.assert_allclose(sub, v[1:3, 1:3, 1:4], rtol=1e-6)
    crop = read_model_file(p, lim=(11, 12, 21, 22, 2, 4))[0]
    assert (crop.x_orig, crop.y_orig, crop.z_orig) == (11, 21, 2)
    d = read_model_file(p, return_type='dict')[0]
    assert d['z'][0] == 1


def test_oriented_plot_exact_pixels(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    model, expected, x, y, depth = _shuffled_model()
    p, s, _ = create_model_file(model, path=tmp_path, step=1,
                                proj_name='NONE')
    g = read_model_file(p)[0]
    fig, axes = plot_model_grids([g], section=(2, 1, 3))
    yz = axes[0, 0].images[0]
    xz = axes[0, 1].images[0]
    xy = axes[0, 2].images[0]
    np.testing.assert_allclose(np.asarray(yz.get_array()), expected[2,:,:].T)
    np.testing.assert_allclose(np.asarray(xz.get_array()), expected[:,1,:].T)
    np.testing.assert_allclose(np.asarray(xy.get_array()), expected[:,:,3].T)
    assert axes[0, 0].yaxis_inverted()
    assert axes[0, 1].yaxis_inverted()
    assert not axes[0, 2].yaxis_inverted()
    np.testing.assert_allclose(yz.get_extent(), (19.5, 22.5, 0.5, 5.5))
    np.testing.assert_allclose(xz.get_extent(), (9.5, 13.5, 0.5, 5.5))
    np.testing.assert_allclose(xy.get_extent(), (9.5, 13.5, 19.5, 22.5))
    assert axes[0, 1].get_ylabel() == 'depth (km)'
    plt.close(fig)


def test_read_model_file_plots_internally(tmp_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    model, expected, *_ = _shuffled_model()
    p, *_ = create_model_file(model, path=tmp_path, step=1,
                              proj_name='NONE')
    grids = read_model_file(p, plot=True, section=(None, 1, 2))
    fig = plt.gcf()
    assert len(grids) == 1
    axes = [ax for ax in fig.axes if ax.images]
    assert len(axes) == 2
    np.testing.assert_allclose(np.asarray(axes[0].images[0].get_array()),
                               expected[:,1,:].T)
    plt.close(fig)
