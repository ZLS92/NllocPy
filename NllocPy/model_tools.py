"""Build, inspect and plot NonLinLoc velocity models from unordered XYZ samples.

Input XYZ/Vp/Vs arrays contain one sample per row, all in km / km/s.
The vertical convention is explicit: z_positive='up' converts elevations
into NonLinLoc depth, while z_positive='down' uses depths unchanged.
Native grids are ordered [x, y, depth], with increasing coordinates and
Z varying fastest in the binary file.
"""

from pathlib import Path
from collections.abc import Mapping

import numpy as np

from .grid_io import NLLGrid, velocity_to_slow_len


_TOL = 1e-6


def _positive_vertical(z_positive):
    if z_positive not in ('up', 'down'):
        raise ValueError("z_positive must be 'up' (elevation) or 'down' (depth)")
    return z_positive == 'up'


def _columns(model, z_positive):
    if not isinstance(model, Mapping):
        raise TypeError("model must be a dict-like mapping of 1D coordinate/velocity arrays")
    keys = ('x', 'y', 'z', 'vp')
    if any(k not in model for k in keys):
        raise ValueError("model must have x, y, z and vp arrays")
    cols = {k: np.asarray(model[k], dtype=float) for k in keys}
    if 'vs' in model and model['vs'] is not None:
        cols['vs'] = np.asarray(model['vs'], dtype=float)
    n = cols['x'].size
    if not n:
        raise ValueError('model contains no samples')
    for key, arr in cols.items():
        if arr.ndim != 1 or arr.size != n:
            raise ValueError(f'{key} must be a 1D array of length {n}')
        if not np.isfinite(arr).all():
            raise ValueError(f'{key} contains NaN or infinite values')
        if key in ('vp', 'vs') and np.any(arr <= 0):
            raise ValueError(f'{key} must contain strictly positive velocities in km/s')
    if _positive_vertical(z_positive):
        cols['z'] = -cols['z']
    return cols


def _axis(lo, hi, step):
    count = int(round((hi - lo) / step))
    if count < 1:
        raise ValueError('Each grid axis must contain at least two nodes')
    if not np.isclose(lo + count * step, hi, rtol=0, atol=step * _TOL):
        raise ValueError('Model extent must be divisible by step; choose compatible lim/step')
    return lo + np.arange(count + 1, dtype=float) * step


def build_regular_model(model, step=None, lim=None, *, z_positive='up',
                        interpolate=False, method='nearest'):
    """Return (x, y, depth, vp, vs) with 3D arrays indexed [x, y, z].

    All unordered input rows are kept together during sorting. A complete,
    regularly spaced source grid is mapped exactly, without interpolation.
    Missing/duplicate nodes or regridding require interpolate=True;
    duplicate coordinates are always rejected as ambiguous. When resampling,
    SciPy is imported lazily, and coverage/NaNs are checked.
    """
    c = _columns(model, z_positive)
    axes = [np.unique(c[k]) for k in ('x', 'y', 'z')]
    orig_ni = tuple(len(a) for a in axes)
    if any(n < 2 for n in orig_ni):
        raise ValueError('At least two distinct coordinate values per axis are required')
    src_index = [np.searchsorted(a, c[k]) for k, a in zip(('x', 'y', 'z'), axes)]
    linear = np.ravel_multi_index(src_index, orig_ni)
    if np.unique(linear).size != len(linear):
        raise ValueError('Duplicate (x,y,z) nodes found; velocities would be ambiguous')

    diffs = [np.diff(a) for a in axes]
    intervals = [float(np.median(d)) for d in diffs]
    uniform = all(np.allclose(d, size, rtol=_TOL, atol=_TOL * size)
                  for d, size in zip(diffs, intervals))
    if step is None:
        step = min(intervals)
    step = float(step)
    if not np.isfinite(step) or step <= 0:
        raise ValueError('step must be finite and positive (km)')

    if lim is None:
        bounds = [(float(a[0]), float(a[-1])) for a in axes]
    else:
        if len(lim) != 6:
            raise ValueError('lim must contain xmin,xmax,ymin,ymax,zmin,zmax')
        if not np.isfinite(lim).all():
            raise ValueError('lim must contain finite coordinates (km)')
        lx = sorted((float(lim[0]), float(lim[1])))
        ly = sorted((float(lim[2]), float(lim[3])))
        lz = sorted((float(lim[4]), float(lim[5])))
        if z_positive == 'up':
            lz = sorted((-lz[0], -lz[1]))
        bounds = [lx, ly, lz]

    targets = tuple(_axis(lo, hi, step) for lo, hi in bounds)
    sizes = tuple(len(a) for a in targets)
    if np.prod(sizes, dtype=np.int64) > 20_000_000:
        raise ValueError('Target grid exceeds 20 million nodes; increase step or crop lim')
    direct = (uniform and len(linear) == np.prod(orig_ni, dtype=np.int64)
              and all(np.array_equal(s, t) or np.allclose(s, t, rtol=0, atol=step * _TOL)
                      for s, t in zip(axes, targets) if len(s) == len(t))
              and all(len(s) == len(t) for s, t in zip(axes, targets)))
    if not direct and not interpolate:
        raise ValueError('Input nodes are incomplete/nonuniform or target step differs; '
                         'set interpolate=True to resample')
    values = {}
    for wave in ('vp', 'vs'):
        if wave not in c:
            values[wave] = None
            continue
        if direct:
            flat = np.empty(int(np.prod(sizes)), dtype=float)
            flat[linear] = c[wave]
            values[wave] = flat.reshape(sizes, order='C')
        else:
            from scipy.interpolate import griddata
            if method not in ('nearest', 'linear', 'cubic'):
                raise ValueError("method must be 'nearest' or 'linear' for 3D grids")
            if method == 'cubic':
                raise ValueError('scipy.griddata cubic interpolation is not available in 3D')
            X, Y, Z = np.meshgrid(*targets, indexing='ij')
            values[wave] = griddata(
                np.column_stack([c['x'], c['y'], c['z']]), c[wave],
                (X, Y, Z), method=method)
            if not np.isfinite(values[wave]).all():
                raise ValueError(f'Interpolation of {wave} left undefined nodes; '
                                 'use method="nearest" or extend input coverage')
            if np.any(values[wave] <= 0):
                raise ValueError(f'Interpolation produced non-positive {wave} values')
    return (*targets, values['vp'], values['vs'])


def create_model_file(model=None, step=None, float_type='FLOAT', lim=None,
                      plot=False, name='layer', method='nearest', path='.',
                      proj_name='SIMPLE', orig_lon=0., orig_lat=0.,
                      interpolate=False, transpose=(0, 1, 2),
                      section=(0, 0, 0), *, z_positive='up',
                      map_rot=0., validate_output=True, **kwargs):
    """Create <name>.P.mod.{hdr,buf} and optional S counterpart.

    transpose is accepted for backward API compatibility but is ignored:
    array order is derived unambiguously from XYZ sample coordinates.
    The old tuple return (P .buf, S .buf or None, last NLLGrid) is preserved.
    """
    if model is None:
        raise ValueError('model cannot be None')
    xyz = build_regular_model(model, step, lim, z_positive=z_positive,
                              interpolate=interpolate, method=method)
    x, y, z, vp, vs = xyz
    dx = float(x[1] - x[0])
    result = []
    folder = Path(path)
    folder.mkdir(parents=True, exist_ok=True)
    last = None
    for phase, velocity in (('P', vp), ('S', vs)):
        if velocity is None:
            result.append(None)
            continue
        grid = NLLGrid()
        grid.array = velocity_to_slow_len(velocity, dx)
        grid.x_orig, grid.y_orig, grid.z_orig = x[0], y[0], z[0]
        grid.dx = grid.dy = grid.dz = dx
        grid.type = 'SLOW_LEN'
        grid.float_type = float_type.upper()
        grid.proj_name = proj_name
        grid.orig_lon, grid.orig_lat, grid.map_rot = orig_lon, orig_lat, map_rot
        grid.basename = str(folder / f'{name}.{phase}.mod')
        Path(grid.basename).parent.mkdir(parents=True, exist_ok=True)
        grid.write_hdr_file()
        grid.write_buf_file()
        if validate_output:
            check = NLLGrid(grid.basename)
            if (check.array.shape != grid.array.shape
                    or not np.allclose(check.array, grid.array,
                                       rtol=1e-6, atol=1e-8)):
                raise IOError(f'NonLinLoc round-trip mismatch in {grid.basename}')
            for field in ('x_orig', 'y_orig', 'z_orig', 'dx', 'dy', 'dz'):
                if not np.isclose(getattr(check, field), getattr(grid, field),
                                  rtol=0, atol=6e-7):
                    raise IOError(f'NonLinLoc header mismatch in {field}')
        result.append(grid.basename + '.buf')
        last = grid
    if plot:
        plot_model_grids([NLLGrid(str(folder / f'{name}.{phase}.mod'))
                          for phase, value in (('P', vp), ('S', vs)) if value is not None],
                         section=section, **kwargs)
    return result[0], result[1], last


def _model_roots(model_data):
    if isinstance(model_data, NLLGrid):
        return [model_data]
    if isinstance(model_data, (list, tuple)):
        return list(model_data)
    path = Path(model_data)
    if path.is_dir():
        return [str(p)[:-4] for p in sorted(path.glob('*.hdr'))
                if any(tag in p.name for tag in ('.mod.hdr', '.time.', '.angle.'))]
    if path.suffix in ('.hdr', '.buf'):
        path = path.with_suffix('')
    return [str(path)]


def _crop(grid, lim):
    if lim is None:
        return
    if len(lim) != 6:
        raise ValueError('lim must contain xmin,xmax,ymin,ymax,zmin,zmax')
    bounds = [(lim[0], lim[1]), (lim[2], lim[3]), (lim[4], lim[5])]
    slices = []
    origins = []
    for (a, b), orig, d, n in zip(bounds,
                                  (grid.x_orig, grid.y_orig, grid.z_orig),
                                  (grid.dx, grid.dy, grid.dz),
                                  (grid.nx, grid.ny, grid.nz)):
        if a is None and b is None:
            lo, hi = 0, n - 1
        else:
            if a is None:
                aa, bb = -np.inf, float(b)
            elif b is None:
                aa, bb = float(a), np.inf
            else:
                aa, bb = sorted((float(a), float(b)))
            coords = orig + np.arange(n) * d
            idx = np.flatnonzero((coords >= aa - d * _TOL) & (coords <= bb + d * _TOL))
            if not len(idx):
                raise ValueError('Requested crop is outside grid extent')
            lo, hi = int(idx[0]), int(idx[-1])
        slices.append(slice(lo, hi + 1))
        origins.append(orig + lo * d)
    grid.array = grid.array[tuple(slices)]
    grid.x_orig, grid.y_orig, grid.z_orig = origins


def _values(grid, slo2vel):
    if grid.array is None:
        raise ValueError('Packed angle grids cannot be converted to a scalar array')
    if slo2vel and grid.type == 'SLOW_LEN':
        return grid.dx / grid.array
    return grid.array


def read_model_file(model_data, plot=False, section=(0, 0, 0),
                    vmin=None, vmax=None, return_type='class',
                    slo2vel=True, lim=None):
    """Read NonLinLoc grids; XYZ dictionaries use Z as depth positive down.

    With plot=True, horizontal XY and depth sections XZ/YZ use the true
    coordinate extents, a correctly oriented image, and depth increasing down.
    """
    inputs = _model_roots(model_data)
    grids = []
    for item in inputs:
        g = item if isinstance(item, NLLGrid) else NLLGrid(item)
        if lim is not None:
            if g.array is None:
                raise ValueError('Cannot crop packed angle grid')
            if isinstance(item, NLLGrid):
                import copy
                g = copy.copy(item)
            _crop(g, lim)
        grids.append(g)
    if plot:
        plot_model_grids(grids, section=section, vmin=vmin, vmax=vmax,
                         slo2vel=slo2vel)
    if return_type == 'class':
        return grids
    if return_type == 'array':
        return [_values(g, slo2vel) for g in grids]
    if return_type == 'dict':
        ret = []
        for g in grids:
            x = g.x_orig + np.arange(g.nx) * g.dx
            y = g.y_orig + np.arange(g.ny) * g.dy
            z = g.z_orig + np.arange(g.nz) * g.dz
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            ret.append({'x': X.ravel(order='C'), 'y': Y.ravel(order='C'),
                        'z': Z.ravel(order='C'),
                        'values': _values(g, slo2vel).ravel(order='C')})
        return ret
    raise ValueError("return_type must be 'class', 'array' or 'dict'")


def plot_model_grids(grids, *, section=(0, 0, 0), vmin=None, vmax=None,
                     slo2vel=True):
    """Plot YZ (X index), XZ (Y index) and XY (Z index) without rotation bugs.

    Returns (figure, axes); read_model_file still returns its original data.
    """
    import matplotlib.pyplot as plt
    sections = tuple(section)
    if len(sections) != 3 or all(v is None for v in sections):
        raise ValueError('section requires at least one index among (ix,iy,iz)')
    active = [(i, int(idx)) for i, idx in enumerate(sections) if idx is not None]
    if not grids:
        raise ValueError('No grids to plot')
    fig, axs = plt.subplots(len(grids), len(active), squeeze=False,
                            constrained_layout=True, figsize=(5 * len(active), 4 * len(grids)))
    for row, g in enumerate(grids):
        vals = _values(g, slo2vel)
        axes = [g.x_orig + np.arange(g.nx) * g.dx,
                g.y_orig + np.arange(g.ny) * g.dy,
                g.z_orig + np.arange(g.nz) * g.dz]
        half = [g.dx / 2, g.dy / 2, g.dz / 2]
        for col, (dim, idx) in enumerate(active):
            if idx < 0 or idx >= vals.shape[dim]:
                raise IndexError(f'Section index {idx} outside axis {dim} extent')
            ax = axs[row, col]
            if dim == 0:
                img = vals[idx, :, :].T
                u, v = 1, 2
                caption = f'YZ at X = {axes[0][idx]:g} km'
            elif dim == 1:
                img = vals[:, idx, :].T
                u, v = 0, 2
                caption = f'XZ at Y = {axes[1][idx]:g} km'
            else:
                img = vals[:, :, idx].T
                u, v = 0, 1
                caption = f'XY at depth = {axes[2][idx]:g} km'
            extent = (axes[u][0] - half[u], axes[u][-1] + half[u],
                      axes[v][0] - half[v], axes[v][-1] + half[v])
            image = ax.imshow(img, extent=extent, origin='lower',
                              interpolation='nearest', aspect='auto',
                              vmin=vmin, vmax=vmax)
            ax.set_xlabel(('X', 'Y', 'depth')[u] + ' (km)')
            ax.set_ylabel(('X', 'Y', 'depth')[v] + ' (km)')
            if v == 2:
                ax.invert_yaxis()
            else:
                ax.set_aspect('equal')
            ax.set_title(f'{g.type}: {caption}')
            fig.colorbar(image, ax=ax, label='km/s' if slo2vel and g.type == 'SLOW_LEN'
                         else g.type)
    return fig, axs
