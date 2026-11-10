# NllocPy

**NllocPy** is a Python wrapper for **NonLinLoc**, an open-source software developed by Anthony Lomax for probabilistic earthquake location.  
This package provides Python tools to automate input file preparation, manage control parameters, run NonLinLoc, and analyze or visualize results.

---

## 🚧 Project Status

**This project is currently under development.**  
Features, interfaces, and file structures may change. Use it at your own risk and feel free to contribute!


---

## ⚙️ Features

- **Input File Automation**: Create NonLinLoc control files, velocity models and travel-time grids.
- **Native Grid I/O**: Read/write NonLinLoc `.hdr`/`.buf` grids via `NllocPy.grid_io`, with no dependency on `nllgrid`.
- **Utility Scripts**: Python tools for handling raster DEMs and geospatial shapefiles.
- **Integrated NonLinLoc Core**: Includes the **original NonLinLoc C source code** (`src/`) for easy compilation.
- **Result Visualization**: Python functions for plotting and interactive analysis.

---

## 📦 Installation

Create a dedicated Python environment, then clone the repository:

```bash
git clone https://github.com/ZLS92/NllocPy.git
cd NllocPy
```

### Install with pip

```bash
pip install .
```

### Install with conda

```bash
conda install conda-build
conda build conda-recipe
conda install --use-local nllocpy
```

---

## Native 3D model grids

`NllocPy.grid_io.NLLGrid` writes a regular 3D grid as an ASCII `.hdr` and binary `.buf` file. Arrays are indexed `[x, y, z]` (Z fastest in the binary file); local coordinates are in km. For velocity models used by `Grid2Time`, the grid type must be `SLOW_LEN` and values must be `step_km / velocity_km_s`. NaN, infinite and non-positive velocities are rejected instead of producing unusable model files.

```python
import numpy as np
from NllocPy.grid_io import NLLGrid, velocity_to_slow_len

vp = np.full((20, 15, 12), 5.5)  # km/s
grid = NLLGrid()
grid.array = velocity_to_slow_len(vp, step_km=0.5)
grid.dx = grid.dy = grid.dz = 0.5
grid.x_orig, grid.y_orig, grid.z_orig = -5, -3, 0
grid.type = "SLOW_LEN"
grid.proj_name = "NONE"  # or SIMPLE with orig_lat/orig_lon
grid.basename = "model/layer.P.mod"  # directory must exist
grid.write_hdr_file()
grid.write_buf_file()

restored = NLLGrid("model/layer.P.mod")
assert restored.array.shape == vp.shape
```

The native implementation is intentionally limited to regular grids. Reading `FLOAT`/`DOUBLE` model and travel-time buffers and packed take-off angle buffers is supported; writing packed angle buffers is not. The header stores projection parameters but does not calculate coordinate transformations. For cross-platform files, the optional `byteorder` argument controls binary endianness. Test grid I/O with `python -m pytest tests/test_grid_io.py` (requires pytest).

## Unordered XYZ velocity models and correctly oriented plots

Input is a **dictionary of 1D arrays**, with one row per velocity sample:
`{'x': x, 'y': y, 'z': z, 'vp': vp, 'vs': vs}`.
Rows may be in any order, even randomly shuffled. The importer always sorts
`(x,y,z,vp,vs)` together, checks for duplicate or missing nodes, and writes a
regular 3D array indexed **[X, Y, depth]**. The NonLinLoc binary buffer stores
**Z fastest, followed by Y, then X**. Distances are in **km** and velocities
are in **km/s**.

**Vertical sign must be explicit**. By default `z_positive='up'` interprets
input Z as elevation and converts to NonLinLoc depth positive down. For
already-positive depths, specify `z_positive='down'`. Sorting is automatic;
vertical sign is *not* guessed from the values.

```python
import numpy as np
from NllocPy.nlloc import create_model_file, read_model_file

# x, y, z, vp and vs are equally-sized 1D arrays; rows may be unordered.
model = {'x': x, 'y': y, 'z': z, 'vp': vp, 'vs': vs}

p_file, s_file, grid = create_model_file(
    model=model, path="model", name="layer",
    step=0.5, z_positive="down",
    proj_name="NONE", validate_output=True
)
# [YZ at X index 5, XZ at Y index 10, XY at depth index 3]
grids = read_model_file(p_file, plot=True, section=(5, 10, 3))
# Returned dictionary uses NonLinLoc depth-positive-down.
points = read_model_file(p_file, return_type="dict")[0]
```

Complete regular input grids are **reordered without interpolation**. If nodes
are missing or you request a different grid resolution, set
`interpolate=True` and choose `method="nearest"` or `"linear"`; SciPy is
only needed for this optional resampling step. Duplicate coordinates are
rejected rather than silently averaged. All generated velocity grids have
uniform `dx=dy=dz`; non-isotropic grids require resampling.

The XY plot shows X increasing rightward and Y upward. The XZ and YZ
sections put **depth increasing downward**, using actual grid-node extents
and the corresponding unrotated/transposed NumPy slices. `SLOW_LEN`
values are displayed as velocity (km/s), and time grids stay in seconds.
Model read/write validation is enabled by default.

Tests for random input order, P/S round trips, C-compatible buffer offsets,
vertical sign, crop extents and plot pixel orientation run under GitHub
Actions. Execute locally with:

```bash
MPLBACKEND=Agg python -m pytest tests/test_grid_io.py tests/test_model_tools.py -q
```

## Python coding conventions and tabular data

NllocPy follows the maintainability rules in [CONTRIBUTING.md](CONTRIBUTING.md).
New grid, model and table modules contain documented standalone functions,
state-only classes and explicit NumPy-based array handling.

For dictionary/DataFrame conversion use:
```python
from NllocPy.table_tools import to_field_dict, restore_field_table

fields = to_field_dict(input_data)  # Mapping or optional pandas DataFrame.
# ... process the aligned 1D fields in a dictionary ...
result = restore_field_table(fields, template=input_data)
```

Pandas is **optional**; numerical methods keep arrays inside dictionaries.
The native 3D model builder accepts either the 1D field dictionary or a
DataFrame with identical columns, preserving the association of XYZ and
velocity samples. Legacy modules are being migrated incrementally to avoid
changes to established scientific calculations.

## 👤 Attribution

- **NllocPy** was developed and is maintained by **Dr. Luigi Sante Zampa**.
- The **C code in `src/`** is the original **NonLinLoc** software developed by Anthony Lomax et al.,.  
  NonLinLoc is available at [ut-beg-texnet/NonLinLoc](https://github.com/ut-beg-texnet/NonLinLoc).
- **NllocPy** only provides a Python automation layer and utility scripts for easier integration.
- The [NLLGrid](https://github.com/claudiodsf/nllgrid) library by Claudio Satriano and collaborators is a useful reference for the NonLinLoc grid format, but NllocPy does not import or depend on it.
- If you use this wrapper, please **also cite NonLinLoc** in your scientific work:

  > Lomax, A., Virieux, J., Volant, P., & Berge, C. (2000).  
  > *Probabilistic earthquake location in 3D and layered models: Introduction of a Metropolis-Gibbs method and comparison with linear locations.*  
  > In **Advances in Seismic Event Location**, Springer, pp. 101–134.  
  > DOI: [10.1007/978-94-015-9536-0_5](https://doi.org/10.1007/978-94-015-9536-0_5)