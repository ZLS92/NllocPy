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

Clone the repository and install locally with pip:

```bash
git clone https://github.com/ZLS92/NllocPy.git
cd NllocPy
pip install .
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