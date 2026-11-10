# NllocPy Python coding conventions

**Maintainer:** Luigi Sante Zampa  
**Project affiliation:** NllocPy (no institutional affiliation asserted)  
**Established:** 2026-11-04

This repository uses consistent structure for scientific Python modules. The goal
is readability, reproducibility and stable behavior during numerical refactors.

## Editable development installation

Until NllocPy has a stable release, **all installations must be editable**.
Use a dedicated environment (conda is recommended for GDAL/pyproj), clone
the repository and run `python -m pip install -e .` at its root.

Do not use `pip install .`, `conda build conda-recipe` or
`conda install --use-local nllocpy` for current development.
The conda recipe is reserved for a possible future non-editable release.
See [README.md](README.md#-installation-editable-mode-only) for
the recommended environment and import instructions.

## Module layout

Each new Python file follows this order:

1. File header with author, affiliation, creation date, last modification date,
   purpose and expected usage. Mark unknown historical dates as unknown instead
   of inventing dates.
2. A single, clearly marked external-import block at the top. Prefer standard
   library imports first, then NumPy/SciPy/Matplotlib/GDAL/pyproj, then local imports.
   Optional dependencies should be checked in this block, not imported inside
   processing functions.
3. Shared constants, physical parameters, geodetic configuration and
   immutable settings immediately below imports.
4. Independent functional groups separated by \`# ===...\` comments and
   every function/class preceded by \`# ---...\`.
5. Complete NumPy-style function docstrings: description, parameters, returns,
   relevant exceptions and an example of usage. Classes document their
   initialization parameters and examples.
6. Explicit names and short, well-commented logical steps, rather than long
   chains of expressions.
7. Multiple returned objects assigned to one tuple/list first:
   \`output = (x, y, z)\` followed by \`return output\`.
8. No nested function declarations, except trivial, necessary lambda callbacks.
9. Classes should initialize state only. General algorithms live in separate
   module functions (preferably above the class). Where compatibility with
   \`NLLGrid.write_buf_file()\` is needed, the class aliases those functions.
10. Related functions are placed together, with group descriptions.
11. Use established scientific libraries rather than reinventing their
    algorithms.
12. Before adding a helper, inspect existing NllocPy modules and imported
    scientific packages for an equivalent operation.
13. Internal table-like data use dictionaries of aligned 1D NumPy arrays.
    Use \`NllocPy.table_tools.to_field_dict\` at input boundaries and
    \`restore_field_table\` at output boundaries if a caller requests pandas.
    Pandas is optional, for analysis and presentation in notebooks.
14. Prefer NumPy array expressions over nested Python loops whenever practical.
    Retain loops when required for readability, plotting multiple independent
    slices, or stateful I/O operations.

## Data conversion example

\`\`\`python
from NllocPy.table_tools import to_field_dict, restore_field_table

# Input can be a dictionary or a pandas DataFrame.
model_fields = to_field_dict(model_input)
model_fields["vp"] = model_fields["vp"] * 1.01

# Return the same table family as the input when requested.
result = restore_field_table(model_fields, template=model_input)
\`\`\`

## Grid conventions to preserve during any refactor

- Grid axis order is strictly \`[X, Y, Z]\`.
- Z is the fastest varying index in \`.buf\` files.
- The caller explicitly selects \`z_positive="up"\` or \`"down"\`.
- Native NonLinLoc grids use kilometers and seconds.
- \`SLOW_LEN\` is grid-cell length divided by velocity, not velocity itself.
- YZ/XZ plots show depth increasing down; XY plots show Y increasing up.
- Exact regular input grids are reordered without interpolation.
- Invalid and duplicate source nodes are rejected rather than silently averaged.

## Validation and migration scope

Run the strict checker for modern modules:

\`\`\`bash
python tools/check_python_style.py --strict
python tools/check_python_style.py --all  # Advisory full-package report.
MPLBACKEND=Agg python -m pytest tests/test_grid_io.py tests/test_model_tools.py tests/test_table_tools.py
\`\`\`

GitHub Actions additionally compiles \`Grid2Time\` and checks an end-to-end
3D travel-time calculation.

The long-standing legacy files \`nlloc.py\`, \`utils.py\`, \`plot.py\`,
\`raster_tools.py\` and \`shp_tools.py\` have standardized file headers and
import section labels, but their thousands of lines of algorithms and
interactive callback functions have **not** yet been rewritten wholesale.
A full migration of these files should follow comprehensive regression tests
for each public method; style refactors must not change scientific results.

The C source in \`NllocPy/src/\` belongs to upstream NonLinLoc and is **outside**
these Python-specific conventions.
