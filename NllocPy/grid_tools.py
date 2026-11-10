# =============================================================================
# Author: Luigi Sante Zampa
# Affiliation: NllocPy project (institutional affiliation not specified)
# Created: 2026-11-03
# Last modified: 2026-11-10
#
# Purpose:
#   Native input/output for regular 3D NonLinLoc .hdr/.buf grid pairs.
#   Arrays use [X, Y, Z] indexing, with Z changing fastest on disk.
#   File origins and spacings are expressed in km. An intentionally thin
#   NLLGrid data object exposes the legacy instance API while independent
#   module-level functions carry out all computations and file operations.
# =============================================================================

# -----------------------------------------------------------------------------
# Import external libraries
# -----------------------------------------------------------------------------
from pathlib import Path

import numpy as np


# -----------------------------------------------------------------------------
# Shared grid constants
# -----------------------------------------------------------------------------
_GRID_DTYPES = {"FLOAT": "f4", "DOUBLE": "f8"}
_GRID_SUFFIXES = (".hdr", ".buf")
_PACKED_ANGLE_TYPES = ("ANGLE", "ANGLE2D")
_VELOCITY_TYPES = ("VELOCITY", "VELOCITY_METERS", "SLOWNESS", "SLOW_LEN")
_SUPPORTED_TRANSFORMS = ("NONE", "SIMPLE", "LAMBERT", "TRANS_MERC",
                         "AZIMUTHAL_EQUIDIST")


# =============================================================================
# Grid geometry, data properties and file naming
# =============================================================================

# -----------------------------------------------------------------------------
def grid_array_get(grid):
    """Return the grid's 3D NumPy array.

    Parameters
    ----------
    grid : NLLGrid
        Grid object.

    Returns
    -------
    numpy.ndarray or None
        Sample values, or None when no scalar grid is loaded.

    Examples
    --------
    >>> grid = NLLGrid()
    >>> grid_array_get(grid) is None
    True
    """
    output = grid._array
    return output


# -----------------------------------------------------------------------------
def grid_array_set(grid, values):
    """Assign and validate a grid array, updating its dimensions.

    Parameters
    ----------
    grid : NLLGrid
        Grid whose values will be replaced.
    values : array-like or None
        Numeric values with shape (nx, ny, nz).

    Returns
    -------
    None
        Array and geometry are updated in-place.

    Examples
    --------
    >>> grid = NLLGrid()
    >>> grid_array_set(grid, np.ones((2, 3, 4)))
    >>> grid.nx, grid.ny, grid.nz
    (2, 3, 4)
    """
    if values is None:
        grid._array = None
        return

    array = np.asarray(values)

    if array.ndim != 3 or 0 in array.shape:
        raise ValueError("A grid must be a nonempty 3D array [nx, ny, nz]")

    if not np.issubdtype(array.dtype, np.number):
        raise TypeError("Grid values must be numeric")

    grid._array = array
    grid._shape = array.shape


# -----------------------------------------------------------------------------
def grid_nx(grid):
    """Return the X-axis node count.

    Parameters
    ----------
    grid : NLLGrid
        Grid descriptor.

    Returns
    -------
    int
        Number of X nodes.

    Examples
    --------
    >>> grid_nx(NLLGrid())
    1
    """
    output = grid._shape[0]
    return output


# -----------------------------------------------------------------------------
def grid_ny(grid):
    """Return the Y-axis node count.

    Parameters
    ----------
    grid : NLLGrid
        Grid descriptor.

    Returns
    -------
    int
        Number of Y nodes.

    Examples
    --------
    >>> grid_ny(NLLGrid())
    1
    """
    output = grid._shape[1]
    return output


# -----------------------------------------------------------------------------
def grid_nz(grid):
    """Return the Z-axis node count.

    Parameters
    ----------
    grid : NLLGrid
        Grid descriptor.

    Returns
    -------
    int
        Number of Z nodes.

    Examples
    --------
    >>> grid_nz(NLLGrid())
    1
    """
    output = grid._shape[2]
    return output


# -----------------------------------------------------------------------------
def grid_root(name):
    """Strip a NonLinLoc header/buffer extension from a file root.

    Parameters
    ----------
    name : str or pathlib.Path
        Full grid root, header or binary buffer path.

    Returns
    -------
    str
        Grid root without .hdr or .buf extension.

    Examples
    --------
    >>> grid_root("layer.P.mod.buf")
    'layer.P.mod'
    """
    root = str(name)

    for extension in _GRID_SUFFIXES:
        if root.endswith(extension):
            root = root[:-len(extension)]
            break

    output = root
    return output


# -----------------------------------------------------------------------------
def grid_filename(grid, extension, basename=None):
    """Resolve a file path and optionally update the grid root.

    Parameters
    ----------
    grid : NLLGrid
        Grid descriptor.
    extension : str
        File extension including its leading dot.
    basename : str or pathlib.Path or None, optional
        Optional new grid root.

    Returns
    -------
    pathlib.Path
        Fully resolved filename using the supplied root.

    Examples
    --------
    >>> grid = NLLGrid()
    >>> grid_filename(grid, ".hdr", "layer.P.mod").name
    'layer.P.mod.hdr'
    """
    if basename is not None:
        grid.basename = grid_root(basename)

    if grid.basename is None:
        raise ValueError("Set basename or pass a file name")

    output = Path(grid.basename + extension)
    return output


# -----------------------------------------------------------------------------
def grid_dtype(grid):
    """Return the NumPy binary dtype selected by the NonLinLoc header.

    Parameters
    ----------
    grid : NLLGrid
        Grid with float_type and byteorder metadata.

    Returns
    -------
    numpy.dtype
        FLOAT as f4 or DOUBLE as f8, with the requested byte order.

    Examples
    --------
    >>> grid_dtype(NLLGrid()).itemsize
    4
    """
    try:
        dtype_code = _GRID_DTYPES[grid.float_type.upper()]
    except KeyError as exc:
        raise ValueError("float_type must be FLOAT or DOUBLE") from exc

    output = np.dtype(grid.byteorder + dtype_code)
    return output


# -----------------------------------------------------------------------------
def grid_get_extent(grid):
    """Return the six cell-edge coordinates for plotting.

    Parameters
    ----------
    grid : NLLGrid
        Grid descriptor with origin, spacing and node counts.

    Returns
    -------
    tuple of float
        (xmin, xmax, ymin, ymax, zmin, zmax), in km.

    Examples
    --------
    >>> grid_get_extent(NLLGrid())
    (-0.5, 0.5, -0.5, 0.5, -0.5, 0.5)
    """
    output = (
        grid.x_orig - grid.dx / 2,
        grid.x_orig + (grid.nx - 0.5) * grid.dx,
        grid.y_orig - grid.dy / 2,
        grid.y_orig + (grid.ny - 0.5) * grid.dy,
        grid.z_orig - grid.dz / 2,
        grid.z_orig + (grid.nz - 0.5) * grid.dz,
    )

    return output


# -----------------------------------------------------------------------------
def grid_get_xyz(grid, i, j, k):
    """Return the spatial coordinates of a grid node.

    Parameters
    ----------
    grid : NLLGrid
        Grid descriptor.
    i, j, k : int
        X, Y and Z node indices.

    Returns
    -------
    tuple of float
        (x, y, z) coordinates in km.

    Examples
    --------
    >>> grid_get_xyz(NLLGrid(), 1, 2, 3)
    (1.0, 2.0, 3.0)
    """
    output = (
        grid.x_orig + i * grid.dx,
        grid.y_orig + j * grid.dy,
        grid.z_orig + k * grid.dz,
    )

    return output


# =============================================================================
# Header metadata: projection and station information
# =============================================================================

# -----------------------------------------------------------------------------
def grid_read_transform(grid, parts):
    """Parse a NonLinLoc TRANSFORM statement into a grid descriptor.

    Parameters
    ----------
    grid : NLLGrid
        Descriptor to update.
    parts : list of str
        Whitespace-separated tokens of a TRANSFORM header statement.

    Returns
    -------
    None
        Projection fields are updated in-place.

    Examples
    --------
    >>> grid = NLLGrid()
    >>> grid_read_transform(grid, ["TRANSFORM", "SIMPLE", "LatOrig", "46"])
    >>> grid.orig_lat
    46.0
    """
    if len(parts) < 2:
        raise ValueError("Invalid TRANSFORM line")

    grid.proj_name = parts[1].upper()

    recognized = (
        "refellipsoid", "latorig", "longorig", "rotcw",
        "firststdparal", "secondstdparal",
    )
    keys = {}

    # The projection header stores keyword/value pairs after the name.
    for index, key in enumerate(parts[:-1]):
        normalized = key.lower()

        if normalized in recognized:
            keys[normalized] = parts[index + 1]

    grid.proj_ellipsoid = keys.get("refellipsoid", grid.proj_ellipsoid)

    field_names = (
        ("latorig", "orig_lat"),
        ("longorig", "orig_lon"),
        ("rotcw", "map_rot"),
        ("firststdparal", "first_std_paral"),
        ("secondstdparal", "second_std_paral"),
    )

    for keyword, attribute in field_names:
        if keyword in keys:
            setattr(grid, attribute, float(keys[keyword]))


# -----------------------------------------------------------------------------
def grid_get_transform_line(grid):
    """Format projection metadata into a valid NonLinLoc header line.

    Parameters
    ----------
    grid : NLLGrid
        Grid with projection fields and optional ellipsoid parameters.

    Returns
    -------
    str
        Complete TRANSFORM header statement.

    Examples
    --------
    >>> grid_get_transform_line(NLLGrid())
    'TRANSFORM NONE'
    """
    name = str(grid.proj_name).upper()

    if name not in _SUPPORTED_TRANSFORMS:
        raise ValueError(f"Unsupported transform for NonLinLoc header: {name}")

    if name == "NONE":
        output = "TRANSFORM NONE"
        return output

    tokens = ["TRANSFORM", name]

    if name != "SIMPLE":
        tokens.extend(("RefEllipsoid", str(grid.proj_ellipsoid)))

    tokens.extend((
        "LatOrig", f"{grid.orig_lat:.6f}",
        "LongOrig", f"{grid.orig_lon:.6f}",
    ))

    if name == "LAMBERT":
        tokens.extend((
            "FirstStdParal", f"{grid.first_std_paral:.6f}",
            "SecondStdParal", f"{grid.second_std_paral:.6f}",
        ))

    tokens.extend(("RotCW", f"{grid.map_rot:.6f}"))

    output = " ".join(tokens)
    return output


# -----------------------------------------------------------------------------
def grid_read_hdr_file(grid, basename=None):
    """Read a NonLinLoc ASCII grid header.

    Parameters
    ----------
    grid : NLLGrid
        Grid descriptor to populate.
    basename : str or pathlib.Path or None, optional
        File root, header filename or buffer filename.

    Returns
    -------
    NLLGrid
        Updated grid object, with any previously loaded array cleared.

    Examples
    --------
    >>> # grid.read_hdr_file("model/layer.P.mod.hdr")
    """
    path = grid_filename(grid, ".hdr", basename)
    lines = path.read_text(encoding="utf-8").splitlines()

    if not lines or len(lines[0].split()) < 10:
        raise ValueError(f"Invalid NonLinLoc header: {path}")

    values = lines[0].split()

    try:
        grid._shape = tuple(map(int, values[:3]))
        grid.x_orig, grid.y_orig, grid.z_orig = map(float, values[3:6])
        grid.dx, grid.dy, grid.dz = map(float, values[6:9])
    except ValueError as exc:
        raise ValueError(f"Invalid grid dimensions or geometry: {path}") from exc

    spacings = (grid.dx, grid.dy, grid.dz)
    invalid_shape = any(count <= 0 for count in grid._shape)
    invalid_spacing = any(not np.isfinite(size) or size <= 0 for size in spacings)

    if invalid_shape or invalid_spacing:
        raise ValueError(f"Non-positive grid dimensions or spacing: {path}")

    grid.type = values[9].upper()
    grid.float_type = values[10].upper() if len(values) > 10 else "FLOAT"

    # Check precision and clear stale in-memory data before parsing metadata.
    _ = grid_dtype(grid)
    grid._array = None
    grid.station = None
    grid.proj_name = "NONE"

    for line in lines[1:]:
        parts = line.split()

        if not parts or parts[0].startswith("#"):
            continue

        if parts[0].upper() in ("TRANS", "TRANSFORM"):
            grid_read_transform(grid, parts)
        elif len(parts) >= 4:
            grid.station = parts[0]
            grid.sta_x, grid.sta_y, grid.sta_z = map(float, parts[1:4])

    output = grid
    return output


# -----------------------------------------------------------------------------
def grid_write_hdr_file(grid, basename=None):
    """Write the ASCII header of a regular NonLinLoc 3D grid.

    Parameters
    ----------
    grid : NLLGrid
        Grid containing geometry, data type and projection metadata.
    basename : str or pathlib.Path or None, optional
        New grid root; otherwise the object's basename is used.

    Returns
    -------
    pathlib.Path
        Path to the generated .hdr file.

    Examples
    --------
    >>> # grid.write_hdr_file("model/layer.P.mod")
    """
    path = grid_filename(grid, ".hdr", basename)
    geometry = (
        grid.x_orig, grid.y_orig, grid.z_orig,
        grid.dx, grid.dy, grid.dz,
    )

    if not np.all(np.isfinite(geometry)):
        raise ValueError("Grid origin/spacing must be finite")

    if min(grid.dx, grid.dy, grid.dz) <= 0:
        raise ValueError("Grid spacing must be positive")

    _ = grid_dtype(grid)

    content = (
        f"{grid.nx} {grid.ny} {grid.nz} "
        f"{grid.x_orig:.6f} {grid.y_orig:.6f} {grid.z_orig:.6f} "
        f"{grid.dx:.6f} {grid.dy:.6f} {grid.dz:.6f} "
        f"{grid.type} {grid.float_type.upper()}\n"
    )

    if grid.station is not None:
        content += (
            f"{grid.station} {grid.sta_x:.6f} "
            f"{grid.sta_y:.6f} {grid.sta_z:.6f}\n"
        )

    content += grid_get_transform_line(grid) + "\n"
    path.write_text(content, encoding="utf-8")

    output = path
    return output


# =============================================================================
# Binary I/O and conversion between velocity and slowness
# =============================================================================

# -----------------------------------------------------------------------------
def grid_write_buf_file(grid, basename=None):
    """Write the scalar grid buffer using NonLinLoc XYZ index order.

    Parameters
    ----------
    grid : NLLGrid
        Grid with a finite numeric array indexed [X, Y, Z].
    basename : str or pathlib.Path or None, optional
        New grid root.

    Returns
    -------
    pathlib.Path
        Path to the generated .buf file.

    Raises
    ------
    ValueError
        Missing/non-finite array or nonpositive velocity/slowness fields.
    NotImplementedError
        Packed angle-grid writing is intentionally not supported.

    Examples
    --------
    >>> # grid.write_buf_file("model/layer.P.mod")
    """
    if grid.type in _PACKED_ANGLE_TYPES:
        raise NotImplementedError("Writing packed angle grids is not supported")

    if grid.array is None:
        raise ValueError("No grid array to write")

    if not np.all(np.isfinite(grid.array)):
        raise ValueError("Grid contains NaN/Inf; cannot write a valid model")

    if grid.type in _VELOCITY_TYPES and np.any(grid.array <= 0):
        raise ValueError("Velocity/slowness values must be positive")

    path = grid_filename(grid, ".buf", basename)
    array = np.asarray(grid.array, dtype=grid_dtype(grid), order="C")
    array.tofile(path)

    output = path
    return output


# -----------------------------------------------------------------------------
def grid_read_buf_file(grid, basename=None):
    """Read a binary scalar or packed-angle NonLinLoc grid buffer.

    Parameters
    ----------
    grid : NLLGrid
        Descriptor previously loaded from .hdr.
    basename : str or pathlib.Path or None, optional
        File root to load.

    Returns
    -------
    NLLGrid
        Grid containing either the scalar array or the angle/dip/quality
        attributes for packed angle grids.

    Examples
    --------
    >>> # grid.read_hdr_file("model/layer.P.mod")
    >>> # grid.read_buf_file()
    """
    path = grid_filename(grid, ".buf", basename)
    size = grid.nx * grid.ny * grid.nz

    if grid.type in _PACKED_ANGLE_TYPES:
        if grid.float_type != "FLOAT":
            raise ValueError("Packed angle grids require FLOAT precision")

        raw = np.fromfile(path, dtype=np.dtype(grid.byteorder + "u2"))

        if raw.size != 2 * size:
            raise ValueError(f"Unexpected size of packed angle grid: {path}")

        # Each cell stores two packed uint16 values (dip/quality and azimuth).
        packed = raw.reshape((-1, 2))
        grid.quality = (packed[:, 0] % 16).reshape(grid._shape)
        grid.dip = (packed[:, 0] // 16 / 10.0).reshape(grid._shape)
        grid.azimuth = (packed[:, 1] / 10.0).reshape(grid._shape)
        grid.dip = np.where(grid.quality == 0, np.nan, grid.dip)
        grid.azimuth = np.where(grid.quality == 0, np.nan, grid.azimuth)
        grid._array = None
    else:
        raw = np.fromfile(path, dtype=grid_dtype(grid))

        if raw.size != size:
            raise ValueError(
                f"Grid contains {raw.size} values; expected {size}: {path}"
            )

        grid.array = raw.reshape(grid._shape, order="C")

    output = grid
    return output


# -----------------------------------------------------------------------------
def velocity_to_slow_len(velocity, step_km):
    """Convert a velocity grid into NonLinLoc SLOW_LEN (seconds per cell).

    Parameters
    ----------
    velocity : numpy.ndarray
        Three-dimensional velocity grid in km/s, indexed [X, Y, Z].
    step_km : float
        Uniform grid cell size in km.

    Returns
    -------
    numpy.ndarray
        SLOW_LEN grid with the same XYZ layout.

    Raises
    ------
    ValueError
        Velocity/spacing is zero, negative, non-finite or invalid in shape.

    Examples
    --------
    >>> velocity_to_slow_len(np.array([[[5.0]]]), 1.0).item()
    0.2
    """
    values = np.asarray(velocity, dtype=float)

    if not np.isfinite(step_km) or step_km <= 0:
        raise ValueError("step_km must be finite and positive")

    if values.ndim != 3:
        raise ValueError("Velocity must be a 3D array of finite positive km/s values")

    if not np.all(np.isfinite(values)) or np.any(values <= 0):
        raise ValueError("Velocity must be a 3D array of finite positive km/s values")

    output = step_km / values
    return output


# =============================================================================
# Grid container: state only, with operations defined above the class
# =============================================================================

# -----------------------------------------------------------------------------
class NLLGrid:
    """Regular 3D NonLinLoc grid, with a backward-compatible instance API.

    This class only initializes persistent state. Its data properties and
    operations reference independent, testable module functions above.

    Parameters
    ----------
    basename : str or pathlib.Path or None, optional
        Existing .hdr/.buf file root, or None for an empty grid.
    byteorder : {'=', '<', '>'}, default '='
        Native, little-endian or big-endian binary representation.

    Examples
    --------
    >>> grid = NLLGrid()
    >>> grid.array = np.ones((2, 3, 4))
    >>> grid.nx, grid.ny, grid.nz
    (2, 3, 4)
    """

    # No computational methods are defined inside the class.
    array = property(grid_array_get, grid_array_set)
    nx = property(grid_nx)
    ny = property(grid_ny)
    nz = property(grid_nz)
    _root = staticmethod(grid_root)
    _filename = grid_filename
    _dtype = property(grid_dtype)
    _read_transform = grid_read_transform
    get_transform_line = grid_get_transform_line
    read_hdr_file = grid_read_hdr_file
    write_hdr_file = grid_write_hdr_file
    read_buf_file = grid_read_buf_file
    write_buf_file = grid_write_buf_file
    get_extent = grid_get_extent
    get_xyz = grid_get_xyz

    def __init__(self, basename=None, *, byteorder="="):
        """Initialize metadata and optionally load a grid from disk.

        Parameters
        ----------
        basename : str or pathlib.Path or None
            Existing grid root, .hdr or .buf filename.
        byteorder : str
            Binary byte order. Use '=' for files written on this machine.

        Examples
        --------
        >>> NLLGrid().float_type
        'FLOAT'
        """
        if byteorder not in ("=", "<", ">"):
            raise ValueError("byteorder must be '=', '<' or '>'")

        # Initialize grid shape, geometry, binary representation and metadata.
        self.byteorder = byteorder
        self.basename = None
        self._shape = (1, 1, 1)
        self._array = None
        self.x_orig = self.y_orig = self.z_orig = 0.0
        self.dx = self.dy = self.dz = 1.0
        self.type = "VELOCITY"
        self.float_type = "FLOAT"
        self.proj_name = "NONE"
        self.orig_lat = self.orig_lon = self.map_rot = 0.0
        self.proj_ellipsoid = "WGS-84"
        self.first_std_paral = self.second_std_paral = 0.0
        self.station = None
        self.sta_x = self.sta_y = self.sta_z = None
        self.azimuth = self.dip = self.quality = None

        # Keep legacy constructor behavior for existing NllocPy callers.
        if basename is not None:
            grid_read_hdr_file(self, basename)
            grid_read_buf_file(self)
