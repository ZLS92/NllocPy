"""Read and write regular NonLinLoc grids without third-party grid libraries.

Each grid consists of an ASCII .hdr file and a binary .buf file.
Data are indexed [x, y, z]; Z varies fastest in the binary buffer.
Coordinates and spacings are in km for local NonLinLoc grids.

This is an independent, limited interface tailored to NllocPy's existing
create_model_file and read_model_file workflows. It does not perform
geographical coordinate transformations; it preserves their parameters.
"""

from pathlib import Path

import numpy as np


class NLLGrid:
    """Regular 3D NonLinLoc grid, optionally initialized from a file root.

    basename can end in .hdr or .buf. The byteorder argument
    is '=' (native), '<' (little-endian) or '>' (big-endian).
    NonLinLoc generally writes buffers in the machine's native byte order.
    """

    def __init__(self, basename=None, *, byteorder="="):
        if byteorder not in ("=", "<", ">"):
            raise ValueError("byteorder must be '=', '<' or '>'")
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
        if basename is not None:
            self.read_hdr_file(basename)
            self.read_buf_file()

    @property
    def array(self):
        return self._array

    @array.setter
    def array(self, values):
        if values is None:
            self._array = None
            return
        a = np.asarray(values)
        if a.ndim != 3 or 0 in a.shape:
            raise ValueError("A grid must be a nonempty 3D array [nx, ny, nz]")
        if not np.issubdtype(a.dtype, np.number):
            raise TypeError("Grid values must be numeric")
        self._array = a
        self._shape = a.shape

    @property
    def nx(self):
        return self._shape[0]

    @property
    def ny(self):
        return self._shape[1]

    @property
    def nz(self):
        return self._shape[2]

    @staticmethod
    def _root(name):
        name = str(name)
        for ext in (".hdr", ".buf"):
            if name.endswith(ext):
                return name[:-len(ext)]
        return name

    def _filename(self, ext, basename=None):
        if basename is not None:
            self.basename = self._root(basename)
        if self.basename is None:
            raise ValueError("Set basename or pass a file name")
        return Path(self.basename + ext)

    @property
    def _dtype(self):
        types = {"FLOAT": "f4", "DOUBLE": "f8"}
        try:
            return np.dtype(self.byteorder + types[self.float_type.upper()])
        except KeyError as exc:
            raise ValueError("float_type must be FLOAT or DOUBLE") from exc

    def read_hdr_file(self, basename=None):
        path = self._filename(".hdr", basename)
        lines = path.read_text(encoding="utf-8").splitlines()
        if not lines or len(lines[0].split()) < 10:
            raise ValueError(f"Invalid NonLinLoc header: {path}")
        v = lines[0].split()
        try:
            self._shape = tuple(map(int, v[:3]))
            self.x_orig, self.y_orig, self.z_orig = map(float, v[3:6])
            self.dx, self.dy, self.dz = map(float, v[6:9])
        except ValueError as exc:
            raise ValueError(f"Invalid grid dimensions or geometry: {path}") from exc
        if any(n <= 0 for n in self._shape) or any(
            not np.isfinite(h) or h <= 0 for h in (self.dx, self.dy, self.dz)
        ):
            raise ValueError(f"Non-positive grid dimensions or spacing: {path}")
        self.type = v[9].upper()
        self.float_type = v[10].upper() if len(v) > 10 else "FLOAT"
        _ = self._dtype
        self._array = None
        self.station = None
        self.proj_name = "NONE"
        for line in lines[1:]:
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            if parts[0].upper() in ("TRANS", "TRANSFORM"):
                self._read_transform(parts)
            elif len(parts) >= 4:
                self.station = parts[0]
                self.sta_x, self.sta_y, self.sta_z = map(float, parts[1:4])
        return self

    def _read_transform(self, parts):
        if len(parts) < 2:
            raise ValueError("Invalid TRANSFORM line")
        self.proj_name = parts[1].upper()
        keys = {
            key.lower(): parts[i + 1]
            for i, key in enumerate(parts[:-1])
            if key.lower() in (
                "refellipsoid", "latorig", "longorig", "rotcw",
                "firststdparal", "secondstdparal"
            )
        }
        self.proj_ellipsoid = keys.get("refellipsoid", self.proj_ellipsoid)
        for key, attr in (
            ("latorig", "orig_lat"), ("longorig", "orig_lon"),
            ("rotcw", "map_rot"), ("firststdparal", "first_std_paral"),
            ("secondstdparal", "second_std_paral")
        ):
            if key in keys:
                setattr(self, attr, float(keys[key]))

    def get_transform_line(self):
        name = str(self.proj_name).upper()
        if name == "NONE":
            return "TRANSFORM NONE"
        if name not in ("SIMPLE", "LAMBERT", "TRANS_MERC", "AZIMUTHAL_EQUIDIST"):
            raise ValueError(f"Unsupported transform for NonLinLoc header: {name}")
        p = ["TRANSFORM", name]
        if name != "SIMPLE":
            p += ["RefEllipsoid", str(self.proj_ellipsoid)]
        p += ["LatOrig", f"{self.orig_lat:.6f}",
              "LongOrig", f"{self.orig_lon:.6f}"]
        if name == "LAMBERT":
            p += ["FirstStdParal", f"{self.first_std_paral:.6f}",
                  "SecondStdParal", f"{self.second_std_paral:.6f}"]
        p += ["RotCW", f"{self.map_rot:.6f}"]
        return " ".join(p)

    def write_hdr_file(self, basename=None):
        path = self._filename(".hdr", basename)
        if any(not np.isfinite(v) for v in (
            self.x_orig, self.y_orig, self.z_orig, self.dx, self.dy, self.dz
        )) or min(self.dx, self.dy, self.dz) <= 0:
            raise ValueError("Grid origin/spacing must be finite, with positive spacing")
        _ = self._dtype
        first = (
            f"{self.nx} {self.ny} {self.nz} "
            f"{self.x_orig:.6f} {self.y_orig:.6f} {self.z_orig:.6f} "
            f"{self.dx:.6f} {self.dy:.6f} {self.dz:.6f} "
            f"{self.type} {self.float_type.upper()}\n"
        )
        if self.station is not None:
            first += (f"{self.station} {self.sta_x:.6f} "
                      f"{self.sta_y:.6f} {self.sta_z:.6f}\n")
        first += self.get_transform_line() + "\n"
        path.write_text(first, encoding="utf-8")
        return path

    def write_buf_file(self, basename=None):
        if self.type in ("ANGLE", "ANGLE2D"):
            raise NotImplementedError("Writing packed angle grids is not supported")
        if self.array is None:
            raise ValueError("No grid array to write")
        if not np.all(np.isfinite(self.array)):
            raise ValueError("Grid contains NaN/Inf; cannot write a valid model")
        if self.type in ("VELOCITY", "VELOCITY_METERS", "SLOWNESS", "SLOW_LEN"):
            if np.any(self.array <= 0):
                raise ValueError("Velocity/slowness values must be positive")
        path = self._filename(".buf", basename)
        np.asarray(self.array, dtype=self._dtype, order="C").tofile(path)
        return path

    def read_buf_file(self, basename=None):
        path = self._filename(".buf", basename)
        size = self.nx * self.ny * self.nz
        if self.type in ("ANGLE", "ANGLE2D"):
            if self.float_type != "FLOAT":
                raise ValueError("Packed angle grids require FLOAT precision")
            raw = np.fromfile(path, dtype=np.dtype(self.byteorder + "u2"))
            if raw.size != 2 * size:
                raise ValueError(f"Unexpected size of packed angle grid: {path}")
            packed = raw.reshape((-1, 2))
            self.quality = (packed[:, 0] % 16).reshape(self._shape)
            self.dip = (packed[:, 0] // 16 / 10.0).reshape(self._shape)
            self.azimuth = (packed[:, 1] / 10.0).reshape(self._shape)
            self.dip = np.where(self.quality == 0, np.nan, self.dip)
            self.azimuth = np.where(self.quality == 0, np.nan, self.azimuth)
            self._array = None
        else:
            raw = np.fromfile(path, dtype=self._dtype)
            if raw.size != size:
                raise ValueError(
                    f"Grid contains {raw.size} values; expected {size}: {path}")
            self.array = raw.reshape(self._shape)
        return self

    def get_extent(self):
        """Extent compatible with NllocPy's existing plotting convention."""
        return (self.x_orig - self.dx / 2,
                self.x_orig + (self.nx - 0.5) * self.dx,
                self.y_orig - self.dy / 2,
                self.y_orig + (self.ny - 0.5) * self.dy,
                self.z_orig - self.dz / 2,
                self.z_orig + (self.nz - 0.5) * self.dz)

    def get_xyz(self, i, j, k):
        return (self.x_orig + i * self.dx,
                self.y_orig + j * self.dy,
                self.z_orig + k * self.dz)


def velocity_to_slow_len(velocity, step_km):
    """Convert a velocity grid (km/s) to slowness-times-cell-size (s)."""
    v = np.asarray(velocity, dtype=float)
    if not np.isfinite(step_km) or step_km <= 0:
        raise ValueError("step_km must be finite and positive")
    if v.ndim != 3 or not np.all(np.isfinite(v)) or np.any(v <= 0):
        raise ValueError("Velocity must be a 3D array of finite positive km/s values")
    return step_km / v
