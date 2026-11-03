"""End-to-end: feed a shuffled native velocity model to the NonLinLoc C engine.

Run with NLL_GRID2TIME=<absolute path to compiled Grid2Time> pytest ...
The test is skipped if no Grid2Time executable is available.
"""
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest

from NllocPy.model_tools import create_model_file
from NllocPy.grid_io import NLLGrid


def test_grid2time_on_native_shuffled_velocity_model(tmp_path):
    exe = os.environ.get('NLL_GRID2TIME')
    if not exe:
        exe = Path(__file__).resolve().parents[1] / 'NllocPy/src/bin/Grid2Time'
    exe = Path(exe).resolve()
    if not exe.is_file():
        pytest.skip('Compile NonLinLoc Grid2Time or set NLL_GRID2TIME')
    x, y, z = [np.arange(n, dtype=float) for n in (9, 9, 7)]
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    rng = np.random.default_rng(91)
    ix = rng.permutation(X.size)
    model = {'x': X.ravel()[ix], 'y': Y.ravel()[ix],
             'z': Z.ravel()[ix],
             'vp': np.full(X.size, 5.0)[ix]}
    create_model_file(model, step=1, z_positive='down',
                      proj_name='NONE', path=tmp_path / 'model')
    (tmp_path / 'time').mkdir()
    params = tmp_path / 'grid2time.in'
    params.write_text('''CONTROL 1 12345
TRANS NONE
GTFILES model/layer time/layer P 0
GTMODE GRID3D ANGLES_NO
GTSRCE STA XYZ 4.0 4.0 2.0 0.0
GT_PLFD 1.0e-3 0
''')
    run = subprocess.run([str(exe), str(params)], cwd=tmp_path,
                         capture_output=True, text=True, timeout=120)
    assert run.returncode == 0, (run.stdout + '\n' + run.stderr)[-5000:]
    root = tmp_path / 'time/layer.P.time.STA'
    assert Path(str(root) + '.buf').exists(), (
        'Grid2Time did not generate travel-time file; output files:\n'
        + '\n'.join(str(p.relative_to(tmp_path)) for p in tmp_path.rglob('*'))
        + '\nSTDOUT:\n' + run.stdout + '\nSTDERR:\n' + run.stderr
    )
    t = NLLGrid(str(root))
    assert t.type == 'TIME'
    assert t.array.shape == (9, 9, 7)
    assert np.isfinite(t.array).all()
    assert t.array[4, 4, 2] == pytest.approx(0, abs=0.005)
    assert t.array[5, 4, 2] == pytest.approx(0.2, abs=0.10)
    assert t.array[8, 4, 2] > t.array[5, 4, 2]
