from pathlib import Path
import shutil

import numpy as np
from numpy.testing import assert_allclose
import pytest

from PyECLOUD import __version__
from PyECLOUD.buildup_simulation import BuildupSimulation


@pytest.mark.parametrize('track_method', ['Boris', 'BorisMultipole'])
def test_short_buildup(tmp_path, monkeypatch, track_method):
    inputs = Path(__file__).parent / 'fixtures' / 'buildup'
    for source in inputs.iterdir():
        shutil.copy2(source, tmp_path / source.name)
    monkeypatch.chdir(tmp_path)
    np.random.seed(12345)
    sim = BuildupSimulation(pyecl_input_folder=str(tmp_path), track_method=track_method)
    particles = sim.cloud_list[0].MP_e
    initial_count = particles.N_mp
    assert initial_count > 0
    initial_charge = particles.nel_mp[:initial_count].sum()
    initial_x = particles.x_mp[:initial_count].copy()
    sim.run(t_end_sim=5e-11)
    assert sim.beamtim.tt_curr >= 5e-11
    assert particles.N_mp == initial_count
    assert_allclose(particles.nel_mp[:particles.N_mp].sum(), initial_charge)
    assert np.any(particles.x_mp[:particles.N_mp] != initial_x)
    for component in ('x_mp', 'y_mp', 'vx_mp', 'vy_mp', 'vz_mp'):
        assert np.all(np.isfinite(getattr(particles, component)[:particles.N_mp]))
    assert np.all(np.isfinite(sim.spacech_ele.phi))
    log = (tmp_path / 'logfile.txt').read_text()
    assert f'PyECLOUD Version {__version__}' in log
