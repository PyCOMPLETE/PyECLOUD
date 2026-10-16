import sys
sys.path.append('../../../')

import PyPARIS_sim_class.Simulation as sim_mod
import PyPARIS.communication_helpers as ch

ring = sim_mod.get_serial_CPUring(init_sim_objects_auto=False)

ring.sim_content.pp.enable_arc_quad = False
ring.sim_content.init_all()
slices = ring.sim_content.init_master()

bunch = ring.sim_content.bunch

import h5py
with h5py.File('generated_bunch.h5', 'w') as fid:
    fid['bunch'] = ch.beam_2_buffer(bunch)
