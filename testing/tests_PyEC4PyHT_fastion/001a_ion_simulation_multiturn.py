import time
import numpy as np
from scipy.constants import c, e, m_p
import pickle
import sys
sys.path.append("../../../")
sys.path.append("../../../PyHEADTAIL")
from PyHEADTAIL.particles.slicing import UniformBinSlicer

print 'Start initialization'
t_start = time.mktime(time.localtime())


# define pyecloud parameters
filename_chm = None
B_multip_per_eV = [0.]
B_multip_per_eV = np.array(B_multip_per_eV)

init_unif_edens_flag = -1
init_unif_edens = None

# rest gas parameters
gas_ion_flag = 1
unif_frac = 0.
P_nTorr = 5.
sigma_ion_MBarn = 1.5
Temp_K = 300.
A = 44.
E_init_ion = 0.157e-27 * A
ion_mass = A * m_p
ion_charge = e

# MPs and MP size (set to match with FASTION code)
N_MP_ele_init = 501 * 156		#500 ion MPs/bunch
N_mp_max = N_MP_ele_init
nel_mp_ref_0 = P_nTorr * sigma_ion_MBarn / 37.89	#number of real ions/MP

# time step (set to bunch spacing)
Dt_ref = 1.0e-9


# define the machine
from CLIC_DR import CLIC_DR
machine_configuration = 'CLIC_DR_1GHz'
optics_mode = 'smooth'

n_segments = 26
n_turns = 10

if optics_mode == 'smooth':
    machine = CLIC_DR(machine_configuration=machine_configuration, n_segments=n_segments)

elif optics_mode == 'non-smooth':
    with open('CLIC_DR_n260_optics.pkl') as fid:
        optics = pickle.load(fid)
    optics.pop('circumference')
    machine = CLIC_DR(machine_configuration=machine_configuration, optics_mode = 'non-smooth',  **optics)


# define bunch
epsn_y = 0.0048e-6
if machine_configuration == 'CLIC_DR_1GHz':
    sigma_z = 1.8e-3
    epsn_x = 0.456e-6
elif machine_configuration == 'CLIC_DR_2GHz':
    sigma_z = 1.6e-3
    epsn_x = 0.472e-6

intensity = 4.1e9
# number of macroparticles per bunch
n_macroparticles = 10000


# make beam
n_bunches = 156
print 'Initializing', n_bunches, 'bunches, of', n_macroparticles, 'macroparticles'

bunches = []
for i_bun in xrange(n_bunches):
    print 'Bunch', i_bun
    bunch = machine.generate_6D_Gaussian_bunch(n_macroparticles=n_macroparticles, intensity=intensity,
                                                   epsn_x=epsn_x, epsn_y=epsn_y, sigma_z=sigma_z)
    bunch.z -= machine.circumference / machine.longitudinal_map.harmonics[0] * i_bun

    print 'Bunch centroid at', bunch.mean_x(), bunch.mean_y(), bunch.mean_z()
    bunches.append(bunch)

beam = sum(bunches)


# compute sigma x and y
inj_optics = machine.transverse_map.get_injection_optics()
sigma_x = np.sqrt(inj_optics['beta_x'] * epsn_x / machine.betagamma)
sigma_y = np.sqrt(inj_optics['beta_y'] * epsn_y / machine.betagamma)
print 'sigma_x = %.2e, sigma_y = %.2e'%(sigma_x, sigma_y)

# define apertures and PIC grid size
chamb_type = 'rect'
x_aper  = 20 * sigma_x
y_aper  = 20 * sigma_y

Dh_x = x_aper / 125.
Dh_y = y_aper / 125.
Dh_sc = [Dh_x, Dh_y]


# find bunch slots and make bunch spacing wide slicer
bunch_slots = beam.get_slices(machine.buncher)
machine.buncher.add_statistics(sliceset=bunch_slots, beam=beam, statistics=True)

mask_filled_slots = bunch_slots.n_macroparticles_per_slice > 0
n_filled_slots = np.sum(mask_filled_slots)
z_cuts_filled_slots = (np.min(bunch_slots.z_bins[mask_filled_slots]),
                       np.max(bunch_slots.z_bins[np.append(mask_filled_slots, True)]))

bunch_slicer = UniformBinSlicer(n_filled_slots, z_cuts=z_cuts_filled_slots)
bunch_slices = beam.get_slices(bunch_slicer)
bunch_slicer.add_statistics(sliceset=bunch_slices, beam=beam, statistics=True)


# define a beam monitor
from PyHEADTAIL.monitors.monitors import SliceMonitor
beam_monitor = SliceMonitor(filename='bunch_evolution_A%d_%db_%dips_%dturns_%.2fnTorr'%(A, n_bunches, n_segments, n_turns, P_nTorr),
                            n_steps = n_turns + 1, slicer=bunch_slicer, write_buffer_every=5)


# initialize ion cloud with single kick per bunch
import PyECLOUD.PyEC4PyHT as PyEC4PyHT
ecloud_sk = PyEC4PyHT.Ecloud(L_ecloud=machine.circumference / n_segments, slicer=bunch_slicer,
            Dt_ref=Dt_ref, pyecl_input_folder='./pyecloud_config', beam_monitor=None,
            chamb_type = chamb_type, PyPICmode = 'FFT_OpenBoundary',
            x_aper=x_aper, y_aper=y_aper,
            filename_chm=filename_chm, Dh_sc=Dh_sc,
            init_unif_edens_flag=init_unif_edens_flag,
            init_unif_edens=init_unif_edens,
            cloud_mass=ion_mass, cloud_charge=ion_charge,
            gas_ion_flag=gas_ion_flag, unif_frac=unif_frac,
            P_nTorr=P_nTorr, sigma_ion_MBarn=sigma_ion_MBarn,
            Temp_K=Temp_K, E_init_ion=E_init_ion,
            N_mp_max=N_mp_max,
            nel_mp_ref_0=nel_mp_ref_0,
            B_multip =B_multip_per_eV * machine.p0 / e * c,
            switch_model='perfect_absorber',
            kick_mode_for_beam_field=True)


#print grid size
nx, ny = ecloud_sk.spacech_ele.PyPICobj.nx, ecloud_sk.spacech_ele.PyPICobj.ny
print 'nx = %d, ny = %d'%(nx, ny)


# install ion clouds in the machine
machine.install_after_each_transverse_segment(ecloud_sk)


# run simulation
beam_monitor.dump(beam)
print 'Start track...'
t_start_sw = time.mktime(time.localtime())
print 'Time for initialization ',(t_start_sw - t_start), 's'
for i_turn in xrange(n_turns):
    print 'Turn %d'%(i_turn + 1)
    machine.track(beam)
    beam_monitor.dump(beam)
t_stop_sw = time.mktime(time.localtime())
print 'Done track in ', (t_stop_sw - t_start_sw), 's'
