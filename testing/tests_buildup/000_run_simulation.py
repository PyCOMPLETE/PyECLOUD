import sys, os
BIN = os.path.expanduser("../../../")
sys.path.append(BIN)

# sim_folder = 'LHC_ArcDipReal_450GeV_sey1.70_2.5e11ppb_bl_1.00ns'
# sim_folder = 'LHC_ArcQuadReal_450GeV_sey1.65_2.5e11ppb_bl_1.00ns'
sim_folder = 'LHC_ArcDipReal_450GeV_sey1.00_2.5e11ppb_bl_1.00ns_gas_ionization'


from PyECLOUD.buildup_simulation import BuildupSimulation


sim = BuildupSimulation(pyecl_input_folder = sim_folder, filen_main_outp = sim_folder+'/Pyecltest.mat')
sim.run()


print ''
print 'Test simulation done!'
print 'To instpect the results you can run:'
print '001_comparison_against_reference.py'
print ''
