all:local

cern:
	/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/f2py -m rhocompute -c compute_rho.f
	/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/f2py -m int_field_for -c interp_field_for.f
	/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/f2py -m hist_for -c compute_hist.f
	/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/f2py -m seg_impact -c update_seg_impact.f
	/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/f2py -m errffor -c errfff.f
	/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/f2py -m boris_step -c boris_step.f
	/afs/cern.ch/project/uslarp/opt/lxplus64/Python-2.7.2/bin/f2py -m vectsum -c vectsum.f
local:
	f2py -m rhocompute -c compute_rho.f
	f2py -m int_field_for -c interp_field_for.f
	f2py -m hist_for -c compute_hist.f
	f2py -m seg_impact -c update_seg_impact.f
	f2py -m errffor -c errfff.f
	f2py -m boris_step -c boris_step.f
	f2py -m vectsum -c vectsum.f
