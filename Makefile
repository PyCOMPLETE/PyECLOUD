
.PHONY: local cern all

# # On arch linux, the command f2py2 has to be used for python2 programs such as PyECLOUD.
# # On ubuntu, there is no command called f2py2
# F2PY2_EXIST := $(shell command -v f2py2 2> /dev/null)
# 
# ifdef F2PY2_EXIST
# 		F2PY = f2py2
# else
# 		F2PY = f2py
# endif

F2PY = f2py

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
	$(F2PY) -m rhocompute -c compute_rho.f
	$(F2PY) -m int_field_for -c interp_field_for.f
	$(F2PY) -m hist_for -c compute_hist.f
	$(F2PY) -m seg_impact -c update_seg_impact.f
	$(F2PY) -m errffor -c errfff.f
	$(F2PY) -m boris_step -c boris_step.f
	$(F2PY) -m vectsum -c vectsum.f
