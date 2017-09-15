
.PHONY: all python2 python3

# On arch linux, the command f2py2 has to be used for python2 programs such as PyECLOUD.
# On ubuntu, there is no command called f2py2
ARCH_LINUX := $(shell command -v pacman &2>/dev/null)

ifdef ARCH_LINUX
		F2PY_py2 = f2py2
		F2PY_py3 = f2py
else
		F2PY_py2 = f2py
		F2PY_py3 = f2py3
endif


all: python2

python2:
	$(F2PY_py2) -m rhocompute -c compute_rho.f
	$(F2PY_py2) -m int_field_for -c interp_field_for.f
	$(F2PY_py2) -m hist_for -c compute_hist.f
	$(F2PY_py2) -m seg_impact -c update_seg_impact.f
	$(F2PY_py2) -m errffor -c errfff.f
	$(F2PY_py2) -m boris_step -c boris_step.f
	$(F2PY_py2) -m vectsum -c vectsum.f
	$@ setup.py build_ext -i

python3:
	$(F2PY_py3) -m rhocompute -c compute_rho.f
	$(F2PY_py3) -m int_field_for -c interp_field_for.f
	$(F2PY_py3) -m hist_for -c compute_hist.f
	$(F2PY_py3) -m seg_impact -c update_seg_impact.f
	$(F2PY_py3) -m errffor -c errfff.f
	$(F2PY_py3) -m boris_step -c boris_step.f
	$(F2PY_py3) -m vectsum -c vectsum.f
	$@ setup.py build_ext -i
