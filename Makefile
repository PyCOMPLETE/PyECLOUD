PYTHON ?= python

.PHONY: all local
all: local

local:
	$(PYTHON) -m pip install -e .
