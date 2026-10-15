# PyECLOUD

PyECLOUD simulates electron cloud effects in particle accelerators. It supports
electron cloud buildup simulations and, together with PyHEADTAIL, simulations
of the interaction between the cloud and the beam.

Start with installation and a buildup simulation, then explore beam dynamics,
input parameters, outputs, and the underlying physical models.

## Simulation modes and packages

- **Buildup:** a prescribed beam drives the electron cloud at an accelerator section.
- **Beam dynamics:** PyHEADTAIL and PyECLOUD exchange beam/cloud kicks for
  single-bunch and coupled-bunch studies.
- **Ion-cloud studies:** additional cloud species and cross-ionization models
  support studies involving ions.

| Package | Role |
| --- | --- |
| `PyECLOUD` | Cloud generation, tracking, emission, and beam coupling |
| `PyPIC` (`pypic-poisson` on PyPI) | Particle-in-cell field solvers; installed as a dependency |
| `PyHEADTAIL` | Beam dynamics for coupled simulations |
| `PyPARIS` | Serial, multiprocessing, and MPI launchers, included with PyECLOUD |
| `PyPARIS_sim_class`, `PyPARIS_CoupledBunch_sim_class` | Configurable simulation classes, included with PyECLOUD |
| `PyKLU` | Default sparse solver; installed as a dependency |
| `nafflib` | Optional frequency analysis, available on PyPI |

Values in the parameter examples illustrate configurations; they are not
complete input files or a list of universal defaults.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
```

```{toctree}
:maxdepth: 1
:caption: Tutorials

tutorials/buildup
tutorials/buildup-outputs
tutorials/single-bunch
tutorials/parallel
tutorials/restart
```

```{toctree}
:maxdepth: 2
:caption: Batch jobs

batch/index
```

```{toctree}
:maxdepth: 2
:caption: Reference

parameters/index
parameters/single-bunch
reference/outputs
physics
citation
```

```{toctree}
:maxdepth: 1
:caption: Contributing

documentation
contributors
```

The [original wiki](https://github.com/PyCOMPLETE/PyECLOUD/wiki) remains available
as a historical source. See [Citing PyECLOUD](citation.md) when using the code
in a publication.
