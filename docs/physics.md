# Physical models and numerical algorithms

These references describe the physics and numerical methods underlying
PyECLOUD. Older presentations document the implementation at the time they were
written; use the current parameter reference for configuration syntax.

## General references

- G. Iadarola, [Electron Cloud Studies for CERN Particle Accelerators and Simulation Code Development](https://repository.cern/records/pb0hw-d9008), doctoral thesis (2014).
- G. Iadarola et al., [Evolution of Python Tools for the Simulation of Electron Cloud Effects](https://cds.cern.ch/record/2289165), IPAC2017.
- [Overview of PyECLOUD's history and features](https://indico.cern.ch/event/580885/).

See [Citing PyECLOUD](citation.md) for the project's recommended citation.

## Particle tracking and fields

| Topic | Reference |
| --- | --- |
| Boris tracking in quadrupole fields | [Tracking presentation](https://indico.cern.ch/event/310641/contributions/1681940/attachments/594402/818138/quad_ecloud_meeting_3.pptx) |
| Shortley–Weller treatment of curved boundaries | [Field-solver presentation](https://indico.cern.ch/event/320287/contributions/741211/attachments/617493/849678/sc_ecloud_meeting_5.pdf) |
| Nested-grid particle-in-cell calculations | [Nested-grid presentation](https://indico.cern.ch/event/547910/contributions/2221660/) |
| Magnetic multipoles | [Implementation discussion](https://indico.cern.ch/event/637703/) |
| Electromagnetic potentials for field calculations | [Technical presentation](https://indico.cern.ch/event/840676/contributions/3532754/) |

## Emission and ionization

| Topic | Reference |
| --- | --- |
| Cosine distribution of secondary-electron emission angles | [Emission-angle discussion](https://indico.cern.ch/event/673160/) |
| Furman–Pivi secondary emission model | [CERN-ACC-NOTE-2019-0029](https://cds.cern.ch/record/2683285) and [presentation](https://indico.cern.ch/event/830024/contributions/3488306/) |
| Electron-induced residual-gas ionization and multiple clouds | [Cross-ionization presentation](https://indico.cern.ch/event/835473/contributions/3525107/) |

## Beam coupling and computational methods

- [PyEC4PyHT coupling and Cython performance work](https://indico.cern.ch/event/394530/contributions/937467/attachments/789851/1082606/005_pyec4pyht_gi.pptx).
- [Parallel beam-dynamics simulation with PyPARIS](https://indico.cern.ch/event/547910/).
- [PyPARIS technical wiki](https://github.com/PyCOMPLETE/PyPARIS/wiki).

For execution instructions, start with the
[single-bunch tutorial](tutorials/single-bunch.md) and
[parallel execution guide](tutorials/parallel.md).
