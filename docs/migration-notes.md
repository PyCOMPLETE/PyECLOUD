---
orphan: true
---

# Migration and audit notes

This page records how the documentation was migrated from the original wiki and
reference manual, and how far its content has been checked against the code.

## Audit and validation notes

A clean documentation build checks markup and references; it does not verify
simulation behaviour. The tutorials still need a complete numerical validation.

The buildup parameter pages have been checked against
`PyECLOUD/default_input_parameters.py` and the code that reads each parameter:

- Every documented parameter name, and every `furman_pivi_surface` key, exists
  in the code. Misspelled names inherited from the manual were corrected
  (`secondary_beams_file_list`, `B_zero_thrhld`, `phem_resc_fac`,
  `theta_r_max`), and `flag_counter_clockwise_chamb`, which does not exist,
  was removed.
- Defaults were checked for parameters marked *Optional*, and accepted option
  values were checked for `switch_model`, `secondary_angle_distribution`,
  `energy_distribution`, `beam_field_file`, `sparse_solver`, and `track_method`.
- Behaviour described in the descriptions was checked where it was changed or
  added, for example `t_ion`, `fact_clean`, the uniform initial density, and
  the photoemission energy distributions.

Not yet checked: the physical meaning and units of the SEY model parameters,
most units elsewhere, and descriptions carried over unchanged from the manual.
The input parameters `stopfile`, `filename_init_MP_state`,
`flag_continuous_emission`, and `cloud_name` are not yet documented.

## Wiki migration record

The initial migration used the public PyECLOUD wiki at revision
`5d30c514549031a5d5baef3eba706a42fa27feed`, reviewed on 7 October 2026.
New pages link to their original sources so authorship and revision history
remain traceable.

The reference manual was converted to Markdown during the initial migration.
The converter is not part of the documentation build; the Markdown pages are
now the sources.

The restart guide was checked against `PyECLOUD/buildup_simulation.py`,
`PyECLOUD/pyecloud_saver.py`, and `PyPARIS_sim_class/Save_Load_Status.py`.
Output field names were checked against `PyECLOUD/pyecloud_saver.py`, including
`xg_hist_det`, which was duplicated as `nel_hist_det` in the original manual.
The output reference also preserves the manual's array dimensions, sampling
definitions, and enabling flags. The manual describes lifetime histograms as
being saved per bunch passage; the current implementation uses
`Dt_lifetime_hist` and saves `t_lifetime_hist`, so the reference documents that
timing instead. Energy normalization and the transverse sampling region for
`cen_density` follow the implementation rather than the manual's older wording.
The parallel execution guide retains the wiki's MPI environment advice;
the old Python and MPI source-compilation recipes remain historical.

The buildup tutorial plots the output of the reader's own simulation. The older
plotting script in `docs/legacy_tex/example/` reads a stored regression
reference file. The plotting guide presents all 20 legacy output examples in
five grouped SVG figures.

| Original material | Destination |
| --- | --- |
| Home: simulation modes, libraries, citation | Homepage, installation, and citation pages |
| Buildup tutorial | Updated buildup tutorial and plotting example |
| PyECLOUD–PyHEADTAIL tutorial | Single-bunch tutorial and parallel execution |
| Reference manual | Input parameter pages, output reference, and downloadable PDF |
| Physical models and algorithms | Physics references |
| Python 3 setup: MPI advice | Parallel execution guide |
| Miniforge installation | Link to the official installation guide |
| Checkpointing presentation link | Restart guide checked against the current code |

Python 2, Ubuntu 14.04, the Python 3 transition announcement, and the obsolete
index remain historical. The old macOS recipe has not been promoted to current
instructions without platform validation. The separate Furman–Pivi and fast
beam-ion wiki pages were placeholders; the useful Furman–Pivi references were
retained from the physics page. The old PDF compilation page is superseded by
the Sphinx build instructions; its LaTeX sources remain in
`docs/legacy_tex/reference/src/`.
