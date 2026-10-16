"""Plot the buildup observables shown in the PyECLOUD documentation."""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PyECLOUD.myloadmat_to_obj import myloadmat_to_obj


def plot_outputs(filename, output_dir):
    data = myloadmat_to_obj(str(filename))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "svg.fonttype": "none"})

    def panel(name, entries):
        fig, axes = plt.subplots(2, 2, figsize=(10, 7), layout="constrained")
        for ax, (x, y, title, xlabel, ylabel) in zip(axes.flat, entries):
            ax.plot(x, y, color="#006b91", linewidth=1.4)
            ax.set(title=title, xlabel=xlabel, ylabel=ylabel)
            ax.grid(alpha=0.2)
            ax.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4))
        fig.savefig(output_dir / f"{name}.svg")
        plt.close(fig)

    time = np.atleast_1d(data.t) * 1e9
    groups = [
        ("population", [
            ("lam_t_array", "Beam line density", "Beam particles / m"),
            ("Nel_timep", "Cloud population", "Electrons / m"),
            ("Nel_imp_time", "Wall impacts", "Electrons / m / saved interval"),
            ("Nel_emit_time", "Secondary emission", "Electrons / m / saved interval"),
        ]),
        ("energy-density", [
            ("En_imp_eV_time", "Energy deposited at the wall", "eV / m / saved interval"),
            ("En_emit_eV_time", "Energy of emitted electrons", "eV / m / saved interval"),
            ("En_kin_eV_time", "Cloud kinetic energy", "eV / m"),
            ("cen_density", "Central electron density", "Electrons / m³"),
        ]),
    ]
    for name, fields in groups:
        panel(name, [(time, np.atleast_1d(getattr(data, key)), title,
                      "Time [ns]", unit) for key, title, unit in fields])

    passage = np.arange(np.size(data.t_hist))
    fields = [
        ("N_mp_pass", "Macroparticles in the chamber", "Macroparticles"),
        ("N_mp_impact_pass", "Macroparticle impacts", "Impacts / passage"),
        ("N_mp_corrected_pass", "Tracking corrections", "Corrections / passage"),
        ("N_mp_ref_pass", "Reference macroparticle weight", "Electrons / m / macroparticle"),
    ]
    panel("macroparticles", [(passage, np.atleast_1d(getattr(data, key)), title,
                             "Recorded passage index", unit)
                            for key, title, unit in fields])

    x = np.atleast_1d(data.xg_hist) * 1e3
    fields = [
        ("nel_hist", "Cloud population", "Electrons / m"),
        ("nel_impact_hist_scrub", "Impacts above the scrubbing threshold", "Electrons / m"),
        ("nel_impact_hist_tot", "All wall impacts", "Electrons / m"),
        ("energ_eV_impact_hist", "Energy deposited at the wall", "eV / m"),
    ]
    histograms = [np.asarray(getattr(data, key)).reshape(len(passage), len(x))
                  for key, _, _ in fields]
    panel("per-passage", [(passage, hist.sum(axis=1), title,
                           "Recorded passage index", unit)
                          for hist, (_, title, unit) in zip(histograms, fields)])
    panel("horizontal", [(x, hist.sum(axis=0), title,
                          "Horizontal position [mm]", unit + " / bin (summed)")
                         for hist, (_, title, unit) in zip(histograms, fields)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("filename", nargs="?", default="Pyecltest.mat")
    parser.add_argument("--output-dir", default="buildup-plots")
    args = parser.parse_args()
    plot_outputs(args.filename, args.output_dir)
