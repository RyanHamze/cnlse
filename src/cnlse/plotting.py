"""Plotting helpers for examples and tutorials."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def plot_power_traces(
    time,
    traces,
    output: str | Path,
    title: str = "Fiber Optics Modulation Simulation",
    subtitle: str | None = None,
) -> Path:
    """Plot one or more optical power traces.

    Parameters
    ----------
    time:
        Time/sample axis.
    traces:
        Sequence of `(label, values)` pairs.
    output:
        File path for the saved plot.
    title:
        Main plot title.
    subtitle:
        Optional smaller subtitle.
    """

    time = np.asarray(time)
    if not traces:
        raise ValueError("at least one trace is required")

    fig, axes = plt.subplots(len(traces), 1, figsize=(12, 2.6 * len(traces)), sharex=True)
    if len(traces) == 1:
        axes = [axes]

    for axis, (label, values) in zip(axes, traces):
        values = np.asarray(values)
        axis.plot(time, values, linewidth=2)
        axis.fill_between(time, 0, values, alpha=0.13)
        axis.set_title(label, loc="left")
        axis.set_ylabel("Power")
        axis.grid(True, alpha=0.25)

    axes[-1].set_xlabel("Time / symbol slots")

    if subtitle:
        fig.suptitle(f"{title}\n{subtitle}", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.92])
    else:
        fig.suptitle(title, fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.94])

    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160)
    plt.close(fig)
    return output
