from __future__ import annotations

from pathlib import Path

import numpy as np


def _configure_matplotlib():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    mpl.rcParams["svg.fonttype"] = "none"
    return plt


def save_heatmap(matrix: np.ndarray, title: str, out_path: Path) -> None:
    plt = _configure_matplotlib()
    plt.figure()
    plt.imshow(matrix, aspect="auto")
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Position")
    plt.colorbar(label="Score")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def save_projection_plot(signal: np.ndarray, title: str, out_path: Path) -> None:
    plt = _configure_matplotlib()
    plt.figure()
    plt.plot(signal)
    plt.title(title)
    plt.xlabel("Position")
    plt.ylabel("Signal")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


def save_fft_plot(
    power_spectrum: np.ndarray,
    title: str,
    out_path: Path,
    minimum_frequency: int,
    maximum_frequency: int,
    highlight_frequency: int | None = None,
) -> None:
    plt = _configure_matplotlib()
    lower_bound = max(1, minimum_frequency)
    upper_bound = min(maximum_frequency, len(power_spectrum) - 1)
    trimmed_power = power_spectrum[lower_bound:upper_bound + 1]
    frequencies = np.arange(lower_bound, lower_bound + len(trimmed_power))

    plt.figure()
    plt.plot(frequencies, trimmed_power)
    if highlight_frequency is not None and highlight_frequency > 0:
        plt.axvline(highlight_frequency, linestyle="--")
    plt.title(title)
    plt.xlabel("Frequency index")
    plt.ylabel("Power")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()
