"""One figure per exposure showing every step from input to final products."""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import AsinhNorm

BIN = 4


def _binned(image: np.ndarray, invalid: np.ndarray) -> np.ndarray:
    data = np.where(invalid, np.nan, image).astype(np.float32)
    ny, nx = data.shape
    data = data[:ny // BIN * BIN, :nx // BIN * BIN]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(data.reshape(ny // BIN, BIN, nx // BIN, BIN), axis=(1, 3))


def plot_steps(path: Path, title: str, panels: list[tuple[str, np.ndarray | None, str]],
               invalid: np.ndarray, sky: float, noise: float) -> None:
    """Draw panels of (label, array, kind) on shared stretches.

    kind is "sky" (sky retained, centred on the sky level), "skyless"
    (centred on zero), "model" (diverging, own symmetric scale) or "mask".
    """
    noise = max(noise, 1e-4)
    ncols = (len(panels) + 1) // 2
    fig, axes = plt.subplots(2, ncols, figsize=(4 * ncols, 8.6), constrained_layout=True)
    for ax in axes.flat:
        ax.set_axis_off()
    for ax, (label, array, kind) in zip(axes.flat, panels):
        ax.set_title(label, fontsize=10)
        if array is None:
            ax.text(0.5, 0.5, "not applied", ha="center", va="center",
                    transform=ax.transAxes, color="0.4")
            continue
        if kind == "mask":
            ax.imshow(array[::BIN, ::BIN], origin="lower", cmap="gray_r",
                      vmin=0, vmax=2, interpolation="nearest")
            continue
        if kind == "model":
            scale = max(float(np.nanpercentile(np.abs(array[::8, ::8]), 99.5)), 0.1 * noise)
            norm, cmap = plt.Normalize(-scale, scale), "coolwarm"
        else:
            # AsinhNorm is linear around zero, so display relative to the sky.
            if kind == "sky":
                array = array - np.float32(sky)
            norm = AsinhNorm(linear_width=noise, vmin=-3 * noise, vmax=25 * noise)
            cmap = "gray"
        image = ax.imshow(_binned(array, invalid), origin="lower", cmap=cmap,
                          norm=norm, interpolation="nearest")
        if kind == "model":
            fig.colorbar(image, ax=ax, shrink=0.7, label="MJy/sr")
    fig.suptitle(title, fontsize=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=80)
    plt.close(fig)
