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

    kind is "sky" (sky retained, displayed relative to the measured median),
    "skyless" (centred on zero), "absolute" (raw MJy/sr values on a shared
    pair stretch), "model" (diverging) or "mask".
    """
    noise = max(noise, 1e-4)
    absolute_norm = AsinhNorm(
        linear_width=noise,
        vmin=min(-3 * noise, sky - 3 * noise),
        vmax=max(25 * noise, sky + 25 * noise))
    ncols = min(5, len(panels))
    nrows = (len(panels) + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(4 * ncols, 4.3 * nrows),
                             squeeze=False, constrained_layout=True)
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
        if kind in ("model", "model_binned"):
            model_image = (array if kind == "model_binned"
                           else _binned(array, invalid))
            scale = max(float(np.nanpercentile(np.abs(model_image), 99.5)), 0.1 * noise)
            norm, cmap = plt.Normalize(-scale, scale), "coolwarm"
        else:
            # AsinhNorm is linear around zero, so display relative to the sky.
            if kind == "sky":
                array = array - np.float32(sky)
            norm = (absolute_norm if kind == "absolute" else
                    AsinhNorm(linear_width=noise, vmin=-3 * noise,
                              vmax=25 * noise))
            cmap = "gray"
            model_image = _binned(array, invalid)
            if kind == "absolute":
                # Final cal and cal_sky_renorm share this exact absolute
                # stretch; do not subtract SKYMED from cal for display.
                ax.text(0.02, 0.02, "shared absolute stretch",
                        transform=ax.transAxes, fontsize=7, color="yellow",
                        ha="left", va="bottom")
                fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                             shrink=0.7, label="MJy/sr")
        image = ax.imshow(model_image, origin="lower", cmap=cmap,
                          norm=norm, interpolation="nearest")
        if kind == "model":
            fig.colorbar(image, ax=ax, shrink=0.7, label="MJy/sr")
    fig.suptitle(title, fontsize=12)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=80)
    plt.close(fig)
