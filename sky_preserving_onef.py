#!/usr/bin/env python3
"""Fit residual NIRCam 1/f from registered dithers with fixed background level.

This accepts either a sky-retaining cal or a sky-subtracted cal_sky_renorm
exposure. The other exposures estimate the sky-coordinate scene; their scalar
levels are matched only in the temporary template. The target's masked
background median is retained, whether that level is sky or near zero.
"""

from __future__ import annotations

import argparse
import json
import shutil
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from matplotlib.colors import AsinhNorm, TwoSlopeNorm
from scipy import ndimage as ndi

from onef_masks import double_source_mask
from skysurfir_clean import base_invalid, robust_sigma


def read_cal(path: Path, science_ext: str = "SCI") -> tuple[np.ndarray, np.ndarray, WCS, str, str, bool]:
    with fits.open(path, memmap=False) as hdul:
        sky_subtracted = bool(hdul["SCI"].header.get("NMFSKY", False))
        image = np.array(hdul[science_ext].data, dtype=np.float32, copy=True)
        dq = np.array(hdul["DQ"].data, dtype=np.uint32, copy=True)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FITSFixedWarning)
            wcs = WCS(hdul["SCI"].header)
        detector = str(hdul[0].header.get("DETECTOR", ""))
        band = str(hdul[0].header.get("FILTER", ""))
    return image, dq, wcs, detector, band, sky_subtracted


def shift_to_target(target_wcs: WCS, peer_wcs: WCS, shape: tuple[int, int]) -> tuple[float, float]:
    ny, nx = shape
    xy = np.array([[nx / 2, ny / 2], [100, 100], [nx - 100, 100],
                   [100, ny - 100], [nx - 100, ny - 100]], dtype=float)
    peer_xy = peer_wcs.all_world2pix(target_wcs.all_pix2world(xy, 0), 0)
    delta = peer_xy - xy
    if np.max(np.linalg.norm(delta - delta[0], axis=1)) > 0.1:
        raise ValueError("WCS differs by more than a shift; use full reprojection")
    # scipy.shift output[p] samples input[p - shift].
    return -float(delta[0, 1]), -float(delta[0, 0])


def aligned_peer(path: Path, target: np.ndarray, good: np.ndarray,
                 target_wcs: WCS, detector: str, band: str,
                 sky_subtracted: bool, science_ext: str = "SCI") -> tuple[np.ndarray, dict]:
    image, dq, wcs, peer_detector, peer_band, peer_sky_subtracted = read_cal(
        path, science_ext=science_ext)
    if (peer_detector, peer_band, image.shape, peer_sky_subtracted) != (
            detector, band, target.shape, sky_subtracted):
        raise ValueError(f"Peer detector/filter/shape/sky status differs: {path}")
    shift = shift_to_target(target_wcs, wcs, target.shape)
    image[base_invalid(dq, image)] = np.nan
    del dq
    aligned = ndi.shift(image, shift, order=1, mode="constant", cval=np.nan,
                        prefilter=False, output=np.float32)
    del image
    sample = np.s_[::8, ::8]
    valid = good[sample] & np.isfinite(aligned[sample])
    if valid.sum() < 10000:
        raise ValueError(f"Too little sky overlap with {path}")
    level = float(np.median((aligned[sample] - target[sample])[valid]))
    aligned -= np.float32(level)  # Template only; no change to the target sky.
    return aligned, {"peer": path.name, "shift_yx": shift, "template_level_shift": level}


def fit_stripes(target: np.ndarray, template: np.ndarray,
                good: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    residual = target - template
    valid = good & np.isfinite(residual)
    center = float(np.median(residual[valid][::8]))
    sigma = robust_sigma(residual[valid][::8])
    valid &= np.abs(residual - center) < 6 * sigma
    if valid.mean() < 0.5:
        raise ValueError("Too few sky pixels to fit residual 1/f")
    residual -= np.float32(center)
    residual[~valid] = np.nan
    ny, nx = target.shape
    row_map = np.empty_like(target)
    row_profiles = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for lo in range(0, nx, nx // 4):
            hi = min(lo + nx // 4, nx)
            profile = np.nanmedian(residual[:, lo:hi], axis=1)
            if np.isfinite(profile).sum() < ny * 0.9:
                raise ValueError("Too many masked rows in an amplifier")
            profile = np.interp(np.arange(ny), np.flatnonzero(np.isfinite(profile)),
                                profile[np.isfinite(profile)]).astype(np.float32)
            row_map[:, lo:hi] = profile[:, None]
            row_profiles.append(profile)
        column = np.nanmedian(residual - row_map, axis=0)
    column = np.nan_to_num(column, nan=0).astype(np.float32)
    model = row_map + column[None, :]
    # Fix the one unconstrained global mode to the incoming product's level.
    model -= np.float32(np.median(model[good]))
    corrected = target - model
    before_sky = float(np.median(target[good]))
    after_sky = float(np.median(corrected[good]))
    model += np.float32(after_sky - before_sky)
    corrected = target - model
    model[~np.isfinite(target)] = 0
    corrected[~np.isfinite(target)] = target[~np.isfinite(target)]
    return corrected, model, {"blank_median_before": before_sky,
                              "blank_median_after": float(np.median(corrected[good])),
                              "residual_sigma": sigma,
                              "model_rms": float(np.sqrt(np.mean(model[good] ** 2))),
                              "fit_fraction": float(valid.mean())}


def row_rms(difference: np.ndarray, good: np.ndarray) -> float:
    valid = good & np.isfinite(difference)
    scalar = float(np.median(difference[valid][::8]))
    rows = []
    width = difference.shape[1] // 4
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for lo in range(0, difference.shape[1], width):
            values = np.where(valid[:, lo:lo + width], difference[:, lo:lo + width], np.nan)
            rows.append(np.nanmedian(values, axis=1) - scalar)
    return float(np.sqrt(np.nanmean(np.asarray(rows, np.float32) ** 2)))


def binned(image: np.ndarray, factor: int) -> np.ndarray:
    ny, nx = image.shape
    view = image[:ny // factor * factor, :nx // factor * factor]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmedian(view.reshape(ny // factor, factor, nx // factor, factor), axis=(1, 3))


def tile_rms(difference: np.ndarray, good: np.ndarray) -> float:
    valid = good & np.isfinite(difference)
    scalar = float(np.median(difference[valid][::8]))
    tiles = binned(np.where(valid, difference - scalar, np.nan), 64)
    return float(np.sqrt(np.nanmean(tiles ** 2)))


def render(path: Path, before: np.ndarray, model: np.ndarray, after: np.ndarray,
           good: np.ndarray, metrics: dict) -> None:
    sky = metrics["blank_median_before"]
    product = "cal_sky_renorm: sky subtracted" if metrics["sky_subtracted"] else "cal: sky retained"
    sigma = max(robust_sigma(before[good][::16]), 0.002)
    sci_norm = AsinhNorm(linear_width=sigma, vmin=sky - 2.5 * sigma,
                         vmax=sky + 6 * sigma)
    scale = max(0.005, float(np.nanpercentile(abs(model[good]), 99.5)))
    mod_norm = TwoSlopeNorm(vmin=-scale, vcenter=0, vmax=scale)
    fig, axes = plt.subplots(3, 3, figsize=(16, 14), constrained_layout=True)
    arrays = (before, model, after)
    labels = (f"Input {product}", "Additional 1/f model subtracted",
              f"Corrected {product}")
    crop = (slice(660, 1320), slice(660, 1320))
    for col, (array, label) in enumerate(zip(arrays, labels)):
        norm = mod_norm if col == 1 else sci_norm
        cmap = "coolwarm" if col == 1 else "gray"
        axes[0, col].imshow(binned(array, 2), cmap=cmap, norm=norm,
                            origin="lower", interpolation="nearest")
        axes[1, col].imshow(array[crop], cmap=cmap, norm=norm,
                            origin="lower", interpolation="nearest")
        background = binned(np.where(good, array, np.nan), 64)
        bg_norm = (mod_norm if col == 1 else
                   TwoSlopeNorm(vmin=sky - 0.02, vcenter=sky, vmax=sky + 0.02))
        axes[2, col].imshow(background, cmap="coolwarm", norm=bg_norm,
                            origin="lower", interpolation="nearest")
        axes[0, col].set_title(label)
        axes[1, col].set_title("Central 660 × 660 pixels, native sampling")
        axes[2, col].set_title("Masked 64 × 64 pixel background medians" +
                               (f" (±{scale:.3f})" if col == 1 else " (center ±0.020 MJy/sr)"))
        for ax in axes[:, col]:
            ax.set_xticks([])
            ax.set_yticks([])
    holdout = (f"held-out stripe RMS {metrics['holdout_row_rms_before']:.5f} → "
               f"{metrics['holdout_row_rms_after']:.5f} MJy/sr"
               if metrics['holdout_row_rms_before'] is not None else "no held-out dithers")
    status = "accepted" if metrics["accepted"] else "diagnostic only: held-out check failed"
    fig.suptitle(f"{metrics['target']} | {status} | masked median {sky:.6f} → "
                 f"{metrics['blank_median_after']:.6f} MJy/sr | {holdout}", fontsize=13)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_cal(target: Path, output: Path, science: np.ndarray, peers: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".writing")
    shutil.copyfile(target, temporary)
    try:
        with fits.open(temporary, mode="update", memmap=True) as hdul:
            hdul["SCI"].data[:] = science
            hdul["SCI"].header["DITH1F"] = (True, "residual 1/f from registered dithers")
            hdul["SCI"].header["DITHN"] = (peers, "number of dither peers in scene template")
            hdul["SCI"].add_checksum()
            hdul.flush()
        temporary.replace(output)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=Path, required=True)
    parser.add_argument("--peers", type=Path, nargs="+", required=True)
    parser.add_argument("--mask", type=Path, required=True,
                        help="diagnostic FITS containing SOURCE_MASK")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--fit-peers", type=int, default=5)
    args = parser.parse_args()
    if not 2 <= args.fit_peers <= len(args.peers):
        parser.error("--fit-peers must be between 2 and the number of peers")
    target, dq, wcs, detector, band, sky_subtracted = read_cal(args.target)
    with fits.open(args.mask, memmap=True) as hdul:
        source = np.asarray(hdul["SOURCE_MASK"].data, bool)
        objects = np.asarray(hdul["OBJECT_MASK"].data, bool) if sky_subtracted and "OBJECT_MASK" in hdul else None
    grown, dilation = double_source_mask(source)
    if objects is not None:
        grown |= objects
    good = ~(base_invalid(dq, target) | grown)
    del dq, source, grown
    if good.mean() < 0.5:
        raise ValueError("Too much of target masked for reliable 1/f fitting")
    fit_paths = args.peers[:args.fit_peers]
    holdout_paths = args.peers[args.fit_peers:]
    aligned = []
    records = []
    for path in fit_paths:
        peer, record = aligned_peer(path, target, good, wcs, detector, band, sky_subtracted)
        aligned.append(peer)
        records.append(record)
    stack = np.stack(aligned)
    del aligned
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        template = np.nanmedian(stack, axis=0, overwrite_input=True).astype(np.float32)
    del stack
    corrected, model, metrics = fit_stripes(target, template, good)
    del template
    validation = []
    for path in holdout_paths:
        peer, record = aligned_peer(path, target, good, wcs, detector, band, sky_subtracted)
        validation.append({**record,
                           "row_rms_before": row_rms(target - peer, good),
                           "row_rms_after": row_rms(corrected - peer, good),
                           "tile_rms_before": tile_rms(target - peer, good),
                           "tile_rms_after": tile_rms(corrected - peer, good)})
        del peer
    metrics.update({"target": args.target.name, "sky_subtracted": sky_subtracted,
                    "fit_peers": records,
                    "holdout_peers": validation, "source_mask_dilation": dilation,
                    "holdout_row_rms_before": float(np.mean([v["row_rms_before"] for v in validation]))
                        if validation else None,
                    "holdout_row_rms_after": float(np.mean([v["row_rms_after"] for v in validation]))
                        if validation else None,
                    "holdout_tile_rms_before": float(np.mean([v["tile_rms_before"] for v in validation]))
                        if validation else None,
                    "holdout_tile_rms_after": float(np.mean([v["tile_rms_after"] for v in validation]))
                        if validation else None})
    metrics["accepted"] = bool(validation) and all(
        v["row_rms_after"] < 0.95 * v["row_rms_before"] and
        v["tile_rms_after"] < 0.95 * v["tile_rms_before"]
        for v in validation)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_cal = args.output_dir / args.target.name
    if metrics["accepted"]:
        write_cal(args.target, out_cal, corrected, len(records))
    render(args.output_dir / (args.target.stem + "_sky_preserving_onef.png"),
           target, model, corrected, good, metrics)
    (args.output_dir / (args.target.stem + "_metrics.json")).write_text(
        json.dumps(metrics, indent=2) + "\n")
    print(json.dumps({"cal": str(out_cal) if metrics["accepted"] else None,
                      "png": str(args.output_dir / (args.target.stem + "_sky_preserving_onef.png")),
                      "accepted": metrics["accepted"],
                      "masked_median_before": metrics["blank_median_before"],
                      "masked_median_after": metrics["blank_median_after"],
                      "heldout_row_rms_before": metrics["holdout_row_rms_before"],
                      "heldout_row_rms_after": metrics["holdout_row_rms_after"],
                      "heldout_tile_rms_before": metrics["holdout_tile_rms_before"],
                      "heldout_tile_rms_after": metrics["holdout_tile_rms_after"]}, indent=2))


if __name__ == "__main__":
    main()
