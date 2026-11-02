#!/usr/bin/env python3
"""Standalone, per-exposure SkySurfIR NIRCam cleanup pilot.

This runner is intentionally outside ``JUMPROPE/``. It reads stage-2 ``*_cal.fits``
files, never edits its inputs, and writes independent ``cal`` and
``cal_sky_renorm`` products. It ports the ProFound scan-style 1/f estimate,
uses Photutils segmentation for object masks and a robust quadratic sky model,
and compares an LW-referenced dynamic wisp model with the installed NMF basis.

Example:
    python skysurfir_clean.py --visit-id 1813005 --cores 64 \
        --input-dir /data2/.../nircam/stage2/cal \
        --output-dir /tmp/skysurfir_1813005_pilot

The pilot chooses the LW-referenced dynamic estimate for subtraction, using
the NMF fit only to fill LW-source masks and areas without LW coverage. The
independent NMF model is saved beside it for inspection; the two models are
never subtracted sequentially.
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import re
import sys
import time
import traceback
import warnings
from functools import lru_cache
from pathlib import Path
from typing import Any

# Each outer process owns one exposure. Prevent nested BLAS/OpenMP thread pools.
for _thread_env in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ[_thread_env] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from photutils.segmentation import SourceFinder
from scipy import ndimage as ndi
from scipy.ndimage import map_coordinates
from scipy.optimize import nnls


DEFAULT_INPUT = Path(
    "/data2/rzi/skysurfir/w200-v1-stage2-pid01813-test-v1/nircam/stage2/cal"
)
DEFAULT_TEMPLATES = Path(__file__).resolve().parent / "templates"
KEEP_HDU_NAMES = {
    "PRIMARY",
    "SCI",
    "ERR",
    "DQ",
    "AREA",
    "VAR_POISSON",
    "VAR_RNOISE",
    "VAR_FLAT",
    "ASDF",
}
WISP_DETECTORS = {"NRCA3", "NRCA4", "NRCB3", "NRCB4"}
WISP_FILTERS = {"F150W", "F200W"}


def robust_sigma(values: np.ndarray) -> float:
    vals = np.asarray(values)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan")
    med = np.median(vals)
    sig = 1.482602218505602 * np.median(np.abs(vals - med))
    return float(sig)


def odd_window(window: int, length: int) -> int:
    out = min(int(window), int(length))
    if out % 2 == 0:
        out -= 1
    return max(1, out)


def scan_1f(
    image: np.ndarray,
    invalid: np.ndarray,
    block_rows: int = 512,
    trend_block: int = 501,
    clip_hi: float = 0.90,
    sky_quantile: float = 0.40,
    keep_trend: bool = True,
) -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    """Python port of JUMPROPE's ProFound ``profoundSkyScan`` XY mode.

    It estimates per-column quantiles in 512-row scan blocks and per-row
    quantiles across the detector width, centers each map, then removes the
    running-median trend used by the configured JUMPROPE ``ow_large`` mode
    only when ``keep_trend`` is true. The ordinary NIRCam mode leaves the
    centered row and column estimates intact.
    Pixels flagged by DQ, zero-valued, nonfinite, or in an optional source mask
    are preserved in the corrected science array.
    """
    orig = np.asarray(image, dtype=np.float32)
    bad = np.asarray(invalid, dtype=bool) | ~np.isfinite(orig) | (orig == 0)
    good = ~bad
    if good.sum() < 1000:
        raise ValueError("too few usable pixels for 1/f scan")
    clip_value = float(np.quantile(orig[good], clip_hi)) if clip_hi < 1 else float("inf")
    scan = orig.astype(np.float64, copy=True)
    scan[bad | (scan > clip_value)] = np.nan

    ny, nx = scan.shape
    qrow = []
    row_slices = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for y0 in range(0, ny, block_rows):
            y1 = min(ny, y0 + block_rows)
            row_slices.append((y0, y1))
            qrow.append(np.nanquantile(scan[y0:y1, :], sky_quantile, axis=0))
    qrow_arr = np.asarray(qrow, dtype=np.float64)
    row_center = float(np.nanmedian(qrow_arr))
    row_map = np.zeros((ny, nx), dtype=np.float32)
    wx = odd_window(trend_block, nx)
    for ib, (y0, y1) in enumerate(row_slices):
        vals = qrow_arr[ib] - row_center
        finite = np.isfinite(vals)
        if finite.any():
            # R's map cleanup turns wholly masked scan columns into zero
            # correction. Fill only for the running-median calculation, then
            # restore those map cells to zero so edge DQ gaps do not poison a
            # 501-pixel trend neighborhood.
            filled = np.interp(np.arange(nx), np.flatnonzero(finite), vals[finite])
            trend = ndi.median_filter(filled, size=wx, mode="nearest") if keep_trend and wx > 1 else 0.0
            piece = np.where(finite, vals - trend, 0.0)
            row_map[y0:y1, :] = piece[None, :]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        qcol = np.nanquantile(scan, sky_quantile, axis=1)
    col_center = float(np.nanmedian(qcol))
    wy = odd_window(trend_block, ny)
    centered_col = qcol - col_center
    finite_col = np.isfinite(centered_col)
    if finite_col.any():
        filled_col = np.interp(np.arange(ny), np.flatnonzero(finite_col), centered_col[finite_col])
        trend_col = ndi.median_filter(filled_col, size=wy, mode="nearest") if keep_trend and wy > 1 else 0.0
        col_map = np.where(finite_col, centered_col - trend_col, 0.0).astype(np.float32)[:, None]
    else:
        col_map = np.zeros((ny, 1), dtype=np.float32)
    total = row_map + col_map
    corrected = orig - total
    corrected[bad] = orig[bad]
    stats = {
        "clip_value_90": clip_value,
        "row_center": row_center,
        "column_center": col_center,
        "rms_map": float(np.sqrt(np.mean(total[good].astype(np.float64) ** 2))),
    }
    return corrected.astype(np.float32), total.astype(np.float32), stats


def base_invalid(dq: np.ndarray, image: np.ndarray) -> np.ndarray:
    dq = np.asarray(dq, dtype=np.uint32)
    bad = (dq & np.uint32(1)) != 0  # JWST DO_NOT_USE, as in JUMPROPE.
    bad = ndi.maximum_filter(bad.astype(np.uint8), size=3, mode="nearest") != 0
    bad |= ~np.isfinite(image) | (image == 0)
    return bad


def source_segmentation(
    image: np.ndarray,
    err: np.ndarray,
    invalid: np.ndarray,
    nsigma: float,
    background: np.ndarray | float | None = None,
) -> tuple[np.ndarray, int, float, float]:
    """Create an undilated Photutils segmentation map for one exposure."""
    valid = ~invalid & np.isfinite(image) & np.isfinite(err) & (err > 0)
    if valid.sum() < 1000:
        return np.zeros(image.shape, dtype=bool), 0, float("nan"), 0.0
    sample = np.asarray(image[valid][:: max(1, int(valid.sum() // 250_000))], dtype=float)
    med = float(np.median(sample))
    sigma = robust_sigma(sample)
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.median(err[valid]))
    floor = max(sigma * 0.5, 1e-12)
    threshold = np.maximum(np.asarray(err, dtype=np.float32), floor) * float(nsigma)
    if background is not None:
        threshold = threshold + np.asarray(background, dtype=np.float32)
    else:
        threshold = threshold + med
    finder = SourceFinder(
        n_pixels=5,
        connectivity=8,
        deblend=True,
        n_levels=24,
        contrast=0.001,
        progress_bar=False,
        n_processes=1,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        seg = finder(np.asarray(image, dtype=np.float32), threshold, mask=invalid)
    if seg is None:
        return np.zeros(image.shape, dtype=bool), 0, sigma, med
    mask = np.asarray(seg.data != 0)
    return mask, int(seg.n_labels), sigma, med


def sky_polynomial(
    image: np.ndarray,
    invalid: np.ndarray,
    object_mask: np.ndarray,
    max_points: int = 300_000,
    degree: int = 2,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Fit a robust 2-D linear or quadratic sky to unmasked pixels."""
    if degree not in (1, 2):
        raise ValueError("sky degree must be 1 or 2")
    ny, nx = image.shape
    valid = ~invalid & ~object_mask & np.isfinite(image)
    ids = np.flatnonzero(valid.ravel())
    if ids.size < 2000:
        raise ValueError(f"too few sky pixels for polynomial fit ({ids.size})")
    step = max(1, int(math.ceil(ids.size / max_points)))
    ids = ids[::step]
    y = (ids // nx).astype(np.float64)
    x = (ids % nx).astype(np.float64)
    xn = (x - (nx - 1) / 2) / max((nx - 1) / 2, 1)
    yn = (y - (ny - 1) / 2) / max((ny - 1) / 2, 1)
    columns = (np.ones_like(xn), xn, yn)
    if degree == 2:
        columns += (xn * xn, xn * yn, yn * yn)
    design = np.column_stack(columns)
    values = np.asarray(image.ravel()[ids], dtype=np.float64)
    use = np.isfinite(values)
    coef = np.zeros(design.shape[1], dtype=np.float64)
    for _ in range(5):
        if use.sum() < 1000:
            break
        coef, *_ = np.linalg.lstsq(design[use], values[use], rcond=None)
        resid = values - design @ coef
        sig = robust_sigma(resid[use])
        if not np.isfinite(sig) or sig <= 0:
            break
        center = float(np.median(resid[use]))
        new_use = np.isfinite(resid) & (np.abs(resid - center) <= 4.0 * sig)
        if np.array_equal(new_use, use):
            break
        use = new_use

    yy = (np.arange(ny, dtype=np.float64) - (ny - 1) / 2) / max((ny - 1) / 2, 1)
    xx = (np.arange(nx, dtype=np.float64) - (nx - 1) / 2) / max((nx - 1) / 2, 1)
    sky = coef[0] + coef[1] * xx[None, :] + coef[2] * yy[:, None]
    if degree == 2:
        sky = (sky + coef[3] * xx[None, :] ** 2
               + coef[4] * yy[:, None] * xx[None, :]
               + coef[5] * yy[:, None] ** 2)
    sky = sky.astype(np.float32)
    rms = robust_sigma(values[use] - design[use] @ coef) if use.any() else float("nan")
    meta = {
        "degree": degree,
        "coefficients_normalized_xy": [float(v) for v in coef],
        "fit_pixels": int(use.sum()),
        "available_sky_pixels": int(valid.sum()),
        "fit_residual_sigma": float(rms),
    }
    return sky, meta


@lru_cache(maxsize=16)
def load_nmf_basis(template_dir: str, detector: str, filt: str):
    """Load the independent NMF basis; upsample the distributed 512px assets."""
    root = Path(template_dir)
    det = detector.lower()
    filt_lower = filt.lower()
    candidates = [
        root / f"{det}_org" / f"{det}_{filt_lower}_wisp.fits.gz",
        root / f"{det}_org" / f"{det}_{filt_lower}_wisp.fits",
        root / det / f"{det}_{filt_lower}_wisp.fits.gz",
        root / det / f"{det}_{filt_lower}_wisp.fits",
    ]
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        return None
    with fits.open(path, memmap=False) as hdul:
        bases = []
        for i in range(10):
            key = f"WISP_{i}"
            if key not in hdul:
                break
            bases.append(np.asarray(hdul[key].data, dtype=np.float32))
        if not bases:
            return None
        support = np.asarray(hdul["MASK"].data, dtype=np.float32)
        high = np.asarray(hdul["MASK_HSNR"].data, dtype=np.float32)
    target = (2048, 2048)
    if bases[0].shape != target:
        zoom = (target[0] / bases[0].shape[0], target[1] / bases[0].shape[1])
        bases = [ndi.zoom(b, zoom, order=1).astype(np.float32) for b in bases]
        support = ndi.zoom(support, zoom, order=1)
        high = ndi.zoom(high, zoom, order=1)
    return np.asarray(bases, dtype=np.float32), support > 0.5, high > 0.5, str(path)


def fit_nmf(
    data: np.ndarray,
    err: np.ndarray,
    invalid: np.ndarray,
    lw_source_mask: np.ndarray,
    basis_info,
) -> tuple[np.ndarray, dict[str, Any]]:
    zeros = np.zeros(data.shape, dtype=np.float32)
    if basis_info is None:
        return zeros, {"available": False, "amplitudes": [], "fit_pixels": 0}
    bases, support, high, basis_path = basis_info
    fit = high & ~invalid & ~lw_source_mask & np.isfinite(data) & np.isfinite(err) & (err > 0)
    sky_pixels = ~invalid & ~lw_source_mask & ~support & np.isfinite(data)
    sky_level = float(np.median(data[sky_pixels])) if sky_pixels.sum() >= 1000 else 0.0
    ids = np.flatnonzero(fit.ravel())
    if ids.size < 500:
        return zeros, {"available": True, "basis": basis_path, "amplitudes": [], "fit_pixels": int(ids.size), "error": "too few fit pixels"}
    ids = ids[:: max(1, int(math.ceil(ids.size / 120_000)))]
    matrix = bases.reshape(bases.shape[0], -1)[:, ids].T.astype(np.float64, copy=False)
    target = np.asarray(data.ravel()[ids], dtype=np.float64) - sky_level
    sigma = np.asarray(err.ravel()[ids], dtype=np.float64)
    weights = 1.0 / np.maximum(sigma, np.nanmedian(sigma) * 0.1)
    keep = np.isfinite(target) & np.isfinite(weights) & (weights > 0) & np.isfinite(matrix).all(axis=1)
    coef = np.zeros(bases.shape[0], dtype=np.float64)
    for _ in range(4):
        if keep.sum() < 500:
            break
        coef, _ = nnls(matrix[keep] * weights[keep, None], target[keep] * weights[keep], maxiter=10 * bases.shape[0])
        residual = target - matrix @ coef
        sig = robust_sigma(residual[keep])
        if not np.isfinite(sig) or sig <= 0:
            break
        center = float(np.median(residual[keep]))
        new_keep = np.isfinite(residual) & (np.abs(residual - center) < 4 * sig)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    model = np.einsum("k,kyx->yx", coef.astype(np.float32), bases, optimize=True)
    return model.astype(np.float32), {
        "available": True,
        "basis": basis_path,
        "amplitudes": [float(v) for v in coef],
        "fit_pixels": int(keep.sum()),
        "support_fraction": float(support.mean()),
        "sky_level": sky_level,
    }


def project_reference(
    ref_data: np.ndarray,
    ref_mask: np.ndarray,
    ref_header,
    target_header,
    target_shape: tuple[int, int],
    chunk_rows: int = 128,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Reproject one matching LW exposure and its Photutils mask onto one SW frame."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        src = WCS(ref_header).celestial
        dst = WCS(target_header).celestial
    if not src.has_celestial or not dst.has_celestial:
        raise ValueError("celestial WCS missing on LW reference or SW target")
    ny, nx = target_shape
    warped = np.full(target_shape, np.nan, dtype=np.float32)
    mask_out = np.zeros(target_shape, dtype=bool)
    covered = np.zeros(target_shape, dtype=bool)
    sy_n, sx_n = ref_data.shape
    for y0 in range(0, ny, chunk_rows):
        y1 = min(ny, y0 + chunk_rows)
        yy, xx = np.mgrid[y0:y1, 0:nx]
        ra, dec = dst.pixel_to_world_values(xx, yy)
        sx, sy = src.world_to_pixel_values(ra, dec)
        valid = np.isfinite(sx) & np.isfinite(sy) & (sx >= 0) & (sx <= sx_n - 1) & (sy >= 0) & (sy <= sy_n - 1)
        sampled = map_coordinates(
            np.nan_to_num(ref_data, nan=0.0),
            [np.nan_to_num(sy, nan=-1), np.nan_to_num(sx, nan=-1)],
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=False,
        )
        sampled_mask = map_coordinates(
            ref_mask.astype(np.uint8),
            [np.nan_to_num(sy, nan=-1), np.nan_to_num(sx, nan=-1)],
            order=0,
            mode="constant",
            cval=0,
            prefilter=False,
        ) != 0
        warped[y0:y1] = np.where(valid, sampled, np.nan)
        mask_out[y0:y1] = valid & sampled_mask
        covered[y0:y1] = valid
    return warped, mask_out, covered


def dynamic_wisp(
    short_sky_sub: np.ndarray,
    lw_warped: np.ndarray,
    lw_mask: np.ndarray,
    covered: np.ndarray,
    base_invalid_mask: np.ndarray,
    support: np.ndarray,
    nmf_model: np.ndarray,
    source_threshold: float = 0.995,
    scale_threshold: float = 0.995,
    clip_threshold: float = 0.9999,
    sigma_hi: float = 2.0,
    sigma_lo: float = 20.0,
) -> tuple[np.ndarray, dict[str, Any]]:
    """LW-referenced residual model with JUMPROPE's quantile/blur settings."""
    usable_ref = covered & np.isfinite(lw_warped) & support
    if usable_ref.sum() < 500:
        return nmf_model.copy(), {"available": False, "reason": "too little LW overlap"}
    ref_values = lw_warped[usable_ref]
    real_source = float(np.quantile(ref_values, source_threshold))
    sky_use = usable_ref & (lw_warped < real_source) & np.isfinite(short_sky_sub)
    sky_level = float(np.quantile(short_sky_sub[sky_use], 0.5)) if sky_use.sum() else 0.0
    ratio_pixels = usable_ref & (lw_warped > real_source) & np.isfinite(short_sky_sub)
    ratios = (short_sky_sub[ratio_pixels] - sky_level) / lw_warped[ratio_pixels]
    ratios = ratios[np.isfinite(ratios) & (ratios > 0)]
    if ratios.size < 20:
        return nmf_model.copy(), {"available": False, "reason": "too few positive high-reference ratios"}
    scale = float(np.quantile(ratios, scale_threshold))

    sample_ok = usable_ref & ~lw_mask & ~base_invalid_mask & np.isfinite(short_sky_sub)
    residual = short_sky_sub - sky_level - scale * np.nan_to_num(lw_warped, nan=0.0)
    vals = residual[sample_ok]
    if vals.size < 500:
        return nmf_model.copy(), {"available": False, "reason": "too few unmasked residual pixels", "scale": scale}
    lo = float(np.quantile(vals, 1.0 - clip_threshold))
    hi = float(np.quantile(vals, clip_threshold))
    sample_ok &= (residual >= lo) & (residual <= hi) & (residual > 0)
    # A normalized Gaussian avoids pulling masked galaxy flux into the dynamic map.
    numerator = ndi.gaussian_filter(np.where(sample_ok, residual, 0.0).astype(np.float32), sigma_hi)
    denominator = ndi.gaussian_filter(sample_ok.astype(np.float32), sigma_hi)
    dynamic = numerator / np.maximum(denominator, 1e-5)
    dynamic[denominator < 0.03] = 0.0
    dynamic[dynamic < 0] = 0.0
    broad = ndi.gaussian_filter(dynamic, sigma_lo)
    dynamic = np.maximum(dynamic, broad)
    # sigma_lo deliberately broadens the dynamic residual. Allow that low-pass
    # component to fade beyond the binary template support instead of cutting
    # it off at a hard mask edge.
    taper_margin = max(32, int(math.ceil(4 * sigma_lo)))
    soft_support = ndi.binary_dilation(support, iterations=taper_margin)
    dynamic[~soft_support] = 0.0

    # Fill source pixels and missing LW footprint using the independent NMF fit.
    projected_source = lw_mask & covered & support
    smooth_source_weight = ndi.gaussian_filter(projected_source.astype(np.float32), sigma_hi)
    smooth_source_weight[projected_source] = 1.0
    nmf_weight = np.clip(smooth_source_weight, 0.0, 1.0)
    dynamic_valid = covered & ~projected_source & soft_support & ~base_invalid_mask
    merged = np.where(dynamic_valid, dynamic, nmf_model)
    merged = (1.0 - nmf_weight) * merged + nmf_weight * nmf_model
    distance = ndi.distance_transform_edt(soft_support)
    taper = np.clip(distance / float(taper_margin), 0.0, 1.0)
    taper = taper * taper * (3.0 - 2.0 * taper)
    merged *= taper.astype(np.float32)
    merged[~soft_support] = 0.0
    return merged.astype(np.float32), {
        "available": True,
        "sky_level": sky_level,
        "real_source_threshold": real_source,
        "scale": scale,
        "ratio_pixels": int(ratios.size),
        "clip_lo": lo,
        "clip_hi": hi,
        "dynamic_fit_fraction": float(dynamic_valid.sum() / max(support.sum(), 1)),
        "source_mask_fraction_in_support": float(projected_source.sum() / max(support.sum(), 1)),
    }


def read_science(path: str):
    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdul:
        primary = hdul[0].header.copy()
        sci_name = "SCI_ORIG" if "SCI_ORIG" in hdul else "SCI"
        sci = np.asarray(hdul[sci_name].data, dtype=np.float32).copy()
        err = np.asarray(hdul["ERR"].data, dtype=np.float32).copy() if "ERR" in hdul else np.ones_like(sci)
        dq = np.asarray(hdul["DQ"].data, dtype=np.uint32).copy() if "DQ" in hdul else np.zeros_like(sci, dtype=np.uint32)
        sci_header = hdul["SCI"].header.copy()
        meta = {
            "source_extension": sci_name,
            "visit_id": str(primary.get("VISIT_ID", "")),
            "obs_id": str(primary.get("OBS_ID", "")),
            "detector": str(primary.get("DETECTOR", sci_header.get("DETECTOR", ""))).upper(),
            "filter": str(primary.get("FILTER", sci_header.get("FILTER", ""))).upper(),
            "module": "A" if str(primary.get("DETECTOR", "")).upper().startswith("NRCA") else "B",
        }
    return primary, sci_header, sci, err, dq, meta


def write_reference(path: Path, data: np.ndarray, mask: np.ndarray, header) -> None:
    hdus = fits.HDUList(
        [
            fits.PrimaryHDU(),
            fits.ImageHDU(np.asarray(data, dtype=np.float32), header=header.copy(), name="SCI"),
            fits.ImageHDU(np.asarray(mask, dtype=np.uint8), header=header.copy(), name="OBJMASK"),
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".writing")
    hdus.writeto(tmp, overwrite=True, output_verify="silentfix")
    os.replace(tmp, path)


def output_hdul(input_path: str, science: np.ndarray, diag: dict[str, np.ndarray], metadata: dict[str, Any]):
    with fits.open(input_path, memmap=False, lazy_load_hdus=True) as src:
        hdus = fits.HDUList()
        for hdu in src:
            name = hdu.name.upper()
            if name not in KEEP_HDU_NAMES:
                continue
            copied = hdu.copy()
            if name == "SCI":
                copied.data = np.asarray(science, dtype=np.float32)
                copied.header["JRPVERS"] = ("0.1-pilot", "standalone SkySurfIR cleanup")
                copied.header["JRPSRC"] = (metadata["source_extension"], "input SCI extension used")
                copied.header["JRP1F"] = (True, "ProFound scan-style 1/f correction")
                copied.header["JRPWISP"] = (metadata["wisp_method"][:20], "selected wisp estimator")
                copied.header["JRPSKY"] = (bool(metadata["sky_subtracted"]), "quadratic sky removed")
            elif name == "PRIMARY":
                copied.header["JRPVERS"] = "0.1-pilot"
                copied.header["JRPVISIT"] = str(metadata.get("visit_id", ""))[-16:]
            hdus.append(copied)
        map_header = src["SCI"].header.copy()

    for name, data in diag.items():
        dtype = np.uint8 if name in {"JROBJMSK", "JRLWMASK"} else np.float32
        h = map_header.copy()
        h["BUNIT"] = "1" if dtype is np.uint8 else str(map_header.get("BUNIT", ""))
        hdus.append(fits.ImageHDU(np.asarray(data, dtype=dtype), header=h, name=name))
    return hdus


def write_products(
    input_path: str,
    cal_path: Path,
    sky_path: Path,
    cal_science: np.ndarray,
    sky_science: np.ndarray,
    diag_cal: dict[str, np.ndarray],
    diag_sky: dict[str, np.ndarray],
    metadata: dict[str, Any],
    overwrite: bool,
) -> None:
    cal_path.parent.mkdir(parents=True, exist_ok=True)
    sky_path.parent.mkdir(parents=True, exist_ok=True)
    for path, science, diag, is_sky in (
        (cal_path, cal_science, diag_cal, False),
        (sky_path, sky_science, diag_sky, True),
    ):
        if path.exists() and not overwrite:
            continue
        file_meta = dict(metadata)
        file_meta["sky_subtracted"] = is_sky
        hdul = output_hdul(input_path, science, diag, file_meta)
        tmp = path.with_name(path.name + ".writing")
        hdul.writeto(tmp, overwrite=True, output_verify="silentfix")
        hdul.close()
        os.replace(tmp, path)


def process_longwave(path: str, cfg: dict[str, Any]) -> dict[str, Any]:
    started = time.time()
    primary, sci_header, image, err, dq, meta = read_science(path)
    invalid = base_invalid(dq, image)
    corrected, model1f, fstats = scan_1f(
        image,
        invalid,
        cfg["block_rows"],
        cfg["trend_block"],
        cfg["clip_hi"],
        cfg["scan_quantile"],
    )
    mask1, n1, sig1, _ = source_segmentation(corrected, err, invalid, cfg["detect_sigma_first"])
    sky0, skymeta0 = sky_polynomial(corrected, invalid, mask1)
    mask2, n2, sig2, _ = source_segmentation(
        corrected - sky0, err, invalid, cfg["detect_sigma_faint"], background=0.0
    )
    objmask = mask1 | mask2
    sky, skymeta = sky_polynomial(corrected, invalid, objmask)
    residual = corrected - sky
    good_sky = ~invalid & ~objmask & np.isfinite(residual)
    pedestal = float(np.median(residual[good_sky])) if good_sky.any() else 0.0
    sky_total = sky + np.float32(pedestal)
    final = corrected - sky_total
    basename = Path(path).name
    out_root = Path(cfg["output_dir"])
    cal_path = out_root / "cal" / basename
    sky_path = out_root / "cal_sky_renorm" / basename
    diag_cal = {"JR1FMOD": model1f, "JROBJMSK": objmask}
    diag_sky = {**diag_cal, "JRSKYMOD": sky_total}
    metadata = {
        **meta,
        "source_extension": meta["source_extension"],
        "wisp_method": "none_longwave_reference",
        "sky_subtracted": True,
        "source_mask_fraction": float(objmask.mean()),
        "source_count_first": n1,
        "source_count_faint": n2,
        "noise_first": sig1,
        "noise_faint": sig2,
        "sky_first_fit": skymeta0,
        "sky_final_fit": skymeta,
        "pedestal_removed": pedestal,
        "one_over_f": fstats,
        "seconds": time.time() - started,
    }
    write_products(path, cal_path, sky_path, corrected, final, diag_cal, diag_sky, metadata, cfg["overwrite"])
    ref_path = Path(cfg["ref_dir"]) / basename
    write_reference(ref_path, final, objmask, sci_header)
    return {"file": basename, "kind": "LW", "status": "ok", "seconds": metadata["seconds"], **metadata}


def choose_reference(meta: dict[str, Any], refs: list[dict[str, Any]], target_header) -> dict[str, Any] | None:
    # Exact matched observation + exposure has the cleanest common footprint.
    stem = Path(meta["basename"]).name
    module_long = "nrcalong" if meta["module"] == "A" else "nrcblong"
    exact_stem = re.sub(r"_nrc[a-z0-9]+_cal\.fits$", f"_{module_long}_cal.fits", stem, flags=re.IGNORECASE)
    for ref in refs:
        if Path(ref["path"]).name == exact_stem:
            return ref
    same_obs = [r for r in refs if r["obs_id"] == meta["obs_id"] and r["module"] == meta["module"]]
    candidates = same_obs or [
        r for r in refs if r["visit_id"] == meta["visit_id"] and r["module"] == meta["module"]
    ]
    if not candidates:
        return None
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        dst = WCS(target_header).celestial
        yy, xx = np.mgrid[0:2048:17j, 0:2048:17j]
        ra, dec = dst.pixel_to_world_values(xx, yy)
        best = None
        best_fraction = -1.0
        for ref in candidates:
            with fits.open(ref["path"], memmap=False, lazy_load_hdus=True) as hdul:
                src = WCS(hdul["SCI"].header).celestial
                shape = hdul["SCI"].shape
            sx, sy = src.world_to_pixel_values(ra, dec)
            inside = np.isfinite(sx) & np.isfinite(sy) & (sx >= 0) & (sx < shape[1]) & (sy >= 0) & (sy < shape[0])
            fraction = float(inside.mean())
            if fraction > best_fraction:
                best = ref
                best_fraction = fraction
    return best if best_fraction > 0.01 else None


def process_shortwave(path: str, cfg: dict[str, Any], refs: list[dict[str, Any]]) -> dict[str, Any]:
    started = time.time()
    primary, sci_header, image, err, dq, meta = read_science(path)
    invalid = base_invalid(dq, image)
    corrected1, map1, fstats1 = scan_1f(
        image,
        invalid,
        cfg["block_rows"],
        cfg["trend_block"],
        cfg["clip_hi"],
        cfg["scan_quantile"],
    )
    mask1, n1, sig1, _ = source_segmentation(corrected1, err, invalid, cfg["detect_sigma_first"])
    sky0, skymeta0 = sky_polynomial(corrected1, invalid, mask1)
    sky_sub0 = corrected1 - sky0
    mask2, n2, sig2, _ = source_segmentation(
        sky_sub0, err, invalid, cfg["detect_sigma_faint"], background=0.0
    )
    objmask = mask1 | mask2
    sky_pre, skymeta_pre = sky_polynomial(corrected1, invalid, objmask)
    sky_sub = corrected1 - sky_pre

    detector = meta["detector"]
    filt = meta["filter"]
    is_wisp_detector = detector in WISP_DETECTORS and filt in WISP_FILTERS
    basis_info = load_nmf_basis(cfg["templates_dir"], detector, filt) if is_wisp_detector else None
    support = basis_info[1] if basis_info is not None else np.zeros(image.shape, dtype=bool)
    selected_ref = choose_reference({**meta, "basename": Path(path).name}, refs, sci_header) if is_wisp_detector else None
    lw_source = np.zeros(image.shape, dtype=bool)
    lw_coverage = np.zeros(image.shape, dtype=bool)
    lw_ref = np.full(image.shape, np.nan, dtype=np.float32)
    ref_name = None
    if selected_ref is not None and is_wisp_detector:
        with fits.open(selected_ref["path"], memmap=False) as rh:
            ref_data = np.asarray(rh["SCI"].data, dtype=np.float32)
            ref_mask = np.asarray(rh["OBJMASK"].data != 0)
            ref_header = rh["SCI"].header.copy()
        lw_ref, lw_source, lw_coverage = project_reference(
            ref_data, ref_mask, ref_header, sci_header, image.shape
        )
        ref_name = Path(selected_ref["path"]).name

    if basis_info is None:
        nmf_model = np.zeros(image.shape, dtype=np.float32)
        nmf_meta = {"available": False, "amplitudes": [], "fit_pixels": 0}
    else:
        # Match the NMF method's source protection: use the actual LW
        # segmentation projected from this exposure, not an SW-derived mask
        # that could have segmented the wisp itself as a source.
        fit_protect = lw_source if selected_ref is not None else objmask
        nmf_model, nmf_meta = fit_nmf(sky_sub, err, invalid, fit_protect, basis_info)

    if is_wisp_detector and basis_info is not None and selected_ref is not None:
        dynamic_model, dynamic_meta = dynamic_wisp(
            sky_sub,
            lw_ref,
            lw_source,
            lw_coverage,
            invalid,
            support,
            nmf_model,
            cfg["source_threshold"],
            cfg["scale_threshold"],
            cfg["clip_threshold"],
            cfg["sigma_hi"],
            cfg["sigma_lo"],
        )
        wisp_method = "lw_dynamic_with_nmf_source_fill" if dynamic_meta.get("available") else "nmf_fallback"
    elif is_wisp_detector and basis_info is not None:
        dynamic_model = nmf_model.copy()
        dynamic_meta = {"available": False, "reason": "no long-wave reference; NMF fallback"}
        wisp_method = "nmf_fallback"
    else:
        dynamic_model = np.zeros(image.shape, dtype=np.float32)
        dynamic_meta = {"available": False, "reason": "detector/filter has no enabled wisp model"}
        wisp_method = "none"

    # The second scan pass follows JUMPROPE's post-wisp correction. Preserve
    # Photutils objects and actual LW objects while estimating residual stripes.
    protected = invalid | objmask | lw_source
    corrected2, map2, fstats2 = scan_1f(
        corrected1 - dynamic_model,
        protected,
        cfg["block_rows"],
        cfg["trend_block"],
        cfg["clip_hi"],
        cfg["scan_quantile"],
    )
    final_resid = corrected2 - sky_pre
    mask3, n3, sig3, _ = source_segmentation(
        final_resid, err, invalid, cfg["detect_sigma_final"], background=0.0
    )
    final_mask = objmask | mask3
    final_sky, skymeta_final = sky_polynomial(corrected2, invalid, final_mask)
    sky_resid = corrected2 - final_sky
    sky_good = ~invalid & ~final_mask & np.isfinite(sky_resid)
    pedestal = float(np.median(sky_resid[sky_good])) if sky_good.any() else 0.0
    sky_total = final_sky + np.float32(pedestal)
    final = corrected2 - sky_total
    map_total = map1 + map2

    basename = Path(path).name
    out_root = Path(cfg["output_dir"])
    cal_path = out_root / "cal" / basename
    sky_path = out_root / "cal_sky_renorm" / basename
    diag_cal = {
        "JR1FMOD": map_total,
        "JRWISPMD": dynamic_model,
        "JRNMFMOD": nmf_model,
        "JROBJMSK": final_mask,
        "JRLWMASK": lw_source,
    }
    diag_sky = {**diag_cal, "JRSKYMOD": sky_total}
    metadata = {
        **meta,
        "source_extension": meta["source_extension"],
        "wisp_method": wisp_method,
        "sky_subtracted": True,
        "longwave_reference": ref_name,
        "longwave_coverage_fraction": float(lw_coverage.mean()),
        "longwave_source_mask_fraction": float(lw_source.mean()),
        "source_mask_fraction": float(final_mask.mean()),
        "source_count_first": n1,
        "source_count_faint": n2,
        "source_count_final": n3,
        "noise_first": sig1,
        "noise_faint": sig2,
        "noise_final": sig3,
        "sky_initial_fit": skymeta0,
        "sky_masked_fit": skymeta_pre,
        "sky_final_fit": skymeta_final,
        "pedestal_removed": pedestal,
        "one_over_f_first": fstats1,
        "one_over_f_post_wisp": fstats2,
        "nmf_candidate": nmf_meta,
        "dynamic_candidate": dynamic_meta,
        "wisp_model_abs_p99": float(np.quantile(np.abs(dynamic_model), 0.99)),
        "wisp_model_peak": float(np.nanmax(np.abs(dynamic_model))),
        "seconds": time.time() - started,
    }
    write_products(path, cal_path, sky_path, corrected2, final, diag_cal, diag_sky, metadata, cfg["overwrite"])
    return {"file": basename, "kind": "SW", "status": "ok", "seconds": metadata["seconds"], **metadata}


def select_inputs(args) -> list[str]:
    inputs = []
    selector = str(args.visit_id or args.obs_id or "").strip()
    selector_digits = re.sub(r"\D", "", selector)
    for path in sorted(Path(args.input_dir).glob("*_cal.fits")):
        try:
            hdr = fits.getheader(path, 0)
        except Exception:
            continue
        visit = str(hdr.get("VISIT_ID", ""))
        obs = str(hdr.get("OBS_ID", ""))
        if args.visit_id:
            digits = re.sub(r"\D", "", visit)
            if selector_digits not in digits:
                continue
        elif args.obs_id:
            if selector.lower() not in obs.lower() and selector_digits not in re.sub(r"\D", "", obs):
                continue
        inputs.append(str(path))
    return inputs


def write_json(path: Path, value: Any):
    def clean(item):
        if isinstance(item, dict):
            return {str(k): clean(v) for k, v in item.items()}
        if isinstance(item, (list, tuple)):
            return [clean(v) for v in item]
        if isinstance(item, (np.bool_, bool)):
            return bool(item)
        if isinstance(item, (np.integer,)):
            return int(item)
        if isinstance(item, (np.floating, float)):
            val = float(item)
            return val if np.isfinite(val) else None
        if isinstance(item, np.ndarray):
            return clean(item.tolist())
        return item

    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".writing")
    tmp.write_text(json.dumps(clean(value), indent=2, allow_nan=False, default=str))
    os.replace(tmp, path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    selector = parser.add_mutually_exclusive_group(required=True)
    selector.add_argument("--visit-id", help="visit/observation prefix, e.g. 1813005")
    selector.add_argument("--obs-id", help="OBS_ID or its digit string, e.g. V01813005001P0000000002101")
    parser.add_argument("--input-dir", default=str(DEFAULT_INPUT))
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--templates-dir", default=str(DEFAULT_TEMPLATES))
    parser.add_argument("--cores", type=int, default=64)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--block-rows", type=int, default=512)
    parser.add_argument("--trend-block", type=int, default=501)
    parser.add_argument("--clip-hi", type=float, default=0.90)
    parser.add_argument("--scan-quantile", type=float, default=0.40)
    parser.add_argument("--detect-sigma-first", type=float, default=2.8)
    parser.add_argument("--detect-sigma-faint", type=float, default=1.8)
    parser.add_argument("--detect-sigma-final", type=float, default=1.8)
    parser.add_argument("--source-threshold", type=float, default=0.995)
    parser.add_argument("--scale-threshold", type=float, default=0.995)
    parser.add_argument("--clip-threshold", type=float, default=0.9999)
    parser.add_argument("--sigma-hi", type=float, default=2.0)
    parser.add_argument("--sigma-lo", type=float, default=20.0)
    args = parser.parse_args(argv)
    if args.cores < 1:
        parser.error("--cores must be positive")
    if not Path(args.input_dir).is_dir():
        parser.error(f"input directory not found: {args.input_dir}")
    if not Path(args.templates_dir).is_dir():
        parser.error(f"wisp templates directory not found: {args.templates_dir}")

    files = select_inputs(args)
    if not files:
        print("No matching stage2 *_cal.fits inputs", file=sys.stderr)
        return 2
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ref_dir = output_dir / "_lw_reference_cache"
    ref_dir.mkdir(parents=True, exist_ok=True)

    lw_files: list[str] = []
    sw_files: list[str] = []
    metadata_by_file = []
    for f in files:
        h = fits.getheader(f, 0)
        det = str(h.get("DETECTOR", "")).upper()
        row = {
            "path": f,
            "basename": Path(f).name,
            "visit_id": str(h.get("VISIT_ID", "")),
            "obs_id": str(h.get("OBS_ID", "")),
            "detector": det,
            "module": "A" if det.startswith("NRCA") else "B",
        }
        metadata_by_file.append(row)
        if det.endswith("LONG"):
            lw_files.append(f)
        else:
            sw_files.append(f)

    cfg = {
        "output_dir": str(output_dir),
        "ref_dir": str(ref_dir),
        "templates_dir": str(Path(args.templates_dir).resolve()),
        "overwrite": args.overwrite,
        "block_rows": args.block_rows,
        "trend_block": args.trend_block,
        "clip_hi": args.clip_hi,
        "scan_quantile": args.scan_quantile,
        "detect_sigma_first": args.detect_sigma_first,
        "detect_sigma_faint": args.detect_sigma_faint,
        "detect_sigma_final": args.detect_sigma_final,
        "source_threshold": args.source_threshold,
        "scale_threshold": args.scale_threshold,
        "clip_threshold": args.clip_threshold,
        "sigma_hi": args.sigma_hi,
        "sigma_lo": args.sigma_lo,
    }
    write_json(output_dir / "run_config.json", {**cfg, "selector": args.visit_id or args.obs_id, "input_count": len(files)})
    print(
        f"Selected {len(files)} stage2 cal files ({len(lw_files)} LW, {len(sw_files)} SW); "
        f"using {args.cores} processes; output={output_dir}",
        flush=True,
    )
    print(
        f"Wisp settings: source={args.source_threshold}, scale={args.scale_threshold}, "
        f"clip={args.clip_threshold}, sigma_hi={args.sigma_hi}, sigma_lo={args.sigma_lo}; "
        "1/f: q90 clipping, q0.4 scan, 512-row blocks, trend=501",
        flush=True,
    )
    print("Phase 1/2: calibrating long-wave references and outputs", flush=True)
    results: list[dict[str, Any]] = []
    failures: list[dict[str, str]] = []
    references: list[dict[str, Any]] = []
    with cf.ProcessPoolExecutor(max_workers=args.cores) as pool:
        future_to_file = {pool.submit(process_longwave, f, cfg): f for f in lw_files}
        done = 0
        for fut in cf.as_completed(future_to_file):
            f = future_to_file[fut]
            done += 1
            try:
                row = fut.result()
                results.append(row)
                references.append(
                    {
                        "path": str(ref_dir / Path(f).name),
                        "visit_id": row["visit_id"],
                        "obs_id": row["obs_id"],
                        "module": row["module"],
                    }
                )
                if done % 8 == 0 or done == len(lw_files):
                    print(f"  longwave {done}/{len(lw_files)} complete", flush=True)
            except Exception as exc:
                failures.append({"file": f, "error": repr(exc), "traceback": traceback.format_exc()})
                print(f"  FAILED longwave {Path(f).name}: {exc}", flush=True)
    print("Phase 2/2: processing short-wave exposures with individual LW matches", flush=True)
    with cf.ProcessPoolExecutor(max_workers=args.cores) as pool:
        future_to_file = {pool.submit(process_shortwave, f, cfg, references): f for f in sw_files}
        done = 0
        for fut in cf.as_completed(future_to_file):
            f = future_to_file[fut]
            done += 1
            try:
                row = fut.result()
                results.append(row)
                if done % 16 == 0 or done == len(sw_files):
                    print(f"  shortwave {done}/{len(sw_files)} complete", flush=True)
            except Exception as exc:
                failures.append({"file": f, "error": repr(exc), "traceback": traceback.format_exc()})
                print(f"  FAILED shortwave {Path(f).name}: {exc}", flush=True)

    results.sort(key=lambda x: x["file"])
    write_json(output_dir / "metrics.json", {"results": results, "failures": failures})
    if failures:
        write_json(output_dir / "failures.json", failures)
    status = "complete" if not failures else "partial"
    print(f"Run {status}: {len(results)} files succeeded, {len(failures)} failed", flush=True)
    print(f"Products: {output_dir / 'cal'} and {output_dir / 'cal_sky_renorm'}", flush=True)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
