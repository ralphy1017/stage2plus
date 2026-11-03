"""Science steps: object masks, dither 1/f, sky, wisp and pedestal.

Each function works on in-memory arrays of one dither group or exposure.
"""

from __future__ import annotations

import math
import warnings
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.wcs import FITSFixedWarning, WCS
from photutils.background import Background2D, BiweightLocationBackground
from photutils.segmentation import detect_sources
from scipy import ndimage as ndi
from scipy.optimize import nnls

# HDUs copied from each input into the output FITS.
KEEP_HDUS = ("PRIMARY", "SCI", "ERR", "DQ", "AREA", "VAR_POISSON",
             "VAR_RNOISE", "VAR_FLAT", "ASDF")

SEGMENT_NSIGMA = 2.8
MASK_GROWTH_CAP = 32       # pixels
WISP_SUPPORT_GROWTH = 8    # pixels excluded around the wisp template support
MAX_HELDOUT = 2            # dithers withheld from each 1/f scene for QC
MIN_DITHER_SCENE_PEERS = 2
SCAN_BLOCK_ROWS = 512
SCAN_TREND_BLOCK = 501

# JADES per-exposure background: photutils.Background2D with these settings.
BKG_BOX = (256, 256)
BKG_FILTER = (3, 3)


@dataclass
class Exposure:
    path: Path
    science: np.ndarray
    err: np.ndarray
    wcs: WCS
    sci_header: fits.Header
    detector: str
    band: str
    invalid: np.ndarray


def read_exposure(path: Path) -> Exposure:
    """Read SCI_ORIG when present (the uncorrected science), else SCI."""
    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdul:
        science = np.array(hdul["SCI_ORIG" if "SCI_ORIG" in hdul else "SCI"].data,
                           dtype=np.float32)
        err = np.array(hdul["ERR"].data, dtype=np.float32)
        dq = np.asarray(hdul["DQ"].data, dtype=np.uint32)
        sci_header = hdul["SCI"].header.copy()
        detector = str(hdul[0].header["DETECTOR"]).upper()
        band = str(hdul[0].header["FILTER"]).upper()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(sci_header)
    # JWST DO_NOT_USE, grown by one pixel, plus unusable values.
    bad = ndi.maximum_filter((dq & 1) != 0, size=3, mode="nearest")
    invalid = bad | ~np.isfinite(science) | (science == 0) | ~np.isfinite(err) | (err <= 0)
    return Exposure(path, science, err, wcs, sci_header, detector, band, invalid)


def robust_sigma(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(values - np.median(values))))


# Object masks ---------------------------------------------------------------

def segment(exp: Exposure, nsigma: float = SEGMENT_NSIGMA) -> np.ndarray:
    """Photutils segmentation footprint above nsigma * ERR over the median."""
    valid = ~exp.invalid
    if valid.sum() < 1000:
        return np.zeros(exp.science.shape, bool)
    sample = exp.science[valid][::max(1, int(valid.sum() // 250_000))]
    sigma = robust_sigma(sample)
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.median(exp.err[valid]))
    threshold = np.maximum(exp.err, max(0.5 * sigma, 1e-12)) * nsigma + np.median(sample)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        # Minimum of 5 pixels, passed positionally (renamed in photutils 3).
        seg = detect_sources(exp.science, threshold, 5, connectivity=8, mask=exp.invalid)
    return np.zeros(exp.science.shape, bool) if seg is None else seg.data != 0


@lru_cache(maxsize=32)
def _disk(radius: int) -> np.ndarray:
    yy, xx = np.ogrid[-radius:radius + 1, -radius:radius + 1]
    return xx * xx + yy * yy <= radius * radius


def grow_mask(mask: np.ndarray, cap: int = MASK_GROWTH_CAP) -> np.ndarray:
    """Grow each object by its equivalent radius (at most cap pixels)."""
    labels, _ = ndi.label(mask)
    areas = np.bincount(labels.ravel())
    grown = mask.copy()
    for index, box in enumerate(ndi.find_objects(labels), start=1):
        if box is None:
            continue
        radius = min(cap, max(1, math.ceil(math.sqrt(areas[index] / math.pi))))
        window = tuple(slice(max(0, s.start - radius), min(n, s.stop + radius))
                       for s, n in zip(box, mask.shape))
        grown[window] |= ndi.binary_dilation(labels[window] == index,
                                             structure=_disk(radius))
    return grown


def project_mask(mask: np.ndarray, mask_header: fits.Header, target: Exposure) -> np.ndarray:
    """Nearest-neighbour reprojection of a mask onto the target's pixels."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        source = WCS(mask_header).celestial
    destination = target.wcs.celestial
    ny, nx = target.science.shape
    out = np.zeros((ny, nx), bool)
    for y0 in range(0, ny, 128):
        yy, xx = np.mgrid[y0:min(y0 + 128, ny), :nx]
        sx, sy = source.world_to_pixel_values(*destination.pixel_to_world_values(xx, yy))
        inside = (np.isfinite(sx) & np.isfinite(sy) & (sx >= 0) & (sy >= 0) &
                  (sx < mask.shape[1]) & (sy < mask.shape[0]))
        sampled = ndi.map_coordinates(mask.astype(np.uint8),
                                      [np.nan_to_num(sy, nan=-1), np.nan_to_num(sx, nan=-1)],
                                      order=0, mode="constant", cval=0, prefilter=False)
        out[y0:y0 + 128] = inside & (sampled != 0)
    return out


def protect_wisp(mask: np.ndarray, high_wisp: np.ndarray, min_row_pixels: int = 64) -> np.ndarray:
    """Add the bright-wisp region unless it empties an amplifier row."""
    out = mask | high_wisp
    width = mask.shape[1] // 4
    for lo in range(0, mask.shape[1], width):
        sparse = np.count_nonzero(~out[:, lo:lo + width], axis=1) < min_row_pixels
        out[sparse, lo:lo + width] = mask[sparse, lo:lo + width]
    return out


# Registered-dither 1/f ------------------------------------------------------

def _aligned(peer: Exposure, target: Exposure, good: np.ndarray) -> np.ndarray:
    """Peer shifted onto the target's pixels, its scalar level matched."""
    ny, nx = target.science.shape
    xy = np.array([[nx / 2, ny / 2], [100, 100], [nx - 100, 100],
                   [100, ny - 100], [nx - 100, ny - 100]], float)
    delta = peer.wcs.all_world2pix(target.wcs.all_pix2world(xy, 0), 0) - xy
    if np.max(np.linalg.norm(delta - delta[0], axis=1)) > 0.1:
        raise ValueError(f"WCS of {peer.path.name} differs by more than a shift")
    # scipy's shift samples input[p - shift] at output p.
    aligned = ndi.shift(np.where(peer.invalid, np.nan, peer.science),
                        (-delta[0, 1], -delta[0, 0]), order=1, mode="constant",
                        cval=np.nan, prefilter=False, output=np.float32)
    sample = np.s_[::8, ::8]
    valid = good[sample] & np.isfinite(aligned[sample])
    if valid.sum() < 10000:
        raise ValueError(f"Too little dither overlap with {peer.path.name}")
    aligned -= np.float32(np.median((aligned[sample] - target.science[sample])[valid]))
    return aligned


def fit_stripes(target: np.ndarray, scene: np.ndarray,
                good: np.ndarray) -> tuple[np.ndarray, dict]:
    """Amplifier-row plus column medians of target - scene, zero mean on sky."""
    residual = target - scene
    valid = good & np.isfinite(residual)
    centre = float(np.median(residual[valid][::8]))
    sigma = robust_sigma(residual[valid][::8])
    valid &= np.abs(residual - centre) < 6 * sigma
    if valid.mean() < 0.5:
        raise ValueError("Too few sky pixels to fit 1/f")
    residual -= np.float32(centre)
    residual[~valid] = np.nan
    ny, nx = target.shape
    rows = np.empty_like(target)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for lo in range(0, nx, nx // 4):
            profile = np.nanmedian(residual[:, lo:lo + nx // 4], axis=1)
            finite = np.isfinite(profile)
            if finite.sum() < 0.9 * ny:
                raise ValueError("Too many masked rows in an amplifier")
            rows[:, lo:lo + nx // 4] = np.interp(np.arange(ny), np.flatnonzero(finite),
                                                 profile[finite])[:, None]
        columns = np.nan_to_num(np.nanmedian(residual - rows, axis=0), nan=0)
    model = (rows + columns[None, :]).astype(np.float32)
    # The fit cannot see a global offset: keep the input's masked sky median.
    model -= np.float32(np.median(model[good]))
    model += np.float32(np.median((target - model)[good]) - np.median(target[good]))
    model[~np.isfinite(target)] = 0
    return model, {"residual_sigma": sigma, "fit_fraction": float(valid.mean()),
                   "model_rms": float(np.sqrt(np.mean(model[good] ** 2)))}


def _row_rms(diff: np.ndarray, good: np.ndarray) -> float:
    valid = good & np.isfinite(diff)
    scalar = float(np.median(diff[valid][::8]))
    width = diff.shape[1] // 4
    rows = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for lo in range(0, diff.shape[1], width):
            part = np.where(valid[:, lo:lo + width], diff[:, lo:lo + width], np.nan)
            rows.append(np.nanmedian(part, axis=1) - scalar)
    return float(np.sqrt(np.nanmean(np.square(rows))))


def _tile_rms(diff: np.ndarray, good: np.ndarray, size: int = 64) -> float:
    valid = good & np.isfinite(diff)
    scalar = float(np.median(diff[valid][::8]))
    ny, nx = diff.shape
    tiles = np.where(valid, diff - scalar, np.nan).reshape(ny // size, size, nx // size, size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return float(np.sqrt(np.nanmean(np.nanmedian(tiles, axis=(1, 3)) ** 2)))


def dither_onef(target: Exposure, peers: list[Exposure],
                good: np.ndarray) -> tuple[np.ndarray, dict]:
    """1/f model of target from a registered median scene of its other dithers.

    Up to MAX_HELDOUT peers (keeping at least two in the scene) are withheld.
    QC passes when their mean amplifier-row and 64-pixel tile RMS against the
    target both improve by more than 5% and no withheld peer worsens by more
    than 5%. A QC failure returns its candidate and diagnostics; the caller
    must select the single-frame fallback instead of applying that candidate.
    """
    if len(peers) < MIN_DITHER_SCENE_PEERS + 1:
        raise ValueError("Dither 1/f requires at least two scene peers and one held-out peer")
    n_fit = max(2, len(peers) - MAX_HELDOUT)
    fit_peers, heldout = peers[:n_fit], peers[n_fit:]
    try:
        stack = np.stack([_aligned(p, target, good) for p in fit_peers])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            scene = np.nanmedian(stack, axis=0, overwrite_input=True)
        del stack
        model, meta = fit_stripes(target.science, scene, good)
        del scene
        corrected = target.science - model
        checks = []
        for peer in heldout:
            aligned = _aligned(peer, target, good)
            checks.append({"peer": peer.path.name,
                           "row_before": _row_rms(target.science - aligned, good),
                           "row_after": _row_rms(corrected - aligned, good),
                           "tile_before": _tile_rms(target.science - aligned, good),
                           "tile_after": _tile_rms(corrected - aligned, good)})
    except (ValueError, FloatingPointError) as exc:
        raise ValueError(f"Dither 1/f fit failed for {target.path.name}: {exc}") from exc
    mean = {k: float(np.mean([c[k] for c in checks]))
            for k in ("row_before", "row_after", "tile_before", "tile_after")}
    accepted = (mean["row_after"] < 0.95 * mean["row_before"] and
                mean["tile_after"] < 0.95 * mean["tile_before"] and
                all(c["row_after"] < 1.05 * c["row_before"] and
                    c["tile_after"] < 1.05 * c["tile_before"] for c in checks))
    return model, {"scene_dithers": [p.path.name for p in fit_peers], **meta,
                   "qc_accepted": bool(accepted), **mean, "heldout": checks}


def _scan_1f_xy(image: np.ndarray, invalid: np.ndarray,
                block_rows: int = SCAN_BLOCK_ROWS,
                trend_block: int = SCAN_TREND_BLOCK,
                clip_hi: float = 0.90,
                sky_quantile: float = 0.40) -> tuple[np.ndarray, dict]:
    """Vectorized ProFound-style scan model in detector XY orientation.

    Estimate column profiles in row blocks and row profiles across the full
    width. Keep their centered 40th-percentile structure after removing a
    501-pixel running-median trend, matching the Jumprope keep_trend=True
    configuration. The upper 10% of unmasked values are excluded from the
    estimates. Masked science pixels are not changed by this model.
    """
    original = np.asarray(image, dtype=np.float32)
    bad = np.asarray(invalid, dtype=bool) | ~np.isfinite(original) | (original == 0)
    good = ~bad
    if good.sum() < 1000:
        raise ValueError("Too few usable pixels for single-frame 1/f scan")
    clip_value = float(np.quantile(original[good], clip_hi))
    # Input SCI is float32; retain that precision and halve the full-frame
    # working copy compared with float64.
    scan = original.copy()
    scan[bad | (scan > clip_value)] = np.nan
    ny, nx = scan.shape

    row_slices = []
    row_profiles = []
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for y0 in range(0, ny, block_rows):
            y1 = min(ny, y0 + block_rows)
            row_slices.append((y0, y1))
            row_profiles.append(np.nanquantile(scan[y0:y1], sky_quantile, axis=0))
        row_profiles = np.asarray(row_profiles, dtype=np.float64)
        if not np.isfinite(row_profiles).any():
            raise ValueError("No valid row profiles for single-frame 1/f scan")
        row_center = float(np.nanmedian(row_profiles))

    row_map = np.zeros((ny, nx), dtype=np.float32)
    wx = min(trend_block, nx)
    if wx % 2 == 0:
        wx -= 1
    for ib, (y0, y1) in enumerate(row_slices):
        values = row_profiles[ib] - row_center
        finite = np.isfinite(values)
        if finite.any():
            filled = np.interp(np.arange(nx), np.flatnonzero(finite), values[finite])
            trend = ndi.median_filter(filled, size=wx, mode="nearest") if wx > 1 else 0.0
            row_map[y0:y1] = np.where(finite, values - trend, 0.0)[None, :]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        column_profile = np.nanquantile(scan, sky_quantile, axis=1)
        if not np.isfinite(column_profile).any():
            raise ValueError("No valid column profile for single-frame 1/f scan")
        column_center = float(np.nanmedian(column_profile))
    centered = column_profile - column_center
    finite_column = np.isfinite(centered)
    wy = min(trend_block, ny)
    if wy % 2 == 0:
        wy -= 1
    if finite_column.any():
        filled_column = np.interp(np.arange(ny), np.flatnonzero(finite_column),
                                  centered[finite_column])
        trend_column = (ndi.median_filter(filled_column, size=wy, mode="nearest")
                        if wy > 1 else 0.0)
        column_map = np.where(finite_column, centered - trend_column, 0.0).astype(np.float32)
    else:
        column_map = np.zeros(ny, dtype=np.float32)

    model = row_map + column_map[:, None]
    model[bad] = 0.0
    stats = {"algorithm": "PROFOUND_SCAN_KEEP_TREND",
             "clip_hi": clip_hi, "clip_value": clip_value,
             "sky_quantile": sky_quantile, "block_rows": block_rows,
             "trend_block": trend_block, "keep_trend": True,
             "row_profile_coverage": float(np.isfinite(row_profiles).mean()),
             "column_profile_coverage": float(finite_column.mean()),
             "fit_fraction": float(np.isfinite(scan).mean()),
             "model_rms": float(np.sqrt(np.mean(model[good].astype(np.float64) ** 2)))}
    return model, stats


def scan_onef(image: np.ndarray, invalid: np.ndarray) -> tuple[np.ndarray, dict]:
    """Single-frame backup using Jumprope's keep_trend=True scan orientation."""
    model_xy, meta = _scan_1f_xy(np.asarray(image).T, np.asarray(invalid).T)
    return np.asarray(model_xy.T, dtype=np.float32), meta


# Sky, wisp and pedestal -----------------------------------------------------

def masked_median(image: np.ndarray, good: np.ndarray) -> float:
    if good.sum() < 10000:
        raise ValueError("Too few unmasked pixels for a sky median")
    return float(np.median(image[good]))


def background2d(image: np.ndarray, mask: np.ndarray,
                 exclude_percentile: float) -> tuple[np.ndarray, dict]:
    """JADES-style residual background: Background2D, subtracted once."""
    bkg = Background2D(np.where(mask, 0, image), BKG_BOX, mask=mask,
                       filter_size=BKG_FILTER, sigma_clip=SigmaClip(),
                       bkg_estimator=BiweightLocationBackground(),
                       exclude_percentile=exclude_percentile)
    model = np.asarray(bkg.background, np.float32)
    used = getattr(bkg, "n_pixels_mesh", None)  # photutils >= 3
    used = bkg.npixels_mesh if used is None else used
    minimum = (1 - exclude_percentile / 100) * np.prod(BKG_BOX)
    return model, {"exclude_percentile": exclude_percentile,
                   "boxes": int(used.size),
                   "boxes_excluded": int(np.count_nonzero(used <= minimum)),
                   "model_p1_p99": [float(v) for v in np.percentile(model[::8, ::8], [1, 99])]}


def basis_path(root: Path, detector: str, band: str) -> Path | None:
    det, filt = detector.lower(), band.lower()
    for directory in (root / f"{det}_org", root / det):
        for extension in (".fits", ".fits.gz"):
            path = directory / f"{det}_{filt}_wisp{extension}"
            if path.is_file():
                return path
    return None


@lru_cache(maxsize=1)
def load_basis(path: Path, shape: tuple[int, int]):
    """NMF wisp basis, support mask and bright-wisp mask at detector size."""
    with fits.open(path, memmap=False) as hdul:
        basis = [np.array(hdul[f"WISP_{i}"].data, np.float32)
                 for i in range(10) if f"WISP_{i}" in hdul]
        support = np.asarray(hdul["MASK"].data, bool)
        high = np.asarray(hdul["MASK_hSNR"].data, bool)
    if not basis:
        raise ValueError(f"Empty NMF basis: {path}")
    if basis[0].shape != shape:
        zoom = (shape[0] / basis[0].shape[0], shape[1] / basis[0].shape[1])
        basis = [ndi.zoom(b, zoom, order=1) for b in basis]
        support = ndi.zoom(support.astype(np.uint8), zoom, order=0) != 0
        high = ndi.zoom(high.astype(np.uint8), zoom, order=0) != 0
    return np.stack(basis), support, high


def fit_wisp(image: np.ndarray, err: np.ndarray, excluded: np.ndarray,
             basis: np.ndarray, support: np.ndarray, high: np.ndarray,
             path: Path) -> tuple[np.ndarray, dict]:
    """Error-weighted non-negative fit of the NMF basis on bright-wisp pixels."""
    valid = ~excluded & np.isfinite(image) & np.isfinite(err) & (err > 0)
    valid[:4] = valid[-4:] = False
    valid[:, :4] = valid[:, -4:] = False
    off_wisp = valid & ~support
    if off_wisp.sum() < 1000:
        raise ValueError("Too few off-wisp pixels for the NMF baseline")
    baseline = float(np.median(image[off_wisp]))
    fit = valid & high
    if fit.sum() < 500:
        raise ValueError(f"Too few NMF fit pixels: {int(fit.sum())}")
    sigma = err[fit].astype(np.float64)
    matrix = np.stack([b[fit] / sigma for b in basis], axis=1)
    amplitudes, _ = nnls(matrix, (image[fit] - baseline) / sigma)
    model = np.einsum("k,kyx->yx", amplitudes.astype(np.float32), basis).astype(np.float32)
    return model, {"basis": str(path), "amplitudes": amplitudes.tolist(),
                   "fit_pixels": int(fit.sum()), "baseline": baseline}


# Output -----------------------------------------------------------------------

def write_fits(exp: Exposure, path: Path, science: np.ndarray, cards: dict,
               history: str) -> None:
    """Copy the input's KEEP_HDUS with SCI replaced and annotated."""
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".writing")
    # Unscaled reads copy the uint32 DQ (stored with BZERO) byte for byte.
    with fits.open(exp.path, memmap=True, lazy_load_hdus=True,
                   do_not_scale_image_data=True) as source:
        hdus, seen = [], set()
        for hdu in source:
            # Some reprocessed inputs repeat SCI after ASDF; keep the first.
            if hdu.name.upper() not in KEEP_HDUS or hdu.name.upper() in seen:
                continue
            seen.add(hdu.name.upper())
            if hdu.name.upper() == "SCI":
                header = hdu.header.copy()
                for key in ("DITHQC", "DITHN"):
                    if key in header:
                        del header[key]
                for key, value in cards.items():
                    header[key] = value
                header.add_history(history)
                hdu = fits.ImageHDU(np.asarray(science, np.float32), header, name="SCI")
            hdus.append(hdu)
        try:
            fits.HDUList(hdus).writeto(temp, overwrite=True, checksum=True,
                                       output_verify="silentfix")
            temp.replace(path)
        finally:
            temp.unlink(missing_ok=True)
