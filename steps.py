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
from scipy.spatial import cKDTree

# HDUs copied from each input into the output FITS.
KEEP_HDUS = ("PRIMARY", "SCI", "ERR", "DQ", "AREA", "VAR_POISSON",
             "VAR_RNOISE", "VAR_FLAT", "ASDF")

SEGMENT_NSIGMA = 2.8
MASK_GROWTH_CAP = 32       # pixels
WISP_SUPPORT_GROWTH = 8    # pixels excluded around the wisp template support
MAX_HELDOUT = 2            # dithers withheld from each 1/f scene for QC
MIN_DITHER_SCENE_PEERS = 1
SCAN_BLOCK_ROWS = 512
SCAN_TREND_BLOCK = 501
REPROJECT_ROWS = 128
REPROJECT_VALID_COVERAGE = 1.0 - 1e-6
SOURCE_MATCH_RADIUS_PIXELS = 5.0
SOURCE_ALIGNMENT_TOL_PIXELS = 0.1
MIN_SOURCE_ALIGNMENT_MATCHES = 3

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


# Registered-dither 1/f ------------------------------------------------------

def _source_centroids(mask: np.ndarray) -> np.ndarray:
    """Return x,y centroids of the source segments used for alignment checks."""
    labels, count = ndi.label(mask)
    if count == 0:
        return np.empty((0, 2), dtype=np.float64)
    areas = np.bincount(labels.ravel())
    ids = np.flatnonzero(areas[1:] >= 5) + 1
    if not ids.size:
        return np.empty((0, 2), dtype=np.float64)
    centers_yx = np.asarray(ndi.center_of_mass(mask, labels, ids), dtype=np.float64)
    return centers_yx[:, ::-1]


def _source_alignment(peer: Exposure, target: Exposure,
                      peer_sources: np.ndarray, target_sources: np.ndarray) -> dict:
    """Estimate a robust peer-to-target pixel shift from mutual source matches.

    Full-WCS coordinates provide the initial alignment. A median translation
    from mutual nearest source centroids refines it, like the relative shift
    correction used before PHANGS multi-tile destriping. The old 0.1-pixel
    residual is recorded as a diagnostic and never rejects the peer.
    """
    peer_xy = _source_centroids(peer_sources)
    target_xy = _source_centroids(target_sources)
    if not len(peer_xy) or not len(target_xy):
        return {"status": "WCS_ONLY", "matches": 0,
                "correction_x_px": 0.0, "correction_y_px": 0.0,
                "reason": "no source centroids in one or both exposures"}

    world = peer.wcs.all_pix2world(peer_xy, 0)
    mapped_xy = target.wcs.all_world2pix(world, 0)
    ny, nx = target.science.shape
    finite = (np.isfinite(mapped_xy).all(axis=1) &
              (mapped_xy[:, 0] >= 0) & (mapped_xy[:, 0] < nx) &
              (mapped_xy[:, 1] >= 0) & (mapped_xy[:, 1] < ny))
    mapped_xy = mapped_xy[finite]
    if not len(mapped_xy):
        return {"status": "WCS_ONLY", "matches": 0,
                "correction_x_px": 0.0, "correction_y_px": 0.0,
                "reason": "no peer source centroids project onto target grid"}

    target_tree = cKDTree(target_xy)
    peer_tree = cKDTree(mapped_xy)
    to_target, target_index = target_tree.query(mapped_xy, k=1)
    _, peer_index = peer_tree.query(target_xy, k=1)
    peers = np.arange(len(mapped_xy))
    mutual = ((peer_index[target_index] == peers) &
              (to_target <= SOURCE_MATCH_RADIUS_PIXELS))
    offsets = target_xy[target_index[mutual]] - mapped_xy[mutual]
    if len(offsets) < MIN_SOURCE_ALIGNMENT_MATCHES:
        return {"status": "WCS_ONLY", "matches": int(len(offsets)),
                "correction_x_px": 0.0, "correction_y_px": 0.0,
                "reason": f"fewer than {MIN_SOURCE_ALIGNMENT_MATCHES} mutual source matches"}

    keep = np.ones(len(offsets), dtype=bool)
    for _ in range(4):
        shift = np.median(offsets[keep], axis=0)
        scatter = np.linalg.norm(offsets - shift, axis=1)
        limit = max(1.0, 3.0 * robust_sigma(scatter[keep]))
        updated = scatter <= limit
        if updated.sum() < MIN_SOURCE_ALIGNMENT_MATCHES or np.array_equal(updated, keep):
            break
        keep = updated
    shift = np.median(offsets[keep], axis=0)
    residuals = np.linalg.norm(offsets[keep] - shift, axis=1)
    pre_residuals = np.linalg.norm(offsets[keep], axis=1)
    return {"status": "SHIFT_APPLIED" if np.any(np.abs(shift) > 0.01) else "WCS_OK",
            "matches": int(keep.sum()),
            "correction_x_px": float(shift[0]),
            "correction_y_px": float(shift[1]),
            "pre_correction_median_px": float(np.median(pre_residuals)),
            "post_correction_median_px": float(np.median(residuals)),
            "post_correction_p95_px": float(np.percentile(residuals, 95)),
            "post_correction_max_px": float(np.max(residuals)),
            "old_0p1px_gate_passed": bool(np.max(residuals) <= SOURCE_ALIGNMENT_TOL_PIXELS),
            "diagnostic_tolerance_px": SOURCE_ALIGNMENT_TOL_PIXELS}


def _aligned(peer: Exposure, target: Exposure, good: np.ndarray,
             peer_exclusion: np.ndarray, peer_sources: np.ndarray,
             target_sources: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict, dict]:
    """Full-WCS reproject peer SCI to target grid, conserving valid coverage.

    Bilinear SCI samples are normalized by their valid interpolation weight.
    A sample is usable only when all contributing source pixels are valid and
    unexcluded; partial and missing coverage remain NaN in the aligned image.
    """
    ny, nx = target.science.shape
    py, px = peer.science.shape
    alignment = _source_alignment(peer, target, peer_sources, target_sources)
    dx = alignment["correction_x_px"]
    dy = alignment["correction_y_px"]
    source_valid = (~peer.invalid & ~peer_exclusion &
                    np.isfinite(peer.science))
    source_values = np.where(source_valid, peer.science, 0).astype(np.float32)
    source_weights = source_valid.astype(np.float32)
    aligned = np.full((ny, nx), np.nan, dtype=np.float32)
    valid_coverage = np.zeros((ny, nx), dtype=np.float32)
    footprint = np.zeros((ny, nx), dtype=bool)

    for y0 in range(0, ny, REPROJECT_ROWS):
        y1 = min(y0 + REPROJECT_ROWS, ny)
        yy, xx = np.mgrid[y0:y1, 0:nx]
        # A positive source correction moves peer data toward +x,+y in the
        # target frame, so sample the inverse shift on the output grid.
        xx = xx - dx
        yy = yy - dy
        xy = np.column_stack((xx.ravel(), yy.ravel()))
        world = target.wcs.all_pix2world(xy, 0)
        peer_xy = peer.wcs.all_world2pix(world, 0)
        sx = peer_xy[:, 0].reshape(y1 - y0, nx)
        sy = peer_xy[:, 1].reshape(y1 - y0, nx)
        geometric = (np.isfinite(sx) & np.isfinite(sy) & (sx >= 0) &
                     (sx <= px - 1) & (sy >= 0) & (sy <= py - 1))
        safe_sx = np.where(geometric, sx, -1e6)
        safe_sy = np.where(geometric, sy, -1e6)
        coords = [safe_sy, safe_sx]
        values = ndi.map_coordinates(source_values, coords, order=1,
                                     mode="constant", cval=0, prefilter=False)
        weights = ndi.map_coordinates(source_weights, coords, order=1,
                                      mode="constant", cval=0, prefilter=False)
        usable = geometric & (weights >= REPROJECT_VALID_COVERAGE)
        block_coverage = valid_coverage[y0:y1]
        block_coverage[:] = np.where(geometric, weights, 0)
        footprint[y0:y1] = geometric
        aligned_block = aligned[y0:y1]
        aligned_block[usable] = values[usable] / weights[usable]

    sample = np.s_[::8, ::8]
    common = good[sample] & np.isfinite(aligned[sample])
    if common.sum() < 10000:
        raise ValueError(f"Too little dither overlap with {peer.path.name}: "
                         f"{int(common.sum())} stride-8 samples; requires 10000")
    background_offset = _sigma_clipped_median(
        aligned[sample][common] - target.science[sample][common])
    aligned -= np.float32(background_offset)
    stats = {"geometric_coverage": float(footprint.mean()),
             "full_valid_coverage": float(np.mean(valid_coverage >= REPROJECT_VALID_COVERAGE)),
             "partial_coverage": float(np.mean((valid_coverage > 0) &
                                                (valid_coverage < REPROJECT_VALID_COVERAGE))),
             "overlap_stride8_samples": int(common.sum()),
             "background_offset": background_offset}
    return aligned, valid_coverage, alignment, stats


def _sigma_clipped_median(values: np.ndarray, sigma: float = 3.0) -> float:
    values = np.asarray(values, dtype=np.float32)
    values = values[np.isfinite(values)]
    if not values.size:
        raise ValueError("No finite samples for sigma-clipped median")
    clipped = SigmaClip(sigma=sigma, maxiters=5)(values, masked=True, copy=False)
    retained = clipped.compressed()
    if not retained.size:
        raise ValueError("Sigma clipping removed all samples")
    return float(np.median(retained))


def fit_stripes(target: np.ndarray, scene: np.ndarray,
                good: np.ndarray) -> tuple[np.ndarray, dict]:
    """PHANGS multi-tile model: sigma-clipped row offsets per amplifier."""
    residual = target - scene
    valid = good & np.isfinite(residual)
    if valid.sum() < 1000:
        raise ValueError("Too few common sky pixels to fit dither 1/f")
    ny, nx = target.shape
    model = np.zeros_like(target, dtype=np.float32)
    row_coverage = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        for amp in range(4):
            lo = amp * nx // 4
            hi = (amp + 1) * nx // 4
            data = np.where(valid[:, lo:hi], residual[:, lo:hi], np.nan)
            clipped = SigmaClip(sigma=3.0, maxiters=5)(
                data, axis=1, masked=True, copy=False)
            profile = np.ma.median(clipped, axis=1).filled(np.nan)
            finite = np.isfinite(profile)
            row_coverage[f"A{amp + 1}"] = float(finite.mean())
            if not finite.any():
                raise ValueError(f"No measurable dither row offsets in amplifier {amp + 1}")
            profile = np.interp(np.arange(ny), np.flatnonzero(finite), profile[finite])
            model[:, lo:hi] = profile[:, None]
    # The fit cannot see a global offset: keep the input's masked sky median.
    model -= np.float32(np.median(model[good]))
    model += np.float32(np.median((target - model)[good]) - np.median(target[good]))
    model[~np.isfinite(target)] = 0
    return model, {"residual_sigma": robust_sigma(residual[valid][::8]),
                   "fit_fraction": float(valid.mean()),
                   "amplifier_row_coverage": row_coverage,
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


def dither_onef(target: Exposure, peers: list[Exposure], good: np.ndarray,
                peer_exclusions: list[np.ndarray], target_sources: np.ndarray,
                peer_sources: list[np.ndarray],
                diagnostics: bool = False) -> tuple[np.ndarray, dict]:
    """Fit amplifier-row offsets from a sigma-clipped median dither scene.

    The default scene uses every overlapping peer. Optional held-out checks
    are diagnostics only and never decide whether a computable model is used.
    Missing or excluded interpolation support stays NaN, never zero sky.
    """
    if target.path in {peer.path for peer in peers}:
        raise ValueError("Target exposure must not appear among its dither peers")
    if not (len(peers) == len(peer_exclusions) == len(peer_sources)):
        raise ValueError("Peer list and peer mask lists have different lengths")
    if len(peers) < MIN_DITHER_SCENE_PEERS:
        raise ValueError("Dither 1/f requires at least one other exposure")
    if diagnostics and len(peers) > 1:
        n_fit = max(1, len(peers) - MAX_HELDOUT)
        fit_peers, heldout = peers[:n_fit], peers[n_fit:]
        fit_exclusions, heldout_exclusions = peer_exclusions[:n_fit], peer_exclusions[n_fit:]
        fit_sources, heldout_sources = peer_sources[:n_fit], peer_sources[n_fit:]
    else:
        fit_peers, heldout = peers, []
        fit_exclusions, heldout_exclusions = peer_exclusions, []
        fit_sources, heldout_sources = peer_sources, []
    try:
        ny, nx = target.science.shape
        scene_stack = np.full((len(fit_peers), ny, nx), np.nan, dtype=np.float32)
        peer_diagnostics = []
        scene_names = []
        n_scene = 0
        for peer, exclusion, sources in zip(fit_peers, fit_exclusions, fit_sources):
            try:
                aligned, coverage, alignment, stats = _aligned(
                    peer, target, good, exclusion, sources, target_sources)
            except (ValueError, FloatingPointError) as exc:
                peer_diagnostics.append({"peer": peer.path.name, "role": "scene",
                                         "used": False, "reason": str(exc)})
                continue
            scene_stack[n_scene] = aligned
            n_scene += 1
            scene_names.append(peer.path.name)
            peer_diagnostics.append({"peer": peer.path.name, "role": "scene",
                                     "used": True, "coverage": stats,
                                     "source_alignment": alignment})
            del aligned, coverage
        if n_scene == 0:
            raise ValueError("No overlapping peers remain for the dither scene")
        scene_stack = scene_stack[:n_scene]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            clipped_scene = SigmaClip(sigma=3.0, maxiters=5)(
                scene_stack, axis=0, masked=True, copy=False)
            scene = np.ma.median(clipped_scene, axis=0).filled(np.nan).astype(np.float32)
        del clipped_scene, scene_stack
        model, meta = fit_stripes(target.science, scene, good)
        del scene
        corrected = target.science - model
        checks = []
        for peer, exclusion, sources in zip(heldout, heldout_exclusions, heldout_sources):
            try:
                aligned, coverage, alignment, stats = _aligned(
                    peer, target, good, exclusion, sources, target_sources)
            except (ValueError, FloatingPointError) as exc:
                peer_diagnostics.append({"peer": peer.path.name, "role": "heldout",
                                         "used": False, "reason": str(exc)})
                continue
            checks.append({"peer": peer.path.name,
                           "row_before": _row_rms(target.science - aligned, good),
                           "row_after": _row_rms(corrected - aligned, good),
                           "tile_before": _tile_rms(target.science - aligned, good),
                           "tile_after": _tile_rms(corrected - aligned, good),
                           "source_alignment": alignment,
                           "coverage": stats})
            peer_diagnostics.append({"peer": peer.path.name, "role": "heldout",
                                     "used": False, "coverage": stats,
                                     "source_alignment": alignment})
            del aligned, coverage
    except (ValueError, FloatingPointError) as exc:
        raise ValueError(f"Dither 1/f fit failed for {target.path.name}: {exc}") from exc
    qc_accepted = None
    mean = {}
    if checks:
        mean = {k: float(np.mean([c[k] for c in checks]))
                for k in ("row_before", "row_after", "tile_before", "tile_after")}
        qc_accepted = (mean["row_after"] < 0.95 * mean["row_before"] and
                       mean["tile_after"] < 0.95 * mean["tile_before"] and
                       all(c["row_after"] < 1.05 * c["row_before"] and
                           c["tile_after"] < 1.05 * c["tile_before"] for c in checks))
    return model, {"scene_dithers": scene_names, **meta,
                   "peer_reprojection": peer_diagnostics,
                   "qc_accepted": qc_accepted,
                   "diagnostics_enabled": diagnostics,
                   "diagnostic_only": True, **mean, "heldout": checks}


def _scan_1f_xy(image: np.ndarray, invalid: np.ndarray,
                block_rows: int = SCAN_BLOCK_ROWS,
                trend_block: int = SCAN_TREND_BLOCK,
                clip_hi: float = 0.90,
                sky_quantile: float = 0.40) -> tuple[np.ndarray, dict]:
    """Vectorized ProFound-style scan model in detector XY orientation.

    Estimate column profiles in row blocks and row profiles across the full
    width. Keep their centered 40th-percentile structure after removing a
    501-pixel running-median trend, matching the Jumprope keep_trend=True
    configuration. No segmentation/source mask is applied; the upper 10% of
    valid values are still clipped regardless of whether they belong to
    sources. Masked science pixels are not changed by this model.
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


def _surface_design(x: np.ndarray, y: np.ndarray, degree: int) -> np.ndarray:
    """Low-order 2D polynomial basis, ordered by total degree."""
    return np.column_stack([x ** i * y ** (total - i)
                            for total in range(degree + 1)
                            for i in range(total + 1)])


def _surface_evaluate(coef: np.ndarray, x: np.ndarray, y: np.ndarray,
                      degree: int) -> np.ndarray:
    """Evaluate a polynomial surface without allocating a full pixel design matrix."""
    result = np.zeros((len(y), len(x)), dtype=np.float32)
    index = 0
    for total in range(degree + 1):
        for i in range(total + 1):
            j = total - i
            result += (np.float32(coef[index]) *
                       (y[:, None] ** j) * (x[None, :] ** i))
            index += 1
    return result


def _robust_surface_fit(x: np.ndarray, y: np.ndarray, z: np.ndarray,
                        degree: int = 3) -> tuple[np.ndarray, dict]:
    """Fit a clipped polynomial to robust tile medians, not individual pixels."""
    design = _surface_design(x, y, degree)
    keep = np.isfinite(z) & np.all(np.isfinite(design), axis=1)
    min_points = design.shape[1] + 5
    if keep.sum() < min_points:
        raise ValueError(f"Too few tiles for degree-{degree} 2D surface fit")
    for _ in range(5):
        coef = np.linalg.lstsq(design[keep], z[keep], rcond=None)[0]
        residual = z - design @ coef
        center = float(np.median(residual[keep]))
        sigma = robust_sigma(residual[keep])
        if not np.isfinite(sigma) or sigma <= 1e-12:
            break
        updated = np.isfinite(residual) & (np.abs(residual - center) <= 3.5 * sigma)
        if updated.sum() < min_points or np.array_equal(updated, keep):
            break
        keep = updated
    coef = np.linalg.lstsq(design[keep], z[keep], rcond=None)[0]
    return coef, {"tiles": int(z.size), "used_tiles": int(keep.sum()),
                  "degree": degree,
                  "fit_rms": float(np.sqrt(np.mean((z[keep] - design[keep] @ coef) ** 2)))}


def _tile_medians(image: np.ndarray, good: np.ndarray, tile: int
                  ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return robust low-k samples and their detector coordinates."""
    ny, nx = image.shape
    xs, ys, values, amplifiers = [], [], [], []
    width = nx // 4
    for y0 in range(0, ny, tile):
        y1 = min(ny, y0 + tile)
        for x0 in range(0, nx, tile):
            x1 = min(nx, x0 + tile)
            amp = min(3, x0 // width)
            valid = good[y0:y1, x0:x1] & np.isfinite(image[y0:y1, x0:x1])
            if valid.mean() < 0.25:
                continue
            values.append(float(np.median(image[y0:y1, x0:x1][valid])))
            xs.append((x0 + x1 - 1) / 2.0)
            ys.append((y0 + y1 - 1) / 2.0)
            amplifiers.append(amp)
    return (np.asarray(xs), np.asarray(ys), np.asarray(values),
            np.asarray(amplifiers, dtype=np.int8))


def fit_amp_local_lowk(image: np.ndarray, invalid: np.ndarray, *,
                       tile: int = 128, global_degree: int = 3,
                       local_degree: int = 2) -> tuple[np.ndarray, dict]:
    """Fit only 2D low-k structure localized to one of four amplifiers.

    A single full-detector polynomial is the shared, continuous low-k scene.
    Each amplifier also gets a local 2D surface. The subtracted term is the
    local surface minus the shared surface, so structure that continues
    smoothly across the detector is retained while amplifier-confined
    curvature or offsets become a correction. Tile medians suppress compact
    sources; only detector-invalid pixels are excluded.
    """
    image = np.asarray(image, dtype=np.float32)
    invalid = np.asarray(invalid, dtype=bool)
    good = ~invalid & np.isfinite(image) & (image != 0)
    ny, nx = image.shape
    if nx % 4:
        raise ValueError("Amplifier-local fit requires four equal-width regions")
    if good.sum() < 10000:
        raise ValueError("Too few valid pixels for amplifier-local low-k fit")
    xs, ys, z, amps = _tile_medians(image, good, tile)
    xg = 2.0 * xs / max(nx - 1, 1) - 1.0
    yg = 2.0 * ys / max(ny - 1, 1) - 1.0
    shared_coef, shared_stats = _robust_surface_fit(xg, yg, z, global_degree)
    yy = 2.0 * np.arange(ny, dtype=np.float32) / max(ny - 1, 1) - 1.0
    xx_global = 2.0 * np.arange(nx, dtype=np.float32) / max(nx - 1, 1) - 1.0
    shared = _surface_evaluate(shared_coef, xx_global, yy, global_degree)
    local = np.zeros_like(image, dtype=np.float32)
    local_stats = {}
    amp_width = nx // 4
    for amp in range(4):
        select = amps == amp
        if select.sum() < (local_degree + 1) * (local_degree + 2) // 2 + 5:
            raise ValueError(f"Too few low-k tiles in amplifier {amp + 1}")
        xlocal = 2.0 * (xs[select] - amp * amp_width) / max(amp_width - 1, 1) - 1.0
        local_coef, fit_stats = _robust_surface_fit(
            xlocal, yg[select], z[select], local_degree)
        x0, x1 = amp * amp_width, (amp + 1) * amp_width
        xl = (2.0 * np.arange(amp_width, dtype=np.float32) /
              max(amp_width - 1, 1) - 1.0)
        local_surface = _surface_evaluate(local_coef, xl, yy, local_degree)
        local[:, x0:x1] = local_surface - shared[:, x0:x1]
        local_stats[f"A{amp + 1}"] = {**fit_stats,
            "local_model_rms": float(np.sqrt(np.mean(
                local[:, x0:x1][good[:, x0:x1]].astype(np.float64) ** 2)))}

    # The low-k component is differential by construction; remove its global
    # median so this step cannot change the exposure's absolute sky pedestal.
    local[~good] = 0.0
    local -= np.float32(np.median(local[good]))
    median_before = float(np.median(image[good]))
    median_after = float(np.median((image - local)[good]))
    local += np.float32(median_after - median_before)
    local[~good] = 0.0
    return local, {"algorithm": "AMP_LOCAL_2D_MINUS_SHARED_SURFACE",
                   "tile_px": tile, "global_degree": global_degree,
                   "local_degree": local_degree,
                   "shared_fit": shared_stats,
                   "amplifiers": local_stats,
                   "shared_surface_rms": float(np.sqrt(np.mean(
                       shared[good].astype(np.float64) ** 2))),
                   "local_model_rms": float(np.sqrt(np.mean(
                       local[good].astype(np.float64) ** 2))),
                   "fit_fraction": float(good.mean()),
                   "global_median_before": median_before,
                   "global_median_after": float(np.median((image - local)[good]))}


def scan_onef_joint(image: np.ndarray, invalid: np.ndarray
                    ) -> tuple[np.ndarray, dict]:
    """Jointly model keep-trend 1/f and amplifier-confined 2D low-k residuals.

    The components are estimated on the same frame and returned as one model
    for one subtraction. The shared full-detector low-k surface is a protected
    reference; only the four amplifier-local residual surfaces are subtracted.
    """
    data = np.asarray(image, dtype=np.float32)
    scan_model, scan_stats = scan_onef(data, invalid)
    amp_model, amp_stats = fit_amp_local_lowk(data - scan_model, invalid)
    model = scan_model + amp_model
    bad = np.asarray(invalid, dtype=bool) | ~np.isfinite(data) | (data == 0)
    model[bad] = 0.0
    return model, {"algorithm": "PROFOUND_SCAN_KEEP_TREND_PLUS_AMP_LOCAL_2D",
                   "keep_trend": True,
                   "scan_model_rms": scan_stats["model_rms"],
                   "scan": scan_stats, "amp_lowk": amp_stats,
                   "model_rms": float(np.sqrt(np.mean(
                       model[~bad].astype(np.float64) ** 2))),
                   "fit_fraction": min(scan_stats["fit_fraction"],
                                       amp_stats["fit_fraction"])}


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
