#!/usr/bin/env python3
"""Batch Stage 2 plus: dither 1/f -> masked median sky -> NMF wisp.

One worker owns an eight-dither detector group. It reads each calibrated input
once, keeps its science, error, DQ, WCS, and masks in memory, and writes one
sky-subtracted FITS per exposure. Long-wave groups first create the Photutils
object masks used by their paired short-wave groups. No single-image 1/f,
four-region sky, polynomial sky, or horizontal sky fit is applied.
"""

from __future__ import annotations

import argparse
import concurrent.futures as futures
import json
import os
import re
import shutil
import sys
import traceback
import warnings
from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

import numpy as np
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi
from scipy.optimize import nnls

from onef_masks import double_source_mask
from sky_preserving_onef import fit_stripes, row_rms, shift_to_target, tile_rms
from skysurfir_clean import KEEP_HDU_NAMES, base_invalid, source_segmentation

WISP_DETECTORS = {"NRCA3", "NRCA4", "NRCB3", "NRCB4"}
DETECTOR_PATTERN = re.compile(r"^NRC[AB](?:[1-4]|LONG)$")


@dataclass
class Frame:
    path: Path
    image: np.ndarray
    err: np.ndarray
    dq: np.ndarray
    wcs: WCS
    sci_header: fits.Header
    detector: str
    band: str
    invalid: np.ndarray
    aligned_source: np.ndarray


def read_frame(path: Path) -> Frame:
    with fits.open(path, memmap=False, lazy_load_hdus=True) as hdul:
        image = np.array(hdul["SCI_ORIG" if "SCI_ORIG" in hdul else "SCI"].data,
                         dtype=np.float32, copy=True)
        err = np.array(hdul["ERR"].data, dtype=np.float32, copy=True)
        dq = np.array(hdul["DQ"].data, dtype=np.uint32, copy=True)
        sci_header = hdul["SCI"].header.copy()
        detector = str(hdul[0].header["DETECTOR"]).upper()
        band = str(hdul[0].header["FILTER"]).upper()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(sci_header)
    invalid = base_invalid(dq, image) | ~np.isfinite(err) | (err <= 0)
    aligned_source = image.copy()
    aligned_source[invalid] = np.nan
    return Frame(path, image, err, dq, wcs, sci_header, detector, band,
                 invalid, aligned_source)


def basis_path(root: Path, detector: str, band: str) -> Path | None:
    det, filt = detector.lower(), band.lower()
    for directory in (root / f"{det}_org", root / det):
        for extension in (".fits", ".fits.gz"):
            path = directory / f"{det}_{filt}_wisp{extension}"
            if path.is_file():
                return path
    return None


@lru_cache(maxsize=1)
def load_basis(root: str, detector: str, band: str):
    path = basis_path(Path(root), detector, band)
    if path is None:
        raise FileNotFoundError(f"Missing NMF wisp basis for {detector}/{band}")
    with fits.open(path, memmap=False) as hdul:
        templates = [np.array(hdul[f"WISP_{i}"].data, dtype=np.float32, copy=True)
                     for i in range(10) if f"WISP_{i}" in hdul]
        support = np.asarray(hdul["MASK"].data, bool)
        high = np.asarray(hdul["MASK_hSNR"].data, bool)
    if not templates:
        raise ValueError(f"Empty NMF basis: {path}")
    shape = (2048, 2048)
    if templates[0].shape != shape:
        zoom = (shape[0] / templates[0].shape[0],
                shape[1] / templates[0].shape[1])
        templates = [ndi.zoom(x, zoom, order=1) for x in templates]
        support = ndi.zoom(support.astype(np.uint8), zoom, order=0) != 0
        high = ndi.zoom(high.astype(np.uint8), zoom, order=0) != 0
    return np.stack(templates), support, high, str(path)


def fit_wisp(image: np.ndarray, err: np.ndarray, invalid: np.ndarray,
             sources: np.ndarray, templates: tuple) -> tuple[np.ndarray, dict]:
    basis, support, high, path = templates
    valid = (~invalid & ~sources & np.isfinite(image) & np.isfinite(err) &
             (err > 0))
    valid[:4] = False
    valid[-4:] = False
    valid[:, :4] = False
    valid[:, -4:] = False
    sky_pixels = valid & ~support
    if sky_pixels.sum() < 1000:
        raise ValueError("Too few off-wisp pixels for the NMF baseline")
    baseline = float(np.median(image[sky_pixels]))
    fit = valid & high
    nfit = int(fit.sum())
    if nfit < 500:
        raise ValueError(f"Too few NMF fit pixels: {nfit}")
    sigma = np.asarray(err[fit], dtype=np.float64)
    matrix = np.empty((nfit, basis.shape[0]), dtype=np.float64)
    for index, template in enumerate(basis):
        matrix[:, index] = np.asarray(template[fit], dtype=np.float64) / sigma
    target = (np.asarray(image[fit], dtype=np.float64) - baseline) / sigma
    amplitudes, _ = nnls(matrix, target)
    model = np.einsum("k,kyx->yx", amplitudes.astype(np.float32), basis,
                      optimize=True).astype(np.float32)
    return model, {"basis": path, "amplitudes": amplitudes.tolist(),
                   "fit_pixels": nfit, "baseline": baseline}


def matched_longwave_name(name: str) -> str:
    match = re.search(r"_nrc([ab])[1-4]_cal\.fits$", name, flags=re.IGNORECASE)
    if not match:
        raise ValueError(f"Not a short-wave cal: {name}")
    return name[:match.start()] + f"_nrc{match.group(1).lower()}long_cal.fits"


def longwave_mask(frame: Frame, output_dir: Path) -> np.ndarray:
    path = output_dir / "_lw_masks" / frame.path.name
    mask, _, _, _ = source_segmentation(frame.image, frame.err, frame.invalid, 2.8)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".writing")
    fits.HDUList([fits.PrimaryHDU(),
                  fits.ImageHDU(mask.astype(np.uint8), header=frame.sci_header,
                                name="OBJMASK")]).writeto(temp, overwrite=True)
    os.replace(temp, path)
    return mask


def projected_longwave_mask(frame: Frame, output_dir: Path) -> np.ndarray:
    path = output_dir / "_lw_masks" / matched_longwave_name(frame.path.name)
    with fits.open(path, memmap=False) as hdul:
        reference = np.asarray(hdul["OBJMASK"].data, dtype=np.uint8)
        ref_header = hdul["OBJMASK"].header.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        source_wcs = WCS(ref_header).celestial
        target_wcs = frame.wcs.celestial
    shape = frame.image.shape
    projected = np.zeros(shape, dtype=bool)
    for y0 in range(0, shape[0], 128):
        y1 = min(y0 + 128, shape[0])
        yy, xx = np.mgrid[y0:y1, :shape[1]]
        world = target_wcs.pixel_to_world_values(xx, yy)
        source_x, source_y = source_wcs.world_to_pixel_values(*world)
        valid = (np.isfinite(source_x) & np.isfinite(source_y) &
                 (source_x >= 0) & (source_y >= 0) &
                 (source_x < shape[1]) & (source_y < shape[0]))
        sampled = ndi.map_coordinates(reference,
                                      [np.nan_to_num(source_y, nan=-1),
                                       np.nan_to_num(source_x, nan=-1)],
                                      order=0, mode="constant", cval=0,
                                      prefilter=False)
        projected[y0:y1] = valid & (sampled != 0)
    return projected


def aligned_peer(peer: Frame, target: Frame, good: np.ndarray) -> np.ndarray:
    shift = shift_to_target(target.wcs, peer.wcs, target.image.shape)
    aligned = ndi.shift(peer.aligned_source, shift, order=1, mode="constant",
                        cval=np.nan, prefilter=False, output=np.float32)
    sample = np.s_[::8, ::8]
    valid = good[sample] & np.isfinite(aligned[sample])
    if valid.sum() < 10000:
        raise ValueError(f"Too little dither overlap: {peer.path.name}")
    level = float(np.median((aligned[sample] - target.image[sample])[valid]))
    aligned -= np.float32(level)  # Peers only; target sky is retained.
    return aligned


def dither_onef(target: Frame, peers: list[Frame], good: np.ndarray):
    fit_peers, heldout = peers[:5], peers[5:]
    try:
        aligned = [aligned_peer(peer, target, good) for peer in fit_peers]
        stack = np.stack(aligned)
        del aligned
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            template = np.nanmedian(stack, axis=0, overwrite_input=True).astype(np.float32)
        del stack
        corrected, model, fit_meta = fit_stripes(target.image, template, good)
        del template
        validation = []
        for peer in heldout:
            aligned = aligned_peer(peer, target, good)
            validation.append({"peer": peer.path.name,
                               "row_before": row_rms(target.image - aligned, good),
                               "row_after": row_rms(corrected - aligned, good),
                               "tile_before": tile_rms(target.image - aligned, good),
                               "tile_after": tile_rms(corrected - aligned, good)})
            del aligned
    except (ValueError, FloatingPointError) as exc:
        raise ValueError(f"Dither 1/f fit failed for {target.path.name}: {exc}") from exc
    row_before = float(np.mean([x["row_before"] for x in validation]))
    row_after = float(np.mean([x["row_after"] for x in validation]))
    tile_before = float(np.mean([x["tile_before"] for x in validation]))
    tile_after = float(np.mean([x["tile_after"] for x in validation]))
    passed = bool(validation) and (
        row_after < 0.95 * row_before and tile_after < 0.95 * tile_before and
        all(x["row_after"] < 1.05 * x["row_before"] and
            x["tile_after"] < 1.05 * x["tile_before"] for x in validation))
    scores = {"accepted": passed, "applied": True, "heldout": validation,
              "qc_rule": "mean row and tile RMS improve >5%; neither peer worsens >5%",
              "row_before": row_before, "row_after": row_after,
              "tile_before": tile_before, "tile_after": tile_after,
              **fit_meta}
    return corrected, model, scores


def masked_sky(image: np.ndarray, invalid: np.ndarray,
               sources: np.ndarray, support: np.ndarray | None) -> tuple[np.ndarray, float, np.ndarray]:
    mask = invalid | sources
    if support is not None:
        mask |= ndi.binary_dilation(support, iterations=8)
    good = ~mask & np.isfinite(image)
    if good.sum() < 10000:
        raise ValueError("Too few pixels for the masked sky median")
    pedestal = float(np.median(image[good]))
    return image - np.float32(pedestal), pedestal, good


def wisp_safe_onef_mask(base_mask: np.ndarray, high_wisp: np.ndarray,
                        min_row_pixels: int = 64) -> tuple[np.ndarray, int]:
    """Keep the bright-wisp mask unless it empties an amplifier row."""
    mask = base_mask | high_wisp
    rescued = 0
    width = mask.shape[1] // 4
    for lo in range(0, mask.shape[1], width):
        hi = lo + width
        sparse = np.count_nonzero(~mask[:, lo:hi], axis=1) < min_row_pixels
        rescued += int(sparse.sum())
        mask[sparse, lo:hi] = base_mask[sparse, lo:hi]
    return mask, rescued


def output_fits(frame: Frame, science: np.ndarray, output: Path, record: dict):
    output.parent.mkdir(parents=True, exist_ok=True)
    temp = output.with_name(output.name + ".writing")
    with fits.open(frame.path, memmap=True, do_not_scale_image_data=True,
                   lazy_load_hdus=True) as source:
        keep = []
        original_header = None
        for hdu in source:
            if hdu.name.upper() not in KEEP_HDU_NAMES:
                continue
            if hdu.name.upper() == "SCI":
                original_header = hdu.header.copy()
                continue
            keep.append(hdu)
        if original_header is None:
            raise ValueError(f"Missing SCI: {frame.path}")
        original_header["STGORDER"] = ("DITH1F-MED-WISP-ZERO", "subtraction order")
        original_header["ONEFMETH"] = ("DITHER", "registered dither 1/f method")
        original_header["DITH1F"] = (record["onef"]["applied"],
                                     "dither 1/f correction applied")
        original_header["DITHQC"] = (record["onef"]["accepted"],
                                     "held-out dither checks passed")
        original_header["DITHN"] = (5, "dithers in scene template")
        original_header["NMFSKY"] = (True, "masked median sky removed")
        original_header["SKYMED"] = (record["sky_median"], "initial sky MJy/sr")
        original_header["POSTMED"] = (record["post_wisp_median"],
                                      "post-wisp scalar MJy/sr")
        original_header["WISPMETH"] = ("NMF" if record["wisp"] else "NONE",
                                       "fixed-basis NMF wisp")
        original_header.add_history("SCI_ORIG: dither 1/f, masked median sky, NMF wisp, scalar zero")
        keep.insert(1, fits.ImageHDU(np.asarray(science, np.float32),
                                     header=original_header, name="SCI"))
        hdus = fits.HDUList(keep)
        try:
            hdus.writeto(temp, overwrite=True, checksum=True,
                         output_verify="silentfix")
            os.replace(temp, output)
        finally:
            temp.unlink(missing_ok=True)


def thumbnail(data: np.ndarray, sky: float, final: np.ndarray,
              good: np.ndarray, size: int = 176) -> np.ndarray:
    sample = final[good][::40]
    median = float(np.median(sample))
    sigma = 1.4826 * float(np.median(np.abs(sample - median)))
    sigma = max(sigma, 0.001)
    high = max(float(np.nanpercentile(final[::8, ::8], 99.5)), 10 * sigma)
    low = -3 * sigma
    shifted = data - np.float32(sky)
    reduced = shifted[::8, ::8]
    scaled = np.arcsinh(np.clip((reduced - low) / (1.5 * sigma), 0, None))
    scaled /= np.arcsinh((high - low) / (1.5 * sigma))
    gray = np.uint8(255 * np.nan_to_num(np.clip(scaled, 0, 1), nan=0))
    return np.asarray(Image.fromarray(gray).resize((size, size), Image.Resampling.BILINEAR))


def process_group(paths: list[str], output_dir: str, template_dir: str,
                  mask_dir: str | None = None,
                  target_names: set[str] | None = None) -> dict:
    out = Path(output_dir)
    frames = [read_frame(Path(path)) for path in sorted(paths)]
    detector, band = frames[0].detector, frames[0].band
    if len(frames) != 8 or any((f.detector, f.band) != (detector, band) for f in frames):
        raise ValueError("Each task must contain eight same-detector/filter dithers")
    templates = (load_basis(template_dir, detector, band)
                 if detector in WISP_DETECTORS else None)
    sources = {}
    for frame in frames:
        if target_names is not None and frame.path.name not in target_names:
            continue
        if detector.endswith("LONG"):
            source = longwave_mask(frame, out)
        else:
            local, _, _, _ = source_segmentation(frame.image, frame.err,
                                                  frame.invalid, 2.8)
            reference_dir = Path(mask_dir) if mask_dir else out
            source = local | projected_longwave_mask(frame, reference_dir)
        sources[frame.path.name] = source
    records = []
    thumbs = []
    for target in frames:
        if target_names is not None and target.path.name not in target_names:
            continue
        source = sources.pop(target.path.name)
        grown, dilation = double_source_mask(source)
        base_mask = target.invalid | grown
        if templates is not None:
            onef_mask, rescued_rows = wisp_safe_onef_mask(base_mask, templates[2])
        else:
            onef_mask, rescued_rows = base_mask, 0
        good = ~onef_mask
        if good.mean() < 0.5:
            raise ValueError(f"Too little clean area for 1/f: {target.path.name}")
        peers = [peer for peer in frames if peer.path != target.path]
        after_onef, model, onef_meta = dither_onef(target, peers, good)
        del model
        skyless, sky_pedestal, sky_good = masked_sky(
            after_onef, target.invalid, grown,
            templates[1] if templates is not None else None)
        if templates is not None:
            wisp, wisp_meta = fit_wisp(skyless, target.err, target.invalid,
                                       grown, templates)
            after_wisp = skyless - wisp
            del wisp
        else:
            wisp_meta = None
            after_wisp = skyless
        post_pedestal = float(np.median(after_wisp[sky_good]))
        final = after_wisp - np.float32(post_pedestal)
        record = {"input": str(target.path), "file": target.path.name,
                  "detector": detector, "filter": band,
                  "onef": onef_meta, "source_fraction": float(source.mean()),
                  "mask_fraction": float(onef_mask.mean()),
                  "wisp_mask_rows_rescued": rescued_rows,
                  "mask_dilation": dilation, "sky_median": sky_pedestal,
                  "post_wisp_median": post_pedestal,
                  "final_masked_median": float(np.median(final[sky_good])),
                  "wisp": wisp_meta}
        output_fits(target, final, out / "cal_sky_renorm" / target.path.name, record)
        thumbs.append((target.path.name,
                       thumbnail(target.image, sky_pedestal + post_pedestal,
                                 final, sky_good),
                       thumbnail(final, 0.0, final, sky_good)))
        records.append(record)
        del source, grown, onef_mask, good, skyless, sky_good, after_wisp, final
    name_parts = frames[0].path.name.split("_")
    return {"key": (name_parts[0], name_parts[1], detector, band),
            "records": records, "thumbs": thumbs}


def font(size: int):
    for path in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",):
        try:
            return ImageFont.truetype(path, size)
        except OSError:
            pass
    return ImageFont.load_default()


def render_contacts(results: list[dict], output_dir: Path) -> list[Path]:
    grouped = defaultdict(dict)
    for result in results:
        observation, group, detector, _ = result["key"]
        grouped[(observation, group, detector[3])][detector] = result["thumbs"]
    paths = []
    tile, gap, label_width, top = 176, 8, 88, 86
    for (observation, group, module), detectors in sorted(grouped.items()):
        names = [name for name in (f"NRC{module}{part}"
                 for part in ("1", "2", "3", "4", "LONG")) if name in detectors]
        width = label_width + 2 * len(names) * (tile + gap)
        height = top + 8 * (tile + gap)
        sheet = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(sheet)
        draw.text((12, 8), f"{observation} group {group} | module {module} | "
                  f"{8 * len(names)} input / final pairs",
                  fill="black", font=font(24))
        draw.text((12, 42), "Input sky pedestal removed for display only; each pair shares a grayscale stretch.",
                  fill="black", font=font(15))
        for j, detector in enumerate(names):
            draw.text((label_width + (2*j)*(tile+gap), 67), detector + " input",
                      fill="black", font=font(13))
            draw.text((label_width + (2*j+1)*(tile+gap), 67), "final skyless",
                      fill="black", font=font(13))
            items = sorted(detectors[detector], key=lambda x: x[0].split("_")[2])
            for row, (name, raw, final) in enumerate(items):
                y = top + row * (tile + gap)
                for col, array in enumerate((raw, final)):
                    sheet.paste(Image.fromarray(array).convert("RGB"),
                                (label_width + (2*j+col)*(tile+gap), y))
        for row in range(8):
            draw.text((10, top + row*(tile+gap) + 78), f"{row+1:05d}",
                      fill="black", font=font(14))
        path = (output_dir / "figures" /
                f"{observation}_{group}_module_{module.lower()}_input_final.png")
        path.parent.mkdir(parents=True, exist_ok=True)
        sheet.save(path, compress_level=4)
        paths.append(path)
    return paths


def patch_contacts(results: list[dict], output_dir: Path) -> None:
    """Replace only repaired detector columns in existing overview PNGs."""
    tile, gap, label_width, top = 176, 8, 88, 86
    for result in results:
        observation, group, detector, _ = result["key"]
        module = detector[3]
        parts = ("1", "2", "3", "4", "LONG")
        column = parts.index(detector[4:])
        path = (output_dir / "figures" /
                f"{observation}_{group}_module_{module.lower()}_input_final.png")
        if not path.exists():
            path = (output_dir / "figures" /
                    f"{group}_module_{module.lower()}_input_final.png")
        sheet = Image.open(path).convert("RGB")
        items = sorted(result["thumbs"], key=lambda x: x[0].split("_")[2])
        for name, raw, final in items:
            row = int(name.split("_")[2]) - 1
            y = top + row * (tile + gap)
            for offset, array in enumerate((raw, final)):
                sheet.paste(Image.fromarray(array).convert("RGB"),
                            (label_width + (2 * column + offset) * (tile + gap), y))
        temp = path.with_name(path.name + ".writing")
        sheet.save(temp, format="PNG", compress_level=4)
        os.replace(temp, path)


def run_phase(tasks: list[list[str]], output_dir: Path, template_dir: Path,
              workers: int, label: str,
              selected: dict[tuple[str, str, str], set[str]] | None = None
              ) -> tuple[list[dict], list[dict]]:
    results, errors = [], []
    if not tasks:
        return results, errors
    with futures.ProcessPoolExecutor(max_workers=min(workers, len(tasks))) as pool:
        pending = {}
        for task in tasks:
            name = Path(task[0]).name
            header = fits.getheader(task[0], 0)
            key = (name.split("_")[0], name.split("_")[1],
                   str(header["DETECTOR"]).upper(),
                   str(header["FILTER"]).upper())
            target_names = selected[key] if selected is not None else None
            pending[pool.submit(process_group, task, str(output_dir),
                                str(template_dir), None, target_names)] = task
        for job in futures.as_completed(pending):
            task = pending[job]
            try:
                result = job.result()
                results.append(result)
                print(f"{label} {len(results)}/{len(tasks)} {result['key']} "
                      f"dither QC passed {sum(x['onef']['accepted'] for x in result['records'])}"
                      f"/{len(result['records'])}",
                      flush=True)
            except Exception as exc:
                errors.append({"group": Path(task[0]).name.split("_")[1],
                               "detector": fits.getheader(task[0], 0).get("DETECTOR"),
                               "error": repr(exc), "traceback": traceback.format_exc()})
                print(f"{label} FAILED {Path(task[0]).name}: {exc}",
                      file=sys.stderr, flush=True)
    return results, errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-list", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--template-dir", type=Path,
                        default=Path(__file__).parent / "templates")
    parser.add_argument("--workers", type=int, default=32)
    parser.add_argument("--plan-only", action="store_true")
    parser.add_argument("--repair-mask-failures", action="store_true",
                        help="rerun only groups whose dither fit lacked clean amplifier rows")
    parser.add_argument("--apply-rejected-onef", action="store_true",
                        help="rerun only prior dither QC failures, applying their successful 1/f fits")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")
    if args.repair_mask_failures and args.apply_rejected_onef:
        parser.error("choose only one targeted rerun mode")
    paths = [str(Path(line.strip()).resolve()) for line in args.input_list.read_text().splitlines()
             if line.strip()]
    if not paths or len(paths) != len(set(paths)) or any(not Path(p).is_file() for p in paths):
        parser.error("input list is empty, duplicated, or has missing files")
    output_dir = args.output_dir.resolve()
    if output_dir in {Path(p).parent for p in paths}:
        parser.error("output directory must differ from the raw input directory")
    groups = defaultdict(list)
    for path in paths:
        header = fits.getheader(path, 0)
        detector = str(header.get("DETECTOR", "")).upper()
        band = str(header.get("FILTER", "")).upper()
        if not DETECTOR_PATTERN.fullmatch(detector) or not band:
            parser.error(f"Unsupported detector/filter in {path}")
        name_parts = Path(path).name.split("_")
        groups[(name_parts[0], name_parts[1], detector, band)].append(path)
    if any(len(group) != 8 for group in groups.values()):
        parser.error("each visit/detector/filter must have exactly eight dithers")
    if len(paths) != len({Path(p).name for p in paths}):
        parser.error("duplicate basenames would collide in output")
    missing_bases = sorted((obs, group, det, band)
                           for obs, group, det, band in groups
                           if det in WISP_DETECTORS and
                           basis_path(args.template_dir, det, band) is None)
    if missing_bases:
        parser.error(f"Missing wisp bases: {missing_bases}")
    lw = [sorted(files) for (obs, group, detector, band), files in sorted(groups.items())
          if detector.endswith("LONG")]
    sw = [sorted(files) for (obs, group, detector, band), files in sorted(groups.items())
          if not detector.endswith("LONG")]
    if args.plan_only:
        example = Path(paths[0])
        with fits.open(example, memmap=True, lazy_load_hdus=True) as hdul:
            bytes_per = sum(hdu._data_size + hdu._data_offset - hdu._header_offset
                            for hdu in hdul if hdu.name.upper() in KEEP_HDU_NAMES)
        print(json.dumps({"exposures": len(paths), "groups": len(groups),
                          "longwave_groups": len(lw), "shortwave_groups": len(sw),
                          "workers": args.workers,
                          "estimated_fits_bytes": int(bytes_per * len(paths)),
                          "free_bytes": shutil.disk_usage(output_dir.parent).free}, indent=2))
        return 0
    if args.apply_rejected_onef:
        summary_path = output_dir / "summary.json"
        summary = json.loads(summary_path.read_text())
        selected = defaultdict(set)
        for record in summary["results"]:
            if not record["onef"]["accepted"]:
                name_parts = record["file"].split("_")
                key = (name_parts[0], name_parts[1], record["detector"],
                       record["filter"])
                selected[key].add(record["file"])
        if not selected:
            print("No dither QC failures to rerun", flush=True)
            return 0
        long_tasks = [sorted(groups[key]) for key in sorted(selected)
                      if key[2].endswith("LONG")]
        short_tasks = [sorted(groups[key]) for key in sorted(selected)
                       if not key[2].endswith("LONG")]
        print(f"Applying dither 1/f to {sum(map(len, selected.values()))} frames "
              f"in {len(selected)} groups", flush=True)
        long_results, long_errors = run_phase(long_tasks, output_dir,
                                              args.template_dir, args.workers,
                                              "LW-RERUN", selected)
        if long_errors:
            print(json.dumps(long_errors, indent=2), file=sys.stderr)
            return 1
        short_results, short_errors = run_phase(short_tasks, output_dir,
                                                args.template_dir, args.workers,
                                                "SW-RERUN", selected)
        if short_errors:
            print(json.dumps(short_errors, indent=2), file=sys.stderr)
            return 1
        replacements = {r["file"]: r for result in long_results + short_results
                        for r in result["records"]}
        if set(replacements) != set().union(*selected.values()):
            raise RuntimeError("Targeted rerun returned an incomplete set of frames")
        patch_contacts(long_results + short_results, output_dir)
        summary["results"] = [replacements.get(r["file"], r)
                              for r in summary["results"]]
        for record in summary["results"]:
            record["onef"].setdefault("applied", record["onef"]["accepted"])
        summary["method"] = "always apply dither 1/f -> masked median sky -> NMF wisp -> scalar zero"
        summary["onef_applied"] = sum(r["onef"]["applied"] for r in summary["results"])
        summary["onef_accepted"] = sum(r["onef"]["accepted"] for r in summary["results"])
        summary["onef_rejected"] = summary["completed"] - summary["onef_accepted"]
        summary.setdefault("repairs", []).append(
            {"method": "apply successful dither 1/f fits regardless of held-out QC",
             "frames": len(replacements)})
        temp = summary_path.with_name(summary_path.name + ".writing")
        temp.write_text(json.dumps(summary, indent=2) + "\n")
        os.replace(temp, summary_path)
        return 0
    if args.repair_mask_failures:
        summary_path = output_dir / "summary.json"
        summary = json.loads(summary_path.read_text())
        keys = sorted({(record["file"].split("_")[0],
                        record["file"].split("_")[1], record["detector"],
                        record["filter"])
                       for record in summary["results"]
                       if record["onef"].get("reason") ==
                       "Too many masked rows in an amplifier"})
        tasks = [sorted(groups[key]) for key in keys]
        if not tasks:
            print("No amplifier-row mask failures to repair", flush=True)
            return 0
        print(f"Repairing {len(tasks)} detector groups with {args.workers} workers",
              flush=True)
        repaired, errors = run_phase(tasks, output_dir, args.template_dir,
                                     args.workers, "REPAIR")
        if errors:
            print(json.dumps(errors, indent=2), file=sys.stderr)
            return 1
        patch_contacts(repaired, output_dir)
        replacements = {r["file"]: r for group_result in repaired
                        for r in group_result["records"]}
        summary["results"] = [replacements.get(r["file"], r)
                              for r in summary["results"]]
        summary["onef_accepted"] = sum(r["onef"]["accepted"]
                                        for r in summary["results"])
        summary["onef_rejected"] = summary["completed"] - summary["onef_accepted"]
        summary["onef_applied"] = sum(r["onef"].get("applied", r["onef"]["accepted"])
                                      for r in summary["results"])
        summary.setdefault("repairs", []).append(
            {"method": "unmask bright-wisp pixels on otherwise empty amplifier rows",
             "groups": [list(key) for key in keys]})
        temp = summary_path.with_name(summary_path.name + ".writing")
        temp.write_text(json.dumps(summary, indent=2) + "\n")
        os.replace(temp, summary_path)
        return 0
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Starting {len(paths)} exposures in {len(groups)} dither groups; "
          f"up to {args.workers} workers", flush=True)
    a, ae = run_phase(lw, output_dir, args.template_dir, args.workers, "LW")
    if ae:
        errors = ae
        b = []
    else:
        b, be = run_phase(sw, output_dir, args.template_dir, args.workers, "SW")
        errors = ae + be
    all_results = a + b
    contact_paths = render_contacts(all_results, output_dir) if not errors else []
    records = sorted([x for result in all_results for x in result["records"]],
                     key=lambda x: x["file"])
    summary = {"method": "always apply dither 1/f -> masked median sky -> NMF wisp -> scalar zero",
               "workers": args.workers, "selected": len(paths),
               "completed": len(records), "onef_accepted": sum(x["onef"]["accepted"] for x in records),
               "onef_rejected": sum(not x["onef"]["accepted"] for x in records),
               "onef_applied": sum(x["onef"]["applied"] for x in records),
               "figures": [str(p) for p in contact_paths],
               "results": records, "errors": errors}
    temp = output_dir / "summary.json.writing"
    temp.write_text(json.dumps(summary, indent=2) + "\n")
    os.replace(temp, output_dir / "summary.json")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
