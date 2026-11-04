#!/usr/bin/env python3
"""Stage 2 plus: clean JWST NIRCam dither groups before drizzling.

For each exposure in a dither group (same visit, detector, filter and pupil)
it writes two products from the input science. The default 1/f model uses all
overlapping peers, full-WCS registration, background matching, a sigma-clipped
median scene, and sigma-clipped per-amplifier row offsets. The single-frame
scan is used only when a dither model cannot be computed:

  cal_sky_renorm/  photutils object mask -> registered-dither 1/f -> masked
                   median sky -> Background2D residual sky -> NMF wisp (where
                   a template exists) -> masked pedestal renormalisation
  cal/             the same corrections without the median and Background2D
                   sky, so the measured sky is retained

Models and masks go to models/, one step-by-step figure per exposure to
figures/, and per-exposure measurements to summary.json. Long-wave groups run
first so their object masks can be projected onto matching short-wave frames.
"""

from __future__ import annotations

import os

for _key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[_key] = "1"

import argparse
import concurrent.futures as futures
import json
import multiprocessing
import re
import shutil
import sys
import traceback
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi

import steps
from figures import plot_steps

VERSION = "2.5"
MIN_DITHERS = 1
DITHER_MIN_DITHERS = 2
DETECTOR_PATTERN = re.compile(r"^NRC[AB](?:[1-4]|LONG)$")
VISIT_KEYS = ("PROGRAM", "OBSERVTN", "VISIT", "VISITGRP", "SEQ_ID", "ACT_ID")


def header_info(path: str) -> dict:
    header = fits.getheader(path, 0)
    value = lambda key: str(header.get(key, "")).strip().upper()
    return {"visit": tuple(value(k) for k in VISIT_KEYS), "exposure": value("EXPOSURE"),
            "detector": value("DETECTOR"), "filter": value("FILTER"), "pupil": value("PUPIL")}


def group_label(info: dict) -> str:
    program, obs, visit, visitgrp, seq, act = info["visit"]
    return "_".join((f"jw{program}{obs}{visit}", f"{visitgrp}{seq}{act}", info["detector"],
                     info["filter"], info["pupil"])).lower()


def fit_onef(exp: steps.Exposure, peers: list[steps.Exposure], good: np.ndarray,
             fallback_bad: np.ndarray, dither_object_mask: bool,
             peer_exclusions: list[np.ndarray], target_sources: np.ndarray,
             peer_sources: list[np.ndarray], dither_diagnostics: bool
             ) -> tuple[np.ndarray, dict]:
    """Use every computable dither model; scan only when the fit cannot run."""
    dither_meta = None
    dither_error = None
    dither_model = None
    try:
        dither_model, dither_meta = steps.dither_onef(
            exp, peers, good, peer_exclusions, target_sources, peer_sources,
            diagnostics=dither_diagnostics)
        status = "PASS" if dither_meta["qc_accepted"] is True else "REVIEW"
        return dither_model, {**dither_meta, "method": "DITHER",
                              "status": status, "applied": True,
                              "object_mask": dither_object_mask,
                              "dither_attempted": True,
                              "fallback_attempted": False,
                              "fallback_reason": None}
    except (ValueError, FloatingPointError) as exc:
        dither_error = f"{type(exc).__name__}: {exc}"

    # Do not keep the rejected full-frame candidate alive while allocating the
    # scan fallback; large groups can otherwise have a high peak RSS here.
    del dither_model
    try:
        # The single-frame scan fallback deliberately includes detected and
        # extended sources and wisp-affected pixels in its profiles. Only
        # detector-invalid pixels are excluded by the supplied mask.
        model, fallback_meta = steps.scan_onef(exp.science, fallback_bad)
    except Exception as exc:
        raise ValueError(f"Dither 1/f failed ({dither_error}); single-frame "
                         f"fallback failed ({type(exc).__name__}: {exc})") from exc
    return model, {"method": "SELF", "status": "REVIEW", "applied": True,
                   "object_mask": False,
                   "fit_mask": "INVALID ONLY; SOURCES+WISP INCLUDED",
                   "qc_accepted": (dither_meta["qc_accepted"]
                                   if dither_meta is not None else None),
                   "dither_attempted": True,
                   "dither": dither_meta,
                   "dither_error": dither_error,
                   "fallback_attempted": True,
                   "fallback_reason": dither_error,
                   "fallback": fallback_meta,
                   "model_rms": fallback_meta["model_rms"],
                   "fit_fraction": fallback_meta["fit_fraction"]}


def process_group(label: str, paths: list[str], cfg: dict) -> dict:
    out = Path(cfg["output_dir"])
    exposures = [steps.read_exposure(Path(p)) for p in paths]
    detector, band = exposures[0].detector, exposures[0].band
    template = steps.basis_path(Path(cfg["template_dir"]), detector, band)
    basis = steps.load_basis(template, exposures[0].science.shape) if template else None

    # 1. Object masks: own segmentation plus the matching long-wave mask.
    objects, grown, onef_masks, fallback_bad_masks = [], [], [], []
    for exp in exposures:
        mask = steps.segment(exp)
        lw_models = cfg["lw_models"].get(exp.path.name)
        if lw_models:
            with fits.open(lw_models, memmap=False) as hdul:
                mask |= steps.project_mask(hdul["OBJMASK"].data != 0,
                                           hdul["OBJMASK"].header, exp)
        objects.append(mask)
        grown.append(steps.grow_mask(mask))
        onef = (exp.invalid | grown[-1]) if cfg["onef_object_mask"] else exp.invalid.copy()
        fallback_bad = exp.invalid.copy()
        onef_masks.append(onef)
        fallback_bad_masks.append(fallback_bad)

    records = []
    support = (ndi.binary_dilation(basis[1], iterations=steps.WISP_SUPPORT_GROWTH)
               if basis is not None else None)
    for i, exp in enumerate(exposures):
        # 2. Registered-dither 1/f, with a source-unmasked scan fallback.
        peer_indices = [j for j in range(len(exposures)) if j != i]
        peers = [exposures[j] for j in peer_indices]
        onef, onef_meta = fit_onef(
            exp, peers, ~onef_masks[i], fallback_bad_masks[i],
            cfg["onef_object_mask"],
            [onef_masks[j] for j in peer_indices], objects[i],
            [objects[j] for j in peer_indices], cfg["dither_diagnostics"])
        after_onef = exp.science - onef
        sky_mask = exp.invalid | grown[i] | ~np.isfinite(after_onef)
        if support is not None:
            sky_mask |= support
        sky_good = ~sky_mask

        # 3. Masked median sky, then 4. JADES Background2D residual sky.
        sky_median = steps.masked_median(after_onef, sky_good)
        skyless = after_onef - np.float32(sky_median)
        background, bkg_meta = steps.background2d(skyless, sky_mask,
                                                  cfg["bkg_exclude_percentile"])
        no_background = skyless - background

        # 5. NMF wisp, fitted on the sky-subtracted frame.
        wisp, wisp_meta = None, None
        if basis is not None:
            wisp, wisp_meta = steps.fit_wisp(no_background, exp.err, exp.invalid | grown[i],
                                             *basis, template)
        after_wisp = no_background - wisp if wisp is not None else no_background

        # 6. Pedestal: rezero the masked median; cal keeps the measured sky.
        pedestal = steps.masked_median(after_wisp, sky_good)
        final_sky = after_wisp - np.float32(pedestal)
        final_cal = after_onef - np.float32(pedestal)
        if wisp is not None:
            final_cal -= wisp

        name = exp.path.name
        cards = {"S2PVER": (VERSION, "stage2plus version"),
                 "ONEFMETH": (onef_meta["method"], "selected 1/f method"),
                 "ONEFSTAT": (onef_meta["status"], "selected method QC status"),
                 "ONEFAPPL": (onef_meta["applied"], "1/f model subtracted"),
                 "ONEFVER": ("DITH-PH25" if onef_meta["method"] == "DITHER"
                             else "SCAN-KT3", "selected 1/f implementation"),
                 "DITH1F": (onef_meta["applied"], "1/f correction applied"),
                 "ONEFOBJ": (onef_meta["object_mask"],
                             "objects masked in selected 1/f fit"),
                 "SKYMED": (sky_median, "masked median sky (MJy/sr)"),
                 "WISPMETH": ("NMF" if wisp is not None else "NONE", "wisp model"),
                 "PEDESTAL": (pedestal, "post-wisp masked median removed (MJy/sr)")}
        if onef_meta["qc_accepted"] is not None:
            cards["DITHQC"] = (onef_meta["qc_accepted"],
                               "held-out dither 1/f check passed")
        if onef_meta["method"] == "DITHER":
            cards["DITHN"] = (len(onef_meta["scene_dithers"]),
                              "dithers in selected 1/f scene")
        steps.write_fits(exp, out / "cal_sky_renorm" / name, final_sky,
                         cards={**cards, "STGORDER": ("MASK-1F-SKY-BKG2D-WISP-PED", "steps"),
                                "SKYSUB": (True, "median sky and Background2D removed"),
                                "BKG2D": (True, "Background2D residual sky removed")},
                         history="SCI = input - 1/f - median sky - Background2D - wisp - pedestal")
        steps.write_fits(exp, out / "cal" / name, final_cal,
                         cards={**cards, "STGORDER": ("MASK-1F-WISP-PED", "steps"),
                                "SKYSUB": (False, "sky retained")},
                         history="SCI = input - 1/f - wisp - pedestal (sky retained)")
        write_models(exp, out / "models" / name, objects[i], sky_mask, onef,
                     background, wisp, sky_median, pedestal, onef_meta)

        noise = steps.robust_sigma(final_sky[sky_good][::40])
        if cfg["figures"]:
            plot_steps(out / "figures" / (Path(name).stem + "_steps.png"),
                       f"{name} | {detector} {band} | 1/f {onef_meta['method']} "
                       f"({onef_meta['status']})",
                       [("Input SCI", exp.science, "sky"),
                        (("Source mask shown (not used by SELF); sky mask (grey)"
                          if onef_meta["method"] == "SELF" else
                          "Object mask (black), sky mask (grey)"),
                         sky_mask.astype(np.uint8) + objects[i], "mask"),
                        (f"{onef_meta['method']} 1/f model", onef, "model"),
                        ("After 1/f", after_onef, "sky"),
                        (f"After median sky ({sky_median:.4f})", skyless, "skyless"),
                        ("Background2D model", background, "model"),
                        ("After Background2D", no_background, "skyless"),
                        ("NMF wisp model", wisp, "model"),
                        (f"Final cal_sky_renorm (pedestal {pedestal:+.2e})", final_sky, "skyless"),
                        ("Final cal (sky retained)", final_cal, "sky")],
                       exp.invalid, sky_median, noise)

        records.append({"input": str(exp.path), "file": name, "group": label,
                        "detector": detector, "filter": band, "dithers": len(exposures),
                        "object_fraction": float(objects[i].mean()),
                        "sky_mask_fraction": float(sky_mask.mean()),
                        "longwave_mask": cfg["lw_models"].get(name),
                        "onef": onef_meta,
                        "sky_median": sky_median, "background2d": bkg_meta,
                        "wisp": wisp_meta, "pedestal": pedestal,
                        "final_masked_median": float(np.median(final_sky[sky_good])),
                        "noise": noise})
    return {"group": label, "records": records}


def write_models(exp, path: Path, objects, sky_mask, onef, background, wisp,
                 sky_median: float, pedestal: float, onef_meta: dict) -> None:
    primary = fits.PrimaryHDU()
    primary.header["SKYMED"] = (sky_median, "masked median sky (MJy/sr)")
    primary.header["PEDESTAL"] = (pedestal, "post-wisp masked median (MJy/sr)")
    primary.header["ONEFMETH"] = (onef_meta["method"], "selected 1/f method")
    primary.header["ONEFSTAT"] = (onef_meta["status"], "selected method QC status")
    primary.header["ONEFAPPL"] = (onef_meta["applied"], "1/f model subtracted")
    primary.header["ONEFVER"] = ("DITH-PH25" if onef_meta["method"] == "DITHER"
                                  else "SCAN-KT3", "selected 1/f implementation")
    primary.header["ONEFOBJ"] = (onef_meta["object_mask"],
                                  "objects masked in selected 1/f fit")
    primary.header["DITH1F"] = (onef_meta["applied"], "1/f correction applied")
    if onef_meta["qc_accepted"] is not None:
        primary.header["DITHQC"] = (onef_meta["qc_accepted"], "dither QC passed")
    hdus = [primary,
            fits.ImageHDU(objects.astype(np.uint8), exp.sci_header, name="OBJMASK"),
            fits.ImageHDU(sky_mask.astype(np.uint8), exp.sci_header, name="SKYMASK"),
            fits.ImageHDU(onef.astype(np.float32), exp.sci_header, name="ONEF"),
            fits.ImageHDU(background.astype(np.float32), exp.sci_header, name="SKYBKG")]
    if wisp is not None:
        hdus.append(fits.ImageHDU(wisp.astype(np.float32), exp.sci_header, name="WISP"))
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + ".writing")
    fits.HDUList(hdus).writeto(temp, overwrite=True, output_verify="silentfix")
    temp.replace(path)


def run_groups(groups: dict[str, list[str]], cfg: dict, workers: int, phase: str):
    results, errors = [], []
    if not groups:
        return results, errors
    context = multiprocessing.get_context("spawn")
    with futures.ProcessPoolExecutor(max_workers=min(workers, len(groups)),
                                     mp_context=context) as pool:
        pending = {pool.submit(process_group, label, paths, cfg): label
                   for label, paths in groups.items()}
        for job in futures.as_completed(pending):
            label = pending[job]
            try:
                result = job.result()
            except Exception as exc:
                errors.append({"group": label, "error": repr(exc),
                               "traceback": traceback.format_exc()})
                print(f"{phase} FAILED {label}: {exc}", file=sys.stderr, flush=True)
                continue
            results.append(result)
            qc_checked = sum(r["onef"]["qc_accepted"] is not None
                             for r in result["records"])
            dither_selected = sum(r["onef"]["method"] == "DITHER"
                                  for r in result["records"])
            fallbacks = sum(r["onef"]["method"] == "SELF"
                            for r in result["records"])
            print(f"{phase} {len(results)}/{len(groups)} {label}: dither models "
                  f"{dither_selected}/{len(result['records'])}; scan fallbacks "
                  f"{fallbacks}; optional QC checked {qc_checked}", flush=True)
    return results, errors


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-list", type=Path, required=True,
                        help="text file with one calibrated NIRCam *_cal.fits path per line")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--template-dir", type=Path,
                        default=Path(__file__).resolve().parent / "templates")
    parser.add_argument("--workers", type=int, default=8,
                        help="dither groups processed in parallel (one CPU each)")
    parser.add_argument("--bkg-exclude-percentile", type=float, default=50.0,
                        help="Background2D exclude_percentile: skip boxes with more than "
                             "this percent masked (default 50; photutils default 10)")
    parser.add_argument("--no-onef-object-mask", action="store_true",
                        help="disable source masking for registered-dither fits; the scan fallback is always source-unmasked")
    parser.add_argument("--dither-diagnostics", action="store_true",
                        help="hold out up to two peers for optional RMS diagnostics; diagnostics never gate model use")
    parser.add_argument("--no-figures", action="store_true")
    parser.add_argument("--plan-only", action="store_true")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be positive")

    paths = [str(Path(line.strip()).resolve()) for line in
             args.input_list.read_text().splitlines() if line.strip()]
    if not paths or len(paths) != len(set(paths)) or any(not Path(p).is_file() for p in paths):
        parser.error("input list is empty, duplicated, or has missing files")
    if len({Path(p).name for p in paths}) != len(paths):
        parser.error("duplicate basenames would collide in the output")
    output_dir = args.output_dir.resolve()
    if output_dir in {Path(p).parent for p in paths}:
        parser.error("output directory must differ from the input directories")

    infos = {p: header_info(p) for p in paths}
    unsupported = [p for p, i in infos.items() if not DETECTOR_PATTERN.fullmatch(i["detector"])]
    if unsupported:
        parser.error(f"not NIRCam imaging detectors: {unsupported[:5]}")
    groups = defaultdict(list)
    for path, info in infos.items():
        groups[group_label(info)].append(path)
    small = sorted(label for label, files in groups.items() if len(files) < MIN_DITHERS)
    if small:
        parser.error(f"{len(small)} groups have fewer than {MIN_DITHERS} dithers: {small[:10]}")
    long_groups = {k: sorted(v) for k, v in sorted(groups.items()) if "long" in k}
    short_groups = {k: sorted(v) for k, v in sorted(groups.items()) if "long" not in k}

    # Short-wave exposures use the object mask of the same long-wave exposure.
    longwave = {(i["visit"], i["exposure"], i["detector"][3]): Path(p).name
                for p, i in infos.items() if i["detector"].endswith("LONG")}
    lw_models = {}
    for path, info in infos.items():
        match = longwave.get((info["visit"], info["exposure"], info["detector"][3]))
        if match and not info["detector"].endswith("LONG"):
            lw_models[Path(path).name] = str(output_dir / "models" / match)

    with_wisp = sorted(k for k, v in groups.items()
                       if steps.basis_path(args.template_dir, infos[v[0]]["detector"],
                                           infos[v[0]]["filter"]))
    if args.plan_only:
        with fits.open(paths[0], memmap=True, lazy_load_hdus=True) as hdul:
            kept = sum(h._data_size + h._data_offset - h._header_offset
                       for h in hdul if h.name.upper() in steps.KEEP_HDUS)
            npix = hdul["SCI"].data.size
        print(json.dumps({"exposures": len(paths), "groups": len(groups),
                          "longwave_groups": len(long_groups),
                          "shortwave_groups": len(short_groups),
                          "dithers_per_group": dict(Counter(len(v) for v in groups.values())),
                          "dither_eligible_groups": sum(len(v) >= DITHER_MIN_DITHERS
                                                        for v in groups.values()),
                          "scan_fallback_groups": sum(len(v) < DITHER_MIN_DITHERS
                                                       for v in groups.values()),
                          "wisp_template_groups": len(with_wisp),
                          "shortwave_with_longwave_mask": len(lw_models),
                          "workers": args.workers,
                          "estimated_output_bytes": len(paths) * (2 * kept + 14 * npix),
                          "free_bytes": shutil.disk_usage(output_dir.parent).free}, indent=2))
        return 0

    output_dir.mkdir(parents=True, exist_ok=True)
    cfg = {"output_dir": str(output_dir), "template_dir": str(args.template_dir.resolve()),
           "figures": not args.no_figures,
           "bkg_exclude_percentile": args.bkg_exclude_percentile, "lw_models": lw_models,
           "onef_object_mask": not args.no_onef_object_mask,
           "dither_diagnostics": args.dither_diagnostics}
    print(f"Starting {len(paths)} exposures in {len(groups)} dither groups "
          f"({len(long_groups)} LW, {len(short_groups)} SW); up to {args.workers} workers",
          flush=True)
    results, errors = run_groups(long_groups, cfg, args.workers, "LW")
    if not errors:
        short_results, errors = run_groups(short_groups, cfg, args.workers, "SW")
        results += short_results
    records = sorted((r for result in results for r in result["records"]),
                     key=lambda r: r["file"])
    summary = {"stage2plus": VERSION,
               "method": "object mask -> PHANGS-style multi-exposure dither 1/f "
                         "or Jumprope-style scan fallback "
                         "-> masked median "
                         "sky -> Background2D -> NMF wisp -> pedestal; cal/ omits both sky steps",
               "settings": {"min_dithers": MIN_DITHERS, "workers": args.workers,
                            "dither_min_dithers": DITHER_MIN_DITHERS,
                            "onef_mode": "AUTO",
                            "fallback_method": "PROFOUND_SCAN_KEEP_TREND",
                            "dither_diagnostics": cfg["dither_diagnostics"],
                            "onef_object_mask": cfg["onef_object_mask"],
                            "fallback_object_mask": False,
                            "dither_fit_mask": ("INVALID+SOURCE; WISP INCLUDED"
                                                if cfg["onef_object_mask"]
                                                else "INVALID ONLY; SOURCES+WISP INCLUDED"),
                            "fallback_fit_mask": "INVALID ONLY; SOURCES+WISP INCLUDED",
                            "dither_wisp_mask": False,
                            "fallback_wisp_mask": False,
                            "dither_reprojection": "full WCS plus matched-source translation",
                            "dither_scene": "sigma-clipped median of all overlapping peers by default",
                            "dither_row_model": "sigma-clipped row median per amplifier",
                            "reproject_valid_coverage_min": steps.REPROJECT_VALID_COVERAGE,
                            "source_alignment_tolerance_px": steps.SOURCE_ALIGNMENT_TOL_PIXELS,
                            "source_match_radius_px": steps.SOURCE_MATCH_RADIUS_PIXELS,
                            "min_source_alignment_matches": steps.MIN_SOURCE_ALIGNMENT_MATCHES,
                            "segment_nsigma": steps.SEGMENT_NSIGMA,
                            "mask_growth_cap": steps.MASK_GROWTH_CAP,
                            "wisp_support_growth": steps.WISP_SUPPORT_GROWTH,
                            "max_heldout": steps.MAX_HELDOUT,
                            "template_dir": cfg["template_dir"]},
               "selected": len(paths), "completed": len(records),
               "onef_applied": sum(r["onef"]["applied"] for r in records),
               "onef_dither_selected": sum(r["onef"]["method"] == "DITHER"
                                           for r in records),
               "onef_fallback_selected": sum(r["onef"]["method"] == "SELF"
                                             for r in records),
               "onef_qc_passed": sum(r["onef"]["qc_accepted"] is True for r in records),
               "onef_qc_failed": sum(r["onef"]["qc_accepted"] is False for r in records),
               "onef_qc_unavailable": sum(r["onef"]["qc_accepted"] is None for r in records),
               "results": records, "errors": errors}
    temp = output_dir / "summary.json.writing"
    temp.write_text(json.dumps(summary, indent=2, default=str) + "\n")
    temp.replace(output_dir / "summary.json")
    print(f"Done: {len(records)}/{len(paths)} exposures, {len(errors)} group errors", flush=True)
    return 0 if not errors else 1


if __name__ == "__main__":
    sys.exit(main())
