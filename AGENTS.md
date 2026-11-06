# Agent guide for stage2plus

This file applies to the standalone `stage2plus` repository. Read
[README.md](README.md) for supported CLI use, [CONTEXT.md](CONTEXT.md) for the
scientific decision history, and [DEVELOPMENT.md](DEVELOPMENT.md) for the
2026-11-03 SKYSURFIR pilot findings and remaining validation goals. v2.7
combines the single-frame
`keep_trend=True` scan with an amplifier-local 2D low-k fit; dither-first is
still opt-in via `--onef-mode auto`. The three-frame test figures are at
`diagnostics/amp-scene-decomposition-20261105/`. The test used existing
production outputs read-only and did not modify the active production run.
The user reviewed the figures and requested publishing the candidate to
`ralphy1017/stage2plus` `main`; this review does not constitute broad science
validation or a production rerun.

The earlier v2.6 20-exposure rerun is at
`/data5/rzi/jp/stage2plus-fallback-default-mini-20261103/`; the user reviewed
those figures and requested that scan default be published. Check the run
summary and this worktree's status before acting.

## Identify the code and run before acting

- Check `git status`, the current commit, and the input/output paths first.
  The local v2.7 default is the source-unmasked scan plus amplifier-local 2D
  fit; `--onef-mode auto` opts into the v2.5 dither-first workflow. The
  three-frame test is not a production rerun or broad science validation.
  `--no-onef-object-mask` controls source masking for the dither fit only.
- `stage2plus.py` owns CLI validation, FITS-header grouping, long-wave then
  short-wave scheduling, product assembly, metadata, and summary writing.
  `steps.py` owns masks, WCS projection, registered-dither and single-frame
  scan/2D amplifier 1/f fits, masked sky, `Background2D`, NMF wisp fitting, and FITS
  copying. `figures.py` renders a step figure for each exposure. `templates/`
  contains 16 existing NMF basis FITS assets. `requirements.txt` lists the
  Python packages.
- Distinguish the old v1 products in the parent `diagnostics/` directory from
  v2 products. The 640-exposure v1 batch is **not** a v2 deployment. The
  verified v2 pilot is `/data5/rzi/jp/stage2plus-test-latest/` (16 exposures).

## Preserve the science contract

1. Keep the current order: Photutils source mask, selected dither-or-scan 1/f,
   object-masked median sky, `Background2D` residual sky, NMF wisp where a
   template exists, then masked scalar pedestal. The sky fits are performed
   before the wisp fit. The local v2.7 candidate adds only the per-frame
   amplifier-local 2D low-k term described below; do not add a separate
   background subtraction pass. PHANGS pjpipe destriping was tried and
   dropped at the user's request.
2. Write both products from the **same** 1/f, wisp, and pedestal models:
   `cal_sky_renorm/` subtracts both sky models; `cal/` retains them. Do not
   apply additional sky subtraction to `cal/`. Do not use sky-subtracted SCI
   to estimate absolute sky surface brightness. Retaining the sky in `cal/`
   is not proof that all diffuse sky modes are unbiased.
3. Use `SCI_ORIG` when present, otherwise input `SCI`. Never alter input FITS.
   Keep PRIMARY, ERR, DQ, AREA, VAR_POISSON, VAR_RNOISE, VAR_FLAT, and ASDF
   data in both outputs, replacing only SCI and adding provenance
   cards/history. Preserve FITS checksums.
4. The local v2.7 default `--onef-mode fallback` fits the Jumprope-style
   single-frame scan (`keep_trend=True`) with per-amplifier 2D surfaces from
   the same frame. A full-detector cubic surface is kept as shared low-k
   structure; each amplifier's quadratic surface minus that shared surface
   joins the scan in one combined model. Only detector-invalid pixels are
   excluded, no object mask or dither peers are used, and the global median is
   preserved. `--onef-mode auto` opts into
   dither fitting; it applies every computable dither model even when optional
   held-out QC fails, and uses the scan only when the dither scene/row model
   cannot be computed. `DITHQC` records optional diagnostics separately; it
   does not select the method. The v2.5 dither method uses all overlapping
   peers, source-based
   translation correction over the full WCS, scalar background matching, a
   sigma-clipped median scene, and sigma-clipped per-amplifier row medians.
   The former 0.1-pixel source residual and 5% held-out RMS rules are optional
   diagnostics, not acceptance gates. `--dither-diagnostics` withholds peers
   and changes the fitted scene; it is for validation, not default processing.
   The scan excludes detector-invalid pixels but does not mask sources or
   wisps; its 90th-percentile value clip
   remains. `DITH1F` and `ONEFAPPL` both record that a 1/f correction was
   applied; `ONEFMETH` identifies `DITHER` or `SELF`. A dither model without a
   passing optional QC is marked `ONEFSTAT=REVIEW`. A selected scan failure,
   or failure of both methods in auto mode, must fail the group rather than
   write an uncorrected product.
5. Group by FITS `PROGRAM`, `OBSERVTN`, `VISIT`, `VISITGRP`, `SEQ_ID`,
   `ACT_ID`, detector, filter, and pupil; require at least one exposure.
   In `auto` mode, dither fitting needs one overlapping peer (two group
   exposures); a group with one exposure uses the scan. Run long-wave groups first so
   matched masks can be projected to short-wave frames. If no paired long-wave
   exposure exists, document that the
   short-wave mask used only its own segmentation. Do not mix different
   upstream calibration versions or duplicate basenames in one manifest.
6. The default object threshold is 2.8σ with at least five connected pixels;
   mask growth is capped at 32 pixels. Neither 1/f method excludes the bright
   wisp region by mask. The sky mask still excludes NMF wisp support grown by
   eight pixels. `Background2D`
   uses 256-pixel boxes, a 3×3 filter, biweight location, `SigmaClip()`, and
   `exclude_percentile=50`. The photutils default of 10 failed the tested
   grown-mask fields; change these settings only with a clear test.
7. Treat missing wisp templates as an explicit review item. v2 currently
   records `WISPMETH=NONE` when no basis is found; `--plan-only` counts
   template-bearing groups but does not reject all missing combinations.
   Check that this is intentional for the chosen filters and detectors.

## Run and verify batches

- On the development machine, use the shared CPU lease launcher at
  `/home/rzi/Work/pipeline/tools/checks/lease.py` and set `--workers` equal
  to the CPUs leased. Put new production runs under `/data5/rzi/jp/<run>/`
  with a separate output directory for each code revision and settings set.
  The host-specific command shape is `python /home/rzi/Work/pipeline/tools/checks/lease.py N python stage2plus.py --input-list PATH --output-dir PATH --workers N`.
  The script uses one process per dither group and limits BLAS/OpenMP threads
  to one in each process. A v2 test peaked near 1.4 GB RSS for the measured
  process; allow headroom for concurrent groups, masks, and I/O cache.
- Start with an immutable newline-separated input manifest. Run
  `--plan-only`, record the Git commit and `git diff` state, check group sizes,
  long-wave matches, wisp-basis coverage, estimated output bytes, and free
  space. For a new program or configuration, run a small long-wave plus
  short-wave pilot before scaling. Capture stdout/stderr in `run.log`; the
  script itself writes `summary.json`, not a log.
- After a run, require `summary.json` to report the intended `selected` and
  `completed` counts and an empty `errors` list. Count files in `cal/`,
  `cal_sky_renorm/`, and `models/`, and count figures unless `--no-figures`
  was selected. Inspect representative figures across detectors, filters,
  dither fits, optional QC failures, and scan fallbacks. For `DITHQC=False`,
  check that `ONEFMETH=DITHER`, `ONEFAPPL=True`, and `ONEFSTAT=REVIEW` agree;
  the computed dither model remains applied. With optional diagnostics off,
  DITHQC is null in JSON and omitted from the SCI header.
  A group that fails after writing some exposures can leave partial products
  outside `summary.json`; reconcile or isolate them before retrying.
- Verify SCI headers (`S2PVER`, `STGORDER`, `ONEFMETH`, `ONEFSTAT`,
  `ONEFAPPL`, `ONEFVER`, `DITH1F`, `DITHQC`, `SKYSUB`, `SKYMED`, `BKG2D`
  where applicable, `WISPMETH`, `PEDESTAL`), FITS
  checksums, and preserved non-SCI arrays on representative outputs. Check
  arithmetic closure against `models/`: skyless SCI should match input minus
  `ONEF`, `SKYMED`, `SKYBKG`, optional `WISP`, and `PEDESTAL` to float32
  precision; `cal/SCI` should match input minus `ONEF`, optional `WISP`, and
  `PEDESTAL`. The final masked skyless median should be near zero.
- Do not infer full-database reliability from the 16-exposure v2 pilot. A
  100k-frame campaign needs staged storage/throughput tests, sharded input
  manifests, failure recovery, and sky-signal injection checks. The current
  CLI has no built-in resume or version compatibility check and can replace
  files in an existing output directory.

## Documentation and Git hygiene

- Keep README as the user run guide, CONTEXT as the historical and scientific
  record, and this file as the operational contract. Update all three when a
  change alters order, masks, fitted models, FITS headers, supported inputs,
  or resource assumptions. Record actual run evidence and limits, not just
  visual impressions.
- Do not stage unrelated local experiments or generated FITS/PNGs with a
  documentation or code change. Check the staged diff before committing and
  the remote branch before pushing. Push only when the user has authorized
  it. Keep machine-specific input/output paths as examples; the clone must
  run with user-supplied paths and bundled templates.
