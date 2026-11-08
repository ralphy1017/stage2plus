# Agent guide for stage2plus

This file applies to the standalone `stage2plus` repository. Read
[README.md](README.md) for supported CLI use, [CONTEXT.md](CONTEXT.md) for the
scientific decision history, and [DEVELOPMENT.md](DEVELOPMENT.md) for recent
SKYSURFIR findings and remaining validation goals. Published v2.11 keeps
object masks out of the keep-trend scan, uses an initial
Photutils mask plus a 512-pixel coarse background and 2.0σ/3-pixel residual
segmentation to capture faint outskirts. The default order is scan, primary
NMF subtraction, masked amplifier-local and detector-wide y-only low-k fits,
then residual NMF. The two-frame NMF-order test is under
`diagnostics/wisp-order-mini-20261107/`; the original production outputs were
not modified. `cal/` subtracts only 1/f and both NMF corrections;
`cal_sky_renorm/` also gets median, `Background2D`, and pedestal subtraction.
This method has not been validated broadly.

The earlier v2.6 20-exposure rerun is at
`/data5/rzi/jp/stage2plus-fallback-default-mini-20261103/`; the user reviewed
those figures and requested that scan default be published. Check the run
summary and this worktree's status before acting.

## Identify the code and run before acting

- Check `git status`, the current commit, and the input/output paths first.
  The v2.11 default fits the source-unmasked scan, subtracts primary NMF,
  then fits object-masked amplifier and detector-wide low-k terms followed by
  a residual NMF pass. `--onef-mode auto` opts into dither-first processing.
  The two-frame injection test is not a production rerun or broad science
  validation.
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

1. Keep the v2.11 order for the default single-frame path: Photutils source
   segmentation and grown masks; source-unmasked keep-trend scan; primary NMF
   fit and subtraction; object-masked amplifier-local and detector-wide
   low-k fits; residual NMF fit and subtraction; object-masked median sky;
   `Background2D` residual sky; then masked scalar pedestal. The low-k maps
   are part of the one per-frame 1/f model, not a separate sky subtraction.
   In optional dither mode, fit the dither correction then one NMF model; the
   single-frame low-k fit is used only if dither fitting fails. PHANGS pjpipe
   destriping was tried and dropped at the user's request.
2. Keep the output branches scientifically distinct. `cal/` subtracts only
   the selected 1/f model (including its masked amplifier-local and
   detector-wide low-k terms) and the wisp template. Do not subtract its
   median, `Background2D`, or pedestal; object-free regions in this product
   are used for absolute sky measurements. `cal_sky_renorm/` additionally
   subtracts the object-masked median, `Background2D` residual, and final
   pedestal to produce a near-zero masked background. Do not use its SCI to
   measure absolute sky brightness.
3. Use `SCI_ORIG` when present, otherwise input `SCI`. Never alter input FITS.
   Keep PRIMARY, ERR, DQ, AREA, VAR_POISSON, VAR_RNOISE, VAR_FLAT, and ASDF
   data in both outputs, replacing only SCI and adding provenance
   cards/history. Preserve FITS checksums.
4. The v2.11 default `--onef-mode fallback` fits the Jumprope-style
   single-frame scan (`keep_trend=True`) using detector-invalid pixels only.
   It subtracts the primary NMF wisp fit before measuring robust 128-pixel
   tile medians outside the grown source mask plus four additional dilation
   pixels. A full-frame cubic surface is kept
   as shared low-k structure; each amplifier's quadratic surface minus that
   shared surface is a local correction. A degree-2 polynomial in detector
   y, constant across x, models broad horizontal structure. Apply the fitted
   low-k surfaces under source pixels and subtract the scan and both surfaces
   once as one model. `LOWKMASK`, `LOWKOBJ`, and `LOWKDIL` record the mask and
   settings. `--onef-mode auto` opts into
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
  `ONEFAPPL`, `ONEFVER`, `LOWKOBJ`, `LOWKDIL`, `DITH1F`, `DITHQC`, `SKYSUB`, `SKYMED`, `BKG2D`
  where applicable, `SEGNSIG`, `SEGNPIX`, `SEGBKG`, `WISPMETH`, `PEDESTAL`), FITS
  checksums, and preserved non-SCI arrays on representative outputs. Check
  arithmetic closure against `models/`: `cal_sky_renorm/SCI` should match
  input minus `ONEF`, `SKYMED`, `SKYBKG`, optional `WISP`, and `PEDESTAL` to
  float32 precision; `cal/SCI` should match input minus `ONEF` and optional
  `WISP`. Check that the skyless masked median is near zero. `cal/` must not
  be forced to any renormalized median.
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
