# Agent guide for stage2plus

This file applies to the standalone `stage2plus` repository. Read
[README.md](README.md) for supported CLI use and [CONTEXT.md](CONTEXT.md) for
the scientific decision history and validation record before changing the
method. The upstream repository is `ralphy1017/stage2plus` on `main`.

## Identify the code and run before acting

- Check `git status`, the current commit, and the input/output paths first.
  The pushed v2.0 baseline is commit `bcc1514`. As of 2026-11-02 this
  development checkout also had an **uncommitted** change in `stage2plus.py`
  adding `--no-onef-object-mask`. Keep that user work separate from other
  commits; a clean clone of `main` may not have the option.
- `stage2plus.py` owns CLI validation, FITS-header grouping, long-wave then
  short-wave scheduling, product assembly, metadata, and summary writing.
  `steps.py` owns masks, WCS projection, registered-dither 1/f, masked sky,
  `Background2D`, NMF wisp fitting, and FITS copying. `figures.py` renders a
  step figure for each exposure. `templates/` contains 16 existing NMF basis
  FITS assets. `requirements.txt` lists the Python packages.
- Distinguish the old v1 products in the parent `diagnostics/` directory from
  v2 products. The 640-exposure v1 batch is **not** a v2 deployment. The
  verified v2 pilot is `/data5/rzi/jp/stage2plus-test-latest/` (16 exposures).

## Preserve the science contract

1. Keep the current order: Photutils source mask, registered-dither 1/f,
   object-masked median sky, `Background2D` residual sky, NMF wisp where a
   template exists, then masked scalar pedestal. The sky fits are performed
   before the wisp fit. The extra four-amplifier/horizontal background fits
   from earlier pilots are not in this workflow. PHANGS pjpipe destriping was
   tried and dropped at the user's request.
2. Write both products from the **same** 1/f, wisp, and pedestal models:
   `cal_sky_renorm/` subtracts both sky models; `cal/` retains them. Do not
   apply additional sky subtraction to `cal/`. Do not use sky-subtracted SCI
   to estimate absolute sky surface brightness. Retaining the sky in `cal/`
   is not proof that all diffuse sky modes are unbiased.
3. Use `SCI_ORIG` when present, otherwise input `SCI`. Never alter input FITS.
   Keep PRIMARY, ERR, DQ, AREA, VAR_POISSON, VAR_RNOISE, VAR_FLAT, and ASDF
   data in both outputs, replacing only SCI and adding provenance
   cards/history. Preserve FITS checksums.
4. Apply every successful dither 1/f model. `DITH1F` records application;
   `DITHQC` is the independent held-out check. A failed QC check does not
   suppress the correction. A **fit failure** must surface as a group error;
   do not silently write an uncorrected product.
5. Group by FITS `PROGRAM`, `OBSERVTN`, `VISIT`, `VISITGRP`, `SEQ_ID`,
   `ACT_ID`, detector, filter, and pupil; require at least four dithers. Run
   long-wave groups first so matched masks can be projected to short-wave
   frames. If no paired long-wave exposure exists, document that the
   short-wave mask used only its own segmentation. Do not mix different
   upstream calibration versions or duplicate basenames in one manifest.
6. The default object threshold is 2.8σ with at least five connected pixels;
   mask growth is capped at 32 pixels. Protect bright wisp regions in 1/f
   fits while retaining at least 64 usable pixels per amplifier row. The sky
   mask also excludes NMF wisp support grown by eight pixels. `Background2D`
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
  and QC failures. `DITHQC=False` is a review flag, not a missing 1/f step.
  A group that fails after writing some exposures can leave partial products
  outside `summary.json`; reconcile or isolate them before retrying.
- Verify SCI headers (`S2PVER`, `STGORDER`, `DITH1F`, `DITHQC`, `SKYSUB`,
  `SKYMED`, `BKG2D` where applicable, `WISPMETH`, `PEDESTAL`), FITS
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
