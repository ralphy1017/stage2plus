# stage2plus: project context and validation record

Updated 2026-11-06. This document explains why the current reducer has its
particular order and what has actually been tested. Use [README.md](README.md)
for installation and CLI examples; use [AGENTS.md](AGENTS.md) for operational
rules when changing or running the code.

## Purpose and scope

`stage2plus` prepares JWST NIRCam Stage 2 calibrated exposures for drizzling.
It operates on groups of dithers from one observation, detector, filter, and
pupil. The science goals are to suppress detector-fixed 1/f striping and NIRCam
wisps while keeping object flux out of the fitted backgrounds, and to provide
both a sky-retaining exposure and a sky-subtracted exposure. It is independent
of Jumprope's execution architecture. The bundled FITS bases are the existing
fixed NMF wisp templates; no Jumprope runtime is invoked.

The v2 pipeline has two outputs for each input:

| Output | SCI content | Intended use |
| --- | --- | --- |
| `cal/` | input − selected 1/f − NMF wisp | Exposure for absolute sky measurements; no median-sky, `Background2D`, or pedestal subtraction. |
| `cal_sky_renorm/` | input − selected 1/f − object-masked median sky − `Background2D` residual − NMF wisp − final scalar pedestal | Sky-subtracted drizzle input with near-zero masked background. |

`cal/` receives only the 1/f and wisp models; it does not receive the median,
`Background2D`, or pedestal subtraction, so object-free regions remain usable
for measuring the true sky. The wisp fit is measured on the temporarily
sky-subtracted image, but only its fitted template is subtracted from `cal/`.
“Renormalization” in `cal_sky_renorm/` means additive median and pedestal
corrections, not multiplication or photometric rescaling. Do not infer
absolute sky surface brightness from `cal_sky_renorm/SCI`. Tests do not prove
that the 1/f model cleanly separates detector artifacts from every true
diffuse astronomical component.

## How the method evolved

1. **Early wisp, sky, and single-exposure 1/f pilots.** Program 1813 tests
   compared NMF wisp removal with sky-first and 1/f-first orders, Photutils
   masks, mask dilation, median and polynomial sky fits, and a `keep_trend`
   setting for the single-exposure 1/f scan. Additional four-amplifier and
   broad horizontal background fits reduced visible steps and bands in some
   sky-subtracted frames. Those broad fits can also remove diffuse sky signal,
   so visual flatness alone was not accepted as a sky-retaining calibration
   criterion. The historical scripts and figures remain in the parent
   `jumprope/diagnostics/` tree and its `handoff.md`; they are not production
   dependencies of this repository.
2. **Registered-dither 1/f pilot.** Other dithers were aligned by WCS and
   level-matched only in a temporary scene template. Amplifier-row and
   detector-column medians of target minus scene then isolated detector-fixed
   stripes. An NRCB2 F200W pilot reduced held-out amplifier-row RMS from
   0.00750 to 0.00517 MJy/sr and 64-pixel sky-tile RMS from 0.00719 to
   0.00487 MJy/sr while preserving its global masked sky median. A synthetic
   smooth sky injected into an already registered target and peer scene was
   recovered with 99th-percentile absolute error 0.00000219 MJy/sr; that
   narrow test did not validate the complete calibration or registration
   chain. Subsequent pilots compared this dither method with the earlier
   four-region/horizontal corrections. PHANGS pjpipe destriping was also
   tried and dropped at the user's request. The adopted code is a custom
   registered-dither method, not a claim to reproduce PHANGS or JADES in full.
3. **v1 batch reducer.** The first master processed one dither group per
   worker, cached peer arrays, generated long-wave masks before short-wave
   masks, and ran dither 1/f → masked median sky → NMF wisp → scalar zero.
   The user chose to apply every successful 1/f fit even if held-out QC
   failed, while keeping that QC result visible. This v1 product had only
   `cal_sky_renorm/`; it did not have the current `Background2D`, `cal/`,
   per-exposure model FITS, or step figures.
4. **v2.0 restructure.** Commit `bcc1514` on `ralphy1017/stage2plus` `main`
   split CLI/orchestration (`stage2plus.py`), science steps (`steps.py`), and
   figures (`figures.py`); it added the sky-retaining `cal/`, saved models,
   per-exposure step figures, header-derived groups of at least four dithers,
   and a JADES-style `photutils.Background2D` residual sky. The default
   photutils `exclude_percentile=10` excluded all boxes after grown object
   masks in the tested deep-field frames. The adopted default is 50: a
   256×256 box needs at least half of its pixels unmasked. This is a local
   algorithm choice validated on the pilot, not a universal sky-fit rule.
5. **v2.1–v2.4 fallback and alignment iterations.** The initial backup policy
   selected a dither model only when held-out QC passed, otherwise using the
   single-frame scan. Later changes added a source-unmasked fallback, full-WCS
   reprojection, source-based alignment diagnostics, and unmasked wisp pixels.
   The 20-exposure v2.3 and v2.4 reruns completed their products, but all
   dither attempts failed the custom 0.1-pixel source-centroid gate before row
   fitting; the fallback ran for all exposures. These runs verified execution
   and FITS assembly, not dither science quality. See the [full-WCS
   mini-test](../../skysurfir/stage2plus-tests/20261103T-wcs-reprojection-mini/FULL_WCS_REPROJECTION_MINITEST_REPORT.md)
   and [wisp-unmasked rerun](../../skysurfir/stage2plus-tests/20261103T-wcs-reprojection-mini/WISP_UNMASKED_RERUN_REPORT.md).
6. **v2.5 PHANGS-inspired multi-exposure destriping (main commit `b88d7db`).** The
   0.1-pixel centroid residual and 5% held-out RMS improvement rules no longer
   reject a computable dither model. Full-WCS reprojection is refined by a
   robust source-matched x/y translation when enough matches exist, and each
   peer is scalar background-matched on common unmasked pixels. By default,
   every usable peer contributes to a sigma-clipped median scene; the
   target-minus-scene residual yields sigma-clipped row medians per amplifier.
   This adapts Williams et al. (2024), §3.6.2; it does not include their
   separate single-tile pass or second multi-tile pass, preserving the user's
   requirement to subtract one 1/f model. Optional held-out diagnostics are
   enabled with `--dither-diagnostics` and never gate application. The scan
   backup is reserved for groups with no usable peer or a dither scene/row fit
   that cannot be computed. The all-peer mini-run completed 20/20 exposures
   with all 20 dither models applied, zero fallbacks, and zero group errors.
   A separate held-out diagnostic run computed 14 dither models: 5 passed and
   9 failed the custom row/tile RMS diagnostic; 6 targets used the scan backup
   because the reduced, held-out scene could not be fit. These diagnostics are
   not gates in the default run. See the run table below; execution results do
   not establish sky/source-flux fidelity.
7. **v2.6 scan-default revision.** After
   inspecting the v2.5 figures and seeing residual 1/f structure, the user
   requested the Jumprope-style `keep_trend=True` scan as the default. The
   registered-dither path remains available only through `--onef-mode auto`.
   The repeated 20-exposure run selected `ONEFMETH=SELF` for every exposure,
   applied 20 scan models, and completed with zero group errors. The user
   reviewed its figures and requested this default be published. This run
   validates execution and file assembly, not science quality.
8. **v2.7 per-frame amplifier-surface correction.**
   The user observed residual backgrounds in the four vertical amplifier
   regions and asked to preserve smooth low-k structure shared across regions
   while correcting 2D low-k structure local to individual regions. The
   candidate fits the `keep_trend=True` scan and additive low-k components
   from one exposure, then subtracts their combined model once. A robust
   full-frame cubic surface is retained; each amplifier's quadratic surface
   minus that shared surface is added to the scan model. The global median is
   recentered exactly. Three existing production examples were tested
   read-only; figures are in
   `diagnostics/amp-scene-decomposition-20261105/`, and the user reviewed
   them and requested publishing this candidate. The 48.1%, 25.8%, and
   1.8% values reported for those frames are fractional changes in the
   256-pixel-smoothed image relative to the current products, not validation
   scores. No production FITS files were modified.

9. **v2.8 source-masked low-k fit and detector-row polynomial (local).** The
   user reported that the v2.7 amplifier fit removed real galaxy flux in
   source-rich amplifier regions. The keep-trend scan remains unmasked for
   striping. Only the amplifier-local surface and a detector-wide quadratic
   polynomial in y use the grown Photutils object mask, plus four extra
   dilation pixels. The shared cubic across the four amplifiers remains a
   protected reference; local quadratic departures and the y-only polynomial
   are added to the scan and subtracted as one 1/f model. The y polynomial is
   constant across x, so it fits horizontal structure without rotating. The
   five requested raw Stage 2 exposures completed in four groups with zero
   errors. Results and figures are in
   `diagnostics/masked-lowk-y-poly-mini-20261106/`; full mini-run products are
   under `/tmp/stage2plus-lowk-masked-mini-20261106/`. The active production
   directories were not changed. Residual amp-median ranges improved in four
   frames and were nearly unchanged in one; row-profile variation improved
   in all five. These are execution and background diagnostics, not source
   flux or diffuse-light validation.

## Current published v2.10 mode and recorded products

`steps.read_exposure` reads `SCI_ORIG` when present, otherwise `SCI`, and
identifies unusable pixels from DO_NOT_USE DQ (grown one pixel), nonfinite or
zero science, and invalid ERR. Header fields `PROGRAM`, `OBSERVTN`, `VISIT`,
`VISITGRP`, `SEQ_ID`, `ACT_ID`, `DETECTOR`, `FILTER`, and `PUPIL` define each
group; filenames are not used to decide dither membership. Input basenames
must be unique. A group needs at least one exposure. The default single-frame
scan does not require a dither peer. `--onef-mode auto` attempts dither fitting
when another group exposure is available and uses every overlapping peer.

1. **Object and fit masks.** Photutils first segments clear sources at 2.8σ
   with a 5-pixel minimum area. It masks that seed and its existing grown
   version while estimating a 512-pixel coarse `Background2D`, then segments
   the residual at 2.0σ with a 3-pixel minimum area and unions those faint
   detections with the seed. The matched long-wave mask is projected by WCS
   onto short-wave exposures when present.
   Each connected object is grown by its equivalent radius, capped at 32
   pixels. Bad pixels and grown sources are excluded from the optional auto
   dither fit. In the default scan path, detector-invalid pixels alone are
   excluded from `keep_trend=True`; sources and wisps participate in its
   profiles subject to the 90th-percentile value clip. For low-k fits only,
   grow the object mask by four more pixels and exclude it along with bad
   pixels while measuring tile and row-band medians. Evaluate the fitted
   surfaces across masked source pixels so the detector background correction
   remains continuous under objects. No bright-wisp region mask is applied.
   The sky mask still excludes invalid pixels, grown objects, and the NMF
   wisp support grown by eight pixels during later sky estimation.
2. **Optional PHANGS-inspired dither 1/f.** In `--onef-mode auto`, full-WCS
   reprojection maps each
   peer onto the target detector grid with fractional valid-coverage handling.
   A robust median x/y shift from mutual source-mask centroids refines WCS
   registration where at least three matches exist; otherwise WCS alignment is
   used. The old 0.1-pixel residual is recorded only as a diagnostic. Each
   peer's scalar background is matched to the target using a 3σ-clipped median
   of common unmasked pixels. By default all overlapping peers form a
   3σ-clipped median scene in auto mode. The target-minus-scene residual is sigma-clipped
   row by row in each of the four amplifier regions. The row-only model is
   anchored to preserve the target's masked global sky median. Any computable
   dither model is applied even if optional QC is false or was not run;
   `ONEFSTAT=REVIEW` unless optional QC passes. `--dither-diagnostics`
   holds out up to two peers and reports the prior 5% row/tile RMS test, but
   that test never selects the method. At least 10,000 valid samples on the
   stride-8 common-overlap grid are required per peer, and each amplifier must
   have measurable row offsets; peers without overlap are skipped. The single-
   scan backup is used only when no peer remains or the dither scene/row model
   cannot be computed. Both outputs use the same one selected 1/f model.
3. **Default single-frame joint model.** The scan component ports Jumprope's
   `profoundSkyScan(..., keep_trend=True)` behavior: transpose the image into
   detector scan orientation, exclude detector-invalid pixels, and do not
   apply source or bright-wisp region masks. It still clips values above the
   90th percentile, so the highest source pixels can be excluded by value. It
   measures 40th-percentile row/column profiles, scans in 512-row blocks, and
   keeps the centered profiles after subtracting a 501-pixel running-median
   trend. On the scan residual, robust 128-pixel tile medians outside the
   dilated object mask fit a cubic surface across all four amplifiers and a
   quadratic 2D surface within each amplifier. The shared surface is only a
   reference; subtract the per-amplifier local-minus-shared difference, not
   the shared sky-like component. Then fit a robust degree-2 polynomial in
   detector y to 64-row-band medians across all amplifiers. This y-only map is
   constant across x and models broad horizontal structure without rotating.
   Both low-k components are recentered on the unmasked background and added
   to the scan in one `ONEF` map. They are evaluated under source pixels but
   fitted only outside the object mask. `ONEFOBJ=False`, `LOWKOBJ=True`, and
   `LOWKDIL=4` distinguish scan masking from low-k fit masking. The method is
   local to `steps.py`, has no dither peers, and has no independent held-out
   QC; it is marked `ONEFSTAT=REVIEW`.
4. **Sky and wisps.** This sequence applies to `cal_sky_renorm/` only. The
   `cal/` branch stops after subtracting the selected 1/f and fitted wisp
   models so absolute sky measurements remain available. In the sky-subtracted
   branch, the masked median is subtracted first. A
   `Background2D` residual model uses 256-pixel boxes, a 3×3 filter,
   `BiweightLocationBackground`, `SigmaClip()`, and 50% box exclusion by
   default. On detector/filter combinations with a bundled basis, the fixed
   NMF components are fitted by ERR-weighted nonnegative least squares after
   measuring an off-wisp baseline. The last masked median is the pedestal.
5. **Outputs.** `models/<basename>.fits` contains `OBJMASK`, `LOWKMASK`,
   `SKYMASK`, `ONEF`, `SKYBKG`, and `WISP` when fitted; its primary header stores
   `SKYMED` and `PEDESTAL`. `figures/<stem>_steps.png` shows input, masks,
   models, intermediate frames, and both final products. `summary.json`
   records per-exposure values and any group errors. Both output FITS copy
   the selected non-SCI arrays and preserve the input files.

The optional held-out diagnostic requires mean amplifier-row and 64-pixel tile
RMS to improve by more than 5%, with no held-out peer worsening by more than
5%. It is not an acceptance gate. In default fallback mode, `DITHQC` is null
in JSON and omitted from FITS because no dither check was attempted. The scan
is recorded as `ONEFMETH=SELF`, `DITH1F=True`, `ONEFAPPL=True`, and
`ONEFSTAT=REVIEW`. In optional auto mode, a computable dither model is applied
even if QC fails or was not run; `ONEFAPPL` records subtraction and `DITH1F`
records that some 1/f correction was applied. If auto dither fitting cannot
be computed, a successful scan is selected as `ONEFMETH=SELF`.
`S2PVER`, `ONEFMETH`, `ONEFSTAT`, `ONEFAPPL`, `ONEFVER`, `LOWKOBJ`, `LOWKDIL`,
`DITH1F`, `DITHQC`, `DITHN`, `SKYSUB`, `SKYMED`, `BKG2D`, `WISPMETH`, and `PEDESTAL` distinguish
method and product provenance.

## What was run and verified

| Run | Code and scope | Recorded outcome |
| --- | --- | --- |
| `diagnostics/1813005_dither_stage2plus/` in the parent workspace | v1, 320 observation-005 exposures | 320 outputs; all 1/f fits applied after the 29 QC-failed targets were rerun. This predates the v2 product layout. |
| `diagnostics/1813_dither_stage2plus/` in the parent workspace | v1, all 640 locally available 1813 exposures (observations 001 and 005), 32 workers | 640/640 sky-subtracted FITS, 16 contact sheets, zero group errors, 6:49.50 wall time; 583 held-out QC passes and 57 failures (27 in observation 001, 30 in observation 005), all 1/f models applied. All SCI headers and one SCI checksum per eight-dither group were checked. |
| `/data5/rzi/jp/stage2plus-test-latest/` | v2.0, eight NRCBLONG F277W and eight NRCB4 F150W exposures, two leased CPUs | 16/16 each of `cal/`, `cal_sky_renorm/`, `models/`, and figures; zero errors; 13/16 held-out QC passes; 8 wisp fits; 2:29.60 wall time. Model subtraction closure and preserved ERR/DQ/ASDF data were verified on this test. |
| `/data5/rzi/jp/stage2plus-test-nomask/` | Experimental local change to omit grown objects from the 1/f fit mask, same 16 inputs | 16/16 outputs, zero errors, 13/16 held-out QC passes, 2:32.56 wall time. This option was uncommitted and absent from pushed `main` when this document was written. The small test does not establish a preferred mask policy. |
| `/data5/rzi/jp/stage2plus-phangs-v25-mini-20261103/` | Local v2.5 default all-peer PHANGS-style dither fit; 20 exposures, four groups, four leased CPUs | 20/20 each of `cal/`, `cal_sky_renorm/`, `models/`, and step figures; 20/20 dither models applied, zero fallbacks, zero group errors. Output FITS checksums passed; 280 non-SCI output HDUs matched input arrays; `cal/` model arithmetic closed exactly; masked final sky medians were within 7.3e-12 MJy/sr of zero. Execution validation only. |
| `/data5/rzi/jp/stage2plus-phangs-v25-diag-mini-20261103/` | Separate optional `--dither-diagnostics --no-figures` run on the same 20 exposures | 14 dither models computed (5 passed, 9 failed the custom held-out 5% RMS check); 6 scan fallbacks occurred because withholding peers left an unusable reduced scene. These results do not gate the default all-peer fit. |
| `/data5/rzi/jp/stage2plus-fallback-default-mini-20261103/` | v2.6 default scan mode, same 20 exposures, four leased CPUs | 20/20 each of `cal/`, `cal_sky_renorm/`, `models/`, and step figures; 20/20 `SELF` scan models applied; zero dither attempts and group errors. Output FITS checksums passed; 280 non-SCI HDUs matched input arrays; `cal/` model arithmetic closed exactly and masked sky medians were within 4.6e-13 MJy/sr of zero. User reviewed the figures and requested publishing the default change; execution checks do not establish science quality. |
| `diagnostics/amp-scene-decomposition-20261105/` | Local v2.7 per-frame candidate on the three named examples from `/data4/rzi/skysurfir/stage2plus/skysurfir-combined-v1/production-20261105-v1/` | Median sky preserved to printed precision on all three; inter-amplifier residual medians moved toward the shared surface. Fractional 256-pixel-smoothed image change versus current products: 48.1%, 25.8%, and 1.8%. Figures are for manual review; this is not a production run or source-flux validation. |
| `diagnostics/masked-lowk-y-poly-mini-20261106/` | Local v2.8 on the five exposures named in the 2026-11-06 user request; single-frame fallback; one worker | 5/5 completed in four groups, zero errors. Comparing old and new `cal/` with the same v2.8 fit mask, the centered four-amplifier median range improved in four frames (0.00425→0.00358, 0.00375→0.00387, 0.00051→0.00004, 0.08351→0.05379, 0.00149→0.00058 MJy/sr). The robust scatter of 64-row background medians decreased in all five (0.00173→0.00139, 0.00453→0.00339, 0.00027→0.00025, 0.25326→0.14003, 0.00248→0.00224 MJy/sr). These are residual-background diagnostics, not source-flux or diffuse-light validation. FITS outputs are in `/tmp/stage2plus-lowk-masked-mini-20261106/`; production outputs were untouched. |
| `diagnostics/flat-background-mini-20261106/` | Local v2.9 on the same five examples; single-frame fallback; one worker | 5/5 completed in four groups, zero errors. `cal/` subtracts only 1/f and wisp; its masked median exactly matches the post-1/f sky median. `cal_sky_renorm/` additionally subtracts median, `Background2D`, and pedestal; its masked median is zero. Both branches closed exactly against saved model arrays on valid pixels. Final pairs use the same absolute display stretch. Outputs are in `/tmp/stage2plus-flat-background-mini-20261106/`; production outputs were untouched. |
| `diagnostics/object-mask-v210-five-20261106/` | Local v2.10 on the same five examples; revised faint-outskirt segmentation; original 1/f-then-wisp order; one worker | 5/5 completed in four groups, zero errors. The 2.8σ/5-pixel seed is extended with 2.0σ/3-pixel residual detections after a seed-masked 512-pixel coarse background. `cal/` subtracts only 1/f+wisp and retains its measured masked sky median; `cal_sky_renorm/` additionally subtracts masked median, `Background2D`, wisp, and pedestal, yielding masked median zero. Both branches close exactly against saved model arrays on valid pixels. The paired final-product sheet uses the same absolute stretch within each exposure pair. Outputs are in `/tmp/stage2plus-objectmask-v210-five-20261106/`; production outputs were untouched. This execution and visual check is not a diffuse-source injection or broad science validation. The wisp-before-1/f experiment was rejected. |

These paths are records on the development machine, **not** data shipped in
the Git repository. The v1 input manifest is
`../diagnostics/1813_dither_stage2plus/inputs.txt`, pointing to the calibrated
source directory
`/data2/rzi/skysurfir/w200-v1-stage2-pid01813-test-v1/nircam/stage2/cal/`.
The v1 640-exposure result is not a v2 full-database deployment. The v2.0
verification covered 16 exposures. The v2.1 fallback execution pilot covered
352 exposures but all dither attempts failed at the 0.1-pixel shift-only WCS
check. The local v2.3 mini-test covered 20 exposures with full-WCS
reprojection; all dither candidates then failed the 0.1-pixel post-reprojection
source-centroid gate, so stripe fitting and held-out QC were not reached. The
v2.2 source-unmasked fallback did execute for those 20 exposures, but its
science quality remains unvalidated. The v2.4 rerun used the same 20-frame
manifest with no wisp exclusion in the 1/f fit masks and completed 20/20
products and figures with zero group errors. All dither candidates again
failed at source alignment; the scan fallback was selected for all 20. This
validates execution and mask provenance only. The local v2.5 all-peer run
then selected the dither method for all 20 exposures, with no fallbacks.
Optional held-out diagnostics selected 6 fallbacks only because the test
intentionally reduced the peer scene; 5 of the 14 computable dither fits
passed the custom diagnostic and 9 failed it. These results validate
execution, not science quality.
The local v2.6 rerun then used the scan as the default on all 20 exposures,
with zero dither attempts, zero group errors, and 20 `SELF` models applied.

## Practical limits and open scientific checks

- `--plan-only` checks group sizes and reports the number of groups with an
  available wisp basis and the number of short-wave exposures matched to a
  long-wave mask. The current v2 code does not require a basis for every
  wisp-sensitive input; absent bases lead to `WISPMETH=NONE`. Review coverage
  before a full database run.
- The v2.8 masked amplifier-surface and y-polynomial candidate was run only on
  the five specified frames. Object masks reduce contamination, but any
  unmasked diffuse galaxy or cluster light can still bias a low-k fit, and
  applying the fitted maps under source masks can affect extended-source
  flux. Review all five figures and run source/diffuse-light injections
  before using it on a full database. The existing production products were
  not changed.
- In optional auto mode, a dither scene needs enough common sky after full-WCS peer reprojection.
  Sparse overlap or heavily masked amplifier rows can still make the fit
  impossible and invoke the scan backup. Source-match residuals are recorded
  as diagnostics, not a 0.1-pixel gate. The all-peer run proves the current
  representative groups can be fit; held-out tests and source-flux/diffuse-sky
  injections are still needed to establish scientific quality.
- A `cal/` image keeps the global masked sky level by construction, but
  detector-correlated diffuse structure may still be absorbed by the 1/f
  model. Absolute sky surface brightness also depends on upstream dark,
  flat, and additive offsets. End-to-end injections into the original dither
  exposures and cross-detector sky comparisons are still needed before
  claiming local sky-surface-brightness fidelity.
- The v2.0 verification covered 16 exposures; v2.5 and v2.6 each covered 20.
  None covers the earlier 640-exposure v1 campaign or a
  100k-frame archive. At large scale, use a
  validated input manifest, a storage estimate, CPU and memory limits,
  isolated output directories per configuration, and staged QA. Each worker
  holds a dither group; the verification run peaked near 1.4 GB RSS for the
  measured process, and memory grows with concurrent groups. There is no
  built-in resume or cross-run product compatibility check.
- The `--no-onef-object-mask` flag and its summary setting are included in
  main v2.5, though they were local changes in earlier revisions. A clean
  clone of the pushed v2.0 baseline does not expose that option. Do not assume
  a run with this flag is the default or that `S2PVER` uniquely identifies the
  code revision; record the Git commit and working-tree state.

The operational checklist and machine-specific CPU-lease convention live in
[AGENTS.md](AGENTS.md). The README remains the user-facing run guide.
