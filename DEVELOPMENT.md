# Development brief: default 1/f correction for science imaging

**Audience:** future agents developing the standalone Stage2Plus reducer and
its SKYSURFIR integration.

**Status:** v2.6 defaults to the single-frame Jumprope-style
`keep_trend=True` scan; the PHANGS-inspired all-peer dither remains opt-in via
`--onef-mode auto`. The 20-exposure v2.6 rerun applied the scan to all 20
exposures with zero errors. The user reviewed its step figures and requested
publishing this default. The run verifies execution and FITS assembly, not
sky/source-flux fidelity or production readiness.

## User goal

The user's current instruction makes the single-frame `keep_trend=True` scan
the default for supported science imaging. The registered-dither method is
available only as an opt-in `--onef-mode auto`. Keep the sky-retaining `cal/`
and sky-subtracted `cal_sky_renorm/` branches consistent with the selected
correction model.

The long-term scope is all selected science imaging, including groups with
fewer than four dithers. Instrument-specific adapters and validation are
required: current Stage2Plus supports NIRCam imaging, and its NIRCam masks,
templates, and thresholds must not be applied blindly to NIRISS or MIRI. The
the current default is the single-frame scan following the user's review of
the v2.6 figures. Unsupported or failed exposures must receive an explicit
hold or failure status, never silently pass through as corrected.

This brief does not authorize a production database run. It describes the
future behavior the user wants and the evidence needed before enabling it.

## Current implementation and pilot evidence

The Stage2Plus code at this document's creation was v2.0, commit `c5dcd92`,
under `stage2plus.py` and `steps.py`. It groups exposures by FITS program,
observation, visit, visit group, sequence, activity, detector, filter, and
pupil. That baseline required at least four dithers. Its dither method aligns
peer images using a translation-only check, forms a median scene, fits
amplifier-row and detector-column residuals, and holds out up to two peers
for QC. A fit exception fails its group. The correction is still written
when `DITHQC=False`.

The user's 2026-11-03 pilot covered 352 NIRCam Stage 2 exposures from visits
1176231, 2738005, and 1727001 (78 groups). The detailed report, manifests,
logs, outputs, figures, and partial-product inventory are in:

[`../../skysurfir/stage2plus-tests/20261103T181509Z-visits-1176231-2738005-1727001/STAGE2PLUS_VISIT_PILOT_REPORT.md`](../../skysurfir/stage2plus-tests/20261103T181509Z-visits-1176231-2738005-1727001/STAGE2PLUS_VISIT_PILOT_REPORT.md)

The baseline run used an isolated copy with full-WCS peer reprojection and
group-level error isolation. It kept Stage2Plus's science estimators and
default object mask. Results:

| Visit ID | Selected | Complete in whole groups | Group errors |
|---|---:|---:|---:|
| 1176231 | 152 | 88 | 16 |
| 2738005 | 120 | 18 | 17 |
| 1727001 | 80 | 0 | 20 |
| **Total** | **352** | **106** | **53 of 78** |

After that pilot, the local v2.1 implementation allows groups with one or
more exposures. Groups of four or more first try the registered-dither model;
groups below four skip to the backup. A dither fit exception or failed
held-out QC discards that candidate and selects one single-frame
`profoundSkyScan`-style model with `keep_trend=True`. The backup uses the
Jumprope settings documented below, the same object/DQ/wisp fit mask, and is
subtracted once in both product branches. It is recorded as `ONEFMETH=SELF`,
`ONEFAPPL=True`, and `ONEFSTAT=REVIEW`; a measurable failed dither QC remains
`DITHQC=False`, while an unavailable check is null in JSON and omitted from
FITS. The v2.1 pilot ran the fallback on all 352 selected exposures with zero
group errors. However, the dither attempts all failed before held-out QC
because that run's `_aligned()` implementation rejected peer WCS differences
that were not a pure shift. The run therefore tested fallback execution and
FITS assembly, not a successful dither model or science quality.
Representative checksums, SCI arithmetic, ERR/DQ arrays, and WCS passed. See
the [v2.1 mini-test report](../../skysurfir/stage2plus-tests/20261103T-v21-auto-fallback/STAGE2PLUS_V21_MINITEST_REPORT.md)
for counts, logs, outputs, and figures.

The v2.1 pilot fallback used the source mask in the 1/f fit. The local v2.2
change passes only detector-invalid and protected bright-wisp masks to the
fallback, so detected/grown sources participate in its profiles. The
90th-percentile high-value clip remains and can still exclude the brightest
source pixels by value. The dither fit remains source-masked by default. The
v2.3 mini-test exercised this fallback on 20 exposures, but this execution
test does not validate its source-flux or diffuse-sky performance.

Of the 106 exposures recorded in `summary.json`, 52 passed and 54 failed the
held-out 1/f QC. There were 40 `Too many masked rows in an amplifier` errors,
11 `Too few sky pixels to fit 1/f` errors, and 2 insufficient dither-overlap
errors. Some failed groups left partial FITS and figure files; those are
listed in `partial-products.json` and are not valid Image3 inputs.

The earlier full-WCS adaptation got past the upstream translation-only WCS
rejection for many groups, but it did not solve missing overlap or low-sky
failures.
A follow-up on 24 exposures tested `--no-onef-object-mask`, one held-out
dither, and a test-only 85% valid-row threshold. It completed 10 exposures,
recorded three group errors, and still had no complete group for 1727001 or
short-wave WISP result for 2738005. **Do not treat an 85% row threshold or
disabling the object mask as an established fix.**

The local v2.3 full-WCS change and its exact mini-test outcomes are documented
in the [full-WCS reprojection mini-test report](../../skysurfir/stage2plus-tests/20261103T-wcs-reprojection-mini/FULL_WCS_REPROJECTION_MINITEST_REPORT.md).
It replaced the pure-shift check with target-grid reprojection and source
centroid matching after WCS mapping. All 20 targets passed the common-overlap
check but failed the retained maximum 0.1-pixel source-centroid tolerance
(median matched residuals 0.63–1.38 px). Thus 6σ fit clipping, amplifier-row
coverage, and held-out RMS gates were not reached. The mini-test completed all
fallback products with zero group errors; it is not a successful dither
validation. Source-centroid scatter may include segmentation variation and
needs confirmation against bright isolated stars or independent image
registration before changing the tolerance.

The v2.4 change removes the bright-wisp region mask from both dither and
single-frame 1/f fitting. Dither fits still exclude detector-invalid and
grown-source pixels by default. The single-frame fallback excludes only
detector-invalid pixels by mask and includes sources and wisps, while its
existing upper-10% value clip remains. The same 20 exposures were rerun with
figures enabled; all 20 products completed, all 20 selected the scan fallback,
and the dither attempts failed at the unchanged source-alignment gate. The
full details and figures are in the [wisp-unmasked rerun report](../../skysurfir/stage2plus-tests/20261103T-wcs-reprojection-mini/WISP_UNMASKED_RERUN_REPORT.md).

The v2.5 implementation removes the 0.1-pixel source-residual rejection and
uses mutual source matches to refine a per-peer translation instead. It uses
all usable peers in the default scene, matches peer backgrounds with a
sigma-clipped overlap median, and fits sigma-clipped amplifier-row medians
from the target-minus-scene residual. Optional held-out row/tile diagnostics
run only with `--dither-diagnostics`; their results never select the scan
backup. A dither fit error or no usable peer still invokes the backup. This is
the multi-exposure component of Williams et al. (2024), §3.6.2; the separate
single-exposure pass and extra PHANGS large-scale pass are omitted to retain
one 1/f model in this workflow.

The default v2.5 mini-run used 20 exposures in four groups (two 4-dither and
two 6-dither groups; both modules and short/long wavelength), with four leased
CPUs. All 20 selected `ONEFMETH=DITHER`; none selected the scan backup. The run
wrote 20 files each in `cal/`, `cal_sky_renorm/`, `models/`, and `figures/`,
with zero group errors. Output checksums passed; non-SCI arrays were preserved
in 280 checked output HDUs; `cal/` arithmetic closed exactly against the saved
`ONEF`, `WISP`, and pedestal models; masked final sky medians were within
7.3e-12 MJy/sr of zero. Products and step figures are under
`/data5/rzi/jp/stage2plus-phangs-v25-mini-20261103/`.

A separate 20-exposure `--dither-diagnostics --no-figures` run used up to two
held-out peers and therefore fit the remaining reduced peer set. It computed
14 dither models (5/14 passed and 9/14 failed the custom 5% row/tile RMS
diagnostic) and selected the scan backup for 6 targets whose reduced scene
could not be fit. Those outcomes did not gate any model in the default run;
they show why this optional diagnostic mode must not be confused with the
all-peer default. Its output is at
`/data5/rzi/jp/stage2plus-phangs-v25-diag-mini-20261103/`. Both are execution
pilots only; visual review and sky/source-flux validation remain necessary.

After the user inspected the v2.5 figures and reported residual 1/f structure,
the v2.6 change made the scan method the CLI default. The 20-exposure
rerun is `/data5/rzi/jp/stage2plus-fallback-default-mini-20261103/`: 20/20
`ONEFMETH=SELF`, 20 applied 1/f models, zero dither attempts, zero group
errors, and 20 figures. Output FITS checksums passed; 280 non-SCI HDUs matched
the input arrays, `cal/` arithmetic closed exactly against the saved models,
and masked final sky medians were within 4.6e-13 MJy/sr of zero. The user
reviewed the figures and requested publication of the scan default.

The plan found matching NMF WISP bases in 31 of 78 groups; 47 groups (212
exposures) had no basis. This is a separate correction-coverage gap from 1/f.

The representative output FITS passed arithmetic closure, checksum, WCS, DQ,
ERR, ASDF, and retained-HDU checks. That validates file assembly for those
examples; it does not establish science quality across the database. The
pilot did not run Image3 or perform a complete source-flux injection study.

## Current default mode and optional dither mode

The current default is the single-frame scan. The dither workflow described
below is opt-in through `--onef-mode auto`; the user selected the scan after
reviewing residual striping in the dither-default figures.

For every exposure in the frozen science manifest, record exactly one
per-exposure outcome. Do not require a complete four-dither group to attempt
1/f correction.

1. **Optional registered-dither mode.** When `--onef-mode auto` is requested,
   use every usable peer from
   the same calibration release and compatible visit, detector, filter, and
   pupil. Register with the complete WCS and keep a per-pixel coverage map.
   Build the scene only from peers with valid local coverage. One usable peer
   is sufficient to attempt a fit; if none overlap, use the backup.
2. **Record diagnostics without making them prerequisites.** Record peer
   coverage by amplifier and row, WCS/source-match residuals, valid fit
   fraction, model amplitude, and optional held-out row/tile RMS. A missing
   peer footprint is not a measured zero. Do not reject a computable all-peer
   model based on the historical 0.1-pixel or 5% criteria. The optional
   `--dither-diagnostics` mode withholds peers and changes the scene being fit;
   use it to measure transfer to held-out peers, not as the default.
3. **Default per-exposure scan.** The v2.5 implementation uses an
   invalid-pixel-masked, source- and wisp-unmasked scan correction based on the
   NEP/JumproPe `scan_1f(..., keep_trend=True)` behavior referenced by
   [`../stage2-plus.py`](../stage2-plus.py). The implementation is pinned in
   `steps.py`; it does not rely on an external JumproPe runtime. It is selected
   by default. In `auto` mode it is selected when no usable peer exists or
   dither fitting raises. A failed optional held-out diagnostic does not
   discard a computed dither model. Exactly one selected model is subtracted in both output
   branches. Validate source flux and diffuse-sky recovery before treating
   this fallback as production-ready.
4. **Use a hybrid model only when its parts are explicit.** If the dither fit
   is supported in some amplifier rows and not others, a per-row or
   per-amplifier fallback may be evaluated. Save a method map and provenance
   for each region. Never add a dither model and a fallback model together
   unless that combined estimator has its own derivation and validation.
5. **Keep unsuccessful cases visible.** If neither method meets its tested
   acceptance criteria, label the exposure `HOLD` or `FAILED` and exclude it
   from automatic Image3 input selection. Do not write an unchanged SCI array
   labeled as corrected. The campaign ledger must still account for the
   exposure.

The current default applies the scan to every exposure. Auto dither remains an
explicit opt-in. Do not switch the default back until the user has inspected
the figures and requests that change. Keep correction method, fit success, and
QC outcome separate.

Suggested provenance fields:

| Field | Meaning |
|---|---|
| `ONEFMETH` | `DITHER`, `SELF`, `HYBRID`, or `NONE` |
| `ONEFSTAT` | `PASS`, `REVIEW`, `HOLD`, or `FAILED` for the selected method |
| `ONEFAPPL` | Whether the selected 1/f model was subtracted |
| `DITHQC` | Optional dither-only held-out result; null in JSON and omitted from FITS when diagnostics were not run or unavailable |
| `DITHN` | Peer count used in the selected dither scene |
| `ONEFCOV` | Valid local coverage summary; detailed maps belong in QA/model products |
| `ONEFVER` | Algorithm version or source revision |

Keep `DITH1F` as the 1/f correction-applied flag; `ONEFMETH` identifies the
selected method. `DITHQC=False` means an applied dither model failed the
optional diagnostic; it never means the fallback was selected. Update
`summary.json` with attempted methods,
fallback reason, thresholds, masks, per-amplifier coverage, before/after
residuals, and the final selected method for every exposure.

## Threshold development

For the optional auto dither mode, the historical 0.1-pixel and held-out 5%
tests are custom diagnostics, not PHANGS acceptance requirements. Treat any
additional feasibility thresholds as measured policy, not a way to make the
error count disappear.

- Parameterize minimum valid fit fraction, peer overlap, valid rows per
amplifier, allowable contiguous row gaps, WCS registration residual, minimum
number of scene peers, and held-out QC tolerances.
- Sweep candidate settings on a stratified sample, including strict v2
settings and the tested 85% row setting. Report product coverage and each
science metric at every setting. Do not select a threshold based only on the
number of groups that finish or visual smoothness.
- Separate dither-fit feasibility from dither-model quality. A fit can be
computable and still leave source-aligned residuals or worsen held-out
striping.
- Evaluate row RMS and tile RMS per held-out peer, per amplifier, and over
common valid coverage. Include source flux, diffuse sky, and overlap
background metrics. Calibrate any acceptance tolerance from these tests.
- Keep a `REVIEW` state for valid models whose QC is borderline. `DITHQC=False`
  must never silently become an unqualified production pass.

## Validation gates before making `auto` the database default

Run small, isolated tests first; leave existing NEP notebooks and live
reductions unchanged.

1. **Synthetic/unit tests:** inject known row/column 1/f patterns, smooth sky,
   point sources, extended sources, DQ defects, NaN edges, partial amplifier
   coverage, and wisp masks. Measure model recovery, sky-level bias, and
   source-flux change.
2. **Method comparison:** compare dither, backup, and any hybrid model on the
   same Stage 2 inputs. Use fixed masks and record why the selected method
   changes. Verify the correction is not absorbing celestial sources or
   large-scale sky.
3. **Stratified JWST pilot:** include the three visits above; 4-, 6-, and
   larger-dither groups; both NIRCam modules; WISP-sensitive and non-WISP
   detectors; groups without matched long-wave masks; sparse and crowded
   fields; and frames with low unmasked-sky fraction. Include NIRISS only
   after a NIRISS-specific adapter and tests exist.
4. **Overlap validation:** compare corrected exposure backgrounds and residual
   stripes in overlaps. Check that mixed dither/fallback methods do not
   introduce visit-dependent zero-point steps.
5. **Publication audit:** require one terminal ledger state per selected
   exposure, atomic per-exposure products, and no untracked partial files.
   Image3 manifests may include only products whose selected method and
   quality state meet the declared policy.
6. **Resource test:** benchmark workers by dither group with single-threaded
   BLAS/OpenMP. Measure peak RSS, read/write throughput, CPU utilization, and
   time by method. Preserve the shared CPU lease and bounded I/O policy.

Do not claim the default is ready until the pilot has no unexplained missing
exposures, every fallback has flux/sky validation, all method and QC states
are represented in the ledger, and partial outputs cannot be mistaken for
completed products.

## Implementation order

1. **Implemented in v2.5:** add full-WCS peer SCI reprojection,
   conservative valid/exclusion coverage mapping, robust source-based
   translation refinement, overlap background matching, sigma-clipped
   all-peer scene construction, and sigma-clipped per-amplifier row fits. The
   20-exposure default mini-run selected dither for every target with zero
   group errors; next validate residual striping, source flux, and diffuse-sky
   preservation before broad use.
2. Refactor dither fitting to return structured feasibility and QC results,
   not just a model or exception. Preserve diagnostics for rows and
   amplifiers without peer coverage.
3. **Implemented locally:** the NIRCam code selects dither or scan
   automatically and records method/status provenance. Add a ledger-level
   hold state before connecting it to automatic Image3 selection.
4. Make each exposure write transactional. Reconcile group completion and
   prevent partial files from entering downstream manifests.
5. Run the validation gates above, update `README.md`, `CONTEXT.md`, and this
   brief with actual evidence, then enable `auto` in the science database
   notebooks by default.

## Related records

- Current science method and history: [`CONTEXT.md`](CONTEXT.md)
- Operational rules for future agents: [`AGENTS.md`](AGENTS.md)
- Pilot report and results: [`SKYSURFIR visit pilot`](../../skysurfir/stage2plus-tests/20261103T181509Z-visits-1176231-2738005-1727001/STAGE2PLUS_VISIT_PILOT_REPORT.md)
- v2.1 fallback mini-test: [`Stage2Plus v2.1 mini-test`](../../skysurfir/stage2plus-tests/20261103T-v21-auto-fallback/STAGE2PLUS_V21_MINITEST_REPORT.md)
- v2.4 wisp-unmasked rerun: [`Stage2Plus wisp-unmasked mini-test`](../../skysurfir/stage2plus-tests/20261103T-wcs-reprojection-mini/WISP_UNMASKED_RERUN_REPORT.md)
- The v2.6 scan-default change follows the user's review of the 20-exposure
  figures. The mini-runs establish execution only, not scientific validation.
- Parent JumProPe history: [`../handoff.md`](../handoff.md)
