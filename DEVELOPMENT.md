# Development brief: default 1/f correction for science imaging

**Audience:** future agents developing the standalone Stage2Plus reducer and
its SKYSURFIR integration.

**Status:** the NIRCam single-frame fallback described below is implemented
in the local v2.1 worktree (2026-11-03), but it has not yet had a science
pilot. The broader instrument, threshold, and database validation goals in
this brief remain open. Do not infer full-database reliability from the code
change alone.

## User goal

The user wants 1/f correction attempted by default for every science-grade JWST
imaging exposure in the SKYSURFIR database. Prefer a registered-dither fit
when the data support it. Automatically try a validated per-exposure fallback
when the dither fit cannot be made or fails its quality checks. Keep the
sky-retaining `cal/` and sky-subtracted `cal_sky_renorm/` branches consistent
with the selected correction model.

The long-term scope is all selected science imaging, including groups with
fewer than four dithers. Instrument-specific adapters and validation are
required: current Stage2Plus supports NIRCam imaging, and its NIRCam masks,
templates, and thresholds must not be applied blindly to NIRISS or MIRI. The
default policy should eventually be `ONEF_MODE=auto` for each supported
instrument; unsupported or failed exposures must receive an explicit hold or
failure status, never silently pass through as successfully corrected.

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
FITS. This change was not run on the 352-exposure pilot and has not been
independently validated.

Of the 106 exposures recorded in `summary.json`, 52 passed and 54 failed the
held-out 1/f QC. There were 40 `Too many masked rows in an amplifier` errors,
11 `Too few sky pixels to fit 1/f` errors, and 2 insufficient dither-overlap
errors. Some failed groups left partial FITS and figure files; those are
listed in `partial-products.json` and are not valid Image3 inputs.

The full-WCS adaptation got past the upstream translation-only WCS rejection
for many groups, but it did not solve missing overlap or low-sky failures.
A follow-up on 24 exposures tested `--no-onef-object-mask`, one held-out
dither, and a test-only 85% valid-row threshold. It completed 10 exposures,
recorded three group errors, and still had no complete group for 1727001 or
short-wave WISP result for 2738005. **Do not treat an 85% row threshold or
disabling the object mask as an established fix.**

The plan found matching NMF WISP bases in 31 of 78 groups; 47 groups (212
exposures) had no basis. This is a separate correction-coverage gap from 1/f.

The representative output FITS passed arithmetic closure, checksum, WCS, DQ,
ERR, ASDF, and retained-HDU checks. That validates file assembly for those
examples; it does not establish science quality across the database. The
pilot did not run Image3 or perform a complete source-flux injection study.

## Required default `auto` behavior

For every exposure in the frozen science manifest, record exactly one
per-exposure outcome. Do not require a complete four-dither group to attempt
1/f correction.

1. **Try registered-dither 1/f when eligible.** Use peers from the same
   calibration release and compatible visit, detector, filter, and pupil.
   Register with the complete WCS and keep a per-pixel/per-row coverage map.
   Build the scene only from peers with valid local coverage. Keep one or more
   independent peers for held-out validation when the group size permits.
2. **Assess the dither model before selecting it.** Record peer coverage by
   amplifier and row, internal and edge gaps, WCS residuals, valid fit
   fraction, model amplitude, held-out row RMS, and tile RMS. A missing peer
   footprint is not a measured zero. Do not extrapolate across large
   unobserved detector regions without a validated maximum-gap rule.
3. **Per-exposure backup.** The local v2.1 worktree now implements an
   object- and DQ-masked single-exposure scan correction based on the
   NEP/JumproPe `scan_1f(..., keep_trend=True)` behavior referenced by
   [`../stage2-plus.py`](../stage2-plus.py). The implementation is pinned in
   `steps.py`; it does not rely on an external JumproPe runtime. It is selected
   when the group is too small for dither QC, dither fitting raises, or
   dither QC fails. Exactly one selected model is subtracted in both output
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

The desired default attempts the dither method and then the backup without a
per-run opt-in for supported science imaging. A fallback is a distinct method,
not a relaxed dither result. Keep correction method, fit success, and QC
outcome separate.

Suggested provenance fields:

| Field | Meaning |
|---|---|
| `ONEFMETH` | `DITHER`, `SELF`, `HYBRID`, or `NONE` |
| `ONEFSTAT` | `PASS`, `REVIEW`, `HOLD`, or `FAILED` for the selected method |
| `ONEFAPPL` | Whether the selected 1/f model was subtracted |
| `DITHQC` | Dither-only held-out result; null in JSON and omitted from FITS when the check was unavailable |
| `DITHN` | Peer count used in the selected dither scene |
| `ONEFCOV` | Valid local coverage summary; detailed maps belong in QA/model products |
| `ONEFVER` | Algorithm version or source revision |

Keep `DITH1F` as the existing dither-method indicator. Do not encode a
per-exposure fallback as `DITHQC=False`; that means a dither QC failure, not
“dither method not selected.” Update `summary.json` with attempted methods,
fallback reason, thresholds, masks, per-amplifier coverage, before/after
residuals, and the final selected method for every exposure.

## Threshold development

The user wants dithered 1/f to work on all exposures and is open to more
relaxed thresholds. Treat thresholds as a measured policy, not a way to make
the error count disappear.

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

1. Add full-WCS alignment and coverage diagnostics in an isolated branch; test
   alignment residuals and flux/noise behavior before broad use.
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
- Parent JumProPe history: [`../handoff.md`](../handoff.md)
