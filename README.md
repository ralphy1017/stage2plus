# stage2plus

Exposure-level cleanup of JWST NIRCam Stage 2 `cal` files before drizzling.
It works on any NIRCam imaging program: exposures are grouped from their FITS
headers into dither groups (same program, observation, visit, visit group,
sequence, activity, detector, filter and pupil). The default 1/f method is a
per-frame unmasked `keep_trend=True` scan plus source-masked amplifier-local
2D and detector-wide horizontal low-k fits. Use `--onef-mode auto` to opt
into dither fitting, with the per-frame method as backup if a dither model
cannot be computed.

## What it does

For every exposure, in this order:

1. **Input**: `SCI_ORIG` if the file has one, otherwise `SCI`.
2. **Object mask**: Photutils first detects clear sources at 2.8σ with a
   5-pixel minimum area. It masks those sources while estimating a coarse
   512×512-pixel `Background2D`, then segments its residual at 2.0σ with a
   3-pixel minimum area and unions those faint detections with the first
   mask. This captures more diffuse galaxy outskirts; connected source masks
   are still grown by the existing equivalent-radius rule (capped at 32 px).
   Short-wave frames add the long-wave mask of the same exposure, reprojected
   through the WCS. On wisp detectors, the bright-wisp region is not masked
   from either 1/f fit.
3. **1/f correction (default)**: use the single-frame Jumprope-style
   `profoundSkyScan(keep_trend=True)` scan. It excludes DQ/invalid pixels,
   includes source and wisp pixels in the scan profiles, and clips values
   above the 90th percentile. Then fit low-k components using robust
   128-pixel tile medians outside the grown Photutils object mask plus four
   extra dilation pixels. Keep a shared detector-wide cubic as the reference;
   each amplifier's quadratic 2D fit minus that reference models
   amplifier-local structure. A degree-2 polynomial in detector y, constant
   across x, models broad horizontal structure without rotating to follow
   sources. Evaluate both low-k maps across masked sources, add them to the
   unmasked scan, and subtract the combined model once. `LOWKMASK`, `LOWKOBJ`,
   and `LOWKDIL` record the mask and its dilation. The method is recorded as
   `ONEFMETH=SELF` and `ONEFSTAT=REVIEW`; it has no independent held-out
   validation. The latest five-frame mini-test figures are in
   `diagnostics/object-mask-v210-five-20261106/`. They show the revised
   faint-outskirt segmentation with the standard 1/f-then-wisp order;
   source-flux and diffuse-light injection tests are still needed before
   broad production use.
   `--onef-mode auto` opts into the registered-dither method: use every
   overlapping peer, reproject
   through full WCS, refine the relative x/y shift from mutual source matches
   when enough matches exist, and match each peer's background to the target
   over common unmasked pixels. A sigma-clipped median of the registered peers
   estimates the shared sky and sources. Subtract that scene temporarily, then
   measure sigma-clipped residual row medians separately in each of the four
   amplifiers. The resulting row-only model is anchored to preserve the
   target's masked global median and is subtracted once. This follows the
   multi-tile method in [Williams et al. (2024), §3.6.2](https://arxiv.org/html/2401.15142v3).
   The paper's single-tile pass is not included; this reducer applies one 1/f
   model. The former 0.1-pixel alignment and 5% held-out improvement checks
   do not reject a computable model. Use `--dither-diagnostics` with
   `--onef-mode auto` to run the optional held-out RMS check; its result is
   recorded but never gates the model. If no peer overlaps or the row model
   cannot be computed, auto mode uses the same single-frame scan.
   `--no-onef-object-mask` controls source masking for the dither fit only.
4. **Masked median sky**: the median of pixels outside the grown objects, bad
   pixels and (on wisp detectors) the wisp template support.
5. **Residual sky structure**: a JADES-style `photutils.Background2D` with
   256×256-pixel boxes, a 3×3 filter, `BiweightLocationBackground` and
   `SigmaClip()`, using the step 4 mask. It is interpolated to full
   resolution and subtracted once; no polynomial is fitted.
6. **Wisp**: an error-weighted non-negative fit of the fixed NMF basis in
   `templates/`, wherever a basis exists for the detector and filter.
7. **Pedestal**: the masked median left after the wisp step is removed.

| Product | SCI |
| --- | --- |
| `cal_sky_renorm/` | input − 1/f − median sky − Background2D − wisp − pedestal |
| `cal/` | input − 1/f − wisp (original sky/background retained) |

The `cal/` product is for measuring the true sky in object-free regions: it
subtracts only the fitted 1/f model (including its masked low-k detector
components) and the NMF wisp model. It does not subtract the median sky,
`Background2D`, or the final pedestal. `cal_sky_renorm/` additionally removes
those sky terms to provide a near-zero background for sky-subtracted mosaics.
Never measure absolute sky brightness from `cal_sky_renorm/`.

Both products keep the input's PRIMARY, ERR, DQ, AREA, VAR_* and ASDF
extensions byte for byte; only SCI is replaced. Inputs are never modified.

### 1/f diagnostics and fallback

The default `--onef-mode fallback` fits a source-unmasked single-frame scan,
then object-masked amplifier-local and detector-wide y-only surfaces per
exposure. It subtracts the combined model once.
`--onef-mode auto` opts into dither fitting and uses that per-frame method
when there is no usable peer or the dither scene/row fit cannot be computed.
Unexpected processing errors fail the group.

Full-WCS reprojection carries valid-pixel coverage and keeps incomplete
interpolation support out of the scene. Source matches estimate a relative
translation; the residuals, including whether they meet the former 0.1-pixel
criterion, are diagnostics only. In `auto` mode every overlapping peer
contributes to the sigma-clipped scene. `--dither-diagnostics` instead holds out up to two
peers and records the prior row/tile RMS comparison, including its 5%
thresholds, without using that outcome to select the model. With diagnostics
off, `DITHQC` is unavailable and omitted from FITS. Because diagnostics
withhold peers from the scene, that mode changes the fit inputs; use it only
for validation, not as the default all-peer processing configuration. A dither
model that is computed is applied even when optional QC fails or was not run;
it is marked `REVIEW` unless optional QC passes. Only expected coverage or fit
failures select the scan backup; unexpected processing errors fail the group.
In `auto` mode, the scan backup is selected when the group has no usable peer
or the dither scene/row fit cannot be computed. If the selected scan fails, or
both fits fail in `auto` mode, the group reports an error rather than writing
an uncorrected exposure.

### Background2D box exclusion

photutils' default `exclude_percentile=10` drops any box with more than 10%
of its pixels masked. After masking grown objects, typical deep-field frames
have no such box, so the default fails. stage2plus uses **50**: a box needs at
least half its pixels unmasked. Boxes inside the masked wisp region are still
excluded and filled by interpolation. Set `--bkg-exclude-percentile` to change
this.

## Install

Requires Python 3.11 or later. Tested with numpy 2.2–2.5, astropy 7.2–8.0
and photutils 2.3–3.0.

```bash
git clone https://github.com/ralphy1017/stage2plus.git
python3 -m venv stage2plus-venv
stage2plus-venv/bin/python -m pip install -r stage2plus/requirements.txt
```

## Run

Write one calibrated `*_cal.fits` path per line, with unique basenames, from
a single calibration version. Include the long-wave exposures with their
short-wave partners so the long-wave masks can be used. Choose an output
directory that is not an input directory.

```bash
CAL=/path/to/stage2/cal
RUN=/path/to/output
mkdir -p "$RUN"
find "$CAL" -maxdepth 1 -name '*_nrc*_cal.fits' | sort > "$RUN/inputs.txt"
stage2plus-venv/bin/python stage2plus/stage2plus.py \
  --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 8 --plan-only
stage2plus-venv/bin/python stage2plus/stage2plus.py \
  --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 8 > "$RUN/run.log" 2>&1
```

`--plan-only` checks the groups, counts wisp-template groups and long-wave
mask matches, and estimates output size. Test a few groups first, for example
by filtering the list with `grep`.

| Option | Meaning |
| --- | --- |
| `--input-list PATH` | Required list of input `cal` FITS paths. |
| `--output-dir PATH` | Required output root. |
| `--workers N` | Dither groups processed in parallel, one CPU each (default 8). |
| `--template-dir PATH` | NMF wisp bases (default `templates/`). |
| `--bkg-exclude-percentile P` | `Background2D` `exclude_percentile` (default 50). |
| `--no-onef-object-mask` | Do not mask objects in the registered-dither fit (later steps still do). |
| `--onef-mode MODE` | `fallback` (default) uses the single-frame scan; `auto` tries dither then scan. |
| `--dither-diagnostics` | In `auto` mode, run optional held-out row/tile RMS diagnostics; never gates dither model application. |
| `--no-figures` | Skip the per-exposure step figures. |
| `--plan-only` | Validate and estimate space without processing. |

Long-wave groups run first, so their object masks are ready for the
short-wave groups; if any long-wave group fails, the short-wave groups are not
started. Each worker holds one dither group in memory, about 1.3 GB for
eight dithers. In testing, one eight-dither group took about 1.2 minutes on
one CPU, including figures.

## Outputs

| Path | Contents |
| --- | --- |
| `cal_sky_renorm/` | Sky-subtracted exposures for sky-subtracted mosaics. |
| `cal/` | 1/f and wisp-corrected exposures with the original sky retained. |
| `models/` | Per exposure: `OBJMASK`, `SKYMASK`, `ONEF` (1/f), `SKYBKG` (Background2D) and, on wisp detectors, `WISP`; `SKYMED` and `PEDESTAL` in the primary header. |
| `figures/` | One PNG per exposure showing each step from input to both products. |
| `summary.json` | Settings and, per exposure: masks, 1/f fit and QC, sky, Background2D, wisp and pedestal measurements; group errors. |

In each exposure figure, the final `cal/` and `cal_sky_renorm/` panels share
the same absolute MJy/sr stretch. The `cal/` panel is not median-centered for
display, so its retained sky can make it look brighter than the sky-subtracted
panel.

SCI headers record `S2PVER`, `STGORDER`, `ONEFMETH` (`DITHER` or `SELF`),
`ONEFSTAT`, `ONEFAPPL`, `ONEFVER`, `DITH1F` (whether a 1/f correction was
applied), `DITHQC` (when optional held-out diagnostics were run), `DITHN` (the
number of peers in the selected scene),
`ONEFOBJ`, `SEGNSIG`, `SEGNPIX`, `SEGBKG`, `SKYSUB`, `SKYMED`, `BKG2D`,
`WISPMETH` and `PEDESTAL`.
`summary.json` records the selected method, scan selection reason, dither
attempts when `auto` is used, and any group errors.
