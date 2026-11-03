# stage2plus

Exposure-level cleanup of JWST NIRCam Stage 2 `cal` files before drizzling.
It works on any NIRCam imaging program: exposures are grouped from their FITS
headers into dither groups (same program, observation, visit, visit group,
sequence, activity, detector, filter and pupil). Every group needs at least
**four dithers**.

## What it does

For every exposure, in this order:

1. **Input**: `SCI_ORIG` if the file has one, otherwise `SCI`.
2. **Object mask**: photutils segmentation (2.8σ over the ERR array, 5
   connected pixels). Short-wave frames add the long-wave mask of the same
   exposure, reprojected through the WCS. Each object is grown by its
   equivalent radius (at most 32 px) to cover galaxy outskirts. On wisp
   detectors the bright-wisp region is also masked for the 1/f fit, except
   on amplifier rows where that would leave fewer than 64 pixels.
3. **Dither 1/f**: the group's other dithers are shifted onto the target
   through their WCS, level-matched, and median-combined into a scene. The
   masked residual (target − scene) gives amplifier-row and detector-column
   medians, which form the 1/f model. The model is zero on average over the
   masked sky, so the input's sky median is kept. `--no-onef-object-mask`
   fits with only bad pixels (and bright wisps) masked. On the 16 test
   exposures this left 6–14% more amplifier-row striping, so the object mask
   is on by default.
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
| `cal/` | input − 1/f − wisp − pedestal (sky retained) |

The two products differ only by the saved sky models. Never measure absolute
sky brightness from `cal_sky_renorm/`.

Both products keep the input's PRIMARY, ERR, DQ, AREA, VAR_* and ASDF
extensions byte for byte; only SCI is replaced. Inputs are never modified.

### 1/f quality check

Up to two other dithers are withheld from each target's scene, keeping at
least two in it: eight dithers use five in the scene and hold out two; four
dithers use two and hold out one. The check passes when the withheld
dithers' mean amplifier-row and 64-pixel tile RMS against the target both
improve by more than 5%, and no withheld dither gets more than 5% worse. The
model is applied either way; `DITHQC` records the result. If a 1/f fit cannot
be made (for example, too little dither overlap or too many masked rows), the
whole group is reported as an error instead of being written uncorrected.

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
| `--no-onef-object-mask` | Do not mask objects in the 1/f fit (later steps still do). |
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
| `cal/` | The same corrections with the sky retained. |
| `models/` | Per exposure: `OBJMASK`, `SKYMASK`, `ONEF` (1/f), `SKYBKG` (Background2D) and, on wisp detectors, `WISP`; `SKYMED` and `PEDESTAL` in the primary header. |
| `figures/` | One PNG per exposure showing each step from input to both products. |
| `summary.json` | Settings and, per exposure: masks, 1/f fit and QC, sky, Background2D, wisp and pedestal measurements; group errors. |

SCI headers record `S2PVER`, `STGORDER`, `ONEFMETH`, `DITH1F`, `DITHQC`,
`DITHN` (dithers in the 1/f scene), `ONEFOBJ` (objects masked for 1/f),
`SKYSUB`, `SKYMED`, `BKG2D`, `WISPMETH` and `PEDESTAL`.
