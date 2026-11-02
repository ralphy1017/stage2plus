# stage2plus

Standalone JWST NIRCam exposure cleanup for eight-dither Stage 2 `cal` groups.
It applies registered-dither 1/f removal, an object-masked median sky
subtraction, fixed-basis NMF wisp removal where applicable, and a final masked
scalar zero point. It writes sky-subtracted `cal_sky_renorm` FITS for drizzle.
The fitted 1/f correction is always applied; its held-out quality score is
recorded separately. No Jumprope runtime or directory layout is required.

## Install

Clone this repository, enter it, and install the Python dependencies (tested
with Python 3.13):

```bash
git clone <repository-url> stage2plus
cd stage2plus
python -m venv ../stage2plus-venv
source ../stage2plus-venv/bin/activate
python -m pip install numpy scipy astropy photutils pillow matplotlib
```

Before a remote URL is configured, clone a local copy with
`git clone --no-local /path/to/stage2plus stage2plus`. The `--no-local` option
also works when the source and destination are on different filesystems.

`stage2plus.py` is the entry point. Its three Python helper modules and the
existing NMF basis FITS in `templates/` are included here. The input exposures
and generated FITS/PNGs are **not** part of this Git repository.

## Input requirements

- One absolute or relative path to a calibrated NIRCam `*_cal.fits` file per
  line in a text input list. Basenames must be unique, and files must contain
  `SCI` (or `SCI_ORIG`), `ERR`, and `DQ` with usable WCS. If `SCI_ORIG` exists,
  the reducer uses it as the uncorrected science input.
- Every selected observation/visit/detector/filter group must contain exactly
  **eight dithers**. Groups are kept separate by observation ID and visit ID.
- For each short-wave wisp-sensitive detector/filter pair (NRCA3, NRCA4,
  NRCB3, NRCB4), a matching NMF basis must exist in `templates/` or in a
  directory supplied with `--template-dir`. Missing bases cause a plan error.
  The bundled assets cover the filters available in their subdirectories;
  this includes every wisp-sensitive pair in the local program 1813 run.
- Supply inputs from a single calibration version. Do not mix duplicate
  basenames from separate Stage 2 reductions in one list.

## One observation

Set the calibrated input directory and an **output directory outside it**.
For example, the prefix below selects observation 001 of program 1813:

```bash
CAL=/path/to/nircam/stage2/cal
RUN=/path/to/output/obs001
mkdir -p "$RUN"
find "$CAL" -maxdepth 1 -type f -name 'jw01813001*_nrc*_cal.fits' | sort > "$RUN/inputs.txt"
python stage2plus.py --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 8 --plan-only
MPLCONFIGDIR=/tmp python stage2plus.py --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 8
```

`--plan-only` validates the groups and template availability, estimates FITS
output bytes, and reports free space without processing exposures. Review its
counts before starting the full command.

## Multiple observations

Put all selected exposures in one list. The reducer keeps overlapping visit
IDs from different observations separate:

```bash
CAL=/path/to/nircam/stage2/cal
RUN=/path/to/output/multi
mkdir -p "$RUN"
find "$CAL" -maxdepth 1 -type f \( -name 'jw01813001*_nrc*_cal.fits' -o -name 'jw01813005*_nrc*_cal.fits' \) | sort > "$RUN/inputs.txt"
python stage2plus.py --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 16 --plan-only
MPLCONFIGDIR=/tmp python stage2plus.py --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 16
```

## Full calibrated database

Build a list from the selected Stage 2 NIRCam calibration directory. For a
database spread across directories, concatenate path lists, sort them, and
remove duplicate basenames before running. The current reducer requires eight
dithers in every selected group; filter or complete groups in the manifest if
your database also contains other observing patterns.

```bash
CAL=/path/to/nircam/stage2/cal
RUN=/path/to/output/full
mkdir -p "$RUN"
find "$CAL" -maxdepth 1 -type f -name '*_nrc*_cal.fits' | sort > "$RUN/inputs.txt"
python stage2plus.py --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 32 --plan-only
MPLCONFIGDIR=/tmp /usr/bin/time -v python stage2plus.py \
  --input-list "$RUN/inputs.txt" --output-dir "$RUN" --workers 32 \
  > "$RUN/run.log" 2>&1
```

The local 1813 batch used 640 exposures in 80 groups, `--workers 32`, and
finished in 6 minutes 49.5 seconds on a 128-thread workstation. Run time
depends strongly on storage speed, filter mix, and available memory.

## Configuration and outputs

| Option | Meaning |
| --- | --- |
| `--input-list PATH` | Required newline-separated list of FITS inputs. |
| `--output-dir PATH` | Required output root; must differ from the input directory. |
| `--workers N` | Maximum parallel detector-group processes; defaults to 32. |
| `--template-dir PATH` | Alternate NMF basis root; defaults to the bundled `templates/`. |
| `--plan-only` | Validate and estimate storage without writing products. |

The reducer first processes long-wave groups to make object masks, then
short-wave groups in parallel. Each worker caches the eight dither exposures
for its group, so memory use grows with `--workers`. Lower this value when RAM
or input-storage throughput is limited. Output storage is roughly the size of
the retained input FITS extensions; use `--plan-only` for a dataset estimate.

The output directory contains:

| Path | Contents |
| --- | --- |
| `cal_sky_renorm/` | One corrected, sky-subtracted FITS per input basename. |
| `_lw_masks/` | Long-wave object masks projected onto short-wave frames. |
| `figures/` | Input/final PNG contact sheets by observation, visit, and module. |
| `summary.json` | Per-exposure fits, masks, sky levels, QC, and errors. |

The command above captures progress and resource timing in `run.log`; the
Python script itself does not create that log. In output SCI headers,
`DITH1F=True` means the 1/f model was applied and `DITHQC` records the
held-out validation result. The input image is preserved separately; the
sky-subtracted SCI product is not suitable for absolute sky-brightness
measurements. If a 1/f fit fails, the affected group is reported as an error
instead of writing an uncorrected final frame.

