# Replot the figure set

This folder contains the portable renderer and the aggregate tables needed to
recreate the figures (`Fig1`–`Fig4` and `FigS1`–`FigS6`). It reads
only the TSV files in `source_data/`; it does not need the original analysis
workspace or patient-level data.

## Requirements

- Python 3.10 or later
- `numpy`
- `matplotlib`

Install the two Python packages if needed:

```bash
python -m pip install numpy matplotlib
```

## Render

From the repository root, run:

```bash
python figures/make_figures.py
```

The script writes one PDF and one 300-dpi PNG per figure, plus
`contact_sheet.png`, to `figures/output/` by default. To select figures
or change the input and output folders:

```bash
python figures/make_figures.py --figures Fig2 FigS5
python figures/make_figures.py --data-dir path/to/source_data --output-dir path/to/output
```

Run `python figures/make_figures.py --help` for the available options.

The publication reference PDFs in `results/reference_figures/` use Arial Bold.
The renderer can use an installed Arial font or accept a local file with
`--font-path /path/to/arialbd.ttf`. On systems without Arial it can use a
fallback font. No font files are required to view the reference PDFs.

## Data and license

The bundled tables contain cohort-, method-, and summary-level values used by
the plotted panels. They contain no patient-level rows or direct identifiers.
The aggregate results are released under **CC BY 4.0**. The renderer is
released under **MIT**; see the repository license files for the full terms.
