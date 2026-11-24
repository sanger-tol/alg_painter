# ALG_painter

This replaces the original sanger-tol/busco_painter scripts.

A Python package for the plotting and painting of ALG units to busco results resulting in a ALG assignment plot.

As of Version 2.0.0, this package has now been re-written to be much more generic, modular and configurable from a single JSON file.

Currently, this package generates 1 style of graph. This is a bar graph, the bars being chromosomes with the location of the busco gene marked by a colored bar (different for each ALG unit). This may be expanded upon in the future.

This package contains 3 subcommands:
- paint - A re-write of the original script (by Charlotte Wright), by Karen van Neikerk, for the paining of ALG's onto busco full table.tsv files. This can optionally call the NCBI API for the sequence report.
- plot - A python re-write of the original R script written by Charlotte Wright, it has been updated to be very modulular for input parameters, allowing for easy customization of the plot. [See mapping.json](assets/alg_mapping.json) which provides a single place for configuring the plot, with alg mappings and palettes.

## Input format

The input ALG mapping file is a 2 column TSV file with the busco gene id in the first column and the ALG unit name in the second column.

```tsv name="Merian Example"
7494at7088	M1
14573at7088	M1
5486at7088	M1
15450at7088	M16
12710at7088	M1
4360at7088	M1
6722at7088	M1
1550at7088	M1
```

More examples exist in the `src/assets/alg_assignments` folder [here](src/assets/alg_assignments).

## Usage

### paint
The first command paints the ALG's onto a busco full table, resulting in upto three files being produced.

```
alg --prefix {output directory} paint {clade} {busco full table}

alg --prefix TESSING paint coleoptera Abax_parallelepipedus_full_table.tsv
```

#### Output

The output files are:
- `{prefix}_paint_all_locations.tsv` - a TSV file with the busco gene id and the ALG unit name
- `{prefix}_paint_chrom_lengths.tsv` - a TSV file with the chromosome lengths (if the tool is provided with `--accession {GCA_accession}`, in the above example `--accession GCA_964197645.1`)
- `{prefix}_paint_summary.tsv` - a TSV file with the summary of the painting

### plot
The second command creates a series of plots showing the ALG's along representations of the chromosomes. By default, the tool will only output the `tiles` plot, but additional plots can be requested using the `--plot-style` option ([tiles|bars|compact|all]).

```
alg --prefix {output directory} plot coleoptera {paint_all_locations.tsv} {paint_chrom_lengths.tsv | samtools fai file}

alg --prefix {output directory} plot coleoptera TESTING_paint_all_locations.tsv TESTING_paint_chrom_lengths.tsv
```

#### Output

The output files are:
- `{prefix}_plot_tiles.{png|svg}` - a PNG or SVG image of the tiles plot
- `{prefix}_plot_bars.{png|svg}` - a PNG or SVG image of the bars plot
- `{prefix}_plot_compact.{png|svg}` - a PNG or SVG image of the compact plot


![ALG_painter_v2 tiles plot](examples/abax_parallelepipedus/TESTING_plots_tiles.png)

![ALG_painter_v2 bars plot](examples/abax_parallelepipedus/TESTING_plots_bars.png)

![ALG_painter_v2 compact plot](examples/abax_parallelepipedus/TESTING_plots_compact.png)

### assign
This command simply created an annotated busco full table, including the ALG units (used in the `TreeVal` pipeline).

```
alg assign {all_locations.tsv} {full_table.tsv} --outdir {output directory}

alg assign TESTING_paint_all_locations.tsv Abax_parallelepipedus_full_table.tsv --outdir TESTING

```

#### Output:
- [`alg_assigned_ancestral.tsv`](./example_data/output/abax_parallelepipedus)

```
OZ078038.1	22521295	22521849	C7	260.7	+	NA
OZ078033.1	31431166	31433351	C1	394.9	+	NA
OZ078036.1	17280006	17249113	C3	454.5	-	NA
OZ078025.1	41885560	41887859	C4	765.8	+	NA
OZ078041.1	18604678	18607537	C1	723.4	+	NA
OZ078025.1	29790416	29788315	C3	555.1	-	NA
OZ078025.1	26950026	26979567	C3	500.7	+	NA
OZ078026.1	2026436	2026798	C7	197.8	+	NA
OZ078030.1	34243263	34241266	C7	728.0	-	NA
OZ078040.1	13416498	13420076	CX	1104.2	+	NA
OZ078031.1	30740719	30602288	C2	769.1	-	NA
OZ078026.1	26288575	26278949	C2	449.6	-	NA
```


### list_the_configs
This command lists the data and information available for a user given clade/group id in the mapping file.

```
alg list_the_configs {alg_clade}

alg list_the_configs coleoptera_odb12
```

#### Output:

```
Mapping file found at: /{USERS INSTALL LOCATION}/site-packages/alg_painter/assets/alg_mapping.json
2026-07-29 14:15:48,962 [ERROR] Palette file not found: assets/alg_assignments/coleopteran_odb12.json
2026-07-29 14:15:48,962 [ERROR] Falling back to built-in palette: /{USERS INSTALL LOCATION}/site-packages/alg_painter/assets/alg_assignments/coleopteran_odb12.json
{
    "legend_title": "Coleopteran ALG's",
    "label_window_default_mb": 20,
    "label_window_min_buscos": 5,
    "source": "https://www.biorxiv.org/content/10.64898/2026.07.17.739156v1",
    "odb": "odb12",
    "palette": {
        "options": "default"
    },
    "alg_file": "assets/alg_assignments/coleopteran_odb12.tsv",
    "custom_order": null,
    "notes": "",
    "y_spacing": 20,
    "bar_height": 0.45,
    "row_height": 0.85,
    "column_count": 2,
    "min_plot_height": 12,
    "tile_width_bp": 50000,
    "has_windowed_labels": true
}
```

### paint
The main command of this tool, `paint` creates an annotated file ( [paint_all_locations.tsv](example_data/output/abax_parallelepipedus/TESTING_paint_all_locations.tsv) ) which we can then use for `plot`

The minimal command 


### plot



```
alg plot coleoptera_odb12 ./example_data/output/abax_parallelepipedus/TESTING_paint_all_locations.tsv /Users/dp24/Documents/alg_painter/example_data/output/abax_parallelepipedus/TESTING_paint_chrom_lengths.tsv
```

#### Output:
```
2026-09-29 13:36:27,974 [ERROR] Palette file not found: assets/alg_assignments/coleopteran_odb12.json
2026-09-29 13:36:27,974 [ERROR] Falling back to built-in palette: /Users/dp24/Documents/alg_painter/.venv/lib/python3.12/site-packages/alg_painter/assets/alg_assignments/coleopteran_odb12.json
2026-09-29 13:36:27,978 [INFO] Valid ALGs count: 3485 out of 3580
2026-09-29 13:36:27,978 [INFO] Labelling dominant ALGs in 20 Mb windows (min BUSCOs: 5)
2026-09-29 13:36:27,985 [INFO] Plotting 3580 BUSCOs across 18 chromosomes/scaffolds (tile style)...
2026-09-29 13:36:27,985 [INFO] Layout: 1 column(s), up to 18 chromosomes/scaffolds per column
2026-09-29 13:36:29,069 [INFO] [INFO] Saved: alg/alg_plots_tiles.png
2026-09-29 13:36:29,642 [INFO] [INFO] Saved: alg/alg_plots_tiles.svg
```

## Installation

```
git clone github.com/sanger-tol/alg_painter.git

cd alg_painter/

uv pip install ./

alg -h
```

or

```
pip install alg_painter
```

![ALG_painter_v2 plot](src/tests/data/alg_plotter_v2.png)
