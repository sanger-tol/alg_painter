# ALG_painter

This replaces the original sanger-tol/busco_painter scripts.

A Python package for the plotting and painting of ALG units to busco results resulting in a ALG assignment plot.

Currently, this is rather specific to the Merian Units assigned by Charlotte Wright based on odb10 output.

This package contains 3 subcommands:
- painter - A re-write of the original script (by Charlotte Wright), based on [Karen van Niekerk's merian-busco-painter](https://github.com/Karenvn/merian-busco-painter), for the painting of ALG's onto busco full table.tsv files. This can optionally call the NCBI API
- plotter - A python re-write of the original R script written by Charlotte Wright.
- plotter_2 - A re-write originally intended for GenomeNote production, so includes specific colour palettes and other logic.

## Installation

```
git clone github.com/sanger-tol/alg_painter.git

cd alg_painter/

uv pip install ./

alg -h
```

![ALG_painter_v2 plot](src/tests/data/alg_plotter_v2.png)
