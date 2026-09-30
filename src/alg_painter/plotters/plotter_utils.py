import os

import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt
import polars as pl
import typer


def setup_font() -> None:
    """Set up the font for plotting."""
    font_path = os.path.join(os.path.dirname(__file__), "fonts", "OpenSans-Regular.ttf")

    font_manager.fontManager.addfont(str(font_path))
    plt.rcParams["font.family"] = "Open Sans"
    plt.rcParams["font.style"] = "normal"
    plt.rcParams["font.weight"] = "normal"


def get_chromosome_order(chrom_lengths: pl.DataFrame, ctx: typer.Context) -> list[str]:
    """Get chromosome order from custom order if available, otherwise sort by length."""
    custom_order = ctx.obj.get("chrom_order")

    if custom_order:
        # Filter to only chromosomes that exist in the data
        available_chroms = set(chrom_lengths.select("query_chr").to_series().to_list())
        # Keep custom order but only include chromosomes that are available
        return [chrom for chrom in custom_order if chrom in available_chroms]
    else:
        # Default: sort by length (descending)
        return chrom_lengths.sort("length", descending=True).select("query_chr").to_series().to_list()


def split_balanced(values: list[str], n_groups: int) -> list[list[str]]:
    """Split values into ordered groups whose sizes differ by at most one."""
    n_groups = max(1, n_groups)
    base_size, extra = divmod(len(values), n_groups)
    groups: list[list[str]] = []
    start = 0
    for group_index in range(n_groups):
        group_size = base_size + (1 if group_index < extra else 0)
        end = start + group_size
        groups.append(values[start:end])
        start = end
    return groups
