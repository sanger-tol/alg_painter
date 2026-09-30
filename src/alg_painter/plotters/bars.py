import logging

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import polars as pl
import typer

from alg_painter.plotters.plotter_utils import get_chromosome_order, setup_font

logger = logging.getLogger(__name__)


def main(
    locations: pl.DataFrame,
    chrom_lengths: pl.DataFrame,
    output_prefix: str,
    ctx: typer.Context,
    palette_data: dict,
    alg_order: list,
    window_mb: float = 0.5,
) -> None:
    """Plot windowed BUSCO counts as stacked bars."""
    setup_font()

    window_bp = window_mb * 1_000_000
    alg_column = "assigned_alg"
    bar_height = ctx.obj.get("bar_height", 0.5)

    # Filter to valid ALGs only
    valid_algs = set(alg_order)
    locations = locations.filter(pl.col(alg_column).is_in(list(valid_algs)))

    # Create windows
    windowed = (
        locations.drop_nulls("position")
        .with_columns((pl.col("position") // window_bp).cast(pl.Int32).alias("window"))
        .with_columns(((pl.col("window") + 1) * window_bp / 1e6).alias("window_pos_mb"))
    )

    # Count ALGs per window per chromosome
    counts = (
        windowed.group_by(["query_chr", "window", "window_pos_mb", alg_column])
        .agg(n=pl.col(alg_column).count())
        .sort(["query_chr", "window"])
    )

    if counts.is_empty():
        logger.warning("No data to plot")
        return

    # Filter chromosomes with enough BUSCOs
    chrom_totals = windowed.group_by("query_chr").agg(total=pl.col(alg_column).count()).filter(pl.col("total") >= 1)

    chrom_order = get_chromosome_order(chrom_lengths.join(chrom_totals, on="query_chr"), ctx)

    n_chroms = len(chrom_order)
    logger.info(f"Plotting {len(locations)} BUSCOs across {n_chroms} chromosomes/scaffolds (bar style)...")

    # Calculate grid dimensions based on column_count
    ncols = max(1, int(ctx.obj.get("column_count", 1)))
    nrows = (n_chroms + ncols - 1) // ncols

    # Calculate subplot heights based on bar_height and number of rows
    subplot_height = max(3, bar_height * 4)
    fig_height = min(subplot_height * nrows, 20)  # Cap max height at 20 inches
    fig_width = 8 * ncols

    # Create faceted plot with calculated heights and shared y-axis
    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(fig_width, fig_height),
        sharex=False,
        sharey=True,
    )

    # Always flatten axes array for consistent 2D indexing
    axes = np.atleast_2d(axes).flatten()

    for ax_idx, chrom in enumerate(chrom_order):
        ax = axes[ax_idx]
        chrom_data = counts.filter(pl.col("query_chr") == chrom)

        if chrom_data.is_empty():
            ax.axis("off")
            continue

        # Pivot to get ALG counts per window
        pivot_data = chrom_data.pivot(on=alg_column, values="n", aggregate_function="sum", sort_columns=False).sort(
            "window_pos_mb"
        )

        window_positions = pivot_data.select("window_pos_mb").to_series().to_list()
        pivot_data = pivot_data.drop(["query_chr", "window", "window_pos_mb"])

        # Stack bars for each ALG with gap between bars
        bar_width_with_gap = window_mb * 0.8  # 80% of window size for the bar, 20% for gap
        bottom = None
        for alg in reversed(alg_order):
            if alg in pivot_data.columns:
                values = pivot_data.select(alg).to_series().to_list()
                color = palette_data.get(alg, "#858585")
                ax.bar(
                    window_positions,
                    values,
                    width=bar_width_with_gap,
                    bottom=bottom if bottom is not None else 0,
                    label=alg,
                    color=color,
                    edgecolor="none",
                )
                if bottom is None:
                    bottom = [v for v in values]
                else:
                    bottom = [b + v for b, v in zip(bottom, values)]

        ax.set_title(chrom, fontsize=20)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(True)
        ax.spines["left"].set_visible(False)
        ax.yaxis.tick_right()
        ax.grid(True, axis="y", alpha=0.3, linestyle="-", linewidth=0.5)
        ax.set_axisbelow(True)
        ax.yaxis.set_major_locator(plt.MultipleLocator(5))
        ax.tick_params(labelbottom=True, labelleft=False, labelright=True)
        ax.set_xlim(0, max(window_positions) + window_mb)

    # Turn off extra empty subplots
    for ax_idx in range(len(chrom_order), len(axes)):
        axes[ax_idx].axis("off")

    # Set x-label on bottom-right subplot
    axes[-1].set_xlabel("Position (Mb)", fontsize=10)
    fig.supylabel("Count", fontsize=10, x=0.02)

    # Create legend with larger size
    legend_elements = [patches.Patch(facecolor=palette_data.get(alg, "#858585"), label=alg) for alg in alg_order]
    if legend_elements:
        fig.legend(
            handles=legend_elements,
            title=ctx.obj["legend_title"],
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=12,
            title_fontsize=13,
            ncol=1,
        )

    plt.tight_layout()

    for ext in ("png", "svg"):
        output_file = f"{output_prefix}/{output_prefix}_plots_bars.{ext}"
        dpi = 300 if ext == "png" else None
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        logger.info(f"[INFO] Saved: {output_file}")

    plt.close()
