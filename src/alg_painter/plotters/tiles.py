import logging
import math

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import polars as pl
import typer

from alg_painter.plotters.plotter_utils import get_chromosome_order, setup_font, split_balanced

logger = logging.getLogger(__name__)


def main(
    locations: pl.DataFrame,
    chrom_lengths: pl.DataFrame,
    output_prefix: str,
    ctx: typer.Context,
    palette_data: dict,
    alg_order: list,
    bar_height: float,
    alg_labels: dict[str, str],
) -> None:
    """Plot individual BUSCO tiles."""
    setup_font()

    chrom_order = get_chromosome_order(chrom_lengths, ctx)
    n_chroms = len(chrom_order)
    logger.info(f"Plotting {len(locations)} BUSCOs across {n_chroms} chromosomes/scaffolds (tile style)...")

    panel_size = max(1, int(ctx.obj["panel_size"]))
    max_columns = max(1, int(ctx.obj["max_panel_columns"]))
    ncols = min(max_columns, max(1, math.ceil(n_chroms / panel_size)))
    panel_chroms_list = split_balanced(chrom_order, ncols)
    chroms_per_panel = max((len(panel_chroms) for panel_chroms in panel_chroms_list), default=0)
    logger.info(f"Layout: {ncols} column(s), up to {chroms_per_panel} chromosomes/scaffolds per column")

    compact_layout = n_chroms > ctx.obj["compact_label_treshold"]
    label_fontsize = 8 if compact_layout else 10
    panel_height = max(
        ctx.obj["min_plot_height"], ctx.obj["row_height"] * min(ctx.obj["max_panel_columns"], max(1, n_chroms))
    )

    panel_limits: list[float] = []
    for panel_chroms in panel_chroms_list:
        if panel_chroms:
            panel_max_length = (
                chrom_lengths.filter(pl.col("query_chr").is_in(panel_chroms)).select("length").max().item()
            )
            panel_limits.append(panel_max_length * (1 + ctx.obj["label_padding_factor"]))
        else:
            panel_limits.append(1.0)

    fig_width = ctx.obj["plot_body_width"] / 2.54
    panel_widths = [1] * ncols if compact_layout else panel_limits
    fig, axes = plt.subplots(
        nrows=1,
        ncols=ncols,
        figsize=(fig_width, panel_height / 2.54),
        squeeze=False,
        gridspec_kw={"width_ratios": panel_widths},
    )
    axes = axes[0]

    for col_idx in range(ncols):
        ax = axes[col_idx]
        panel_chroms = panel_chroms_list[col_idx]
        if not panel_chroms:
            ax.axis("off")
            continue

        y_spacing = ctx.obj.get("y_spacing", 1.5)
        y_positions = {chrom: i * y_spacing for i, chrom in enumerate(reversed(panel_chroms))}
        panel_limit = panel_limits[col_idx]

        for chrom in panel_chroms:
            y = y_positions[chrom]
            length = chrom_lengths.filter(pl.col("query_chr") == chrom).select("length").item()
            bar_bottom = y - bar_height / 2

            rect = patches.Rectangle(
                (0, bar_bottom),
                length,
                bar_height,
                facecolor="white",
                edgecolor="black",
                linewidth=0.5,
            )
            ax.add_patch(rect)

            chrom_buscos = locations.filter(
                (pl.col("query_chr") == chrom) & pl.col("assigned_alg").is_in(list(alg_order))
            )

            for row in chrom_buscos.iter_rows(named=True):
                if row["position"] is not None:
                    color = palette_data.get(row["assigned_alg"], "#858585")
                    tile = patches.Rectangle(
                        (row["position"] - ctx.obj["tile_width_bp"] / 2, bar_bottom),
                        ctx.obj["tile_width_bp"],
                        bar_height,
                        facecolor=color,
                        edgecolor="none",
                    )
                    ax.add_patch(tile)

            if chrom in alg_labels:
                ax.text(
                    length * (1 + ctx.obj["label_offset_factor"]),
                    y,
                    alg_labels[chrom],
                    va="center",
                    ha="left",
                    fontsize=label_fontsize,
                    color="#333333",
                    linespacing=0.85,
                )

        ax.set_xlim(0, panel_limit)
        ax.set_ylim(-(bar_height / 2 + 0.5), (len(panel_chroms) - 1) * y_spacing + bar_height / 2 + 0.5)
        ax.set_xlabel("")
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f"{x / 1e6:.0f}"))
        ax.set_yticks([y_positions[chrom] for chrom in panel_chroms])
        ax.set_yticklabels(panel_chroms, fontsize=label_fontsize)
        ax.set_ylabel("")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)

    fig.supxlabel("Position (Mb)", fontsize=11)
    plt.tight_layout()

    # Create legend
    legend_elements = [patches.Patch(facecolor=palette_data.get(alg, "#858585"), label=alg) for alg in alg_order]
    if legend_elements:
        fig.legend(
            handles=legend_elements,
            title=ctx.obj["legend_title"],
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=10,
            title_fontsize=11,
            ncol=ctx.obj["column_count"],
            columnspacing=1.2,
        )

    for ext in ("png", "svg"):
        output_file = f"{output_prefix}/{output_prefix}_plots_tiles.{ext}"
        dpi = 300 if ext == "png" else None
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        logger.info(f"[INFO] Saved: {output_file}")

    plt.close()
