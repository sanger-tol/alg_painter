import logging
from pathlib import Path

import polars as pl
import typer

from alg_painter.plotters.bars import main as plot_bars
from alg_painter.plotters.compact import main as plot_compact
from alg_painter.plotters.tiles import main as plot_tiles

logger = logging.getLogger(__name__)


def normalize_location_columns(locations: pl.DataFrame) -> pl.DataFrame:
    """Normalize column names and data types for consistency across styles."""
    if "assigned_alg" not in locations.columns:
        raise ValueError("Location table must contain assigned_alg column")

    locations = locations.with_columns(
        pl.col("query_chr").str.replace(r":.*", ""),
        pl.col("position").cast(pl.Float64, strict=False),
        pl.col("assigned_alg").cast(pl.Utf8),
    )
    return locations


def detect_lengths_format(lengths_file: Path) -> str:
    """Return 'final' for Chrom/Length_Mb tables, or 'draft' for .fai files."""
    with open(lengths_file) as fh:
        first_line = fh.readline().strip().split("\t")

    if first_line[:2] == ["Chrom", "Length_Mb"]:
        return "final"
    if len(first_line) >= 2:
        return "draft"
    raise ValueError(f"Could not detect lengths file format: {lengths_file}")


def load_lengths(lengths_file: Path, assembly_mode: str = "auto") -> pl.DataFrame:
    """Load either final chrom_lengths.tsv or draft .fai lengths into bp."""
    detected = detect_lengths_format(lengths_file)
    if assembly_mode != "auto" and assembly_mode != detected:
        expected = "Chrom/Length_Mb TSV" if assembly_mode == "final" else ".fai"
        raise ValueError(
            f"{lengths_file} looks like {detected!r} lengths, but --assembly-mode {assembly_mode} expects {expected}"
        )

    if detected == "final":
        chrom_lengths = pl.read_csv(lengths_file, separator="\t")
        chrom_lengths = (
            chrom_lengths.with_columns(pl.col("Length_Mb").mul(1e6).alias("length"))
            .select("Chrom", "length")
            .rename({"Chrom": "query_chr"})
        )
        return chrom_lengths

    chrom_lengths = pl.read_csv(
        lengths_file,
        separator="\t",
        has_header=False,
        new_columns=["query_chr", "length"],
    ).select("query_chr", "length")
    return chrom_lengths


def load_data(
    location_file: Path, lengths_file: Path | None = None, assembly_mode: str = "auto"
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Load BUSCO locations and chromosome/scaffold lengths."""
    locations = pl.read_csv(location_file, separator="\t")
    locations = normalize_location_columns(locations)

    if lengths_file:
        chrom_lengths = load_lengths(lengths_file, assembly_mode=assembly_mode)
    else:
        if assembly_mode in {"final", "draft"}:
            logging.warning(
                "[WARN] No lengths file supplied; estimating lengths from BUSCO "
                f"positions despite --assembly-mode {assembly_mode}"
            )
        chrom_lengths = locations.groupby("query_chr")["position"].max().reset_index()  # type: ignore
        chrom_lengths["length"] = chrom_lengths["position"] * 1.05

    return locations, chrom_lengths


def format_label(labels_list: list[str], wrap: int = 4) -> str:
    """Return a label, optionally wrapped by element count."""
    if wrap <= 0:
        return "; ".join(labels_list)

    wrapped_lines = ["; ".join(labels_list[index : index + wrap]) for index in range(0, len(labels_list), wrap)]
    return "\n".join(wrapped_lines)


def collapse_consecutive(values: list[str]) -> list[str]:
    """Collapse consecutive duplicate values."""
    collapsed: list[str] = []
    for value in values:
        if not collapsed or collapsed[-1] != value:
            collapsed.append(value)
    return collapsed


def calculate_windowed_alg_labels(
    locations: pl.DataFrame,
    alg_column: str,
    alg_order: list[str],
    window_mb: float,
    min_buscos: int,
    min_fraction: float,
    wrap: int,
) -> dict[str, str]:
    """Return labels from dominant ALGs in non-overlapping genomic windows."""
    if window_mb <= 0:
        raise ValueError("window_mb must be greater than 0")
    if not 0 < min_fraction <= 1:
        raise ValueError("min_fraction must be > 0 and <= 1")

    alg_order_index = {label: index for index, label in enumerate(alg_order)}
    window_bp = window_mb * 1_000_000

    windowed = locations.drop_nulls("position").with_columns(
        (pl.col("position") // window_bp).cast(pl.Int32).alias("window")
    )

    counts = (
        windowed.group_by(["query_chr", "window", alg_column])
        .agg(n=pl.col(alg_column).count())
        .sort(["query_chr", "window", "n"], descending=[False, False, True])
    )

    if counts.is_empty():
        return {}

    counts = counts.with_columns(
        pl.col(alg_column)
        .map_elements(lambda alg: alg_order_index.get(alg, len(alg_order)), return_dtype=pl.UInt32)
        .alias("alg_sort")
    )
    counts = counts.sort(["query_chr", "window", "n", "alg_sort"], descending=[False, False, True, False])

    dominant = counts.group_by(["query_chr", "window"]).head(1)
    totals = windowed.group_by(["query_chr", "window"]).agg(window_n=pl.col(alg_column).count())
    dominant = dominant.join(totals, on=["query_chr", "window"])
    dominant = dominant.with_columns((pl.col("n") / pl.col("window_n")).alias("fraction")).filter(
        (pl.col("n") >= min_buscos) & (pl.col("fraction") >= min_fraction)
    )

    labels: dict[str, str] = {}
    for chrom in dominant.select("query_chr").unique().to_series():
        chrom_windows = dominant.filter(pl.col("query_chr") == chrom).sort("window")
        algs = collapse_consecutive(chrom_windows.select(alg_column).to_series().to_list())
        if algs:
            labels[str(chrom)] = format_label(algs, wrap=wrap)

    return labels


def calculate_alg_labels(
    locations: pl.DataFrame,
    alg_column: str,
    alg_order: list[str],
    threshold: int = 5,
    wrap: int = 4,
) -> dict[str, str]:
    """Return ALG labels ordered by median BUSCO position per chromosome."""
    alg_order_index = {label: index for index, label in enumerate(alg_order)}

    counts = (
        locations.group_by(["query_chr", alg_column])
        .agg(
            n=pl.col(alg_column).count(),
            median_position=pl.col("position").median(),
        )
        .filter(pl.col("n") >= threshold)
    )

    labels: dict[str, str] = {}
    for chrom in locations.select("query_chr").unique().to_series():
        chrom_counts = counts.filter(pl.col("query_chr") == chrom).clone()
        chrom_counts = chrom_counts.with_columns(
            pl.col(alg_column)
            .map_elements(lambda alg: alg_order_index.get(alg, len(alg_order)), return_dtype=pl.UInt32)
            .alias("alg_sort")
        )
        chrom_counts = chrom_counts.sort(["median_position", "alg_sort"])
        algs = chrom_counts.select(alg_column).to_series().to_list()
        if algs:
            labels[str(chrom)] = format_label(algs, wrap=wrap)

    return labels


def paint_to_plot(
    ctx: typer.Context,
    assembly_mode: str,
    prefix: str,
    locations: Path,
    palette_data: dict,
    index_file: Path,
    minimum: int,
    differences: bool,
    bar_height: int,
    bar_width: float,
    plot_style: str = "tiles",
):
    """Main entry point - routes to appropriate plotting function."""
    logger.debug(ctx.obj)

    output_dir = Path(prefix)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get the order from the palette data (preserves order from palette)
    alg_order = list(palette_data.keys())

    locations_df, chrom_lengths = load_data(
        location_file=locations,
        lengths_file=index_file,
        assembly_mode=assembly_mode,
    )

    # Validate ALGs
    is_valid_alg = locations_df.select("assigned_alg").to_series().is_in(list(alg_order))
    logger.info(f"Valid ALGs count: {is_valid_alg.sum()} out of {len(locations_df)}")

    plotted_chroms = set(chrom_lengths.select("query_chr").to_series().to_list())
    chrom_lengths = chrom_lengths.filter(pl.col("query_chr").is_in(list(plotted_chroms)))

    if chrom_lengths.is_empty():
        raise ValueError("No plotted chromosomes/scaffolds have matching lengths")

    label_locations = locations_df.filter((pl.col("buscoID") != "NA") & is_valid_alg)

    # Calculate labels (with windowing if specified)
    if ctx.obj["has_windowed_labels"] and ctx.obj["label_window_mb"] and ctx.obj["label_window_mb"] > 0:
        logger.info(
            f"Labelling dominant ALGs in {ctx.obj['label_window_mb']:g} Mb windows (min BUSCOs: {ctx.obj['label_window_min_buscos']})"
        )
        alg_labels = calculate_windowed_alg_labels(
            label_locations,
            alg_column="assigned_alg",
            alg_order=alg_order,
            window_mb=ctx.obj["label_window_mb"],
            min_buscos=ctx.obj["label_window_min_buscos"],
            min_fraction=ctx.obj["label_window_min_fraction"],
            wrap=ctx.obj["label_wrap"],
        )
    else:
        alg_labels = calculate_alg_labels(
            label_locations,
            alg_column="assigned_alg",
            alg_order=alg_order,
            threshold=ctx.obj["label_window_min_buscos"],
            wrap=ctx.obj["label_wrap"],
        )

    # Route to appropriate plotting function
    match plot_style:
        case "tiles":
            plot_tiles(
                locations=locations_df,
                chrom_lengths=chrom_lengths,
                output_prefix=prefix,
                ctx=ctx,
                palette_data=palette_data,
                alg_order=alg_order,
                bar_height=bar_height,
                alg_labels=alg_labels,
            )
        case "bars":
            plot_bars(
                locations=locations_df,
                chrom_lengths=chrom_lengths,
                output_prefix=prefix,
                ctx=ctx,
                palette_data=palette_data,
                alg_order=alg_order,
                window_mb=0.5,
            )
        case "compact":
            plot_compact(
                locations=locations_df,
                chrom_lengths=chrom_lengths,
                output_prefix=prefix,
                ctx=ctx,
                palette_data=palette_data,
                alg_order=alg_order,
                window_mb=0.5,
            )
        case "all":
            plot_tiles(
                locations=locations_df,
                chrom_lengths=chrom_lengths,
                output_prefix=prefix,
                ctx=ctx,
                palette_data=palette_data,
                alg_order=alg_order,
                bar_height=bar_height,
                alg_labels=alg_labels,
            )
            plot_bars(
                locations=locations_df,
                chrom_lengths=chrom_lengths,
                output_prefix=prefix,
                ctx=ctx,
                palette_data=palette_data,
                alg_order=alg_order,
                window_mb=0.5,
            )
            plot_compact(
                locations=locations_df,
                chrom_lengths=chrom_lengths,
                output_prefix=prefix,
                ctx=ctx,
                palette_data=palette_data,
                alg_order=alg_order,
                window_mb=0.5,
            )
        case _:
            raise ValueError(f"Unknown plot style: {plot_style}. Choose tiles, bars, or compact.")
