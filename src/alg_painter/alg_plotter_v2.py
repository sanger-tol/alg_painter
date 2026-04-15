import logging
from pathlib import Path

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib import font_manager

from alg_painter.plotter_colours import MERIAN_COLOURS

logger = logging.getLogger(__name__)


def _load_locations_file(location_file: Path) -> pd.DataFrame:
    """
    Load BUSCO locations from a TSV file.
    """
    # Read locations
    locations = pd.read_csv(location_file, sep="\t")

    # Clean chromosome names (remove any :scaffolds)
    locations["query_chr"] = locations["query_chr"].str.replace(":.*", "", regex=True)

    return locations


def load_data(location_file: Path, lengths_file=None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load BUSCO locations and chromosome lengths"""
    locations = _load_locations_file(location_file)

    # Read or estimate chromosome lengths
    if lengths_file:
        chrom_lengths = pd.read_csv(
            lengths_file,
            sep="\t",
            header=None,
            names=["query_chr", "length", "offset", "linebases", "linewidth"],
        )
    else:
        # Fallback: estimate from max position + 5% padding
        chrom_lengths = locations.groupby("query_chr")["position"].max().reset_index()
        chrom_lengths["length"] = chrom_lengths["position"] * 1.05

    return locations, chrom_lengths


def get_palette(palette_name: str) -> dict:
    return MERIAN_COLOURS[palette_name]


def filter_chromosomes(locations, minimum_buscos=3):
    """
    Filter chromosomes by minimum BUSCO count
    """
    busco_counts = locations.groupby("query_chr").size().reset_index(name="n_busco")
    valid_chroms = busco_counts[busco_counts["n_busco"] >= minimum_buscos]["query_chr"]
    return locations[locations["query_chr"].isin(valid_chroms)]


def setup_font(font_path: Path) -> None:
    """
    Set up the font for plotting
    """
    font_manager.fontManager.addfont(str(font_path))
    plt.rcParams["font.family"] = "Open Sans"
    plt.rcParams["font.style"] = "normal"
    plt.rcParams["font.weight"] = "normal"


def calculate_merian_labels(locations, threshold=5):
    """
    Calculate Merian element labels for each chromosome.
    Returns {chromosome: 'M1; M3; M18'} for elements with >= threshold BUSCOs.
    """
    # Count BUSCOs per chromosome-Merian combination
    counts = locations.groupby(["query_chr", "assigned_chr"]).size().reset_index(name="n")

    # Filter by threshold
    counts = counts[counts["n"] >= threshold]

    # Group by chromosome and join Merian elements
    labels = {}
    for chrom in counts["query_chr"].unique():
        chrom_merians = counts[counts["query_chr"] == chrom]["assigned_chr"].tolist()
        # Sort Merian elements (MZ first, then M1-M31 numerically)
        sorted_merians = sorted(
            chrom_merians, key=lambda x: (x != "MZ", int(x[1:]) if x != "MZ" else 0)
        )
        labels[chrom] = "; ".join(sorted_merians)

    return labels


def plot_merian_chromosomes(
    locations: pd.DataFrame,
    chrom_lengths: pd.DataFrame,
    output_prefix: str,
    minimum_buscos: int,
    palette: dict,
    label_threshold: int,
):
    """
    Create the main Merian plot with chromosome labels
    """

    # Keep track of placeholder rows in the filtered locations (where buscoID is NA)
    is_placeholder = locations["buscoID"] == "NA"

    # Filter only valid Merian assignments for actual BUSCOs
    valid_merians = ["MZ"] + [f"M{i}" for i in range(1, 32)]
    locations_filtered = locations[
        is_placeholder | locations["assigned_chr"].str.upper().isin(valid_merians)
    ]
    locations_filtered["assigned_chr"] = locations_filtered["assigned_chr"].str.upper()

    # Calculate Merian labels for each chromosome
    merian_labels = calculate_merian_labels(locations, threshold=label_threshold)

    # Order chromosomes by length (descending) - USE ALL CHROMOSOMES FROM chrom_lengths
    chrom_order = chrom_lengths.sort_values("length", ascending=False)["query_chr"].tolist()

    # DO NOT FILTER - keep all chromosomes even if they have no BUSCOs

    # Prepare y-positions (reversed so longest is at top)
    n_chroms = len(chrom_order)
    y_positions = {chrom: i for i, chrom in enumerate(reversed(chrom_order))}

    # Calculate plot height dynamically
    plot_height = max(15, 0.9 * n_chroms)

    # Create figure with extra space for labels
    fig, ax = plt.subplots(figsize=(30 / 2.54, plot_height / 2.54))

    # Plot chromosome bars
    for chrom in chrom_order:
        y = y_positions[chrom]
        length = chrom_lengths[chrom_lengths["query_chr"] == chrom]["length"].values[0]
        # Draw chromosome backbone (white rectangle with black border)
        rect = patches.Rectangle(
            (0, y - 0.4), length, 0.8, facecolor="white", edgecolor="black", linewidth=0.5
        )
        ax.add_patch(rect)

        # Plot BUSCO positions
        chrom_buscos = locations[locations["query_chr"] == chrom]
        for _, busco in chrom_buscos.iterrows():
            if pd.notna(busco["position"]):
                color = palette.get(busco["assigned_chr"], (0.85, 0.85, 0.85))
                tile = patches.Rectangle(
                    (busco["position"] - 25000, y - 0.4),
                    50000,
                    0.8,
                    facecolor=color,
                    edgecolor="none",
                )
                ax.add_patch(tile)

        # Add Merian label to the right of the chromosome
        if chrom in merian_labels:
            label_text = merian_labels[chrom]
            ax.text(
                length * 1.02, y, label_text, va="center", ha="left", fontsize=10, color="#333333"
            )

    # Set axis limits and labels
    max_length = chrom_lengths["length"].max()
    ax.set_xlim(0, max_length * 1.25)
    ax.set_ylim(-0.6, n_chroms - 0.4)

    # X-axis: position in Mb
    ax.set_xlabel("Position (Mb)", fontsize=11)
    ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f"{x / 1e6:.0f}"))

    # Y-axis: chromosome names
    ax.set_yticks([y_positions[c] for c in chrom_order])
    ax.set_yticklabels(chrom_order, fontsize=10)
    ax.set_ylabel("")

    # Remove top and right spines
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)

    # Create legend
    legend_elements = []
    for merian in ["MZ"] + [f"M{i}" for i in range(1, 32)]:
        legend_elements.append(patches.Patch(facecolor=palette[merian], label=merian))

    ax.legend(
        handles=legend_elements,
        title="Merian elements",
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        fontsize=10,
        title_fontsize=11,
    )

    plt.tight_layout(rect=(0, 0, 0.82, 1))

    # Save outputs
    for ext in ["png", "svg"]:
        output_file = f"{output_prefix}.{ext}"
        dpi = 300 if ext == "png" else None
        plt.savefig(output_file, dpi=dpi, bbox_inches="tight")
        logger.info(f"[Plotter V2] Saved: {output_file}")

    plt.close()


def plotter_v2_main(arguments):
    logger.info("[Plotter V2] Loading data from %s", arguments.file)
    locations, chrom_lengths = load_data(arguments.file, arguments.lengths_file)

    logger.info(
        f"[Plotter V2] Plotting {len(locations)} BUSCOs across {len(chrom_lengths)} chromosomes"
    )
    setup_font(arguments.font)

    # Filter by minimum BUSCOs, but keep track of all chromosomes from chrom_lengths
    # This ensures sex chromosomes (like W) are always plotted even if empty
    # all_chromosomes = set(chrom_lengths["query_chr"].tolist())
    locations_filtered = filter_chromosomes(locations, arguments.minimum)

    final_palette: dict = get_palette(arguments.palette)

    plot_merian_chromosomes(
        locations_filtered,
        chrom_lengths,
        arguments.prefix,
        arguments.minimum,
        final_palette,
        arguments.label_threshold,
    )
