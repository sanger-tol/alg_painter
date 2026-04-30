import argparse
import logging
import os
from pathlib import Path

from alg_painter.alg_painter import painting_main as alg_painter
from alg_painter.alg_plotter_v1 import plotter_v1_main as alg_plotter_v1
from alg_painter.alg_plotter_v2 import plotter_v2_main as alg_plotter_v2

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("alg_painter.log"),  # logs to file
        logging.StreamHandler(),  # logs to console
    ],
)
logger = logging.getLogger("alg_painter")

VERSION = "0.1.0"


def file_validator(in_file: str, file_type) -> Path | None:
    """
    Validate files exist and are infact files
    """
    file_checks = {"font": {"format": [".ttf"], "validated": False}}

    file_path = Path(in_file)

    if not file_path.exists() or not file_path.is_file():
        raise argparse.ArgumentTypeError(f"{in_file} is not a valid file")

    if (
        file_type == "font"
        and file_path.suffix not in file_checks["font"]["format"]
    ):
        raise argparse.ArgumentTypeError(
            f"{in_file} has unsupported extension '{file_path.suffix}', expected one of {file_checks['font']['format']}"
        )

    return file_path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="alg_painter", description="Ancestral Linkage Group Painter"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version=VERSION,
    )

    subparsers = parser.add_subparsers(
        help="Subcommand help", dest="command", required=True
    )

    painter = subparsers.add_parser(
        "painter",
        help="Run the Ancestral Linkage Group Painting command",
        description="Paint the BUSCOs colored by their putative ALG assignments",
    )
    painter.add_argument(
        "--alg_table",
        "-a",
        type=Path,
        required=True,
        help="The reference ALG table for mapping",
    )
    painter.add_argument(
        "--query_table",
        "-q",
        type=Path,
        required=True,
        help="The Query table to paint (e.g. your assembly's ALG table)",
    )
    painter.add_argument(
        "--prefix",
        "-p",
        type=str,
        default="alg_painter",
        help="Output prefix of files",
    )
    painter.add_argument("--accession", help="Assembly Accession")
    painter.add_argument(
        "--write_summary",
        action="store_true",
        help="Write a summary of the results",
    )
    painter.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"ALG-{VERSION}-PAINTER-2.1.0",
    )

    plotter = subparsers.add_parser(
        name="plotter",
        help="Original plotter for generic plotting of BUSCO locations",
        description=(
            "Generate a plot showing locations of putative Merian elements "
            "(or species-level orthologs) in a genomic assembly."
        ),
    )
    plotter.add_argument(
        "-f",
        "--file",
        type=Path,
        # required=True,
        help="The output file from the painter command, should be a locations.tsv file",
    )
    plotter.add_argument(
        "-i",
        "--index",
        type=Path,
        default=None,
        help="Genome index file (.fai)",
    )
    plotter.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="alg_plotter_v2",
        help="Output prefix of files",
    )
    plotter.add_argument(
        "-m",
        "--merians",
        action="store_true",
        default=False,
        help="Are you comparing a genome to Merian elements?",
    )
    plotter.add_argument(
        "-d",
        "--differences",
        action="store_true",
        default=False,
        help="Colour only BUSCOs moved from the dominant chromosome",
    )
    plotter.add_argument(
        "-n",
        "--minimum",
        type=int,
        default=3,
        help="Minimum number of BUSCOs per contig to retain (default: 3)",
    )
    plotter.add_argument(
        "--bar_height",
        type=int,
        default=12,
        help="Height of each BUSCO bar in bp (default: 12)",
    )
    plotter.add_argument(
        "--bar_width",
        type=float,
        default=2e4,
        help="Half-width of each BUSCO bar in bp (default: 2e4)",
    )
    plotter.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"ALG-{VERSION}-PLOTTER-2.0.0",
    )

    plotter_2 = subparsers.add_parser(
        "plotter2",
        help="Newer Plotting process, currently in use for genomenotes, Plot BUSCO locations colored by ALG assignments",
        description="Plot the BUSCO labelled by the ALG painter onto a per chromosome plot",
    )
    plotter_2.add_argument(
        "-f",
        "--file",
        type=Path,
        required=True,
        help="The output file from the painter command, should be a locations.tsv file",
    )
    plotter_2.add_argument(
        "-l",
        "--lengths_file",
        default=None,
        type=Path,
        help="Path to a chromosome lengths TSV or genome index file (.fai, optional)",
    )
    plotter_2.add_argument(
        "-p",
        "--prefix",
        type=str,
        default="alg_plotter_v2",
        help="Output prefix of files",
    )
    plotter_2.add_argument(
        "-m",
        "--minimum",
        type=int,
        default=3,
        help="Minimum number of BUSCOs to plot per chromosome (default: 3)",
    )
    plotter_2.add_argument(
        "--palette",
        choices=["categorical", "spectrum", "merianbow", "merianbow4"],
        default="categorical",
        help="Palette to use for coloring BUSCOs (default: categorical), spectrum (turbo), merianbow (original for lepidoptera), merianbow4 (CVD-optimised for lepidoptera)",
    )
    plotter_2.add_argument(
        "--label-threshold",
        type=int,
        default=5,
        help="Minimum number of BUSCOs for an ALG unit appear on the plot (default: 5)",
    )

    # Font Arguments
    font_path = os.path.join(
        os.path.dirname(__file__), "fonts", "OpenSans-Regular.ttf"
    )
    plotter_2.add_argument(
        "--font",
        type=lambda s: file_validator(s, "font"),
        default=font_path,
        help="Path to the font file to use for plotting (default: OpenSans-Regular.ttf)",
    )
    plotter_2.add_argument(
        "-v",
        "--version",
        action="version",
        version=f"ALG-{VERSION}-PLOTTER_V2-2.1.0",
    )

    return parser.parse_args()


def validate_locations(file_path) -> None:
    with open(file_path, "r") as f:
        first_line = f.readline()
        print(first_line.split("\t"))
        if first_line.split("\t") == [
            "buscoID",
            "query_chr",
            "position",
            "assigned_chr",
            "status\n",
        ]:
            raise ValueError(f"File ({file_path}): first line is empty")


def main():
    args = parse_args()

    if args.command == "plotter" and args.file:
        validate_locations(args.file)

    match args.command:
        case "painter":
            alg_painter(args)
        case "plotter":
            alg_plotter_v1(args)
        case "plotter2":
            alg_plotter_v2(args)


if __name__ == "__main__":
    main()
