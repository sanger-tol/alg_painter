import json
import logging
import os
import sys
from pathlib import Path

# Suppress matplotlib font debug messages
import matplotlib
import typer

from alg_painter.assign import assign_alg_to_full_table
from alg_painter.cli_options import AlgClade, IndexFile, LocationFile, QueryFile
from alg_painter.enums import AssemblyMode, GraphType, LogLevel
from alg_painter.painter import paint_genomes
from alg_painter.plotter import paint_to_plot

matplotlib.set_loglevel("WARNING")
logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)

app = typer.Typer()


def setup_logging(log_level: LogLevel, log_file: str) -> None:
    if log_level not in LogLevel:
        raise ValueError(f"Invalid log level: {log_level}\nSelect from {LogLevel}")

    logging.basicConfig(
        level=log_level.value,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),  # logs to file
            logging.StreamHandler(),  # logs to console
        ],
    )

    logging.getLogger(__name__)


def load_palette_file(file_path: Path) -> dict[str, dict[str, str]]:
    """Load palette JSON, with fallback to built-in if not found."""
    try:
        with open(file_path, "r") as f:
            return json.load(f)
    except FileNotFoundError:
        # Get just the filename from the path and look for it in the package assets
        file_name = Path(file_path).name
        fall_back = Path(__file__).parent / "assets" / "alg_assignments" / file_name
        logging.error(f"Palette file not found: {file_path}")
        logging.error(f"Falling back to built-in palette: {fall_back}")
        try:
            with open(fall_back, "r") as f:
                return json.load(f)
        except FileNotFoundError:
            sys.exit(f"Palette file not found at {file_path} or {fall_back}")
        except json.JSONDecodeError:
            sys.exit(f"Invalid JSON in fallback palette: {fall_back}")
    except json.JSONDecodeError:
        sys.exit(f"Invalid JSON in palette file: {file_path}")


def get_palette_data(file_path: Path, name: str = "") -> dict[str, str]:
    valid_palette_names = []

    palette_data = load_palette_file(file_path)

    for key, value in palette_data.items():
        valid_palette_names.append(key)

    if name == "":
        return {"options": ", ".join(valid_palette_names)}

    if name not in valid_palette_names:
        sys.exit(f"No palette found for {name} -- valid options are: {', '.join(valid_palette_names)}")

    return palette_data[name]


def check_clades(alg_clade: str, mapping_data: dict) -> None:
    """
    Check if the given clade is valid against the mapping data.
    """
    clades = list(mapping_data.keys())
    clades.remove("CONSTANTS")

    if alg_clade not in clades:
        sys.exit(f"No mapping found for clade: {alg_clade} - try one of: {', '.join(clades)}")


@app.callback()
def app_callback(
    ctx: typer.Context,
    mapping_file: Path = typer.Option(
        default=os.path.join(os.path.dirname(__file__), "assets", "alg_mapping.json"),
        help="Path to ALG mapping file",
        show_default=True,
        dir_okay=False,
        exists=True,
    ),
    prefix: str = typer.Option("alg", help="Prefix for output files", show_default=True),
    log_level: LogLevel = typer.Option(
        LogLevel.INFO,
        help="Logging level",
        show_default=True,
        case_sensitive=False,
    ),
    log_file: str = typer.Option("alg_painter.log", help="Path to log file", show_default=True),
):
    setup_logging(log_level, log_file)
    # If no mapping file provided, use built-in fallback

    ctx.obj = {"mapping_data": mapping_file, "prefix": prefix}


@app.command("paint", help="Paint ALG locations into a genome.")
def painter(
    ctx: typer.Context,
    alg_clade: AlgClade,
    alg_query: QueryFile,
    busco_version: str = typer.Option(None, help="BUSCO version", show_default=True),
    accession: str = typer.Option(default="", help="GCA number", show_default=True),
    ncbi_api_key: str = typer.Option(None, envvar="NCBI_API_KEY", help="NCBI API key", show_default=True),
):
    full_prefix = f"{ctx.obj['prefix']}"

    try:
        mapping_file = Path(ctx.obj.get("mapping_data"))
        if mapping_file.is_file() and alg_clade:
            with open(mapping_file) as f:
                mapping_data = json.load(f)
                check_clades(alg_clade, mapping_data)
                clade_config = mapping_data.get(alg_clade, {})

                ctx.obj.update(
                    {
                        "alg_clade": alg_clade,
                        "alg_file": clade_config["alg_file"],
                        "alg_version": busco_version or clade_config.get("odb"),
                    }
                )
    except (FileNotFoundError, json.JSONDecodeError, ValueError):
        sys.exit(1)

    paint_genomes(
        ctx,
        accession,
        full_prefix,
        alg_query,
        ncbi_api_key,
    )


@app.command("plot", help="Plot ALG locations into a graph.")
def plotter(
    ctx: typer.Context,
    alg_clade: AlgClade,
    location_file: LocationFile,
    index_file: IndexFile,
    palette_name: str = typer.Option("default", help="Name of palette to use"),
    assembly_mode: AssemblyMode = typer.Option(AssemblyMode.AUTO, case_sensitive=False, show_default=True),
    bar_height: int = 12,
    plot_style: GraphType = typer.Option(GraphType.TILES, help="Plot style (tiles, bars, or compact)"),
    alg_file: Path = typer.Option(None, help="Path to ALG/BUSCO mappings", show_default=True, dir_okay=False),
    palette: Path = typer.Option(None, help="Path to palette file", show_default=True, dir_okay=False),
    label_threshold: int = typer.Option(None, help="Lower threshold of units for labeling", show_default=True),
    panel_size: int = typer.Option(None, help="Size of panels in plots", show_default=True),
    has_windowed_labels: bool = typer.Option(None, help="Whether to use windowed labels", show_default=True),
    label_window_min_buscos: int = typer.Option(
        None, help="Minimum number of BUSCOs for windowed labels", show_default=True
    ),
    min_plot_height: int = typer.Option(None, help="Minimum height of the plot", show_default=True),
    column_count: int = typer.Option(None, help="Number of columns in the plot", show_default=True),
    row_height: int = typer.Option(None, help="Height of rows in the plot", show_default=True),
    tile_width_bp: int = typer.Option(None, help="Width of tiles in the plot", show_default=True),
    legend_title: str = typer.Option(None, help="Title of the legend", show_default=True),
    y_spacing: float = typer.Option(None, help="Spacing between rows in the plot", show_default=True),
    label_window_mb: int = typer.Option(None, help="Size of window for windowed labels in Mb", show_default=True),
):
    clade_config = {}
    full_prefix = f"{ctx.obj['prefix']}"

    try:
        mapping_file = Path(ctx.obj.get("mapping_data"))
        if mapping_file.is_file() and alg_clade:
            with open(mapping_file) as f:
                mapping_data = json.load(f)
                check_clades(alg_clade, mapping_data)
                clade_config = mapping_data.get(alg_clade, {})
                constants = mapping_data.get("CONSTANTS", {})
                clade_config.update(constants)
    except (FileNotFoundError, json.JSONDecodeError, ValueError):
        sys.exit(1)

    if clade_config == {}:
        sys.exit(f"No configuration found for alg_clade: {alg_clade}")

    ctx.obj = {
        "alg_clade": alg_clade,
        "palette": palette or clade_config.get("palette", ""),
        "alg_file": alg_file or clade_config.get("alg_file", ""),
        "alg_unit_full": clade_config.get("alg_unit", {}).get("full", ""),
        "label_threshold": label_threshold or clade_config.get("label_window_min_buscos", 5),
        "panel_size": panel_size or clade_config.get("panel_size"),
        "has_windowed_labels": has_windowed_labels or clade_config.get("has_windowed_labels", False),
        "label_window_mb": label_window_mb or clade_config.get("label_window_default_mb"),
        "label_window_min_buscos": label_window_min_buscos or clade_config.get("label_window_min_buscos", 5),
        "column_count": column_count or clade_config.get("column_count", 1),
        "min_plot_height": min_plot_height or clade_config.get("min_plot_height"),
        "row_height": row_height or clade_config.get("row_height"),
        "bar_height": bar_height or clade_config.get("bar_height"),
        "tile_width_bp": tile_width_bp or clade_config.get("tile_width_bp"),
        "legend_title": legend_title or clade_config.get("legend_title"),
        "y_spacing": y_spacing or clade_config.get("y_spacing"),
        "chrom_order": clade_config.get("custom_order"),
    }

    ctx.obj.update({k: v for k, v in constants.items() if k not in ctx.obj})

    palette_data = get_palette_data(ctx.obj["palette"], palette_name)

    logging.debug(ctx.obj)

    paint_to_plot(
        ctx,
        assembly_mode.value,
        full_prefix,
        location_file,
        palette_data,
        index_file,
        3,
        False,
        bar_height,
        2e4,
        plot_style.value,
    )


@app.command("assign")
def assign_ancestral_algs(
    ctx: typer.Context,
    location_file: LocationFile,
    full_table: Path = typer.Option(..., help="Path to busco full table", show_default=True),
    outdir: Path = typer.Option(Path("."), help="Output directory", show_default=True),
):
    full_prefix = f"{ctx.obj['prefix']}_assigned"
    assign_alg_to_full_table(location_file, full_table, full_prefix, outdir)


@app.command("list_the_configs")
def list_the_configs(ctx: typer.Context, alg_clade: AlgClade):
    try:
        mapping_file = Path(ctx.obj.get("mapping_data"))
        print(f"Mapping file found at: {mapping_file}")
        if mapping_file.is_file() and alg_clade:
            with open(mapping_file) as f:
                mapping_data = json.load(f)
                clade_config = mapping_data.get(alg_clade, {})
                if clade_config == {}:
                    print(f"Config for {alg_clade} doesn't exist!")
                else:
                    palette_options = get_palette_data(clade_config.get("palette"))
                    clade_config.update({"palette": palette_options})
                    print(json.dumps(clade_config, indent=4))
    except (FileNotFoundError, json.JSONDecodeError, ValueError):
        sys.exit(1)


if __name__ == "__main__":
    app()
