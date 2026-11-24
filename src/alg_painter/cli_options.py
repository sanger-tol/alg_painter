from pathlib import Path
from typing import Annotated

import typer

AlgClade = Annotated[
    str,
    typer.Argument(
        ...,
        help="Clade associated with Ancestral Linkage Group",
        show_default=True,
        case_sensitive=True,
    ),
]

QueryFile = Annotated[
    Path,
    typer.Argument(
        ...,
        help="Path to the query file",
        show_default=True,
        dir_okay=False,
        exists=True,
    ),
]

LocationFile = Annotated[
    Path,
    typer.Argument(
        ...,
        help="Path to the all_location.tsv file",
        show_default=True,
        dir_okay=False,
        exists=True,
    ),
]

IndexFile = Annotated[
    Path,
    typer.Argument(
        ...,
        help="Index file for assembly mode",
        show_default=True,
        dir_okay=False,
        exists=True,
    ),
]
