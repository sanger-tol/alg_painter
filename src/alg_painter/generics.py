from pathlib import Path


def write_tsv(lines: list[str], path: Path) -> None:
    """
    Writes a list of strings to a TSV file at the given path.
    """
    path.write_text("\n".join(lines) + "\n")
