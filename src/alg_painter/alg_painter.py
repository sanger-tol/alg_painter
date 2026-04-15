import csv
import logging
from collections import Counter
from pathlib import Path

from alg_painter.generics import write_tsv
from alg_painter.ncbi_api import (
    chrom_lengths_with_unloc,
    fetch_sequence_report,
)

logger = logging.getLogger(__name__)


def parse_reference_map(ref_path: Path) -> dict[str, str]:
    """
    Return {BUSCO-ID → ALG element} from the reference full_table.
    """

    # IF LEPIDOPTERAN... SHOULD THIS SIMPLY BE ALG_Z, ALG_1
    alg_set = {"MZ"} | {f"M{i}" for i in range(1, 32)}
    ref_map: dict[str, str] = {}
    with ref_path.open(newline="") as reference_file:
        reader = csv.reader(reference_file, delimiter="\t")
        for row in reader:
            if (
                not row
                or len(row) < 3
                or row[0].startswith("#")
                or row[0].lower().startswith("busco")
            ):
                continue
            alg_unit = row[2].upper().strip()
            if alg_unit in alg_set:
                ref_map[row[0]] = alg_unit
    return ref_map


def parse_busco_table(path: Path) -> tuple[list[tuple], list[str]]:
    """
    Return list of (busco_id, chr, start, stop) tuples and sorted unique chr list.
    This does not filter out any duplicates
    """
    tbl: list[tuple] = []
    chroms: set[str] = set()
    keep_status = {"Complete", "Duplicated"}

    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or len(row) < 5:
                continue
            bid, status, chrom, start, stop = row[:5]
            if status not in keep_status:
                continue
            try:
                start_coord, end_coord = int(start), int(stop)
            except ValueError:
                continue
            tbl.append((bid, chrom, start_coord, end_coord))
            chroms.add(chrom)
    return tbl, sorted(chroms)


def build_location_rows(
    reference_map: dict[str, str], query_table: list[tuple]
) -> list[str]:
    rows = ["buscoID\tquery_chr\tposition\tassigned_chr\tstatus"]
    for busco_id, query_chromosome, start, end in query_table:
        position = (start + end) / 2
        assigned = reference_map.get(busco_id, "NA")
        rows.append(
            f"{busco_id}\t{query_chromosome}\t{position}\t{assigned}\t{assigned}"
        )
    return rows


def painting_main(arguments):
    logger.info(f"Running PAINTER with: {arguments}")

    prefix_location = Path(arguments.prefix)
    logging.info(f"Making sure output dir exists at: {prefix_location}")
    prefix_location.mkdir(parents=True, exist_ok=True)

    output_file_dict = {
        "all_locs": Path(
            f"{prefix_location}/{arguments.prefix}_all_locations.tsv"
        ),
        "chrom_lengths": Path(
            f"{prefix_location}/{arguments.prefix}_chrom_lengths.tsv"
        ),
        "summary": Path(f"{prefix_location}/{arguments.prefix}_summary.tsv"),
    }

    # Parse the BUSCO tables into the needed formats
    reference_map = parse_reference_map(arguments.alg_table)

    query_table, query_chromosomes = parse_busco_table(arguments.query_table)

    # Build the rows of ALL buscos
    all_busco_data = build_location_rows(reference_map, query_table)

    # If ACCESSION is given then call the NCBI API to get chromosome lenths
    if arguments.accession:
        logger.info(
            f"[Painter] Fetching sequence report for accession: {arguments.accession}"
        )
        seq_report = fetch_sequence_report(arguments.accession)
        pairs = chrom_lengths_with_unloc(seq_report)
        length_lines = ["Chrom\tLength_Mb"] + [
            f"{chromosome}\t{basepairs / 1e6:.3f}"
            for chromosome, basepairs in pairs
        ]
        write_tsv(length_lines, output_file_dict["chrom_lengths"])
        chrom_order = [c for c, _ in pairs]
    else:
        chrom_order = query_chromosomes.copy()

    query_chroms_set = {chrom for _, chrom, _, _ in query_table}
    missing = [c for c in chrom_order if c not in query_chroms_set]
    for c in missing:
        all_busco_data.append(f"NA\t{c}\tNA\tNA\tNA")

    logger.info(f"[Painter] Writing tsv to: {output_file_dict['all_locs']}")
    write_tsv(all_busco_data, output_file_dict["all_locs"])

    # Optional summary generation
    if arguments.write_summary:
        logger.info(
            f"[Painter] Writing summary to: {output_file_dict['summary']}"
        )
        counts = Counter(chrom for _, chrom, _, _ in query_table)
        counts.update({c: 0 for c in missing})
        sum_lines = ["query_chr\tbusco_hits"] + [
            f"{c}\t{counts[c]}" for c in chrom_order
        ]
        write_tsv(sum_lines, output_file_dict["summary"])
