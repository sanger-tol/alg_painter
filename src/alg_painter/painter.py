import csv
import logging
import os
import sys
from collections import Counter
from pathlib import Path

import requests
import typer

logger = logging.getLogger(__name__)


def setup_output_dir(prefix: str) -> Path:
    output_dir = Path(prefix)
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


def parse_reference_map(ref_path: Path) -> dict[str, str]:
    """
    Return {BUSCO-ID → ALG element} from the reference full_table.
    """
    ref_map: dict[str, str] = {}
    with open(os.path.join(os.path.dirname(__file__), ref_path), newline="") as reference_file:
        reader = csv.reader(reference_file, delimiter="\t")
        for row in reader:
            if not row or len(row) < 2 or row[0].startswith("#"):
                continue
            alg_unit = row[1].upper().strip()
            ref_map[row[0]] = alg_unit
    return ref_map


def parse_busco_table(path: Path) -> tuple[list[tuple], list[str], list[list[str]]]:
    """
    Return list of (busco_id, chr, start, stop) tuples and sorted unique chr list.
    This does not filter out any duplicates
    """
    tbl: list[tuple] = []
    chroms: set[str] = set()
    header: list[list[str]] = []
    keep_status = {"Complete", "Duplicated"}

    with path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if row[0].startswith("#"):
                header.append(row)
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
    return tbl, sorted(chroms), header


def build_location_rows(reference_map: dict[str, str], query_table: list[tuple]) -> list[str]:
    rows = ["buscoID\tquery_chr\tposition\tassigned_alg\tstatus"]
    for busco_id, query_chromosome, start, end in query_table:
        position = (start + end) / 2
        assigned = reference_map.get(busco_id, "NA")
        rows.append(f"{busco_id}\t{query_chromosome}\t{position}\t{assigned}\t{assigned}")
    return rows


def fetch_sequence_report(accession: str, ncbi_api_key: str) -> list[dict]:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/sequence_reports"
    headers = {"accept": "application/json", "User-Agent": "buscopainter"}
    params = {}
    if ncbi_api_key:
        params["api_key"] = ncbi_api_key
    logger.info(f"[Fetch seq report] Fetching chromosome info from NCBI for {accession}...")
    response = requests.get(url, headers=headers, params=params, timeout=60)
    response.raise_for_status()
    return response.json().get("sequence_report", {}).get("records") or response.json().get("reports", [])


def write_tsv(lines: list[str], path: Path) -> None:
    """
    Writes a list of strings to a TSV file at the given path.
    """
    path.write_text("\n".join(lines) + "\n")


def chrom_lengths_with_unloc(records):
    """
    Return (accession, total_bp) where accession is the main assembled-molecule
    GenBank accession for each chromosome.  Length sums main + its unloc scaffolds.
    """
    # first map chr_name → main accession
    main_acc = {}
    for rec in records:
        if rec.get("role") == "assembled-molecule" and rec.get("assigned_molecule_location_type") == "Chromosome":
            main_acc[rec["chr_name"]] = rec["genbank_accession"]

    bp_tot: dict[str, int] = {acc: 0 for acc in main_acc.values()}

    for rec in records:
        role = rec.get("role")
        loc = rec.get("assigned_molecule_location_type", "")
        if role == "assembled-molecule" and loc == "Chromosome":
            acc = rec["genbank_accession"]
            bp_tot[acc] += int(rec.get("length", 0))
        elif role == "unlocalized-scaffold":
            parent = rec.get("chr_name")
            acc = main_acc.get(parent)
            if acc:
                bp_tot[acc] += int(rec.get("length", 0))

    return sorted(bp_tot.items(), key=lambda x: -x[1])


def fetch_from_ncbi(accession: str, chrom_lengths_output: Path, ncbi_api_key: str) -> list[str]:
    """
    When an accesstion is provided, fetch the sequence report from NCBI and
    write the chrom lengths to the output file.
    """
    seq_report = fetch_sequence_report(accession, ncbi_api_key)
    pairs: dict = chrom_lengths_with_unloc(seq_report)
    length_lines: list[str] = ["Chrom\tLength_Mb"] + [
        f"{chromosome}\t{basepairs / 1e6:.3f}" for chromosome, basepairs in pairs
    ]
    write_tsv(length_lines, chrom_lengths_output)
    return [c for c, _ in pairs]


def paint_genomes(
    ctx: typer.Context,
    accession: str,
    prefix: str,
    alg_query: Path,
    ncbi_api_key: str,
):
    output_dir = setup_output_dir(prefix)
    logger.info(f"Output dir: {output_dir}")

    output_file_dict = {
        "all_locs": Path(f"{output_dir}/{prefix}_paint_all_locations.tsv"),
        "chrom_lengths": Path(f"{output_dir}/{prefix}_paint_chrom_lengths.tsv"),
        "summary": Path(f"{output_dir}/{prefix}_paint_summary.tsv"),
    }

    reference_map = parse_reference_map(ctx.obj["alg_file"])

    query_table, query_chromosomes, header = parse_busco_table(alg_query)

    valid = True if ctx.obj["alg_version"] in header[1][0] else False
    logger.info(f"BUSCO version requested: {ctx.obj['alg_version']}")
    logger.info(f"Does this match the full-table? ({valid})")

    if not valid:
        sys.exit("Requested BUSCO version doesn't match the data!")

    all_busco_data = build_location_rows(reference_map, query_table)

    if accession:
        logger.info(f"[Painter] Fetching sequence report for accession: {accession}")
        chrom_order: list[str] = fetch_from_ncbi(accession, output_file_dict["chrom_lengths"], ncbi_api_key)
    else:
        logger.info(f"[Painter] Using query chromosomes as sequence report: {query_chromosomes}")
        chrom_order = query_chromosomes.copy()

    query_chroms_set = {chrom for _, chrom, _, _ in query_table}
    missing = [c for c in chrom_order if c not in query_chroms_set]

    for c in missing:
        all_busco_data.append(f"NA\t{c}\tNA\tNA\tNA")

    logger.info(f"[Painter] Writing tsv to: {output_file_dict['all_locs']}")
    write_tsv(all_busco_data, output_file_dict["all_locs"])

    logger.info(f"[Painter] Writing summary to: {output_file_dict['summary']}")
    counts: Counter[str] = Counter(chrom for _, chrom, _, _ in query_table)
    counts.update({c: 0 for c in missing})
    sum_lines = ["query_chr\tbusco_hits"] + [f"{c}\t{counts[c]}" for c in chrom_order]
    write_tsv(sum_lines, output_file_dict["summary"])
