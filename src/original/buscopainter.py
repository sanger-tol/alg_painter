"""BUSCO → Merian painter

Mandatory output
----------------
* `all_location.tsv` every Complete **and** Duplicated BUSCO hit

Optional outputs (when switches present)
---------------------------------------
* `chrom_lengths.tsv` main‑chromosome + unloc scaffold lengths (`--accession`)
* `summary.tsv`       simple BUSCO count per chromosome (`--write_summary`)

If `--prefix` ends with a path separator **or** is an existing directory it is
used as an output directory and the fixed filenames above are written inside
it. Otherwise the prefix is treated as a stem and the files are written as
`<stem>_all_location.tsv`, etc.
"""

from __future__ import annotations

import argparse
import csv
import os
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple

import requests
from Bio import Entrez

API_KEY = os.getenv("NCBI_API_KEY")


# ───────────────────────────────────── BUSCO full_table parsing ──
# Do not filter out duplicates here, as we want to include them in the plot.
def parse_busco_table(path: Path) -> Tuple[List[Tuple], List[str]]:
    """Return list of (busco_id, chr, start, stop) tuples and sorted unique chr list."""
    tbl: List[Tuple] = []
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
                s_i, e_i = int(start), int(stop)
            except ValueError:
                continue
            tbl.append((bid, chrom, s_i, e_i))  # Append instead of overwrite
            chroms.add(chrom)
    return tbl, sorted(chroms)


# ───────────────────────────── BUSCO → Merian lookup table ──────


def build_ref_map(ref_path: Path) -> Dict[str, str]:
    """Return {BUSCO-ID → Merian element} from the reference full_table."""
    mset = {"MZ"} | {f"M{i}" for i in range(1, 32)}
    ref_map: Dict[str, str] = {}
    with ref_path.open(newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            if not row or row[0].startswith("#") or row[0].lower().startswith("busco"):
                continue
            merian = row[2].upper().strip()
            if merian in mset:
                ref_map[row[0]] = merian
    return ref_map


# ───────────────────────────────── BUSCO → Merian rows ──────────


def build_location_rows(ref_map: Dict[str, str], qry_tbl: List[Tuple]) -> List[str]:
    rows = ["buscoID\tquery_chr\tposition\tassigned_chr\tstatus"]
    for bid, qchr, s, e in qry_tbl:
        pos = (s + e) / 2
        assigned = ref_map.get(bid, "NA")
        rows.append(f"{bid}\t{qchr}\t{pos}\t{assigned}\t{assigned}")
    return rows


# ─────────────────────────────── chromosome lengths ─────────────


def fetch_sequence_report(accession: str) -> List[dict]:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/sequence_reports"
    headers = {"accept": "application/json", "User-Agent": "buscopainter"}
    params = {}
    if API_KEY:
        params["api_key"] = API_KEY
    print(f"  → Fetching chromosome info from NCBI for {accession}...")
    r = requests.get(url, headers=headers, params=params, timeout=60)
    r.raise_for_status()
    return r.json().get("sequence_report", {}).get("records") or r.json().get("reports", [])


def chrom_lengths_with_unloc(records):
    """Return (accession, total_bp) where accession is the main assembled-molecule
    GenBank accession for each chromosome.  Length sums main + its unloc scaffolds.
    """
    # first map chr_name → main accession
    main_acc = {}
    for rec in records:
        if (
            rec.get("role") == "assembled-molecule"
            and rec.get("assigned_molecule_location_type") == "Chromosome"
        ):
            main_acc[rec["chr_name"]] = rec["genbank_accession"]

    bp_tot: Dict[str, int] = {acc: 0 for acc in main_acc.values()}

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


# ───────────────────────────────────────── write helper ─────────


def write_tsv(lines: List[str], path: Path) -> None:
    path.write_text("\n".join(lines) + "\n")


# ───────────────────────────────────────── main ─────────────────


def main() -> None:
    ap = argparse.ArgumentParser(description="BUSCO → Merian mapper + length table")
    ap.add_argument("--reference_table", "-r", type=Path, required=True)
    ap.add_argument("--query_table", "-q", type=Path, required=True)
    ap.add_argument("--prefix", "-p", type=Path, default="buscopainter")
    ap.add_argument("--accession", "-a", help="assembly accession (optional)")
    ap.add_argument("--write_summary", action="store_true")
    args = ap.parse_args()

    # Derive output paths
    pref = Path(args.prefix)
    if str(pref).endswith(("/", "\\")) or pref.is_dir():
        out_dir = pref
        out_all = out_dir / "all_location.tsv"
        out_len = out_dir / "chrom_lengths.tsv"
        out_sum = out_dir / "summary.tsv"
    else:
        out_dir = pref.parent
        stem = pref.name
        out_all = out_dir / f"{stem}_all_location.tsv"
        out_len = out_dir / f"{stem}_chrom_lengths.tsv"
        out_sum = out_dir / f"{stem}_summary.tsv"

    # 1. Parse BUSCO tables
    ref_map = build_ref_map(args.reference_table)
    qry_tbl, qry_chrs = parse_busco_table(args.query_table)

    # 2. Build rows (all BUSCOs)
    all_rows = build_location_rows(ref_map, qry_tbl)

    # 3. Optional chromosome lengths + placeholder logic
    chrom_order = qry_chrs.copy()
    wrote_len = False
    if args.accession:
        pairs = chrom_lengths_with_unloc(fetch_sequence_report(args.accession))
        length_lines = ["Chrom\tLength_Mb"] + [f"{c}\t{bp / 1e6:.3f}" for c, bp in pairs]
        write_tsv(length_lines, out_len)
        wrote_len = True
        chrom_order = [c for c, _ in pairs]

    # add placeholder rows
    missing = [c for c in chrom_order if c not in qry_chrs]
    for c in missing:
        all_rows.append(f"NA\t{c}\tNA\tNA\tNA")

    # 4. Write all_location.tsv
    write_tsv(all_rows, out_all)

    # 5. Optional summary
    wrote_sum = False
    if args.write_summary:
        counts = Counter(ch for ch, *_ in qry_tbl.values())
        counts.update({c: 0 for c in missing})
        sum_lines = ["query_chr\tbusco_hits"] + [f"{c}\t{counts[c]}" for c in chrom_order]
        write_tsv(sum_lines, out_sum)
        wrote_sum = True

    # 6. Log
    print("[+] Outputs written:")
    print("[+]\t" + str(out_all))
    if wrote_len:
        print("[+]\t" + str(out_len))
    if wrote_sum:
        print("[+]\t" + str(out_sum))


if __name__ == "__main__":
    main()
