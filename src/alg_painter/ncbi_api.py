import logging
import os

import requests

API_KEY = os.getenv("NCBI_API_KEY")

logger = logging.getLogger(__name__)


def fetch_sequence_report(accession: str) -> list[dict]:
    url = f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accession}/sequence_reports"
    headers = {"accept": "application/json", "User-Agent": "buscopainter"}
    params = {}
    if API_KEY:
        params["api_key"] = API_KEY
    logger.info(f"[Fetch seq report] Fetching chromosome info from NCBI for {accession}...")
    response = requests.get(url, headers=headers, params=params, timeout=60)
    response.raise_for_status()
    return response.json().get("sequence_report", {}).get("records") or response.json().get(
        "reports", []
    )


def chrom_lengths_with_unloc(records):
    """
    Return (accession, total_bp) where accession is the main assembled-molecule
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
