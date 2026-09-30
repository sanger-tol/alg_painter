from pathlib import Path

import polars as pl


def assign_alg_to_full_table(location_file, full_table, prefix, outdir):
    fulltable_colnames_10 = [
        "buscoID",
        "Status",
        "Sequence",
        "Gene Start",
        "Gene End",
        "Strand",
        "Score",
        "Length",
        "OrthoDB url",
        "Description",
    ]

    fulltable_colnames_8 = ["buscoID", "Status", "Sequence", "Gene Start", "Gene End", "Strand", "Score", "Length"]

    location_data = pl.read_csv(location_file, separator="\t", comment_prefix="#")

    if len(location_data.columns) == 10:
        full_table_colnames = fulltable_colnames_10
    else:
        full_table_colnames = fulltable_colnames_8

    full_table_data = pl.read_csv(
        full_table, separator="\t", comment_prefix="#", has_header=False, new_columns=full_table_colnames
    )

    joint_data = location_data.join(full_table_data, on="buscoID")

    # Build column selection, checking if OrthoDB url exists
    columns_to_select = [
        pl.col("Sequence"),
        pl.col("Gene Start"),
        pl.col("Gene End"),
        pl.col("assigned_alg"),
        pl.col("Score"),
        pl.col("Strand"),
    ]
    if "OrthoDB url" in joint_data.columns:
        columns_to_select.append(pl.col("OrthoDB url"))
    else:
        columns_to_select.append(pl.lit("NA").alias("OrthoDB url"))

    final_data = (
        joint_data.select(*columns_to_select)
        .fill_null("NA")
        .filter(pl.col("Sequence") != "NA")
        .with_columns(
            pl.col("Gene End").cast(pl.Int64),
            pl.col("Gene Start").cast(pl.Int64),
        )
        .with_columns(pl.col("Sequence").str.replace(r":.*", ""))
        .with_columns(
            gene_start_swapped=pl.when(pl.col("Gene Start") < pl.col("Gene End"))
            .then(pl.col("Gene Start"))
            .otherwise(pl.col("Gene End")),
            gene_end_swapped=pl.when(pl.col("Gene Start") < pl.col("Gene End"))
            .then(pl.col("Gene End"))
            .otherwise(pl.col("Gene Start")),
        )
        .select(
            pl.col("Sequence"),
            pl.col("Gene Start"),
            pl.col("Gene End"),
            pl.col("assigned_alg"),
            pl.col("Score"),
            pl.col("Strand"),
            pl.col("OrthoDB url"),
        )
    )

    output = Path(f"{outdir}/{prefix}_assigned_ancestral.tsv")
    output.parent.mkdir(parents=True, exist_ok=True)
    final_data.write_csv(output, include_header=False, separator="\t")
