## Structure of the ALG assets folder

NOTE: `default` for lepidoptera is `merianbow4`

### alg_mapping.json
This file contains the information required for `alg` to generate mappings and plots.

Below is an example taken from the `alg_mapping.json` file, as you can see it contains the data needed to modify output on a per clade basis. A requirement of the Faculty and Genome Note teams at Sanger.

The clade grouping / ALG set ("lepidoptera", in this case), could feasibly be any identifier such as `Psyche` or `ToL`. Clades are simply easier as that's what ALG's are typically based on.

As you can see, we do point to the original article for each ALG set.

Custom ordering of chrom units are also supported, particularly needed for Nematode Nigon Units which follow a roman numeral ordering. This takes the format of a list of the preferred ordering of chrom units ("custom_order": ["I", "II", "III", "IV", "V", "X"]).



```
{
    "CONSTANTS": {
        "label_offset_factor": 0.02,
        "label_padding_factor": 0.06,
        "panel_size": 20,
        "plot_body_width": 22,
        "max_panel_columns": 3,
        "compact_label_treshold": 80,
        "label_wrap": 4,
        "label_window_min_fraction": 0.5
    },
    "lepidoptera": {
        "legend_title": "Merian elements",
        "label_window_default_mb": null,
        "label_window_min_buscos": 5,
        "source": "https://www.nature.com/articles/s41559-024-02329-4",
        "odb": "odb10",
        "palette": "assets/alg_assignments/lepidoptera_merians_odb10.json",
        "alg_file": "assets/alg_assignments/lepidoptera_merians_odb10.tsv",
        "custom_order": null,
        "alg_unit_count": 35,
        "notes": "",
        "y_spacing": 16,
        "bar_height": 0.62,
        "row_height": 2.10,
        "column_count": 2,
        "min_plot_height": 25,
        "tile_width_bp": 50000,
        "has_windowed_labels": false
    }
}
```
