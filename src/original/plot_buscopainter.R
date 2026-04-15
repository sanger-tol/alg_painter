#!/usr/bin/env Rscript

### Original script by  : Charlotte Wright
### Original Repo       : https://github.com/charlottewright/lep_busco_painter/
### Modified by         : Damon-Lee Pointon
### Date of modification:
###   2025-04-24
###     Converting double space to tabs for linting and added
###     version information
###   2025-08-26
###     Conform script to linting standards
###
### Description         :
###   This script generates a plot showing locations of putative
###   merian elements in a genomic assembly


### Load packages
library(optparse)
suppressMessages(library(tidyverse))
suppressMessages(library(scales))

version <- "1.0.1"


### Funcions for making busco paints in R ####

prepare_data <- function(args1) {
    locations <- read_tsv(args1, col_types = cols())

    # format location data
    locations <- locations |> filter(!grepl(":", query_chr))

    # format location data
    locations <- locations |> filter(!grepl(":", assigned_chr))
    locations <- locations |>
        group_by(query_chr) |>
        mutate(length = max(position)) |>
        ungroup()
    locations$start <- 0

    locations
}

prepare_data_with_index <- function(args1, args2) {
    locations <- read_tsv(args1, col_types = cols())
    contig_lengths <- read_tsv(args2, col_names = FALSE, col_types = cols())
    colnames(contig_lengths) <- c("Seq", "length", "offset",
        "linebases", "linewidth"
    )

    # format location data
    locations <- locations |> filter(!grepl(":", query_chr))

    # format location data
    locations <- locations |> filter(!grepl(":", assigned_chr))
    locations <- merge(
        locations, contig_lengths, by.x = "query_chr", by.y = "Seq"
    )
    locations$start <- 0

    locations
}

# minimum of buscos to be present
filter_buscos <- function(locations, minimum) {
    locations_filt <- locations  |>

        # filter df to only keep query_chr with >=3 buscos to remove shrapnel
        group_by(query_chr) |>

        # make a new column reporting number buscos per query_chr
        mutate(n_busco = n()) |>
        ungroup() |>
        filter(n_busco >= minimum)

    locations_filt
}

# Set mapping of Merian element to colour when only plot
set_merian_colour_mapping <- function(location_set) {
    merian_order <- c(
        "MZ", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10",
        "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18", "M19", "M20",
        "M21", "M22", "M23", "M24", "M25", "M26", "M27", "M28", "M29", "M30",
        "M31", "self"
    )
    colour_palette <- append(hue_pal()(32), "grey")
    status_merians <- unique(location_set$status)
    subset_merians <- subset(colour_palette, merian_order %in% status_merians)

    subset_merians
}

busco_paint_theme <- theme(
    legend.position = "right",
    strip.text.x = element_text(
        margin = margin(0, 0, 0, 0, "cm")
    ),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 20)
)

bp_no_facet_labels_theme <- theme(
    legend.position = "right",
    strip.text.x = element_blank(),
    panel.background = element_rect(fill = "white", colour = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 15),
    strip.text.y = element_text(angle = 0),
    strip.background = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "italic", size = 20),
    plot.subtitle = element_text(hjust = 0.5, size = 20)
)

# plot only buscos that have moved - paint by Merians
paint_merians_differences_only <- function(
    spp_df, subset_merians, num_col, title, karyotype
) {
    merian_order <- c(
        "MZ", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10",
        "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18", "M19", "M20",
        "M21", "M22", "M23", "M24", "M25", "M26", "M27", "M28", "M29", "M30",
        "M31", "self"
    )
    spp_df$status_f <- factor(spp_df$status, levels = merian_order)
    chr_levels <- subset(
        spp_df,
        select = c(query_chr, length)
    ) |>
        unique() |>
        arrange(length, decreasing = TRUE)

    chr_levels <- chr_levels$query_chr

    # set chr order as order for plotting
    spp_df$query_chr_f <- factor(spp_df$query_chr, levels = chr_levels)
    sub_title <- paste("n contigs =", karyotype)

    the_plot <- ggplot(data = spp_df) +
        scale_colour_manual(
            values = subset_merians,
            aesthetics = c("colour", "fill")
        ) +
        geom_rect(
            aes(xmin = start, xmax = length, ymax = 0, ymin  = 12),
            colour = "black",
            fill = "white"
        ) +
        geom_rect(
            aes(
                xmin = position - 2e4,
                xmax = position + 2e4,
                ymax = 0,
                ymin = 12,
                fill = status_f
            )
        ) +
        facet_wrap(query_chr_f ~ ., ncol = num_col) +
        guides(scale = "none") +
        xlab("Position (Mb)") +
        scale_x_continuous(labels = function(x) x / 1e6, expand = c(0.005, 1)) +
        scale_y_continuous(breaks = NULL) +
        ggtitle(label = title, subtitle = sub_title)  +
        guides(fill = guide_legend("Merian element"), color = "none")

    # busco_paint_theme
    the_plot
}


# plot only buscos that have moved - paint by species
paint_species_differences_only <- function(spp_df, num_col, title, karyotype) {
    chr_levels <- subset(
        spp_df,
        select = c(query_chr, length)
    ) |>
        unique() |>
        arrange(length, decreasing = TRUE)

    chr_levels <- chr_levels$query_chr
    chr_levels <- chr_levels [! chr_levels %in% "self"]

    # set chr order as order for plotting query chr
    spp_df$query_chr_f <- factor(spp_df$query_chr, levels = chr_levels)
    legend_levels <- unique(spp_df$status)

    # remove "self" from list
    legend_levels <- legend_levels[legend_levels != "self"]

    # then put "self" back in to have it in first position as want "self"
    # to always be painted grey.
    legend_levels <- c("self", legend_levels)
    num_colours <- length(legend_levels)
    col_palette <- hue_pal()(num_colours)
    col_palette[1] <- "grey"

    # set chr order as order for plotting
    spp_df$status_f <- factor(spp_df$status, levels = legend_levels)

    sub_title <- paste("n contigs =", karyotype)
    the_plot <- ggplot(data = spp_df) +
        scale_colour_manual(
            values = col_palette,
            aesthetics = c("fill"),
            breaks = legend_levels
        ) +
        geom_rect(
            aes(xmin = start, xmax = length, ymax = 0, ymin = 12),
            colour = "black",
            fill = "white"
        ) +
        geom_rect(
            aes(
                xmin = position - 2e4,
                xmax = position + 2e4,
                ymax = 0,
                ymin = 12,
                fill = status_f
            )
        ) +
        facet_wrap(
            query_chr_f ~ ., ncol = num_col, strip.position = "right"
        ) +
        guides(scale = "none") +
        xlab("Position (Mb)") +
        scale_x_continuous(labels = function(x) x / 1e6, expand = c(0.005, 1)) +
        scale_y_continuous(breaks = NULL) +
        ggtitle(label = title, subtitle = sub_title)  +
        guides(fill = guide_legend("Query chromosome"), color = "none") +
        busco_paint_theme

    the_plot
}

paint_merians_all <- function(spp_df, num_col, title, karyotype) {
    colour_palette <- append(hue_pal()(32), "grey")
    merian_order <- c(
        "MZ", "M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "M9", "M10",
        "M11", "M12", "M13", "M14", "M15", "M16", "M17", "M18", "M19", "M20",
        "M21", "M22", "M23", "M24", "M25", "M26", "M27", "M28", "M29", "M30",
        "M31", "self"
    )
    spp_df$assigned_chr_f <- factor(spp_df$assigned_chr, levels = merian_order)
    chr_levels <- subset(
        spp_df, select = c(query_chr, length)
    ) |>
        unique() |>
        arrange(length, decreasing = TRUE)

    chr_levels <- chr_levels$query_chr

    # set chr order as order for plotting
    spp_df$query_chr_f <- factor(spp_df$query_chr, levels = chr_levels)
    sub_title <- paste("n contigs =", karyotype)
    the_plot <- ggplot(data = spp_df) +
        scale_colour_manual(
            values = colour_palette,
            aesthetics = c("colour", "fill")
        ) +
        geom_rect(
            aes(xmin = start, xmax = length, ymax = 0, ymin  = 12),
            colour = "black",
            fill = "white"
        ) +
        geom_rect(
            aes(
                xmin = position - 2e4,
                xmax = position + 2e4,
                ymax = 0,
                ymin = 12,
                fill = assigned_chr_f
            )
        ) +
        facet_wrap(
            query_chr_f ~ ., ncol = num_col, strip.position = "right"
        ) +
        guides(scale = "none") +
        xlab("Position (Mb)") +
        scale_x_continuous(labels = function(x) x / 1e6, expand = c(0.005, 1)) +
        scale_y_continuous(breaks = NULL) +
        ggtitle(label = title, subtitle = sub_title)  +
        guides(fill = guide_legend("Merian element"), color = "none") +
        busco_paint_theme

    the_plot
}

# paint all buscos by species
paint_species_all <- function(spp_df, num_col, title, karyotype) {
    chr_levels <- subset(
        spp_df, select = c(query_chr, length)
    ) |>
        unique() |>
        arrange(length, decreasing = TRUE)

    chr_levels <- chr_levels$query_chr
    chr_levels <- chr_levels [! chr_levels %in% "self"]

    # set chr order as order for plotting
    spp_df$query_chr_f <- factor(spp_df$query_chr, levels = chr_levels)
    legend_levels <- subset(spp_df, select = c(assigned_chr)) |> unique()
    legend_levels <- legend_levels$assigned_chr
    num_colours <- length(legend_levels)
    col_palette <- hue_pal()(num_colours)

    # set chr order as order for plotting
    spp_df$assigned_chr_f <- factor(spp_df$assigned_chr, levels = legend_levels)

    sub_title <- paste("n contigs =", karyotype)
    the_plot <- ggplot(data = spp_df) +
        scale_colour_manual(
            values = col_palette,
            aesthetics = c("fill"),
            breaks = legend_levels
        ) +
        geom_rect(
            aes(xmin = start, xmax = length, ymax = 0, ymin = 12),
            colour = "black",
            fill = "white"
        ) +
        geom_rect(
            aes(
                xmin = position - 2e4,
                xmax = position + 2e4,
                ymax = 0,
                ymin = 12,
                fill = assigned_chr_f
            )
        ) +
        facet_wrap(
            query_chr_f ~ .,
            ncol = num_col,
            strip.position = "right"
        ) +
        guides(scale = "none") +
        xlab("Position (Mb)") +
        scale_x_continuous(labels = function(x) x / 1e6, expand = c(0.005, 1)) +
        scale_y_continuous(breaks = NULL) +
        ggtitle(label = title, subtitle = sub_title)  +
        guides(fill = guide_legend("Query chromosome"), color = "none") +
        busco_paint_theme

    the_plot
}

### get args
option_list <- list(
    make_option(
        c("-f", "--file"),
        type = "character",
        default = NULL,
        help = "location.tsv file",
        metavar = "character"
    ),
    make_option(
        c("-p", "--prefix"),
        type = "character",
        default = "Query species",
        help = "prefix for plot title",
        metavar = "character"
    ),
    make_option(
        c("-i", "--index"),
        type = "character",
        default = "False",
        help = "genome index file",
        metavar = "character"
    ),
    make_option(
        c("-m", "--merians"),
        type = "character",
        default = "False",
        help = "Are you comparing a genome to Merian elements",
        metavar = "character"
    ),
    make_option(
        c("-d", "--differences"),
        type = "character",
        default = "False",
        help = "Colour buscos moved from the dominant chromosome",
        metavar = "character"
    ),
    make_option(
        c("-n", "--minimum"),
        type = "integer",
        default = 3,
        help = "minimum number of buscos ",
        metavar = "number"
    ),
    make_option(
        c("-v", "--version"),
        action = "store_true",
        default = FALSE,
        help = "show version and exit"
    )
)

opt_parser <- OptionParser(
    usage = "%prog [options] input",
    option_list = option_list
)
opt <- parse_args(opt_parser)

if (isTRUE(opt$version)) {
    cat(paste(version, "\n"))
    quit(status = 0)
}

locations <- opt$file
prefix <- opt$prefix

index <- opt$index
merians <- opt$merians
differences_only <- opt$differences
minimum <- opt$minimum

# if no index supplied
if (index == "False") {
    location_set <- prepare_data(locations)
    locations_filt <- filter_buscos(location_set, minimum)

    # if index supplied
} else {
    location_set <- prepare_data_with_index(locations, index)
    locations_filt <- filter_buscos(location_set, minimum)
}

# total number of query_chr before filtering
total_contigs <- length(unique(location_set$query_chr))

# number of query_chr after filtering
num_contigs <- as.character(length(unique(locations_filt$query_chr)))
num_removed_contigs <- length(
    unique(location_set$query_chr)
) -
    length(
        unique(locations_filt$query_chr)
    )

print(
    paste(
        "Number of contigs before filtering by number of BUSCOs:",
        total_contigs
    )
)

print(
    paste(
        "Number of contigs removed by filtering :",
        num_removed_contigs
    )
)

print(
    paste(
        "Number of contigs post-filtering:",
        num_contigs
    )
)

# if Merian elements are being used as the comparator
if (merians != "False") {
    subset_merians <- set_merian_colour_mapping(locations_filt)
}

# generate the plot - four possible options based on given arguments to script
# plot only buscos that have moved - paint by Merians

# if comparing two species
if (merians == "False") {

    # if colouring all orthologs
    if (differences_only == "False") {
        p <- paint_species_all(locations_filt, 1, prefix, num_contigs)
        # if only colouring orthologs that have moved
    } else {
        p <- paint_species_differences_only(
            locations_filt, 1, prefix, num_contigs
        )
    }

    # comparing one species to Merian elements
} else {
    # if colouring all orthologs
    if (differences_only == "False") {
        p <- paint_merians_all(locations_filt, 1, prefix, num_contigs)

        # if only colouring orthologs that have moved
    } else {
        if (length(locations_filt$query_chr) < 100) {
            p <- paint_merians_differences_only(
                locations_filt, subset_merians, 1, prefix, num_contigs
            )
            p <- p + busco_paint_theme
        } else {
            p <- paint_merians_differences_only(
                locations_filt, subset_merians, 3, prefix, num_contigs
            )
            p <- p + bp_no_facet_labels_theme
        }
    }
}


ggsave(
    paste(
        as.character(opt$file),
        "_buscopainter.png",
        sep = ""
    ),
    plot = p,
    width = 15,
    height = 30,
    units = "cm",
    device = "png"
)

pdf(NULL)
ggsave(
    paste(
        as.character(opt$file),
        "_buscopainter.pdf",
        sep = ""
    ),
    plot = p,
    width = 15,
    height = 30,
    units = "cm",
    device = "pdf"
)

ggsave(
    paste(
        as.character(opt$file),
        "_buscopainter.png",
        sep = ""
    ),
    plot = p,
    width = 15,
    height = 30,
    units = "cm",
    device = "png"
)

pdf(NULL)
ggsave(
    paste(
        as.character(opt$file),
        "_buscopainter.pdf",
        sep = ""
    ),
    plot = p,
    width = 15,
    height = 30,
    units = "cm",
    device = "pdf"
)
