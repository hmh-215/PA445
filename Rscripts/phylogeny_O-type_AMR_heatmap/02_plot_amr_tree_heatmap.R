#!/usr/bin/env Rscript

# 02_plot_amr_tree_heatmap.R  (Strategy 2, step 2)
# -----------------------------------------------------------------------------
# Build a phylogeny-aligned AMR heatmap with geography/serotype annotations
# from a merged table (produced by 01_build_merged_table.R) and a tree.
#
# Usage:
#   Rscript 02_plot_amr_tree_heatmap.R \
#       --merged merged_tables.csv \
#       --tree   PA445.phylogeny.treefile \
#       --out-prefix AMR_tree_heatmap \
#       --edge-scale 3 \
#       --heatmap-width 12 \
#       --AMR-hmap-contrast
# -----------------------------------------------------------------------------

# ---- argument parsing --------------------------------------------------------
parse_args <- function(args) {
  opt <- list(
    merged = NA, 
    tree = "PA445.phylogeny.treefile",
    out_prefix = "AMR_tree_heatmap", 
    edge_scale = 3,
    heatmap_width = NA,              # Manual control over AMR heatmap width
    amr_contrast = FALSE             # Toggle binary high-contrast mode
  )
  i <- 1
  while (i <= length(args)) {
    a <- args[i]; take <- function() { i <<- i + 1; args[i] }
    if      (a %in% c("--merged", "-m"))     opt$merged <- take()
    else if (a %in% c("--tree", "-t"))       opt$tree <- take()
    else if (a %in% c("--out-prefix", "-o")) opt$out_prefix <- take()
    else if (a == "--edge-scale")            opt$edge_scale <- as.numeric(take())
    else if (a == "--heatmap-width")         opt$heatmap_width <- as.numeric(take())
    else if (a == "--AMR-hmap-contrast")     opt$amr_contrast <- TRUE
    else if (a %in% c("--help", "-h"))       { cat("See header of this script for usage.\n"); quit(save="no") }
    else stop("Unknown argument: ", a)
    i <- i + 1
  }
  opt
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

install_if_missing <- function(packages) {
  installed <- rownames(installed.packages())
  missing <- setdiff(packages, installed)
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cran.rstudio.com")
  }
}

required_pkgs <- c("tidyverse", "ggtree", "treeio", "RColorBrewer", "ape", "ggnewscale")
install_if_missing(required_pkgs)

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(RColorBrewer)
  library(ape)
  library(ggnewscale)
})

# Input files
merged_file <- if (!is.na(opt$merged)) {
  if (!file.exists(opt$merged)) stop("Merged CSV not found: ", opt$merged)
  opt$merged
} else if (file.exists("merged_tables.csv")) {
  "merged_tables.csv"
} else if (file.exists("merged_table.csv")) {
  "merged_table.csv"
} else {
  stop("No merged_table CSV file found. Please provide merged_tables.csv or merged_table.csv.")
}
tree_file <- opt$tree

# Output files
out_pdf <- paste0(opt$out_prefix, ".pdf")
out_png <- paste0(opt$out_prefix, ".png")

# Read merged table, using its sample/metadata/AMR columns directly
merged_raw <- read_csv(merged_file, show_col_types = FALSE, na = c("", "NA", "missing", "N/A"))

required_meta_cols <- c("accession", "strain", "geo_location", "O-type")
if (!all(required_meta_cols %in% names(merged_raw))) {
  stop("Merged table is missing required columns: ", paste(setdiff(required_meta_cols, names(merged_raw)), collapse = ", "))
}

merged_raw <- merged_raw %>%
  mutate(
    strain = as.character(strain),
    accession = as.character(accession),
    sample = if_else(!is.na(strain) & strain != "", strain, accession),
    O_type = `O-type`,
    geo_location = na_if(geo_location, "N/A"),
    geo_location = na_if(geo_location, "missing"),
    geo_location = na_if(geo_location, "NA"),
    country = if_else(is.na(geo_location), NA_character_, str_trim(str_extract(geo_location, "^[^:]+"))),
    country = replace_na(country, "unknown"),
    # Standardize explicit variations of "Others" to title case
    country = if_else(tolower(country) == "others", "Others", country),
    O_type = replace_na(O_type, "unknown"),
    host = "unknown"
  )

amr_cols <- setdiff(names(merged_raw), c("accession", "strain", "geo_location", "O-type", "O_type", "sample", "country", "host"))
if (length(amr_cols) == 0) {
  stop("No AMR columns found in merged table.")
}

amr_wide <- merged_raw %>%
  select(sample, all_of(amr_cols)) %>%
  distinct(sample, .keep_all = TRUE) %>%
  column_to_rownames("sample")

norm_name <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  x <- gsub("^GCF_|^GCA_", "", x)
  x <- gsub("[: ,()]+", "_", x)
  x <- gsub("[^A-Za-z0-9_\\-\\.]+", "", x)
  x <- gsub("\\n", "", x)
  x <- trimws(x)
  x
}

rownames(amr_wide) <- norm_name(rownames(amr_wide))
merged_meta <- merged_raw %>%
  distinct(sample, .keep_all = TRUE) %>%
  mutate(sample = norm_name(sample)) %>%
  arrange(sample)

# Keep one unique sample per normalized name
merged_meta <- merged_meta %>% distinct(sample, .keep_all = TRUE)

# Remove duplicate normalized AMR rows if present
amr_wide <- amr_wide[!duplicated(rownames(amr_wide)), , drop = FALSE]

# Read tree and align samples
tree <- read.tree(tree_file)

# --- Dynamic Scale Bar Calculation ---
tree_height_raw <- max(ape::node.depth.edgelength(tree), na.rm = TRUE)
suggested_width <- tree_height_raw * 0.225 
exponent <- floor(log10(suggested_width))
scalebar_unscaled <- round(suggested_width / (10^exponent)) * (10^exponent)
if (scalebar_unscaled == 0) scalebar_unscaled <- 10^exponent

if (!is.null(tree$edge.length) && any(!is.na(tree$edge.length))) {
  tree$edge.length <- tree$edge.length * opt$edge_scale
}
scalebar_width <- scalebar_unscaled * opt$edge_scale
# -------------------------------------

tree$tip.label <- norm_name(tree$tip.label)

matched_tips <- intersect(tree$tip.label, merged_meta$sample)
if (length(matched_tips) == 0) {
  stop("No tree tips could be matched to merged table sample names.")
}

if (length(matched_tips) < length(tree$tip.label)) {
  removed_tips <- setdiff(tree$tip.label, matched_tips)
  message("Removing ", length(removed_tips), " tree tips not present in merged_table: ", paste(head(removed_tips, 10), collapse = ", "))
}

tree <- keep.tip(tree, matched_tips)

missing_in_amr <- setdiff(tree$tip.label, rownames(amr_wide))
if (length(missing_in_amr) > 0) {
  warning("The following tree tips are missing from the merged AMR matrix: ", paste(missing_in_amr, collapse = ", "))
}

amr_wide <- amr_wide[intersect(tree$tip.label, rownames(amr_wide)), , drop = FALSE]

# ---- Process AMR matrices based on --AMR-hmap-contrast state -----------------
if (opt$amr_contrast) {
  # Strict/Perfect hits (2) -> "1" (Yes), Loose/No hits (0,1) -> "0" (No)
  amr_wide[] <- lapply(amr_wide, function(x) {
    factor(if_else(as.character(x) == "2", "1", "0"), levels = c("0", "1"))
  })
} else {
  # Default original mapping
  amr_wide[] <- lapply(amr_wide, function(x) factor(as.character(x), levels = c("0", "1", "2")))
}

# Build metadata for plotting from merged table
meta_plot <- merged_meta %>%
  filter(sample %in% tree$tip.label) %>%
  distinct(sample, .keep_all = TRUE) %>%
  mutate(
    sample = factor(sample, levels = tree$tip.label),
    O_type = factor(O_type),
    host = factor(host)
  )

if (nrow(meta_plot) == 0) {
  stop("No metadata rows could be matched to the tree tips. Check sample naming consistency.")
}

# ---- Order Countries explicitly with "Others" and "unknown" at bottom -------
all_countries <- sort(unique(as.character(meta_plot$country)))
base_countries <- setdiff(all_countries, c("Others", "unknown"))

ordered_countries <- base_countries
if ("Others" %in% all_countries)  ordered_countries <- c(ordered_countries, "Others")
if ("unknown" %in% all_countries) ordered_countries <- c(ordered_countries, "unknown")

meta_plot$country <- factor(meta_plot$country, levels = ordered_countries)
meta_plot$country <- forcats::fct_explicit_na(meta_plot$country, na_level = "unknown")

# High-contrast country palette configuration
contrast_base <- c(brewer.pal(8, "Dark2"), brewer.pal(7, "Set1"))
known_countries <- setdiff(ordered_countries, c("Others", "unknown"))

known_colors <- if (length(known_countries) == 0) {
  character(0)
} else if (length(known_countries) <= length(contrast_base)) {
  contrast_base[seq_along(known_countries)]
} else {
  colorRampPalette(contrast_base)(length(known_countries))
}

country_colors <- set_names(known_colors, known_countries)
if ("Others" %in% ordered_countries)  country_colors <- c(country_colors, "Others" = "#D3D3D3") # Light Grey
if ("unknown" %in% ordered_countries) country_colors <- c(country_colors, "unknown" = "#AAAAAA") # Neutral Grey

# Create O-type annotation dataframe
o_types_unique <- sort(unique(as.character(meta_plot$O_type)))
otype_mapping <- set_names(seq_along(o_types_unique), o_types_unique)

otype_df <- meta_plot %>%
  arrange(factor(sample, levels = tree$tip.label)) %>%
  select(sample, O_type) %>%
  mutate(O_type_num = otype_mapping[as.character(O_type)]) %>%
  select(sample, O_type_num) %>%
  column_to_rownames("sample") %>%
  as.data.frame()

otype_df[] <- lapply(otype_df, function(x) factor(as.character(x), levels = as.character(seq_along(o_types_unique))))

n_tips <- length(tree$tip.label)

tree_plot <- ggtree(tree, size = 0.35) %<+% meta_plot +
  theme_tree() +
  # NOTE: geom_treescale() always *labels* the bar with the raw `width` value
  # it is given. Since tree$edge.length was multiplied by opt$edge_scale above
  # (to spread branches visually), passing scalebar_width there would draw a
  # correctly-sized bar but print a label inflated by edge_scale rather than
  # the true branch length from the .treefile. Draw the segment/label
  # manually so the two can be decoupled: geometric length = scalebar_width
  # (matches the scaled tree coordinates), printed label = scalebar_unscaled
  # (the true, unscaled substitutions/site value).
  annotate("segment",
           x = 0, xend = scalebar_width,
           y = n_tips + 1, yend = n_tips + 1,
           color = "black", linewidth = 0.5) +
  annotate("text",
           x = scalebar_width / 2, y = n_tips + 1 + 0.1,
           label = format(scalebar_unscaled, scientific = FALSE, trim = TRUE),
           size = 3, color = "black", vjust = 0) +
  theme(legend.position = "right",
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8))

otype_palette <- if (length(o_types_unique) <= 12) {
  brewer.pal(max(3, length(o_types_unique)), "Paired")
} else {
  colorRampPalette(brewer.pal(12, "Paired"))(length(o_types_unique))
}

# Use manual override if --heatmap-width parameter is supplied
heatmap_width <- if (!is.na(opt$heatmap_width)) opt$heatmap_width else max(6, min(18, ncol(amr_wide) / 4))

otype_col_width <- (heatmap_width / ncol(amr_wide)) * 1
label_offset <- 0.45

gheat <- gheatmap(tree_plot,
                  otype_df,
                  offset = label_offset,
                  width = otype_col_width,
                  colnames = FALSE,
                  color = NULL) +
  scale_fill_manual(
    name = "O-type",
    values = set_names(otype_palette[seq_along(o_types_unique)], seq_along(o_types_unique)),
    breaks = seq_along(o_types_unique),
    labels = o_types_unique,
    na.value = "#CCCCCC",
    guide = guide_legend(override.aes = list(shape = NA, color = NA))
  )

gheat <- gheat + ggnewscale::new_scale_fill()

tree_x_range <- diff(range(tree_plot$data$x, na.rm = TRUE))
amr_offset   <- label_offset + otype_col_width * tree_x_range

# Create flexible fill details for standard vs high-contrast profiles
amr_scale <- if (opt$amr_contrast) {
  scale_fill_manual(name = "AMR presence",
                    values = c("0" = "#E8E8E8", "1" = "#000000"),
                    labels = c("0" = "no", "1" = "yes"),
                    na.value = "#F5F5F5",
                    guide = guide_legend(override.aes = list(shape = NA, color = NA)))
} else {
  scale_fill_manual(name = "AMR presence",
                    values = c("0" = "#E8E8E8", "1" = "#808080", "2" = "#000000"),
                    labels = c("0" = "no hit", "1" = "loose hit", "2" = "strict hit"),
                    na.value = "#F5F5F5",
                    guide = guide_legend(override.aes = list(shape = NA, color = NA)))
}

gheat <- gheatmap(gheat,
                  amr_wide,
                  offset = amr_offset,
                  width = heatmap_width,
                  colnames_angle = 90,
                  colnames_offset_y = -0.5,
                  font.size = 3.5,
                  hjust = 1,
                  color = NULL) +
  amr_scale +
  theme(legend.position = "right",
        legend.key.height = unit(0.4, "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

tip_coords <- gheat$data %>%
  filter(isTip) %>%
  select(node, x, y, label) %>%
  left_join(meta_plot %>% select(sample, country, strain, O_type),
            by = c("label" = "sample")) %>%
  mutate(tip_label = if_else(!is.na(strain) & strain != "", strain, label))

final_plot <- gheat +
  geom_point(data = tip_coords,
             aes(x = x, y = y, color = country),
             size = 2.8,
             inherit.aes = FALSE,
             show.legend = TRUE) +
  geom_text(data  = tip_coords,
            aes(x = x, y = y, label = tip_label),
            hjust       = 0,
            nudge_x     = 0.02,
            size        = 3.5,
            inherit.aes = FALSE) +
  scale_color_manual(name = "Geography (country)", values = country_colors, na.value = "grey70") +
  guides(color = guide_legend(override.aes = list(size = 4, stroke = 0))) +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(t = 0.5, r = 0.5, b = 10, l = 0.5, unit = "cm"))

ggsave(out_pdf, final_plot, width = 22, height = 20, units = "in", limitsize = FALSE)
ggsave(out_png, final_plot, width = 22, height = 20, units = "in", dpi = 300, limitsize = FALSE)

message("Figure saved to: ", out_pdf, " and ", out_png)