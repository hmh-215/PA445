#!/usr/bin/env Rscript

# parse_merge_utils.R
# -----------------------------------------------------------------------------
# Shared helpers for building the merged AMR/metadata table that
# plot_amr_tree_heatmap.R consumes.
#
# Inputs supported (all optional except at least one of --rgi / --accession):
#   * NCBI accession table  (accession, strain, geo_location, O-type columns)
#   * RGI heatmap matrix    (genes x samples or samples x genes, 0/1/2 values)
#   * Pasty serotyping      (O-type per sample; overrides metadata O-type)
#   * Newick phylogeny tree (adds an `in_tree` TRUE/FALSE column)
#
# Deliberately written in BASE R only (no tidyverse): the parse/merge step is
# lightweight, has no plotting dependencies, runs anywhere, and starts instantly.
# The heavy packages (ggtree/treeio/ggplot2) live only in the plotting script.
#
# These functions are sourced by:
#   * 01_build_merged_table.R   (Strategy 2, step 1)
#   * run_all.R                 (Strategy 1, one-shot convenience)
# -----------------------------------------------------------------------------

# ---- column auto-detection --------------------------------------------------

# Case-/punctuation-insensitive lookup of the first matching candidate name.
# Returns the *actual* column name from `cols`, or NA_character_ if none match.
.detect_col <- function(cols, candidates) {
  norm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  cols_norm <- norm(cols)
  for (cand in candidates) {
    hit <- which(cols_norm == norm(cand))
    if (length(hit) > 0) return(cols[hit[1]])
  }
  NA_character_
}

# Default aliases for the four metadata columns we need. Extend freely.
.DEFAULT_COL_ALIASES <- list(
  accession    = c("accession", "assembly", "assembly_accession", "acc",
                   "genome", "genome_accession", "sample", "id", "isolate_id"),
  strain       = c("strain", "isolate", "isolate_name", "name", "strain_name"),
  geo_location = c("geo_location", "geo_loc_name", "geographic_location",
                   "geo", "location", "country", "geo_loc"),
  `O-type`     = c("o-type", "o_type", "otype", "serotype", "serogroup",
                   "o_serotype", "oantigen", "o_antigen")
)

# ---- value helpers ----------------------------------------------------------

# NCBI geo_location is "<country>: <city/institute>". Take the country part.
# Uninformative tokens collapse to NA (the plotting script maps NA -> "unknown").
parse_country <- function(geo) {
  geo <- as.character(geo)
  bad <- c("", "na", "n/a", "missing", "not collected", "not applicable",
           "unknown", "none", "not provided", "restricted access")
  country <- sub(":.*$", "", geo)          # everything before the first colon
  country <- trimws(country)
  country[tolower(trimws(geo))    %in% bad] <- NA_character_
  country[tolower(trimws(country)) %in% bad] <- NA_character_
  country
}

# Normalise an identifier for cross-table matching.
# Strips common genome-file suffixes, tool-appended suffixes (_rgi, .pasty),
# and non-alphanumerics, then lowercases. This makes matching work across:
#   * RGI heatmap column names  (LYT4_rgi       -> lyt4)
#   * Pasty sample names        (LYT4.pasty     -> lyt4)
#   * NCBI accession strings    (GCF_012971705.1 -> gcf012971705)
#   * Strain names              (Pae1255-NDM1   -> pae1255ndm1)
# NOTE: GCF_/GCA_ prefixes are KEPT so the two assembly versions stay distinct.
# Returns "" for empty/NA input.
normalize_id <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- sub("(\\.fna|\\.fasta|\\.fa|\\.gbk|\\.gbff)$", "", x, ignore.case = TRUE)
  x <- sub("_genomic$", "", x, ignore.case = TRUE)
  x <- sub("_rgi$",     "", x, ignore.case = TRUE)   # RGI heatmap column suffix
  x <- sub("\\.pasty$", "", x, ignore.case = TRUE)   # Pasty output sample suffix
  x <- gsub("[^A-Za-z0-9]+", "", x)                  # punctuation / spaces
  tolower(trimws(x))
}

# Coerce RGI cell values to the "0"/"1"/"2" scheme the plot expects.
#   0 = no hit, 1 = loose hit, 2 = strict hit   (per this project's RGI heatmap)
# If the input is already 0/1/2 it passes straight through. Otherwise a
# best-effort text mapping is applied and anything unrecognised becomes NA
# (with a warning) so problems are visible rather than silently wrong.
coerce_rgi_values <- function(v) {
  s <- trimws(as.character(v))
  s[is.na(s)] <- ""

  already <- s %in% c("", "0", "1", "2")
  if (all(already)) {
    out <- s
    out[out == ""] <- "0"               # blank = absent
    return(out)
  }

  map <- c(
    "0" = "0", "no hit" = "0", "no hits" = "0", "none" = "0",
    "absent" = "0", "false" = "0",
    "1" = "1", "loose" = "1", "present" = "1",
    "2" = "2", "strict" = "2", "perfect" = "2", "confident" = "2",
    "confident hit" = "2", "true" = "2", "nudged" = "2"
  )
  key <- tolower(s)
  out <- unname(map[key])
  out[key == ""] <- "0"                  # blank cell = absent
  unknown <- is.na(out) & key != ""
  if (any(unknown)) {
    warning("coerce_rgi_values(): unrecognised RGI value(s) set to NA: ",
            paste(unique(s[unknown]), collapse = ", "))
  }
  out
}

# ---- readers ----------------------------------------------------------------

# Read the NCBI accession/metadata table (.tsv or .csv) and return a tidy
# data.frame with exactly: accession, strain, geo_location, O-type, country.
# `col_map` optionally overrides auto-detection, e.g. list(accession = "Run").
read_metadata <- function(path, col_map = list()) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  raw <- read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "\"",
                    na.strings = c("", "NA", "N/A", "missing"))
  cols <- names(raw)

  pick <- function(field) {
    if (!is.null(col_map[[field]])) {
      if (!col_map[[field]] %in% cols)
        stop(sprintf("Metadata column '%s' (mapped to '%s') not found.",
                     col_map[[field]], field))
      return(col_map[[field]])
    }
    .detect_col(cols, .DEFAULT_COL_ALIASES[[field]])
  }

  acc_col   <- pick("accession")
  str_col   <- pick("strain")
  geo_col   <- pick("geo_location")
  otype_col <- pick("O-type")

  if (is.na(acc_col))
    stop("Could not find an accession column in metadata. ",
         "Looked for: ", paste(.DEFAULT_COL_ALIASES$accession, collapse = ", "),
         ". Use --map accession=<yourcolumn> to specify it.")

  get <- function(col) if (is.na(col)) NA_character_ else as.character(raw[[col]])

  out <- data.frame(
    accession    = trimws(get(acc_col)),
    strain       = trimws(get(str_col)),
    geo_location = trimws(get(geo_col)),
    `O-type`     = trimws(get(otype_col)),
    check.names  = FALSE,
    stringsAsFactors = FALSE
  )
  out$country <- parse_country(out$geo_location)
  # Keep a row if it has EITHER an accession OR a strain to match on.
  has_id <- (!is.na(out$accession) & out$accession != "") |
            (!is.na(out$strain)    & out$strain    != "")
  out <- out[has_id, , drop = FALSE]

  message(sprintf(
    "Metadata: %d rows | accession='%s' strain='%s' geo='%s' O-type='%s'",
    nrow(out), acc_col,
    ifelse(is.na(str_col),   "<none>", str_col),
    ifelse(is.na(geo_col),   "<none>", geo_col),
    ifelse(is.na(otype_col), "<none>", otype_col)))
  out
}

# Read a raw RGI heatmap matrix (genes x samples OR samples x genes) and return
# a data.frame with one row per sample and one "0/1/2" column per AMR gene,
# plus a leading `sample` column. Orientation is auto-detected by overlap with
# the provided `sample_ids` (normalised).
read_rgi_matrix <- function(path, sample_ids = character(0)) {
  sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
  raw <- read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "\"", row.names = NULL)
  if (ncol(raw) < 2) stop("RGI matrix has fewer than 2 columns; check delimiter.")

  idcol  <- raw[[1]]                 # first column = labels for axis 1
  body   <- raw[, -1, drop = FALSE]  # remaining columns = labels for axis 2
  axis2  <- names(body)

  want <- normalize_id(sample_ids)
  ov_rows <- if (length(want) > 0) mean(normalize_id(idcol) %in% want) else 0
  ov_cols <- if (length(want) > 0) mean(normalize_id(axis2) %in% want) else 0

  if (ov_cols >= ov_rows) {
    # genes are rows, samples are columns (typical `rgi heatmap` output) -> transpose
    mat <- t(as.matrix(body))
    colnames(mat) <- idcol           # genes
    samples <- axis2
    orient <- "genes x samples (transposed to samples x genes)"
  } else {
    # samples are rows, genes are columns already
    mat <- as.matrix(body)
    colnames(mat) <- axis2           # genes
    samples <- idcol
    orient <- "samples x genes (kept as-is)"
  }

  df <- as.data.frame(mat, check.names = FALSE, stringsAsFactors = FALSE)
  df[] <- lapply(df, coerce_rgi_values)
  df <- cbind(sample = as.character(samples), df, stringsAsFactors = FALSE)

  best_ov <- max(ov_rows, ov_cols)
  message(sprintf(
    "RGI matrix: %d samples x %d genes | orientation: %s | sample-id overlap: %.0f%%",
    nrow(df), ncol(df) - 1, orient, 100 * best_ov))
  if (best_ov == 0 && length(want) > 0)
    message("  NOTE: 0% overlap between RGI sample names and metadata IDs after ",
            "normalisation. Check that RGI column names correspond to strains or ",
            "accessions in the metadata file.")
  df
}

# Read Pasty serotyping results (.tsv or .csv) and return a data.frame with
# columns: sample_orig, sample_norm, otype.
#
# Pasty writes <strain>.pasty as the sample name and the O-type in a column
# called "type". Both the sample column and type column are auto-detected;
# override with col_map[["sample"]] or col_map[["O-type"]].
read_pasty_serotypes <- function(path, col_map = list()) {
  sep <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
  raw <- read.delim(path, sep = sep, header = TRUE, check.names = FALSE,
                    stringsAsFactors = FALSE, quote = "\"",
                    na.strings = c("", "NA", "N/A"))
  cols <- names(raw)

  pick_pasty <- function(field, candidates) {
    if (!is.null(col_map[[field]])) {
      if (!col_map[[field]] %in% cols)
        stop(sprintf("Pasty file '%s': column '%s' not found.", path, col_map[[field]]))
      return(col_map[[field]])
    }
    hit <- .detect_col(cols, candidates)
    if (is.na(hit))
      stop(sprintf("Pasty file '%s': cannot find '%s' column. Looked for: %s.\n",
                   path, field, paste(candidates, collapse = ", ")),
           "Use --map ", field, "=<yourcolumn> to override.")
    hit
  }

  s_col <- pick_pasty("sample",
                       c("sample", "id", "name", "isolate", "strain", "accession"))
  t_col <- pick_pasty("O-type",
                       c("type", "o-type", "otype", "serotype", "o_antigen",
                         "predicted_type", "serogroup"))

  out <- data.frame(
    sample_orig = trimws(as.character(raw[[s_col]])),
    otype       = trimws(as.character(raw[[t_col]])),
    stringsAsFactors = FALSE
  )
  # normalize_id already strips the .pasty suffix before collapsing punctuation
  out$sample_norm <- normalize_id(out$sample_orig)
  out <- out[nzchar(out$sample_norm), , drop = FALSE]

  message(sprintf("Pasty serotypes: %d samples | sample='%s' type='%s'",
                  nrow(out), s_col, t_col))
  out
}

# Extract tip labels from a Newick tree file (base R only, no ape dependency).
# Tip labels appear after '(' or ',' and before ':' in the Newick string;
# internal-node bootstrap values appear after ')' and are thus excluded.
# Returns a character vector of the original (un-normalised) tip labels.
read_tree_tips <- function(path) {
  nwk <- paste(readLines(path, warn = FALSE), collapse = "")
  # Lookbehind for '(' or ','; lookahead for ':'; body excludes the five
  # Newick-special characters , ( ) ; :
  tips <- regmatches(nwk,
                     gregexpr("(?<=[,(])[^,();:]+(?=:)", nwk, perl = TRUE))[[1]]
  tips <- trimws(tips)
  tips <- tips[nzchar(tips)]
  # Remove any pure-numeric tokens that slip through on unusual topologies
  tips <- tips[!grepl("^[0-9.]+$", tips)]
  if (length(tips) == 0)
    warning("read_tree_tips(): no tip labels found in '", path, "'.")
  else
    message(sprintf("Phylogeny: %d tip labels read from '%s'.",
                    length(tips), basename(path)))
  tips
}

# ---- the merge --------------------------------------------------------------

# Build the merged table from any combination of:
#   metadata_path  - NCBI accession table (accession, strain, geo, O-type)
#   rgi_path       - RGI heatmap matrix (0/1/2 per gene per sample)
#   serotype_path  - Pasty O-typing results (overrides metadata O-type)
#   tree_path      - Newick phylogeny (adds `in_tree` TRUE/FALSE column)
#
# At least one of metadata_path / rgi_path must be supplied.
#
# Returns a data.frame in the schema:
#   accession, strain, geo_location, O-type [, in_tree], <gene columns...>
build_merged_table <- function(metadata_path = NULL, rgi_path = NULL,
                               serotype_path = NULL, tree_path = NULL,
                               col_map = list(), keep_raw_geo = FALSE) {

  # ---- 1. Read metadata -------------------------------------------------------
  if (!is.null(metadata_path)) {
    meta <- read_metadata(metadata_path, col_map = col_map)
  } else {
    meta <- NULL
  }

  # ---- 2. Read RGI matrix -----------------------------------------------------
  if (!is.null(rgi_path)) {
    ids  <- if (!is.null(meta))
              unique(c(meta$accession, meta$strain[!is.na(meta$strain)]))
            else character(0)
    rgi  <- read_rgi_matrix(rgi_path, sample_ids = ids)
    gene_cols <- setdiff(names(rgi), "sample")
  } else {
    rgi       <- NULL
    gene_cols <- character(0)
  }

  if (is.null(meta) && is.null(rgi))
    stop("At least one of metadata_path / rgi_path must be provided.")

  # ---- 3. Build the core merged data.frame ------------------------------------
  blank_na <- function(x) { x[x == ""] <- NA; x }

  if (!is.null(rgi) && !is.null(meta)) {
    # Primary path: merge metadata LEFT JOIN RGI on normalised ID.
    rgi_key  <- blank_na(normalize_id(rgi$sample))
    acc_key  <- blank_na(normalize_id(meta$accession))
    strn_key <- blank_na(normalize_id(meta$strain))

    # For each RGI sample find the matching metadata row (accession first, then strain).
    midx <- match(rgi_key, acc_key)
    need <- is.na(midx)
    midx[need] <- match(rgi_key[need], strn_key)

    matched <- !is.na(midx)
    if (!any(matched))
      stop("Merge produced 0 rows. RGI sample names do not match the metadata ",
           "accession/strain columns even after normalisation. Check that RGI ",
           "sample names correspond to the metadata accession or strain column.")

    meta_part <- meta[midx[matched],
                      c("accession", "strain", "geo_location", "O-type", "country"),
                      drop = FALSE]
    gene_part <- if (length(gene_cols) > 0)
                   rgi[matched, gene_cols, drop = FALSE]
                 else data.frame(row.names = seq_len(sum(matched)))
    merged <- cbind(meta_part, gene_part)

    rgi_only  <- rgi$sample[!matched]
    meta_only <- setdiff(
      meta$accession[!is.na(meta$accession) & meta$accession != ""],
      merged$accession)
    if (length(meta_only) > 0)
      message("  ", length(meta_only), " metadata sample(s) had no RGI match: ",
              paste(head(meta_only, 10), collapse = ", "),
              if (length(meta_only) > 10) " ..." else "")
    if (length(rgi_only) > 0)
      message("  ", length(rgi_only), " RGI sample(s) had no metadata match: ",
              paste(head(rgi_only, 10), collapse = ", "),
              if (length(rgi_only) > 10) " ..." else "")

  } else if (!is.null(rgi)) {
    # RGI-only: derive strain from RGI sample names (strip _rgi suffix).
    message("No metadata provided - using RGI sample names as strain identifiers.")
    merged <- data.frame(
      accession    = NA_character_,
      strain       = sub("_rgi$", "", rgi$sample, ignore.case = TRUE),
      geo_location = NA_character_,
      `O-type`     = NA_character_,
      country      = NA_character_,
      check.names  = FALSE,
      stringsAsFactors = FALSE
    )
    if (length(gene_cols) > 0)
      merged <- cbind(merged, rgi[, gene_cols, drop = FALSE])

  } else {
    # Metadata-only: no AMR gene columns.
    message("No RGI matrix provided - output contains only metadata columns.")
    merged <- data.frame(
      accession    = meta$accession,
      strain       = meta$strain,
      geo_location = meta$geo_location,
      `O-type`     = meta$`O-type`,
      country      = meta$country,
      check.names  = FALSE,
      stringsAsFactors = FALSE
    )
  }
  rownames(merged) <- NULL

  # ---- 4. Merge Pasty O-types (optional) --------------------------------------
  if (!is.null(serotype_path)) {
    pasty <- read_pasty_serotypes(serotype_path, col_map = col_map)

    row_acc  <- blank_na(normalize_id(merged$accession))
    row_str  <- blank_na(normalize_id(merged$strain))

    # Match each merged row to a Pasty entry (accession first, then strain).
    pidx <- match(row_acc, pasty$sample_norm)
    need <- is.na(pidx)
    pidx[need] <- match(row_str[need], pasty$sample_norm)

    has_pasty <- !is.na(pidx)
    if (any(has_pasty))
      merged$`O-type`[has_pasty] <- pasty$otype[pidx[has_pasty]]

    message(sprintf("  Pasty O-types: %d/%d samples matched.",
                    sum(has_pasty), nrow(merged)))
    if (sum(has_pasty) < nrow(merged))
      message("  ", nrow(merged) - sum(has_pasty),
              " sample(s) had no Pasty match (O-type left as-is).")
  }

  # ---- 5. Add in_tree flag (optional) -----------------------------------------
  if (!is.null(tree_path)) {
    tips      <- read_tree_tips(tree_path)
    tips_norm <- normalize_id(tips)

    row_acc <- blank_na(normalize_id(merged$accession))
    row_str <- blank_na(normalize_id(merged$strain))

    merged$in_tree <- (!is.na(row_acc) & row_acc %in% tips_norm) |
                      (!is.na(row_str) & row_str %in% tips_norm)

    message(sprintf("  Tree: %d/%d samples present in the phylogeny.",
                    sum(merged$in_tree), nrow(merged)))
  }

  # ---- 6. Assemble output schema ----------------------------------------------
  if (!keep_raw_geo) merged$geo_location <- merged$country

  fixed_cols <- c("accession", "strain", "geo_location", "O-type")
  if (!is.null(tree_path)) fixed_cols <- c(fixed_cols, "in_tree")
  out <- merged[, c(fixed_cols, gene_cols), drop = FALSE]
  out <- out[order(out$accession, na.last = TRUE), , drop = FALSE]
  rownames(out) <- NULL

  message(sprintf(
    "Merged table: %d samples | %d AMR genes | O-type source: %s | in_tree: %s",
    nrow(out), length(gene_cols),
    if (!is.null(serotype_path)) "Pasty" else "metadata",
    if (!is.null(tree_path)) "included" else "not included"))
  out
}

# Convenience wrapper: build + write CSV in one call.
write_merged_table <- function(metadata_path = NULL, rgi_path = NULL,
                               out_path = "merged_tables.csv",
                               serotype_path = NULL, tree_path = NULL,
                               col_map = list(), keep_raw_geo = FALSE) {
  out <- build_merged_table(
    metadata_path = metadata_path,
    rgi_path      = rgi_path,
    serotype_path = serotype_path,
    tree_path     = tree_path,
    col_map       = col_map,
    keep_raw_geo  = keep_raw_geo
  )
  write.csv(out, out_path, row.names = FALSE, na = "")
  message("Wrote merged table -> ", out_path)
  invisible(out)
}
