#!/usr/bin/env Rscript

# 01_build_merged_table.R  (Strategy 2, step 1)
# -----------------------------------------------------------------------------
# Merge an RGI heatmap with any combination of: an NCBI accession/metadata
# table, Pasty O-typing results, and a Newick phylogeny tree.  Writes the
# merged_tables.csv consumed by 02_plot_amr_tree_heatmap.R.
#
# Usage:
#   Rscript 01_build_merged_table.R \
#       --accession selected_genomes.tsv \
#       --rgi       PA445.AMR_heatmap-72.csv \
#       --serotype  PA_serotypes.tsv \
#       --tree      PA445.phylogeny.treefile \
#       --out       merged_tables.csv
#
# Options:
#   --accession PATH  NCBI accession/metadata table (.tsv or .csv)
#   --metadata  PATH  alias for --accession (legacy name)
#   --rgi       PATH  RGI heatmap matrix (.csv or .tsv)             [required]
#   --serotype  PATH  Pasty serotyping results (.tsv or .csv)       [optional]
#   --tree      PATH  Newick phylogeny tree (.treefile / .nwk)      [optional]
#   --out       PATH  output merged CSV          [default: merged_tables.csv]
#   --keep-raw-geo    keep "Country: city" in geo_location instead of
#                     writing only the parsed country
#   --map field=col   override a column name (repeatable), e.g.
#                       --map accession=Run
#                       --map O-type=serotype
#                       --map sample=isolate_id    (Pasty sample column)
#
# Column auto-detection covers common NCBI/assembly/Pasty aliases; use --map
# only when a column is named unusually.
#
# Output columns:
#   accession    - from the NCBI accession table (blank when no table given)
#   strain       - from the NCBI accession table (or RGI sample name)
#   geo_location - country parsed from geo_loc_name (or full string with --keep-raw-geo)
#   O-type       - from Pasty (if --serotype given) or from the accession table
#   in_tree      - TRUE/FALSE, present only when --tree is given
#   <gene cols>  - one column per AMR gene from the RGI matrix (0=absent,
#                  1=loose hit, 2=strict hit)
# -----------------------------------------------------------------------------

# Locate parse_merge_utils.R next to this script (works regardless of CWD).
.this_file <- (function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f[1]) else "01_build_merged_table.R"
})()
source(file.path(dirname(.this_file), "parse_merge_utils.R"))

# ---- argument parser --------------------------------------------------------
parse_args <- function(args) {
  opt <- list(accession = NULL, rgi = NULL, serotype = NULL, tree = NULL,
              out = "merged_tables.csv", keep_raw_geo = FALSE, map = list())
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    take <- function() { i <<- i + 1; args[i] }
    if      (a %in% c("--accession", "--metadata", "-a", "-m")) opt$accession <- take()
    else if (a %in% c("--rgi",       "-r"))                     opt$rgi       <- take()
    else if (a %in% c("--serotype",  "-s"))                     opt$serotype  <- take()
    else if (a %in% c("--tree",      "-t"))                     opt$tree      <- take()
    else if (a %in% c("--out",       "-o"))                     opt$out       <- take()
    else if (a == "--keep-raw-geo")                              opt$keep_raw_geo <- TRUE
    else if (a == "--map") {
      kv <- strsplit(take(), "=", fixed = TRUE)[[1]]
      if (length(kv) != 2) stop("--map expects field=column, e.g. --map accession=Run")
      opt$map[[kv[1]]] <- kv[2]
    }
    else if (a %in% c("--help", "-h")) { cat_help(); quit(save = "no") }
    else stop("Unknown argument: ", a)
    i <- i + 1
  }
  opt
}

cat_help <- function() {
  lines <- readLines(.this_file)
  cat(lines[grep("^#", lines)], sep = "\n")
  cat("\n")
}

args <- commandArgs(trailingOnly = TRUE)
opt  <- parse_args(args)

if (is.null(opt$rgi))
  stop("--rgi is required. Run with --help for usage.")

# Validate that all supplied paths exist.
paths_to_check <- Filter(Negate(is.null),
                         list(opt$accession, opt$rgi, opt$serotype, opt$tree))
for (p in paths_to_check)
  if (!file.exists(p)) stop("Input file not found: ", p)

write_merged_table(
  metadata_path = opt$accession,
  rgi_path      = opt$rgi,
  serotype_path = opt$serotype,
  tree_path     = opt$tree,
  out_path      = opt$out,
  col_map       = opt$map,
  keep_raw_geo  = opt$keep_raw_geo
)
