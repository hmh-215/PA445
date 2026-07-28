#!/usr/bin/env Rscript

# run_all.R  (Strategy 1: all-in-one convenience runner)
# -----------------------------------------------------------------------------
# Go from raw inputs straight to the final figure in one command.
# Internally calls the same two modular pieces: build the merged table
# (parse_merge_utils.R), then plot (02_plot_amr_tree_heatmap.R).
#
# Usage:
#   Rscript run_all.R \
#       --accession selected_genomes.tsv \
#       --rgi       PA445.AMR_heatmap-72.csv \
#       --serotype  PA_serotypes.tsv \
#       --tree      PA445.phylogeny.treefile \
#       --out-prefix AMR_tree_heatmap
#
# Options:
#   --accession  PATH   NCBI accession/metadata table (.tsv or .csv)
#   --metadata   PATH   alias for --accession (legacy name)
#   --rgi        PATH   RGI heatmap matrix                           [required]
#   --serotype   PATH   Pasty serotyping results (.tsv or .csv)      [optional]
#   --tree       PATH   Newick phylogeny tree    [default: PA445.phylogeny.treefile]
#   --merged     PATH   intermediate merged CSV  [default: merged_tables.csv]
#   --out-prefix STR    figure basename          [default: AMR_tree_heatmap]
#   --edge-scale NUM    branch-length multiplier [default: 3]
#   --keep-raw-geo      keep "Country: city" in geo_location
#   --map field=col     override a column name (repeatable)
# -----------------------------------------------------------------------------

.this_file <- (function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) normalizePath(f[1]) else "run_all.R"
})()
here <- dirname(.this_file)
source(file.path(here, "parse_merge_utils.R"))

parse_args <- function(args) {
  opt <- list(accession = NULL, rgi = NULL, serotype = NULL,
              tree = "PA445.phylogeny.treefile",
              merged = "merged_tables.csv", out_prefix = "AMR_tree_heatmap",
              edge_scale = "3", keep_raw_geo = FALSE, map = list())
  i <- 1
  while (i <= length(args)) {
    a <- args[i]; take <- function() { i <<- i + 1; args[i] }
    if      (a %in% c("--accession", "--metadata", "-a", "-m")) opt$accession  <- take()
    else if (a %in% c("--rgi",      "-r"))                      opt$rgi        <- take()
    else if (a %in% c("--serotype", "-s"))                      opt$serotype   <- take()
    else if (a %in% c("--tree",     "-t"))                      opt$tree       <- take()
    else if (a == "--merged")                                    opt$merged     <- take()
    else if (a %in% c("--out-prefix", "-o"))                    opt$out_prefix <- take()
    else if (a == "--edge-scale")                                opt$edge_scale <- take()
    else if (a == "--keep-raw-geo")                              opt$keep_raw_geo <- TRUE
    else if (a == "--map") {
      kv <- strsplit(take(), "=", fixed = TRUE)[[1]]
      if (length(kv) != 2) stop("--map expects field=column")
      opt$map[[kv[1]]] <- kv[2]
    }
    else if (a %in% c("--help", "-h")) {
      cat("See header of run_all.R for usage.\n"); quit(save = "no")
    }
    else stop("Unknown argument: ", a)
    i <- i + 1
  }
  opt
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(opt$rgi))
  stop("--rgi is required. Run with --help for usage.")

paths_to_check <- Filter(Negate(is.null),
                         list(opt$accession, opt$rgi, opt$serotype, opt$tree))
for (p in paths_to_check)
  if (!file.exists(p)) stop("Input file not found: ", p)

# Step 1: parse + merge -> merged CSV
message("== Step 1/2: building merged table ==")
write_merged_table(
  metadata_path = opt$accession,
  rgi_path      = opt$rgi,
  serotype_path = opt$serotype,
  tree_path     = NULL,             # tree used by plotter, not embedded in CSV
  out_path      = opt$merged,
  col_map       = opt$map,
  keep_raw_geo  = opt$keep_raw_geo
)

# Step 2: plot (delegate to the standalone plotting script)
message("== Step 2/2: plotting ==")
plot_script <- file.path(here, "02_plot_amr_tree_heatmap.R")
status <- system2("Rscript", c(shQuote(plot_script),
                               "--merged",     shQuote(opt$merged),
                               "--tree",       shQuote(opt$tree),
                               "--out-prefix", shQuote(opt$out_prefix),
                               "--edge-scale", shQuote(opt$edge_scale)))
if (status != 0) stop("Plotting step failed (exit ", status, ").")
message("== Done ==")
