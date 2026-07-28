#!/usr/bin/env Rscript
#
# plot_time_tree.R -- Plot a BEAST/MCC time tree with tip labels, posterior
# support values, and 95% HPD node-height bars, optionally rescaled to
# calendar years.
#
# Usage:
#   Rscript plot_time_tree.R --tree <mcc.tree> [--prefix <out_prefix>]
#                             [--rtt <rtt.csv>] [--present <decimal_year>]
#
# Options:
#   -t, --tree      Path to the input BEAST/MCC Nexus tree file. (required)
#   -p, --prefix    Output file prefix; writes <prefix>.pdf and <prefix>.png.
#                   Defaults to the --tree filename (without extension).
#       --rtt       Optional TreeTime rtt.csv used to infer the present time
#                   (the most recent tip's decimal-year date) when --present
#                   is not given.
#       --present   Present time as a decimal year (e.g. 2026.5), used to
#                   anchor the x-axis to real calendar years. If omitted,
#                   the script tries to infer it from --rtt; if neither is
#                   given, the x-axis shows raw branch-length units from the
#                   root instead of calendar years.  The latest *observed*
#                   tip date is the appropriate value; in this analysis it
#                   is PA445 = 2026.5.
#
# Examples:
#   Rscript plot_time_tree.R --tree PA445.lineage.mcc.tree --prefix PA445_time_tree
#   Rscript plot_time_tree.R --tree PA445.lineage.mcc.tree --prefix PA445_time_tree_year --rtt rtt.csv
#   Rscript plot_time_tree.R --tree PA445.lineage.mcc.tree --prefix PA445_time_tree_year --present 2026.5
#
# Requires: optparse, treeio, ggtree, ggplot2

required_packages <- c("optparse", "treeio", "ggtree", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required R package(s): ", paste(missing_packages, collapse = ", "),
       ". Install them, then rerun this script.", call. = FALSE)
}

library(optparse)
library(treeio)
library(ggtree)
library(ggplot2)

# 1. Define command-line options
option_list <- list(
  make_option(c("-t", "--tree"), type = "character", default = NULL,
              help = "Path to the input BEAST/MCC Nexus tree file", metavar = "character"),
  make_option(c("-p", "--prefix"), type = "character", default = NULL,
              help = "Output file prefix. Script will write <prefix>.pdf and <prefix>.png [default derived from --tree]", metavar = "character"),
  make_option(c("--rtt"), type = "character", default = NULL,
              help = "Optional TreeTime rtt.csv (for inferring/validating present time).", metavar = "character"),
  make_option(c("--present"), type = "double", default = NA_real_,
              help = "Present time as decimal year, used to anchor the x-axis to real calendar years. If not set, script attempts to infer it from --rtt.", metavar = "number")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$tree)) {
  print_help(opt_parser)
  stop("Error: You must provide an input tree file using --tree or -t.", call. = FALSE)
}

# 2. Load the tree file and clean up column names
message("Reading tree file: ", opt$tree)
mcc_tree <- read.beast(opt$tree)

# BEAST node heights are the posterior time estimates (in years before the
# most recent sample).  The Newick branch lengths in an MCC tree are often
# independent posterior summaries and need not equal parent-height minus
# child-height.  ggtree draws from those Newick branch lengths, which is why a
# direct plot can put dated tips in the wrong places.  Reconstruct the edge
# lengths from the MCC node heights so the plotted chronology agrees exactly
# with the BEAST height annotations.
# In an MCC summary, `height` is commonly a posterior mean.  Means of node
# ages can be mutually inconsistent after trees are summarized.  The median
# node heights in this file retain a valid ancestor-before-descendant order,
# so prefer them for a drawable dated tree.
time_height_col <- intersect(c("height_median", "height"), colnames(mcc_tree@data))
if (length(time_height_col) == 0) {
  stop("The MCC tree has no 'height' or 'height_median' annotation, so calendar dates cannot be reconstructed.",
       call. = FALSE)
}
time_height_col <- time_height_col[1]
node_heights <- suppressWarnings(as.numeric(mcc_tree@data[[time_height_col]]))
# BEAST may omit height_median for the most-recent tip (its height is exactly
# zero); use its annotated height in that case.
if (time_height_col == "height_median" && "height" %in% colnames(mcc_tree@data)) {
  annotated_height <- suppressWarnings(as.numeric(mcc_tree@data$height))
  node_heights[is.na(node_heights)] <- annotated_height[is.na(node_heights)]
}
names(node_heights) <- as.character(mcc_tree@data$node)
edge <- mcc_tree@phylo$edge
parent_height <- unname(node_heights[as.character(edge[, 1])])
child_height <- unname(node_heights[as.character(edge[, 2])])

if (anyNA(parent_height) || anyNA(child_height)) {
  stop("Could not match every phylogeny node to a BEAST height annotation.", call. = FALSE)
}

height_edge_lengths <- parent_height - child_height
if (any(height_edge_lengths < -sqrt(.Machine$double.eps))) {
  # MCC node ages are marginal summaries.  Consequently, two independently
  # summarized adjacent nodes can differ by a small amount in the impossible
  # direction even when every sampled BEAST tree is chronological.  Project
  # the drawing heights onto a valid tree by raising each parent to its oldest
  # child; this affects only the displayed edge(s), not the MCC file or HPDs.
  original_heights <- node_heights
  for (iteration in seq_len(nrow(edge))) {
    for (i in rev(seq_len(nrow(edge)))) {
      parent_id <- as.character(edge[i, 1])
      child_id <- as.character(edge[i, 2])
      node_heights[parent_id] <- max(node_heights[parent_id], node_heights[child_id])
    }
  }
  parent_height <- unname(node_heights[as.character(edge[, 1])])
  child_height <- unname(node_heights[as.character(edge[, 2])])
  height_edge_lengths <- parent_height - child_height
  adjusted_nodes <- sum(node_heights > original_heights + sqrt(.Machine$double.eps))
  message("Adjusted ", adjusted_nodes,
          " MCC parent-node median age(s) upward to remove a marginal-summary inversion for plotting.")
}
mcc_tree@phylo$edge.length <- pmax(height_edge_lengths, 0)
message("Drawing branches from BEAST '", time_height_col, "' node ages.")

# Rename/sanitize HPD column name variants for geom_range() safety
hpd_candidates <- c(
  "height_95%_HPD", "height_95_HPD",
  "height_0.95_HPD", "height_0.95_hpd",
  "height_0.95_Hpd", "height_0.95_hPD",
  "length_95%_HPD", "length_95_HPD",
  "length_0.95_HPD", "length_0.95_Hpd"
)

existing <- intersect(hpd_candidates, colnames(mcc_tree@data))
if (length(existing) == 0) {
  stop(
    "Could not find an HPD column for the 95% interval in mcc_tree@data.\n" ,
    "Available columns include: ", paste(head(colnames(mcc_tree@data), 30), collapse = ", "),
    "\nLooked for: ", paste(hpd_candidates, collapse = ", ")
  )
}

# If height_95%_HPD exists but height_95_HPD doesn't, sanitize the name
if ("height_95%_HPD" %in% existing && !("height_95_HPD" %in% existing)) {
  colnames(mcc_tree@data)[colnames(mcc_tree@data) == "height_95%_HPD"] <- "height_95_HPD"
}

# Prefer height-based HPD range column if available
preferred <- c("height_95_HPD", "height_0.95_HPD", "height_0.95_Hpd", "height_0.95_hpd")
hpd_range_col <- intersect(preferred, colnames(mcc_tree@data))
if (length(hpd_range_col) == 0) {
  hpd_range_col <- existing[1]
} else {
  hpd_range_col <- hpd_range_col[1]
}

# 3. Generate the time-calibrated plot
message("Generating plot...")
# Optional conversion: if present time is provided/inferred, show calendar year
present <- opt$present

if (is.na(present) && !is.null(opt$rtt)) {
  message("Reading rtt file: ", opt$rtt)
  rtt <- read.csv(opt$rtt, comment.char = "#", stringsAsFactors = FALSE)
  if ("date" %in% colnames(rtt)) {
    # Use the maximum inferred date as “present”. This matches your rtt.csv example where PA445 ~ 2026.5.
    present <- max(as.numeric(rtt$date), na.rm = TRUE)
    message("Inferred present time from rtt.csv: ", present)
  } else {
    warning("rtt.csv does not contain a 'date' column; cannot infer present time.")
  }
}

if (!is.na(present)) {
  # mrsd anchors the most-recent tip to the real calendar date and lets ggtree
  # shift the whole x-axis into calendar time itself (see ggtree:::scaleX_by_time_from_mrsd);
  # decimal2Date/mrsd must be used together rather than relabeling axis ticks by hand,
  # since ggtree's x is measured forward from the root, not backward from the present.
  mrsd <- decimal2Date(present)
  x_label <- paste0("Calendar Year (present=", round(present, 3), ")")
  root_year <- present - max(node_heights, na.rm = TRUE)
  message("MCC root estimate: ", format(root_year, scientific = FALSE, trim = TRUE),
          " calendar years (", round(max(node_heights, na.rm = TRUE)),
          " years before the latest sample).")
} else {
  x_label <- "Time (branch length from root)"
  mrsd <- NULL
}

p <- ggtree(mcc_tree, mrsd = mrsd) +
  theme_tree2() +
  labs(x = x_label)

# 4. Add customized layout elements
p_custom <- p + 
  geom_tiplab(size = 4, color = "black", fontface = "bold",
              offset = max(1, 0.001 * max(node_heights))) + 
  geom_text2(aes(label = ifelse(!is.na(posterior) & as.numeric(posterior) >= 0.5, 
                               round(as.numeric(posterior), 2), "")), 
             hjust = 1.2, vjust = -0.5, size = 3, color = "darkred")

# 5. Add 95% HPD error bars
p_hpd <- p_custom + 
  geom_range(range = hpd_range_col, color = 'blue', alpha = 0.4, size = 1.5)

p_final <- p_hpd + hexpand(0.3, direction = 1)

# Determine output prefix
prefix <- opt$prefix
if (is.null(prefix)) {
  # Derive from tree filename
  prefix <- tools::file_path_sans_ext(basename(opt$tree))
}

pdf_out <- paste0(prefix, ".pdf")
png_out <- paste0(prefix, ".png")

# 6. Save the outputs
message("Saving plot to: ", pdf_out)
ggsave(pdf_out, plot = p_final, width = 10, height = 7)

message("Saving plot to: ", png_out)
ggsave(png_out, plot = p_final, width = 10, height = 7, dpi = 300)

message("Done!")
