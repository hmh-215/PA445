# AMR tree-heatmap pipeline

Refactor of `plot_amr_tree_heatmap.R` so you no longer have to build
`merged_tables.csv` by hand. Two ways to run it (pick whichever suits you —
they share the same code):

| | Strategy 2 (modular) | Strategy 1 (all-in-one) |
|---|---|---|
| Build merged table | `01_build_merged_table.R` | (done internally) |
| Draw figure | `02_plot_amr_tree_heatmap.R` | (done internally) |
| One command | — | `run_all.R` |

The merge/parse logic lives in `parse_merge_utils.R` (base R, no extra
packages). The plotting code is your original script, unchanged except for a
small `--merged/--tree/--out-prefix` argument layer (with the same defaults, so
`Rscript 02_plot_amr_tree_heatmap.R` behaves exactly as before).

## Strategy 2 — two steps

```bash
# 1. parse countries + merge RGI -> merged_tables.csv
Rscript 01_build_merged_table.R \
    --metadata subsampled_metadata.tsv \
    --rgi      rgi_heatmap.csv \
    --out      merged_tables.csv

# 2. draw the figure from the merged table + tree
Rscript 02_plot_amr_tree_heatmap.R \
    --merged merged_tables.csv \
    --tree   PA445.phylogeny.treefile \
    --out-prefix AMR_tree_heatmap
```

Re-plotting (tweaking colours, edge scale, etc.) only re-runs step 2 — no need
to re-merge.

## Strategy 1 — one step

```bash
Rscript run_all.R \
    --metadata subsampled_metadata.tsv \
    --rgi      rgi_heatmap.csv \
    --tree     PA445.phylogeny.treefile \
    --out-prefix AMR_tree_heatmap
```

## What the inputs are assumed to look like

I didn't have your raw `.tsv` / raw RGI files, so the reader is deliberately
forgiving and self-reports what it detected. Check the console messages on the
first run.

**Raw subsampling metadata (`--metadata`, .tsv or .csv)**
- One row per sample.
- Needs columns for *accession*, *strain*, *geo_location*, *O-type*. Names are
  auto-detected case/punctuation-insensitively from common aliases
  (e.g. `assembly_accession`, `geo_loc_name`, `serotype`, `isolate`, …).
- `geo_location` is expected in NCBI form `Country: city/institute`; the
  country (text before the first `:`) is extracted. `missing` / `N/A` / blank
  become unknown.
- Override any column with `--map field=column`, e.g.
  `--map accession=Run --map O-type=O_serotype`.

**Raw RGI heatmap matrix (`--rgi`, .csv or .tsv)**
- A genes × samples (typical `rgi heatmap` output) **or** samples × genes
  matrix. Orientation is auto-detected by which axis overlaps your sample IDs
  and transposed if needed, so the output always has samples as rows.
- Cell values map to the plot's `0/1/2` scheme
  (**0 = no hit, 1 = loose hit, 2 = strict hit**). Values already in `0/1/2`
  pass straight through; common text categories (`Loose`→1, `Strict`→2,
  `Perfect`→2, blank→0) are mapped otherwise, and anything unrecognised becomes
  `NA` with a warning.

**Sample matching.** RGI sample names are matched to metadata by accession
first, then by strain (so the accession-less reference `PA445`, whose RGI
sample is named by its strain, still joins). Matching ignores case,
punctuation, and `.fna`/`_genomic`-style suffixes, but keeps the `GCF_`/`GCA_`
prefix so the two assembly versions of one genome stay separate rows.

## Options worth knowing

- `--keep-raw-geo` — write the full `Country: city` string into
  `geo_location` instead of just the country. (The plot extracts the country
  either way; default matches your existing `merged_tables.csv`, which stored
  just the country.)
- `--edge-scale N` (plotting) — branch-length multiplier; default `3`, matching
  the original hard-coded value.

## Validation

`01_build_merged_table.R` was checked by reconstructing your supplied
`merged_tables.csv` from synthesised raw inputs (non-standard column names + a
transposed genes×samples RGI matrix + a strain-named reference). Result: the
166-gene × 74-sample AMR matrix came back **bit-identical**, columns and order
matched, and metadata matched (the reconstruction additionally normalises
`N/A`/`missing` to blank, which the plotting script already treats as missing).
The plotting script's internals are unchanged from your working version.
