#!/usr/bin/env python3
"""
Reconcile Prokka annotations against RGI calls to identify true-positive
resistance gene hits for downstream clinker synteny analysis.

BACKGROUND
----------
Prokka tends to annotate broadly based on sequence similarity, which can
inflate hit counts (false positives). RGI uses the CARD database with
curated cut-offs (Strict / Perfect), giving higher-confidence calls.
This script treats RGI as the ground truth and flags Prokka annotations
that are supported by a matching RGI call as "true positives".

Matching strategy
-----------------
Gene names in Prokka (--prokka-gene) and RGI (--rgi-gene) are often
non-identical (e.g. "blaNDM-1" vs "NDM-1"). Both values are normalised
(lowercased, hyphens/spaces/underscores removed) before comparison so
partial matches still resolve correctly.  You can also supply alternate
alias tokens with --aliases if further flexibility is needed.

Directory conventions (matching PA_comparative.bash)
-----------------------------------------------------
  prokka inputs : {PROKKA_PA}/<strain>.prokka/<strain>.prokka.tsv
  rgi inputs    : {RGI_PA}/<strain>_rgi.txt   (RGI main --output basename)

Usage examples
--------------
# Basic: compare blaNDM-1 (prokka) vs NDM-1 (rgi) across all strains
python prokka_rgi_reconcile.py \\
    --prokka-inputs /path/to/prokka/ \\
    --rgi-inputs    /path/to/rgi/ \\
    --prokka-gene   blaNDM-1 \\
    --rgi-gene      NDM-1 \\
    --outdir        /path/to/reconcile_out/

# Also write a file of true-positive GBK paths ready to copy to clinker
python prokka_rgi_reconcile.py \\
    --prokka-inputs /path/to/prokka/ \\
    --rgi-inputs    /path/to/rgi/ \\
    --prokka-gene   blaNDM-1 \\
    --rgi-gene      NDM-1 \\
    --outdir        /path/to/reconcile_out/ \\
    --truepos-prokka-paths \\
    --rgi-cutoff Strict

# Then copy true-positive GBKs to clinker input folder:
#   while IFS= read -r p; do cp "$p" /path/to/clinker/gbk_inputs/; done \\
#       < /path/to/reconcile_out/truepos_prokka_gbk_paths.txt
"""

import argparse
import csv
import glob
import os
import re
import sys


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalise(text: str) -> str:
    """Lowercase and strip hyphens, underscores, spaces for fuzzy matching."""
    return re.sub(r"[-_ ]", "", text.lower())


def find_file(directory: str, pattern: str) -> str | None:
    """Return the first file matching a glob pattern inside directory, or None."""
    matches = glob.glob(os.path.join(directory, pattern))
    return matches[0] if matches else None


def strain_from_prokka_dir(path: str) -> str:
    """
    Extract the strain label from a Prokka output directory name.
    e.g. /some/path/PA445.prokka  ->  PA445
         /some/path/LESB58.prokka ->  LESB58
    """
    basename = os.path.basename(path.rstrip("/"))
    # Strip trailing .prokka if present
    return basename.removesuffix(".prokka")


def strain_from_rgi_file(path: str) -> str:
    """
    Extract the strain label from an RGI output filename.
    RGI main writes  <prefix>_rgi.txt  (with --output <prefix>_rgi).
    e.g. PA445_rgi.txt -> PA445
    """
    basename = os.path.basename(path)
    # Strip known RGI suffixes
    for suffix in ("_rgi.txt", "_rgi.tsv", ".txt", ".tsv"):
        if basename.endswith(suffix):
            basename = basename[: -len(suffix)]
            break
    return basename


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_prokka_tsv(tsv_path: str, gene_query: str) -> list[dict]:
    """
    Parse a Prokka .tsv annotation table and return rows where the 'gene'
    or 'product' column contains gene_query (fuzzy match).

    Prokka TSV columns (tab-separated, no quoting):
        locus_tag  ftype  length_bp  gene  EC_number  COG  product
    """
    query_norm = normalise(gene_query)
    hits = []
    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            gene    = row.get("gene", "") or ""
            product = row.get("product", "") or ""
            if query_norm in normalise(gene) or query_norm in normalise(product):
                hits.append(dict(row))
    return hits


def parse_rgi_txt(txt_path: str, gene_query: str, cutoff: str | None) -> list[dict]:
    """
    Parse an RGI main output file (.txt / .tsv) and return rows where
    'Best_Hit_ARO' or 'ARO Term' contains gene_query (fuzzy match).

    Optionally filter by Cut_Off level (Perfect / Strict / Loose).

    Key RGI columns:
        ORF_ID  Contig  Start  Stop  Orientation  Cut_Off
        Best_Hit_ARO  Best_Identities  ARO  Model_type  ...
    """
    query_norm = normalise(gene_query)
    allowed_cutoffs = {normalise(c) for c in cutoff.split(",")} if cutoff else None
    hits = []
    with open(txt_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            aro_term = row.get("Best_Hit_ARO", "") or row.get("ARO Term", "") or ""
            cut      = row.get("Cut_Off", "") or ""
            if allowed_cutoffs and normalise(cut) not in allowed_cutoffs:
                continue
            if query_norm in normalise(aro_term):
                hits.append(dict(row))
    return hits


# ---------------------------------------------------------------------------
# Core reconciliation
# ---------------------------------------------------------------------------

def reconcile(
    prokka_dir: str,
    rgi_dir: str,
    prokka_gene: str,
    rgi_gene: str,
    rgi_cutoff: str | None,
) -> list[dict]:
    """
    Walk every Prokka subdirectory, find the matching RGI file for the same
    strain, and classify each Prokka hit as TRUE_POSITIVE or FALSE_POSITIVE.

    Returns a list of result dicts, one per (strain, prokka_hit) pair.
    """
    results = []

    # Collect all Prokka strain directories
    prokka_dirs = sorted(
        d for d in glob.glob(os.path.join(prokka_dir, "*"))
        if os.path.isdir(d)
    )

    if not prokka_dirs:
        sys.exit(f"ERROR: No subdirectories found in --prokka-inputs: {prokka_dir}")

    for pdir in prokka_dirs:
        strain = strain_from_prokka_dir(pdir)

        # Locate Prokka TSV (try <strain>.prokka.tsv first, then any .tsv)
        tsv_path = (
            find_file(pdir, f"{strain}.prokka.tsv")
            or find_file(pdir, "*.tsv")
        )
        if not tsv_path:
            print(f"  [{strain}] WARNING: No Prokka TSV found in {pdir} - skipping.")
            continue

        # Locate matching RGI file for this strain
        # RGI files are named <strain>_rgi.txt by convention in PA_comparative.bash
        # The RGI header-cleanup step prepends the strain name, so filenames may be
        # e.g. PA445_rgi.txt or LESB58_rgi.txt
        rgi_path = (
            find_file(rgi_dir, f"{strain}_rgi.txt")
            or find_file(rgi_dir, f"{strain}.txt")
            or find_file(rgi_dir, f"{strain}_rgi.tsv")
            or find_file(rgi_dir, f"{strain}.tsv")
        )

        prokka_hits = parse_prokka_tsv(tsv_path, prokka_gene)
        rgi_hits    = parse_rgi_txt(rgi_path, rgi_gene, rgi_cutoff) if rgi_path else []

        if not rgi_path:
            print(
                f"  [{strain}] WARNING: No RGI file found in {rgi_dir}"
                f"all Prokka hits will be marked as UNVERIFIED."
            )

        rgi_hit_count = len(rgi_hits)

        # Locate the Prokka GBK for this strain (needed for --truepos-prokka-paths)
        gbk_path = (
            find_file(pdir, f"{strain}.prokka.gbk")
            or find_file(pdir, "*.gbk")
        )

        for ph in prokka_hits:
            # A Prokka hit is a TRUE_POSITIVE if the strain has at least one
            # matching RGI call (at the requested cut-off level).
            # We do a per-strain verdict rather than per-locus because RGI
            # operates on contigs and its locus IDs differ from Prokka's.
            if rgi_path is None:
                verdict = "UNVERIFIED"
            elif rgi_hit_count > 0:
                verdict = "TRUE_POSITIVE"
            else:
                verdict = "FALSE_POSITIVE"

            results.append({
                "strain":            strain,
                "verdict":           verdict,
                "prokka_locus":      ph.get("locus_tag", ""),
                "prokka_gene":       ph.get("gene", ""),
                "prokka_product":    ph.get("product", ""),
                "prokka_length_bp":  ph.get("length_bp", ""),
                "rgi_hit_count":     rgi_hit_count,
                "rgi_best_aro":      "; ".join(
                                        r.get("Best_Hit_ARO", "") for r in rgi_hits
                                     ) if rgi_hits else "",
                "rgi_cutoff":        "; ".join(
                                        r.get("Cut_Off", "") for r in rgi_hits
                                     ) if rgi_hits else "",
                "rgi_identity":      "; ".join(
                                        r.get("Best_Identities", "") for r in rgi_hits
                                     ) if rgi_hits else "",
                "prokka_tsv":        tsv_path,
                "rgi_file":          rgi_path or "",
                "prokka_gbk":        gbk_path or "",
            })

    return results


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

TSV_COLUMNS = [
    "strain",
    "verdict",
    "prokka_locus",
    "prokka_gene",
    "prokka_product",
    "prokka_length_bp",
    "rgi_hit_count",
    "rgi_best_aro",
    "rgi_cutoff",
    "rgi_identity",
    "prokka_tsv",
    "rgi_file",
    "prokka_gbk",
]


def write_tsv(results: list[dict], outpath: str) -> None:
    with open(outpath, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=TSV_COLUMNS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)
    print(f"  Reconciliation table : {outpath}")


def write_summary(results: list[dict], prokka_gene: str, rgi_gene: str,
                  rgi_cutoff: str | None, outpath: str) -> None:
    tp = [r for r in results if r["verdict"] == "TRUE_POSITIVE"]
    fp = [r for r in results if r["verdict"] == "FALSE_POSITIVE"]
    uv = [r for r in results if r["verdict"] == "UNVERIFIED"]

    tp_strains = sorted({r["strain"] for r in tp})
    fp_strains = sorted({r["strain"] for r in fp})
    uv_strains = sorted({r["strain"] for r in uv})

    lines = [
        "=============================================================",
        " Prokka vs RGI Reconciliation Summary",
        "=============================================================",
        f"  Prokka gene query : {prokka_gene}",
        f"  RGI gene query    : {rgi_gene}",
        f"  RGI cut-off filter: {rgi_cutoff or 'none (all hits included)'}",
        "",
        f"  Total Prokka hits examined : {len(results)}",
        f"  TRUE_POSITIVE              : {len(tp)}",
        f"  FALSE_POSITIVE             : {len(fp)}",
        f"  UNVERIFIED (no RGI file)   : {len(uv)}",
        "",
        f"  Strains with TRUE_POSITIVE hits  ({len(tp_strains)}):",
    ] + [f"    {s}" for s in tp_strains] + [
        "",
        f"  Strains with FALSE_POSITIVE hits ({len(fp_strains)}):",
    ] + [f"    {s}" for s in fp_strains] + [
        "",
        f"  Strains with UNVERIFIED hits     ({len(uv_strains)}):",
    ] + [f"    {s}" for s in uv_strains] + [
        "=============================================================",
    ]

    with open(outpath, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    # Also print to stdout
    print("\n".join(lines))
    print(f"\n  Summary file         : {outpath}")


def write_truepos_paths(results: list[dict], outpath: str) -> None:
    """
    Write a plain-text file listing the Prokka GBK path for every
    TRUE_POSITIVE strain (one path per line, deduplicated).
    Intended for use with a bash while-read loop to bulk-copy GBKs.
    """
    seen = set()
    lines = []
    for r in results:
        if r["verdict"] == "TRUE_POSITIVE" and r["prokka_gbk"]:
            if r["prokka_gbk"] not in seen:
                seen.add(r["prokka_gbk"])
                lines.append(r["prokka_gbk"])

    with open(outpath, "w") as fh:
        fh.write("\n".join(sorted(lines)) + ("\n" if lines else ""))

    print(f"  True-positive GBK paths : {outpath}")
    print( "  Copy to clinker inputs with:")
    print(f"    while IFS= read -r p; do cp \"$p\" /path/to/clinker/gbk_inputs/; done \\")
    print(f"        < {outpath}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Reconcile Prokka annotations against RGI calls to identify "
            "true-positive resistance gene hits for clinker synteny analysis."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # Required
    parser.add_argument(
        "--prokka-inputs", required=True, metavar="DIR",
        help=(
            "Directory containing per-strain Prokka output subdirectories "
            "(each named <strain>.prokka/ and containing a <strain>.prokka.tsv "
            "and <strain>.prokka.gbk)."
        ),
    )
    parser.add_argument(
        "--rgi-inputs", required=True, metavar="DIR",
        help=(
            "Directory containing per-strain RGI output files "
            "(named <strain>_rgi.txt as produced by 'rgi main --output <strain>_rgi')."
        ),
    )
    parser.add_argument(
        "--prokka-gene", required=True, metavar="GENE",
        help=(
            "Gene name to search in Prokka TSV 'gene' and 'product' columns "
            "(case-insensitive, partial/fuzzy match). E.g. blaNDM-1"
        ),
    )
    parser.add_argument(
        "--rgi-gene", required=True, metavar="GENE",
        help=(
            "Gene/ARO term to search in RGI 'Best_Hit_ARO' column "
            "(case-insensitive, partial/fuzzy match). E.g. NDM-1"
        ),
    )
    parser.add_argument(
        "--outdir", required=True, metavar="DIR",
        help="Directory for output files (created if absent).",
    )

    # Optional
    parser.add_argument(
        "--rgi-cutoff", default="Strict,Perfect", metavar="LEVEL(S)",
        help=(
            "Comma-separated RGI Cut_Off level(s) to accept. "
            "Options: Perfect, Strict, Loose. "
            "Default: Strict,Perfect  (excludes Loose hits). "
            "Pass 'all' to include every cut-off level."
        ),
    )
    parser.add_argument(
        "--truepos-prokka-paths", action="store_true",
        help=(
            "If set, write a plain-text file listing the Prokka GBK path for "
            "every true-positive strain. Useful for bulk-copying GBKs to a "
            "clinker input folder with a bash while-read loop."
        ),
    )

    return parser.parse_args()


def main():
    args = parse_args()

    rgi_cutoff = None if args.rgi_cutoff.lower() == "all" else args.rgi_cutoff

    os.makedirs(args.outdir, exist_ok=True)

    print(f"\nProkka inputs : {args.prokka_inputs}")
    print(f"RGI inputs    : {args.rgi_inputs}")
    print(f"Prokka gene   : {args.prokka_gene}")
    print(f"RGI gene      : {args.rgi_gene}")
    print(f"RGI cut-off   : {rgi_cutoff or 'all'}")
    print(f"Output dir    : {args.outdir}\n")

    results = reconcile(
        prokka_dir  = args.prokka_inputs,
        rgi_dir     = args.rgi_inputs,
        prokka_gene = args.prokka_gene,
        rgi_gene    = args.rgi_gene,
        rgi_cutoff  = rgi_cutoff,
    )

    if not results:
        sys.exit(
            f"\nNo Prokka hits found for '{args.prokka_gene}' in any strain. "
            "Check --prokka-gene and the Prokka TSV files."
        )

    # Safe gene label for filenames
    safe_gene = re.sub(r"[-_ /]", "", args.prokka_gene)

    write_tsv(
        results,
        os.path.join(args.outdir, f"{safe_gene}.prokka_vs_rgi.tsv"),
    )

    write_summary(
        results,
        prokka_gene = args.prokka_gene,
        rgi_gene    = args.rgi_gene,
        rgi_cutoff  = rgi_cutoff,
        outpath     = os.path.join(args.outdir, f"{safe_gene}.reconciliation_summary.txt"),
    )

    if args.truepos_prokka_paths:
        write_truepos_paths(
            results,
            os.path.join(args.outdir, f"{safe_gene}.truepos_prokka_gbk_paths.txt"),
        )

    print("\nDone.")


if __name__ == "__main__":
    main()
