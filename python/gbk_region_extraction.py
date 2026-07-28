#!/usr/bin/env python3
"""
Extract flanking regions surrounding any target gene from a GenBank file.

Usage:
    python gbk_region_extraction.py <input.gbk> <output_dir> <strain> \
        [--gene GENE_NAME] [--flank FLANK_BP] [--distance-cap BP] \
        [--coord-table PATH]

Arguments:
    input.gbk       Input GenBank file (single or multi-record, e.g. Prokka output)
    output_dir      Directory to write extracted region GBK(s)
    strain          Strain label used in output filenames and descriptions

    --gene          Gene name to search for (case-insensitive, partial match).
                    Searches 'gene' and 'product' qualifiers of CDS features.
                    Default: blaNDM-1
    --flank         Flanking region size in bp on each side. Default: 20000
    --distance-cap  If set, gene copies on the same contig that are MORE than
                    this many bp apart (edge-to-edge gap between adjacent copies)
                    are split into SEPARATE flanking regions instead of one giant
                    region spanning them all. Copies within the cap stay merged
                    (e.g. tandem duplications). Default: off (single merged region,
                    the original behaviour). Flanks are NOT trimmed, so two split
                    regions may overlap if the cap is smaller than 2x--flank;
                    that's intentional - hide the redundant end groups in clinker.
    --coord-table   Optional combined TSV to APPEND original-coordinate rows to
                    (header written once). Default: a per-strain table
                    <output_dir>/<strain>.<gene>.region_coords.tsv

Coordinate preservation:
    A GenBank record is always numbered 1-based from the first base it contains,
    so the extracted features cannot literally carry their ~2.06 Mbp genomic
    numbering. Instead the original coordinates are preserved as metadata that
    does not interfere with clinker:
      * the DEFINITION line records the genomic window (contig:start-end),
      * every feature gets a /note="orig_coord=contig:start..end",
      * a sidecar TSV lists each target-gene copy's original coordinates.
    Local coordinates are kept in the GBK itself so clinker's synteny view works.
"""

import argparse
import os
import sys
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract flanking regions around a target gene from a GenBank file."
    )
    parser.add_argument("gbk_file",   help="Input GenBank file")
    parser.add_argument("output_dir", help="Output directory for extracted GBK(s)")
    parser.add_argument("strain",     help="Strain label for output filenames")
    parser.add_argument(
        "--gene", default="blaNDM-1",
        help="Gene name to search (case-insensitive, partial match). Default: blaNDM-1")
    parser.add_argument(
        "--flank", type=int, default=20000,
        help="Flanking region size in bp on each side. Default: 20000")
    parser.add_argument(
        "--distance-cap", type=int, default=None,
        help="Split copies farther apart than this (edge-to-edge bp) into separate "
             "regions. Default: off (one merged region).")
    parser.add_argument(
        "--coord-table", default=None,
        help="Combined original-coordinate TSV to append to. "
             "Default: per-strain table in output_dir.")
    return parser.parse_args()


def normalise(text):
    """Lowercase and strip spaces/hyphens for fuzzy matching."""
    return text.lower().replace("-", "").replace(" ", "")


def find_target_features(record, gene_query):
    """Return all CDS features whose 'gene'/'product' contains gene_query."""
    query_norm = normalise(gene_query)
    hits = []
    for feat in record.features:
        if feat.type != "CDS":
            continue
        gene    = ";".join(feat.qualifiers.get("gene",    []))
        product = ";".join(feat.qualifiers.get("product", []))
        if query_norm in normalise(f"{gene} {product}"):
            hits.append(feat)
    return hits


def cluster_features(features, distance_cap):
    """
    Group target features by genomic proximity (single-linkage).

    Sorted by start; a new group begins whenever the edge-to-edge gap to the
    previous copy exceeds distance_cap. With distance_cap=None everything stays
    in one group (original behaviour). Tandem/continuous arrays stay together
    even if their extremes are far apart, as long as each consecutive gap is
    within the cap.
    """
    feats = sorted(features, key=lambda f: int(f.location.start))
    if distance_cap is None or len(feats) <= 1:
        return [feats]
    clusters = []
    current = [feats[0]]
    current_end = int(feats[0].location.end)
    for f in feats[1:]:
        start = int(f.location.start)
        if start - current_end > distance_cap:
            clusters.append(current)
            current = [f]
            current_end = int(f.location.end)
        else:
            current.append(f)
            current_end = max(current_end, int(f.location.end))
    clusters.append(current)
    return clusters


def strand_sym(strand):
    return {1: "+", -1: "-"}.get(strand, ".")


def extract_region(record, cluster, flank, strain, gene_query, outdir, suffix):
    """Extract a GBK sub-record spanning one cluster of target features +/- flank."""
    starts = [int(f.location.start) for f in cluster]
    ends   = [int(f.location.end)   for f in cluster]

    region_start = max(0, min(starts) - flank)
    region_end   = min(len(record.seq), max(ends) + flank)

    region_record = record[region_start:region_end]

    contig_id = record.id
    # Unique LOCUS/record id per region so clinker labels two split regions
    # distinctly; unchanged (= contig_id) when a contig yields a single region.
    label = contig_id if suffix is None else f"{contig_id}_{suffix}"
    region_record.id = label
    region_record.name = label
    region_record.description = (
        f"{strain} {gene_query} region {contig_id}:{region_start + 1}-{region_end}")
    # GenBank writing needs molecule_type; slicing can drop annotations.
    region_record.annotations["molecule_type"] = record.annotations.get(
        "molecule_type", "DNA")

    # Re-clip features to the window, re-offset to local coords, and stamp each
    # with its ORIGINAL genomic coordinate as a /note.
    new_features = []
    for feat in record.features:
        try:
            fstart = int(feat.location.start)
            fend   = int(feat.location.end)
        except (TypeError, ValueError):
            continue
        if fend < region_start or fstart > region_end:
            continue

        new_start = max(fstart, region_start) - region_start
        new_end   = min(fend,   region_end)   - region_start

        quals = {k: list(v) for k, v in feat.qualifiers.items()}  # copy, don't mutate
        if feat.type == "source":
            quals["strain"] = [strain]
            quals.setdefault("note", []).append(
                f"extracted_region={contig_id}:{region_start + 1}-{region_end}")
        else:
            quals.setdefault("note", []).append(
                f"orig_coord={contig_id}:{fstart + 1}..{fend}({strand_sym(feat.location.strand)})")

        new_features.append(SeqFeature(
            FeatureLocation(new_start, new_end, strand=feat.location.strand),
            type=feat.type, qualifiers=quals))

    region_record.features = new_features

    safe_gene = gene_query.replace("/", "_").replace(" ", "_")
    outfile = os.path.join(
        outdir, f"{strain}.{contig_id}.{safe_gene}.{region_start + 1}-{region_end}.region.gbk")
    SeqIO.write(region_record, outfile, "genbank")

    copies = [(int(f.location.start) + 1, int(f.location.end), strand_sym(f.location.strand))
              for f in cluster]
    return outfile, region_start + 1, region_end, copies


def main():
    args = parse_args()

    if not os.path.isfile(args.gbk_file):
        sys.exit(f"ERROR: Input file not found: {args.gbk_file}")
    os.makedirs(args.output_dir, exist_ok=True)

    safe_gene = args.gene.replace("/", "_").replace(" ", "_")
    coord_path = args.coord_table or os.path.join(
        args.output_dir, f"{args.strain}.{safe_gene}.region_coords.tsv")
    append = bool(args.coord_table) and os.path.isfile(coord_path) and os.path.getsize(coord_path) > 0

    records_processed = 0
    records_with_hits = 0
    coord_rows = []

    for record in SeqIO.parse(args.gbk_file, "genbank"):
        records_processed += 1
        hits = find_target_features(record, args.gene)
        if not hits:
            print(f"  [{record.id}] No CDS features matching '{args.gene}' - skipping.")
            continue

        records_with_hits += 1
        clusters = cluster_features(hits, args.distance_cap)
        n_regions = len(clusters)
        for idx, cluster in enumerate(clusters, 1):
            suffix = idx if n_regions > 1 else None
            outfile, reg_start, reg_end, copies = extract_region(
                record, cluster, args.flank, args.strain, args.gene, args.output_dir, suffix)
            print(f"  [{record.id}] region {idx}/{n_regions}: "
                  f"{len(copies)} copy/copies of '{args.gene}', "
                  f"{reg_start}-{reg_end} (+/-{args.flank} bp) -> {os.path.basename(outfile)}")
            for cs, ce, strand in copies:
                coord_rows.append([os.path.basename(outfile), args.strain, record.id,
                                   str(reg_start), str(reg_end), str(len(copies)),
                                   args.gene, str(cs), str(ce), strand])

    if coord_rows:
        mode = "a" if append else "w"
        with open(coord_path, mode, encoding="utf-8") as fh:
            if not append:
                fh.write("region_file\tstrain\tcontig\tregion_start\tregion_end\t"
                         "copies_in_region\tgene\tcopy_orig_start\tcopy_orig_end\tstrand\n")
            for row in coord_rows:
                fh.write("\t".join(row) + "\n")
        print(f"\nOriginal-coordinate table -> {coord_path}")

    print(f"\nDone. Processed {records_processed} record(s); "
          f"extracted from {records_with_hits} record(s); "
          f"wrote {len(coord_rows)} gene-copy row(s).")

    if records_with_hits == 0:
        sys.exit(f"WARNING: '{args.gene}' was not found in any record. "
                 f"Check the gene name or inspect the GBK file manually.")


if __name__ == "__main__":
    main()
