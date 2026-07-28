#!/usr/bin/env python3
"""
subsample_by_st.py
==================
Pick a background sub-sample for phylogeny/heatmap that removes *clonal
redundancy* (the real problem) instead of geographic rarity.

Order of operations (this order matters):
  1. PROTECT a fixed set (NDM-1 strains, references, focal isolate) - always kept,
     bypassing every cap.
  2. ST DEDUP background: keep 1 representative per sequence type (or a few for
     STs that share a lineage with the protected/NDM-1 strains - "context STs").
     The representative is chosen by assembly quality, not at random, so the
     result is deterministic and your tree inputs are the best assemblies.
  3. ST-AWARE COUNTRY CAP: where a country still exceeds --max-per-country, drop
     the strains whose ST is already represented in another country FIRST, so a
     country's unique lineages survive and only redundant ones are trimmed.

Inputs
------
  --metadata   accessions.tsv (columns: dataset, accession, strain, ...,
               assembly_level, contig_n50, geo_location, ...)
  --mlst       all.mlst.tsv from `mlst` (no header; col1=FILE, col3=ST, 4..=alleles)
  --protect    file of strain names to always keep (e.g. NDM-1 strain list)
  --protect-names  comma list of extra names to always keep (e.g. PAO1,PA14,LESB58,PA445)

Key options
-----------
  --max-per-st         background reps per ST (default 1)
  --context-per-st     reps for STs that also contain a protected/NDM-1 strain
                       (default 2; gives focal lineages visible non-NDM context)
  --max-per-country    ST-aware country cap (default 3)

Outputs
-------
  --out-strains  final strain set (protected + chosen background), one per line
  --out-table    audit TSV: strain, ST, country, dataset, assembly_level,
                 contig_n50, role, kept, drop_reason
"""

import argparse
import csv
import os
import sys

_LEVEL_RANK = {  # higher = better assembly
    "complete genome": 4, "complete": 4, "chromosome": 3,
    "scaffold": 2, "contig": 1, "na": 0, "": 0,
}


def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)


def parse_country(geo):
    if not geo or geo.strip().upper() == "NA":
        return "Unknown"
    return geo.split(":")[0].strip() or "Unknown"


def st_key(st, alleles):
    """Stable ST key. Clean integer -> 'ST<n>'. Novel/partial -> allelic profile,
    so different novel profiles stay distinct instead of all collapsing."""
    st = (st or "").strip()
    if st.isdigit():
        return f"ST{st}"
    prof = "_".join(a.strip().replace("~", "").replace("?", "").replace("*", "")
                    for a in alleles)
    prof = prof.strip("_")
    return f"profile:{prof}" if prof else "ST_unknown"


def load_mlst(path):
    """strain -> ST key. mlst output is headerless & positional."""
    out = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            strain = os.path.basename(cols[0])
            for ext in (".fasta", ".fna", ".fa"):
                if strain.endswith(ext):
                    strain = strain[: -len(ext)]
                    break
            out[strain] = st_key(cols[2], cols[3:])
    return out


def safe_strain_name(strain, accession):
    """Match the bash download loop's naming: each of ' /\\:*?"<>|' -> '_'
    (no collapsing), and fall back to the accession when strain is NA/empty.
    This is what the FASTA files and MLST FILE basenames are keyed on."""
    s = "".join("_" if c in ' /\\:*?"<>|' else c for c in str(strain))
    if not s or s == "NA":
        s = accession
    return s


def load_metadata(path):
    rows = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            s = (r.get("strain") or "").strip()
            acc = (r.get("accession") or "").strip()
            key = safe_strain_name(s, acc)
            if not key:
                continue
            try:
                n50 = int((r.get("contig_n50") or "0").strip() or 0)
            except ValueError:
                n50 = 0
            rows[key] = {
                "accession": acc,
                "dataset": (r.get("dataset") or "").strip(),
                "country": parse_country(r.get("geo_location", "")),
                "level": (r.get("assembly_level") or "NA").strip(),
                "n50": n50,
            }
    return rows


def quality(meta):
    return (_LEVEL_RANK.get(meta["level"].lower(), 0), meta["n50"])


def read_list(path):
    if not path or not os.path.isfile(path):
        return []
    with open(path, encoding="utf-8") as fh:
        return [x.strip() for x in fh if x.strip()]


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        prog="subsample_by_st.py",
        description="ST-aware background subsampling (dedup clones, then ST-aware country cap).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--mlst", required=True)
    ap.add_argument("--protect", default=None, help="file of strains to always keep")
    ap.add_argument("--protect-names", default="", help="comma list of extra strains to keep")
    ap.add_argument("--max-per-st", type=int, default=1)
    ap.add_argument("--context-per-st", type=int, default=2)
    ap.add_argument("--max-per-country", type=int, default=3)
    ap.add_argument("--out-strains", required=True)
    ap.add_argument("--out-table", required=True)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    for f in (args.metadata, args.mlst):
        if not os.path.isfile(f):
            eprint(f"ERROR: not found: {f}")
            sys.exit(2)

    mlst = load_mlst(args.mlst)
    meta = load_metadata(args.metadata)

    protect = set(read_list(args.protect))
    protect |= {x.strip() for x in args.protect_names.split(",") if x.strip()}

    # STs that contain a protected/NDM-1 strain -> keep more background context there
    context_sts = {mlst.get(s) for s in protect if mlst.get(s)}

    # Background pool = metadata rows tagged BACKGROUND and not protected
    background = {s: m for s, m in meta.items()
                  if m["dataset"].upper() == "BACKGROUND" and s not in protect}

    # ---- 2. ST dedup -------------------------------------------------------
    by_st = {}
    for s, m in background.items():
        by_st.setdefault(mlst.get(s, "ST_unknown"), []).append(s)

    audit = {}            # strain -> [role, kept, reason]
    selected_bg = []
    for st, members in by_st.items():
        members.sort(key=lambda s: (-quality(meta[s])[0], -quality(meta[s])[1],
                                    meta[s]["accession"]))
        cap = args.context_per_st if st in context_sts else args.max_per_st
        # spread context reps across countries when possible
        chosen, seen_countries = [], set()
        for s in members:
            if len(chosen) >= cap:
                break
            c = meta[s]["country"]
            if cap > 1 and c in seen_countries and len(members) > cap:
                continue  # prefer a new country first
            chosen.append(s)
            seen_countries.add(c)
        for s in members:
            if len(chosen) >= cap:
                break
            if s not in chosen:
                chosen.append(s)
        for s in members:
            if s in chosen:
                selected_bg.append(s)
                audit[s] = ["background", "yes",
                            f"ST rep ({st}; {'context' if st in context_sts else 'dedup'})"]
            else:
                audit[s] = ["background", "no", f"redundant in {st}"]

    # ---- 3. ST-aware country cap ------------------------------------------
    st_of = {s: mlst.get(s, "ST_unknown") for s in selected_bg}
    by_country = {}
    for s in selected_bg:
        by_country.setdefault(meta[s]["country"], []).append(s)

    final_bg = set(selected_bg)
    for country, members in by_country.items():
        if len(members) <= args.max_per_country:
            continue
        # how many other countries also hold each ST (redundancy score)
        def redundancy(s):
            st = st_of[s]
            return sum(1 for o in selected_bg
                       if o != s and st_of[o] == st and meta[o]["country"] != country)
        # drop most-redundant, lowest-quality first; keep unique STs
        ordered = sorted(members,
                         key=lambda s: (redundancy(s), quality(meta[s])[0], quality(meta[s])[1]),
                         reverse=True)
        keep = set(ordered[:args.max_per_country])
        for s in members:
            if s not in keep:
                final_bg.discard(s)
                reason = ("country cap (ST redundant elsewhere)"
                          if redundancy(s) > 0 else "country cap (unique ST dropped)")
                audit[s] = ["background", "no", reason]

    # ---- assemble + write --------------------------------------------------
    for s in protect:
        audit[s] = ["protected", "yes", "always kept"]

    final = sorted(protect) + sorted(final_bg)
    with open(args.out_strains, "w", encoding="utf-8") as fh:
        for s in final:
            fh.write(s + "\n")

    with open(args.out_table, "w", encoding="utf-8") as fh:
        fh.write("strain\tST\tcountry\tdataset\tassembly_level\tcontig_n50\trole\tkept\tdrop_reason\n")
        all_strains = set(meta) | protect
        for s in sorted(all_strains):
            m = meta.get(s, {})
            role, kept, reason = audit.get(s, ["background", "no", "not in background pool"])
            fh.write("\t".join([
                s, mlst.get(s, "NA"), m.get("country", "NA"), m.get("dataset", "NA"),
                m.get("level", "NA"), str(m.get("n50", "NA")), role, kept, reason]) + "\n")

    distinct_sts = len({st_of[s] for s in final_bg})
    eprint("==================== subsample_by_st ====================")
    eprint(f"  protected (always kept) : {len(protect)}")
    eprint(f"  background ST groups     : {len(by_st)}")
    eprint(f"  background after dedup   : {len(selected_bg)}")
    eprint(f"  background after cap     : {len(final_bg)}  ({distinct_sts} distinct STs)")
    eprint(f"  FINAL total              : {len(final)}")
    eprint(f"  context STs (extra reps) : {sorted(context_sts) if context_sts else 'none'}")
    eprint(f"  strains -> {args.out_strains}")
    eprint(f"  audit   -> {args.out_table}")
    eprint("=========================================================")
    return 0


if __name__ == "__main__":
    sys.exit(main())
