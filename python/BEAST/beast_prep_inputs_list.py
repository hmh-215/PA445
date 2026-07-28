#!/usr/bin/env python3
"""
prep_beast_inputs.py
====================
Assemble the taxon set and tip dates for a within-lineage BEAST run.

Time-calibrated dating is only valid WITHIN a clonal lineage, so this selects
the strains sharing PA445's sequence type (or an explicit --focal-st), pulls
their collection dates from the metadata, converts them to decimal years, and
flags anything undated. It does NOT date across the whole species.

Inputs
------
  --mlst        all.mlst.tsv (headerless; col1=FILE, col3=ST, 4..=alleles)
  --metadata    accessions.tsv (needs: strain, accession, collection_date)
  --focal       focal isolate name whose ST defines the lineage (default: PA445)
  --focal-st    override: use this ST key (e.g. ST308) instead of the focal's ST
  --focal-date  collection date for the focal isolate if absent from metadata
  --extra       optional file of extra strain names to force-include
  --min-strains minimum dated taxa to proceed (default 5; dating below this is unreliable)

Outputs
-------
  --out-strains  focal-lineage strains that have a usable date (BEAST taxon set)
  --out-dates    TSV: taxon, decimal_date, raw_date, precision
  --out-report   TSV: every focal-ST strain with date status (incl. undated)

Exit codes: 0 ok; 2 bad input; 3 too few dated taxa.
"""

import argparse
import csv
import os
import re
import sys


def eprint(*a, **k):
    print(*a, file=sys.stderr, **k)


def safe_strain_name(strain, accession):
    s = "".join("_" if c in ' /\\:*?"<>|' else c for c in str(strain))
    if not s or s == "NA":
        s = accession
    return s


def st_key(st, alleles):
    st = (st or "").strip()
    if st.isdigit():
        return f"ST{st}"
    prof = "_".join(a.strip().replace("~", "").replace("?", "").replace("*", "") for a in alleles).strip("_")
    return f"profile:{prof}" if prof else "ST_unknown"


def load_mlst(path):
    out = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            s = os.path.basename(cols[0])
            for ext in (".fasta", ".fna", ".fa"):
                if s.endswith(ext):
                    s = s[:-len(ext)]
                    break
            out[s] = st_key(cols[2], cols[3:])
    return out


_DAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]


def _leap(y):
    return y % 4 == 0 and (y % 100 != 0 or y % 400 == 0)


def to_decimal(year, month=None, day=None):
    if month is None:
        return year + 0.5, "year"          # mid-year point for year-only
    days = _DAYS[:]
    if _leap(year):
        days[1] = 29
    doy = sum(days[:month - 1]) + (day if day else 15)
    total = 366 if _leap(year) else 365
    return round(year + (doy - 1) / total, 4), ("day" if day else "month")


def parse_date(raw):
    """Return (decimal_year, precision) or (None, reason)."""
    if not raw or raw.strip().upper() in {"", "NA", "MISSING", "UNKNOWN", "NOT COLLECTED", "NOT APPLICABLE"}:
        return None, "missing"
    raw = raw.strip()
    # range like 2018/2019 or 2018-2019 (years only) -> midpoint
    m = re.fullmatch(r"(\d{4})\s*[/-]\s*(\d{4})", raw)
    if m:
        y1, y2 = int(m.group(1)), int(m.group(2))
        return round((y1 + y2) / 2 + 0.5, 4), "year-range"  # midpoint of the span
    # ISO-ish: YYYY[-MM[-DD]]
    m = re.fullmatch(r"(\d{4})(?:-(\d{1,2}))?(?:-(\d{1,2}))?", raw)
    if m:
        y = int(m.group(1))
        mo = int(m.group(2)) if m.group(2) else None
        d = int(m.group(3)) if m.group(3) else None
        if mo and not (1 <= mo <= 12):
            mo = None
        dec, prec = to_decimal(y, mo, d)
        return dec, prec
    # bare year embedded
    m = re.search(r"(19|20)\d{2}", raw)
    if m:
        return int(m.group(0)) + 0.5, "year"
    return None, f"unparsed:{raw}"


def parse_args(argv=None):
    ap = argparse.ArgumentParser(prog="prep_beast_inputs.py",
                                 description="Select focal-ST taxa and tip dates for BEAST.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mlst", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--focal", default="PA445")
    ap.add_argument("--focal-st", default=None)
    ap.add_argument("--focal-date", default=None)
    ap.add_argument("--extra", default=None)
    ap.add_argument("--min-strains", type=int, default=5)
    ap.add_argument("--out-strains", required=True)
    ap.add_argument("--out-dates", required=True)
    ap.add_argument("--out-report", required=True)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    for f in (args.mlst, args.metadata):
        if not os.path.isfile(f):
            eprint(f"ERROR: not found: {f}")
            sys.exit(2)

    mlst = load_mlst(args.mlst)
    focal_st = args.focal_st or mlst.get(args.focal)
    if not focal_st:
        eprint(f"ERROR: could not determine focal ST (no MLST entry for '{args.focal}'); "
               f"pass --focal-st explicitly.")
        sys.exit(2)
    eprint(f"Focal lineage: {focal_st} (from {args.focal})")

    # metadata: safe-name -> raw collection_date
    dates_raw = {}
    with open(args.metadata, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = safe_strain_name((r.get("strain") or "").strip(), (r.get("accession") or "").strip())
            dates_raw[key] = (r.get("collection_date") or "NA").strip()

    members = sorted(s for s, st in mlst.items() if st == focal_st)
    if args.extra and os.path.isfile(args.extra):
        for s in open(args.extra, encoding="utf-8"):
            s = s.strip()
            if s and s not in members:
                members.append(s)

    # Normalize the focal date: an empty/whitespace string (e.g. from an unset
    # shell variable) must become None, NOT be treated as a real value and then
    # silently fall through to a missing metadata lookup.
    focal_date = (args.focal_date or "").strip() or None
    if args.focal_date is not None and focal_date is None:
        eprint(f"WARNING: --focal-date was given but is empty for focal '{args.focal}'; "
               f"falling back to metadata (which usually lacks the focal isolate).")

    report, dated = [], []
    for s in members:
        raw = focal_date if (s == args.focal and focal_date) else dates_raw.get(s, "NA")
        dec, prec = parse_date(raw)
        if dec is None:
            report.append((s, focal_st, raw, "NA", prec, "undated-excluded"))
        else:
            report.append((s, focal_st, raw, f"{dec}", prec, "dated"))
            dated.append((s, dec, raw, prec))

    with open(args.out_strains, "w", encoding="utf-8") as fh:
        for s, *_ in dated:
            fh.write(s + "\n")
    with open(args.out_dates, "w", encoding="utf-8") as fh:
        fh.write("taxon\tdecimal_date\traw_date\tprecision\n")
        for s, dec, raw, prec in dated:
            fh.write(f"{s}\t{dec}\t{raw}\t{prec}\n")
    with open(args.out_report, "w", encoding="utf-8") as fh:
        fh.write("strain\tST\traw_date\tdecimal_date\tprecision\tstatus\n")
        for row in report:
            fh.write("\t".join(row) + "\n")

    # The focal isolate is the whole point of the run: if it belongs to the lineage
    # but has no usable date, it would be silently dropped from the BEAST XML. Fail.
    dated_names = {s for s, *_ in dated}
    if args.focal in members and args.focal not in dated_names:
        eprint(f"ERROR: focal isolate '{args.focal}' is in lineage {focal_st} but has no usable date "
               f"(focal-date={focal_date!r}, metadata={dates_raw.get(args.focal, 'NA')!r}). "
               f"It would be dropped from the tree. Pass a valid --focal-date.")
        sys.exit(3)

    span = (max(d for _, d, _, _ in dated) - min(d for _, d, _, _ in dated)) if dated else 0
    eprint("==================== prep_beast_inputs ====================")
    eprint(f"  focal ST                : {focal_st}")
    eprint(f"  ST members found        : {len(members)}")
    eprint(f"  dated (BEAST taxa)      : {len(dated)}")
    eprint(f"  undated (excluded)      : {len(members) - len(dated)}")
    eprint(f"  sampling span (years)   : {span:.2f}")
    eprint(f"  strains -> {args.out_strains}")
    eprint(f"  dates   -> {args.out_dates}")
    eprint("===========================================================")
    if span < 1 and dated:
        eprint("WARNING: sampling span < 1 year - temporal signal will likely be too weak to date.")
    if len(dated) < args.min_strains:
        eprint(f"ERROR: only {len(dated)} dated taxa (< --min-strains {args.min_strains}); "
               f"dating this lineage is not advisable.")
        sys.exit(3)
    return 0


if __name__ == "__main__":
    sys.exit(main())