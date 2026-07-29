#!/usr/bin/env python3
"""
ncbi_accessions_filter.py
=========================
Filter the *background* genome set produced by `datasets summary genome ...`
for a comparative-genomics workflow, while passing the *NDM* set through
untouched.

Design contract (matches PA_comparative_1.bash)
-----------------------------------------------
  * Two NCBI `datasets summary genome` JSON files are taken as input:
        1. NDM_JSON          - candidate blaNDM-1 genomes
        2. BACKGROUND_JSON   - everything else
  * NDM accessions are NEVER filtered. Every accession in NDM_JSON is written
    verbatim to `accessions_NDM.txt`. (True-NDM-1 confirmation happens later,
    via prokka_rgi_reconcile.py - not here.)
  * BACKGROUND accessions are filtered by the criteria below and written to
    `accessions_background.txt`.

Background filter criteria (all optional, opt-in via flags)
-----------------------------------------------------------
  --geography            keep only genomes that report a usable geo_loc_name
  --host HOST            keep only genomes whose biosample host matches HOST
                         ("Homo sapiens" and "human" are treated as synonyms)
  --max-contigs N        drop genomes with more than N contigs (off by default)
  --min-n50 N            drop genomes with contig_N50 below N (off by default)
  --dedup-bioproj        keep only ONE genome per BioProject (reduces clonal
                         over-representation); genomes with no BioProject are
                         each kept

Outputs (written to --out-dir, default ".")
-------------------------------------------
  accessions_NDM.txt          one NDM accession per line (unfiltered)
  accessions_background.txt   one surviving background accession per line
  filter_summary.txt          human-readable funnel of what was kept/dropped

Exit codes
----------
  0  success
  2  bad input (missing file, invalid JSON, no "reports" key)
  3  background result empty AND --fail-on-empty was given

Example
-------
  python3 ncbi_accessions_filter.py \\
      summary_NDM.json summary_background.json \\
      --out-dir download_logs \\
      --geography --host "Homo sapiens" --dedup-bioproj
"""

import argparse
import json
import os
import sys

# Values that mean "this field is effectively empty / unusable".
_MISSING_TOKENS = {
    "", "na", "n/a", "none", "null", "unknown", "missing",
    "not collected", "not applicable", "not provided", "not determined",
    "not available", "missing: control sample", "-",
}

# Host strings that should be treated as "human".
_HUMAN_SYNONYMS = {"homo sapiens", "human", "homo sapiens (human)"}


def eprint(*args, **kwargs):
    """Print to stderr (so it shows up in the bash log without polluting stdout)."""
    print(*args, file=sys.stderr, **kwargs)


def is_missing(value):
    """True if a metadata value should be treated as absent."""
    if value is None:
        return True
    return str(value).strip().lower() in _MISSING_TOKENS


def normalise(value):
    return str(value).strip().lower() if value is not None else ""


def load_json(path):
    """Load a JSON file, exiting with a clear message on any failure."""
    if not os.path.isfile(path):
        eprint(f"ERROR: input file not found: {path}")
        sys.exit(2)
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except json.JSONDecodeError as exc:
        eprint(f"ERROR: {path} is not valid JSON ({exc}).")
        eprint("       Check that `datasets summary genome` actually wrote a report,")
        eprint("       e.g. the query may have returned nothing or hit a network error.")
        sys.exit(2)
    except OSError as exc:
        eprint(f"ERROR: could not read {path}: {exc}")
        sys.exit(2)


def get_reports(data, source_name):
    """Return the list under data['reports'], with helpful diagnostics."""
    if not isinstance(data, dict):
        eprint(f"ERROR: {source_name}: top-level JSON is not an object.")
        sys.exit(2)
    reports = data.get("reports")
    if reports is None:
        # `datasets` emits {"total_count": 0} (no "reports") when nothing matched.
        eprint(f"WARNING: {source_name}: no 'reports' key found "
               f"(total_count={data.get('total_count', 'NA')}). Treating as empty.")
        return []
    if not isinstance(reports, list):
        eprint(f"ERROR: {source_name}: 'reports' is not a list.")
        sys.exit(2)
    return reports


def biosample_attr(biosample, key):
    """
    Robustly fetch a biosample attribute.

    NCBI `datasets` v2 usually nests metadata under biosample['attributes']
    as a list of {"name": ..., "value": ...}. Some versions expose a few of
    these as direct keys. We check both so geography/host detection does not
    silently fail.
    """
    if not isinstance(biosample, dict):
        return None
    # 1. direct key
    if key in biosample and not is_missing(biosample.get(key)):
        return biosample.get(key)
    # 2. attributes list
    for attr in biosample.get("attributes", []) or []:
        if isinstance(attr, dict) and normalise(attr.get("name")) == normalise(key):
            val = attr.get("value")
            if not is_missing(val):
                return val
    return None


def host_matches(host_value, wanted):
    """Match host flexibly; treat 'Homo sapiens' and 'human' as the same thing."""
    if is_missing(host_value):
        return False
    hv = normalise(host_value)
    wanted_n = normalise(wanted)
    if wanted_n in _HUMAN_SYNONYMS:
        return hv in _HUMAN_SYNONYMS or "homo sapiens" in hv or hv == "human"
    return wanted_n in hv  # substring match for any other requested host


def country_of(geo_value):
    """Parse 'Country: city/institute' -> 'Country' (for the summary only)."""
    if is_missing(geo_value):
        return "Unknown"
    return str(geo_value).split(":")[0].strip() or "Unknown"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="ncbi_accessions_filter.py",
        description="Filter background NCBI genomes; pass NDM genomes through unfiltered.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "example:\n"
            "  python3 ncbi_accessions_filter.py summary_NDM.json summary_background.json \\\n"
            "      --out-dir download_logs --geography --host \"Homo sapiens\" --dedup-bioproj\n"
        ),
    )
    parser.add_argument("ndm_json", metavar="NDM_JSON",
                        help="datasets-summary JSON for NDM genomes (kept unfiltered)")
    parser.add_argument("background_json", metavar="BACKGROUND_JSON",
                        help="datasets-summary JSON for background genomes (filtered)")
    parser.add_argument("--out-dir", default=".",
                        help="directory for the output .txt files (default: current dir)")
    parser.add_argument("--geography", action="store_true",
                        help="require a usable geo_loc_name on background genomes")
    parser.add_argument("--host", default=None, metavar="HOST",
                        help='require this biosample host, e.g. "Homo sapiens"')
    parser.add_argument("--max-contigs", type=int, default=None, metavar="N",
                        help="drop background genomes with more than N contigs "
                             "(default: no contig filter)")
    parser.add_argument("--min-n50", type=int, default=None, metavar="N",
                        help="drop background genomes with contig_N50 below N "
                             "(default: no N50 filter)")
    parser.add_argument("--dedup-bioproj", action="store_true",
                        help="keep only one background genome per BioProject")
    parser.add_argument("--ndm-out", default="accessions_NDM.txt", metavar="FILE",
                        help="filename for NDM accessions (default: accessions_NDM.txt)")
    parser.add_argument("--bg-out", default="accessions_background.txt", metavar="FILE",
                        help="filename for filtered background accessions "
                             "(default: accessions_background.txt)")
    parser.add_argument("--fail-on-empty", action="store_true",
                        help="exit non-zero if no background genome survives filtering")
    return parser.parse_args(argv)


def filter_background(reports, args):
    """Apply the background filters, returning (kept_accessions, funnel_counts)."""
    funnel = {
        "total": len(reports),
        "no_accession": 0,
        "drop_geography": 0,
        "drop_host": 0,
        "drop_contigs": 0,
        "drop_n50": 0,
        "drop_dedup": 0,
        "kept": 0,
    }
    seen_bioprojects = set()
    countries = {}
    kept = []

    for r in reports:
        accession = r.get("accession")
        if is_missing(accession):
            funnel["no_accession"] += 1
            continue

        info = r.get("assembly_info", {}) or {}
        stats = r.get("assembly_stats", {}) or {}
        biosample = info.get("biosample", {}) or {}

        geo = biosample_attr(biosample, "geo_loc_name")
        host = biosample_attr(biosample, "host")

        if args.geography and is_missing(geo):
            funnel["drop_geography"] += 1
            continue

        if args.host is not None and not host_matches(host, args.host):
            funnel["drop_host"] += 1
            continue

        if args.max_contigs is not None:
            contigs = stats.get("number_of_component_sequences")
            try:
                if contigs is not None and int(contigs) > args.max_contigs:
                    funnel["drop_contigs"] += 1
                    continue
            except (TypeError, ValueError):
                pass  # unpar. leave it in rather than drop on bad metadata

        if args.min_n50 is not None:
            n50 = stats.get("contig_n50")
            try:
                if n50 is not None and int(n50) < args.min_n50:
                    funnel["drop_n50"] += 1
                    continue
            except (TypeError, ValueError):
                pass

        if args.dedup_bioproj:
            bioproj = info.get("bioproject_accession")
            if not is_missing(bioproj):
                if bioproj in seen_bioprojects:
                    funnel["drop_dedup"] += 1
                    continue
                seen_bioprojects.add(bioproj)

        kept.append(accession)
        funnel["kept"] += 1
        c = country_of(geo)
        countries[c] = countries.get(c, 0) + 1

    return kept, funnel, countries


def write_lines(path, lines):
    with open(path, "w", encoding="utf-8") as fh:
        for line in lines:
            fh.write(f"{line}\n")


def main(argv=None):
    args = parse_args(argv)

    try:
        os.makedirs(args.out_dir, exist_ok=True)
    except OSError as exc:
        eprint(f"ERROR: cannot create --out-dir {args.out_dir}: {exc}")
        sys.exit(2)

    # ---- NDM set: pass through, unfiltered ----
    ndm_data = load_json(args.ndm_json)
    ndm_reports = get_reports(ndm_data, "NDM_JSON")
    ndm_accessions, seen = [], set()
    for r in ndm_reports:
        acc = r.get("accession")
        if is_missing(acc):
            eprint("WARNING: NDM report with no accession - skipped.")
            continue
        if acc not in seen:
            seen.add(acc)
            ndm_accessions.append(acc)

    # ---- Background set: filter ----
    bg_data = load_json(args.background_json)
    bg_reports = get_reports(bg_data, "BACKGROUND_JSON")
    bg_accessions, funnel, countries = filter_background(bg_reports, args)
    # de-duplicate while preserving order
    bg_accessions = list(dict.fromkeys(bg_accessions))

    ndm_path = os.path.join(args.out_dir, args.ndm_out)
    bg_path = os.path.join(args.out_dir, args.bg_out)
    summary_path = os.path.join(args.out_dir, "filter_summary.txt")

    write_lines(ndm_path, ndm_accessions)
    write_lines(bg_path, bg_accessions)

    # ---- Summary (stderr + file) ----
    active = []
    if args.geography:
        active.append("geography=required")
    if args.host:
        active.append(f'host~="{args.host}"')
    if args.max_contigs is not None:
        active.append(f"max_contigs={args.max_contigs}")
    if args.min_n50 is not None:
        active.append(f"min_n50={args.min_n50}")
    if args.dedup_bioproj:
        active.append("dedup_bioproject")
    active_str = ", ".join(active) if active else "none (pass-through)"

    lines = [
        "==================== ncbi_accessions_filter.py ====================",
        f"NDM input        : {args.ndm_json}",
        f"Background input : {args.background_json}",
        f"Active filters   : {active_str}",
        "-------------------------------------------------------------------",
        f"NDM accessions kept (unfiltered) : {len(ndm_accessions)}",
        "-------------------- background funnel -----------------------------",
        f"  candidates                : {funnel['total']}",
        f"  dropped: no accession     : {funnel['no_accession']}",
        f"  dropped: no geography     : {funnel['drop_geography']}",
        f"  dropped: host mismatch    : {funnel['drop_host']}",
        f"  dropped: too many contigs : {funnel['drop_contigs']}",
        f"  dropped: N50 too low      : {funnel['drop_n50']}",
        f"  dropped: bioproject dedup : {funnel['drop_dedup']}",
        f"  KEPT                      : {funnel['kept']}",
        "-------------------- background by country -------------------------",
    ]
    for country, n in sorted(countries.items(), key=lambda kv: (-kv[1], kv[0])):
        lines.append(f"  {country:<30} {n}")
    lines += [
        "-------------------------------------------------------------------",
        f"NDM accessions  -> {ndm_path}",
        f"Background       -> {bg_path}",
        "===================================================================",
    ]
    report_text = "\n".join(lines)
    eprint(report_text)
    try:
        with open(summary_path, "w", encoding="utf-8") as fh:
            fh.write(report_text + "\n")
    except OSError as exc:
        eprint(f"WARNING: could not write summary file {summary_path}: {exc}")

    if funnel["kept"] == 0:
        eprint("WARNING: no background genome passed the filters.")
        eprint("         Common causes: host string mismatch, or geography missing")
        eprint("         in this `datasets` version. Re-run without --host/--geography")
        eprint("         to sanity-check, or inspect filter_summary.txt.")
        if args.fail_on_empty:
            sys.exit(3)

    return 0


if __name__ == "__main__":
    sys.exit(main())
