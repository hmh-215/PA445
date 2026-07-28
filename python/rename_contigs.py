#!/usr/bin/env python3
"""
rename_contigs.py
=================
Rewrite a FASTA into an RGI/Prodigal-safe form, regardless of how the original
NCBI/SPAdes/Unicycler genome was formatted.

Why bare `contigNNNNN` names
----------------------------
RGI runs Prodigal, which names ORFs `{contig}_{ORFnumber}`, then tries to
recover the contig name from the ORF id. Its recovery is not robust to
underscores in the contig name: depending on the code path it may strip the
last `_<field>` or split on the first `_`, so ANY underscore in a contig name
risks:
    "Requested rname <mangled_name> does not exist! Please check your FASTA file."
A contig name containing NO underscores (e.g. `contig00001`) is recovered
correctly by every one of those rules, because the only `_` in the ORF id is
the one Prodigal added before the ORF number.

So contigs are renamed to `contig00001`, `contig00002`, ... per file. The
strain/isolate identity is NOT encoded in the contig name - the FASTA filename
already carries that - and is instead recorded in a single shared mapping TSV.

What it also normalises
-----------------------
  * re-wraps every sequence to a uniform width (faidx-safe).
  * uppercases bases; converts anything that is not a valid IUPAC nucleotide
    code to 'N' (handles stray '-', '*', '.', etc.).

Mapping TSV (one combined file for all samples)
-----------------------------------------------
Columns: sample <TAB> new_name <TAB> original_header <TAB> length_bp
With --append, rows are added to the given --map file (header written only when
the file is new/empty), so a whole batch of genomes lands in one table. This is
what lets you trace a renamed contig back to its original header later - e.g.
to tell chromosome from plasmid in the RGI 'Contig' column.

Usage
-----
  rename_contigs.py INPUT.fasta OUTPUT.fasta \
      [--sample NAME] [--map MAP.tsv] [--append] \
      [--prefix contig] [--pad 5] [--width 80] [--start 1]

Exit codes
----------
  0  success
  2  bad input (missing/empty/not-a-FASTA)
"""

import argparse
import os
import re
import sys

# Valid IUPAC nucleotide codes (everything else -> N).
_VALID = set("ACGTRYSWKMBDHVN")
# A contig prefix must contain NO underscores (and no other odd characters),
# otherwise we reintroduce the very ambiguity we are trying to avoid.
_PREFIX_STRIP = re.compile(r"[^A-Za-z0-9]+")
_FASTA_EXT = re.compile(r"\.(fa|fna|fasta|fas|seq)(\.gz)?$", re.IGNORECASE)


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def clean_prefix(prefix):
    """Strip the prefix down to [A-Za-z0-9] so the contig name has no underscore."""
    p = _PREFIX_STRIP.sub("", prefix)
    return p or "contig"


def default_sample(input_path):
    return _FASTA_EXT.sub("", os.path.basename(input_path))


def clean_seq(seq):
    return "".join(ch if ch in _VALID else "N" for ch in seq.upper())


def wrap(seq, width):
    if width <= 0:
        return [seq]
    return [seq[i:i + width] for i in range(0, len(seq), width)]


def read_fasta(path):
    """Yield (header_without_'>', sequence) pairs, streaming."""
    header, chunks = None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header, chunks = line[1:].strip(), []
            else:
                chunks.append(line)
        if header is not None:
            yield header, "".join(chunks)


def parse_args(argv=None):
    ap = argparse.ArgumentParser(
        prog="rename_contigs.py",
        description="Rename contigs to bare contigNNNNN and normalise a FASTA so RGI/Prodigal cannot choke on it.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "example (single file):\n"
            "  rename_contigs.py in.fasta out.fasta --sample STRAIN --map all.contigmap.tsv --append\n"
        ),
    )
    ap.add_argument("input", metavar="INPUT.fasta", help="source FASTA (any formatting)")
    ap.add_argument("output", metavar="OUTPUT.fasta", help="RGI-safe FASTA to write")
    ap.add_argument("--sample", default=None,
                    help="sample/strain label recorded in the map "
                         "(default: input filename without extension); "
                         "NOT used in the contig name")
    ap.add_argument("--map", default=None, metavar="MAP.tsv",
                    help="combined mapping TSV (default: <OUTPUT>.contigmap.tsv)")
    ap.add_argument("--append", action="store_true",
                    help="append to --map instead of overwriting "
                         "(header written only when the file is new/empty)")
    ap.add_argument("--prefix", default="contig",
                    help="contig name prefix; underscores are stripped (default: contig)")
    ap.add_argument("--pad", type=int, default=5, metavar="N",
                    help="zero-padding width for the contig number (default: 5)")
    ap.add_argument("--width", type=int, default=80, metavar="N",
                    help="sequence line width; 0 = single line (default: 80). "
                         "Use 0 for RGI inputs: one line per contig makes the "
                         "faidx 'inconsistent line length' error impossible.")
    ap.add_argument("--validate", action="store_true",
                    help="after writing, build a faidx index on the OUTPUT with "
                         "pysam/htslib (the same library RGI uses) and fail if it "
                         "cannot be indexed")
    ap.add_argument("--start", type=int, default=1, metavar="N",
                    help="first contig index (default: 1)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if not os.path.isfile(args.input):
        eprint(f"ERROR: input FASTA not found: {args.input}")
        sys.exit(2)
    if os.path.getsize(args.input) == 0:
        eprint(f"ERROR: input FASTA is empty: {args.input}")
        sys.exit(2)

    prefix = clean_prefix(args.prefix)
    sample = args.sample if args.sample else default_sample(args.input)
    map_path = args.map or (args.output + ".contigmap.tsv")

    os.makedirs(os.path.dirname(os.path.abspath(args.output)) or ".", exist_ok=True)
    map_dir = os.path.dirname(os.path.abspath(map_path)) or "."
    os.makedirs(map_dir, exist_ok=True)

    # Header is written when overwriting, or when appending to a new/empty file.
    map_mode = "a" if args.append else "w"
    write_header = True
    if args.append and os.path.isfile(map_path) and os.path.getsize(map_path) > 0:
        write_header = False

    n = 0
    idx = args.start
    try:
        with open(args.output, "w", encoding="utf-8") as out, \
             open(map_path, map_mode, encoding="utf-8") as mp:
            if write_header:
                mp.write("sample\tnew_name\toriginal_header\tlength_bp\n")
            for original_header, seq in read_fasta(args.input):
                if not seq:
                    eprint(f"WARNING: [{sample}] skipping '{original_header}' (no sequence).")
                    continue
                new_name = f"{prefix}{idx:0{args.pad}d}"
                seq = clean_seq(seq)
                out.write(f">{new_name}\n")
                for chunk in wrap(seq, args.width):
                    out.write(chunk + "\n")
                mp.write(f"{sample}\t{new_name}\t{original_header}\t{len(seq)}\n")
                idx += 1
                n += 1
    except OSError as exc:
        eprint(f"ERROR: writing output failed: {exc}")
        sys.exit(2)

    if n == 0:
        eprint(f"ERROR: no contigs written - is {args.input} really a FASTA?")
        sys.exit(2)

    if args.validate:
        try:
            import pysam  # provided by the rgi_env / bioinformatics envs
            pysam.faidx(args.output)  # builds <output>.fai via htslib
            pysam.FastaFile(args.output).close()
        except ImportError:
            eprint("WARNING: --validate requested but pysam is not importable here; "
                   "skipping faidx check.")
        except Exception as exc:  # noqa: BLE001 - surface htslib's reason verbatim
            eprint(f"ERROR: faidx validation FAILED on {args.output}: {exc}")
            eprint("       RGI would hit 'Requested rname ... does not exist' on this file.")
            sys.exit(2)

    eprint(f"rename_contigs.py: [{sample}] {args.input} -> {args.output} "
           f"({n} contigs as {prefix}{'#' * args.pad}; map -> {map_path})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
