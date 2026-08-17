#!/usr/bin/env python3
"""
Clean up a gffread -y protein FASTA before handing it to DIAMOND:
  1. Replace any character outside the standard protein alphabet in
     sequence lines with 'X', leaving headers untouched.
  2. Re-wrap sequences to a fixed line width.

Why (1): gffread inserts non-standard placeholder characters -- most
commonly '.' -- for CDS regions it can't cleanly translate, most often
an incomplete/partial codon at the edge of a CDS that doesn't start or
end on a complete codon boundary. This is common in whole-genome NCBI/
Gnomon annotations (lots of partial gene models flagged "partial=true")
-- it's what produced DIAMOND's hard failure on the dsim
(GCF_016746395.2) proteome:
    Error: Invalid character in sequence: '.'
DIAMOND's protein alphabet accepts the 20 standard amino acids, the
ambiguity/rare codes B/X/Z/J/U/O, and a trailing '*' for stop codons;
anything else gets replaced with 'X' (unknown residue) here -- the same
thing gffread's '.' was already standing in for.

Why (2): re-wrapping guards against (and helps diagnose, via the
longest-sequence report below) any pathological single-very-long-line
behavior in downstream tools. D. melanogaster's sallimus/titin-like
proteins in particular can be tens of thousands of residues on one
unwrapped line -- worth ruling in/out for the dmel-side
`diamond makedb` hang this script was written to investigate (2026-08-16:
dmel_proteins.fasta consistently hung diamond at "Loading sequences..."
across multiple nodes/filesystems/thread counts, while dsim's file
failed fast on the character issue above -- same diamond process,
different content, so the content is the prime suspect for both).

Usage:
    python sanitize_protein_fasta.py in.fasta out.fasta [--width 60]
"""
import argparse
import sys

VALID = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("in_fasta")
    ap.add_argument("out_fasta")
    ap.add_argument("--width", type=int, default=60, help="line-wrap width (default: 60)")
    args = ap.parse_args()

    n_seqs = 0
    n_seqs_changed = 0
    n_chars_changed = 0
    changed_chars = {}
    longest = (0, None)

    def flush(out_fh, header, seq_chars):
        nonlocal n_seqs, n_seqs_changed, n_chars_changed, longest
        if header is None:
            return
        n_seqs += 1
        seq_changed = False
        cleaned = []
        for ch in seq_chars:
            if ch.upper() in VALID:
                cleaned.append(ch)
            else:
                cleaned.append("X")
                n_chars_changed += 1
                seq_changed = True
                changed_chars[ch] = changed_chars.get(ch, 0) + 1
        if seq_changed:
            n_seqs_changed += 1
        if len(cleaned) > longest[0]:
            longest = (len(cleaned), header)
        out_fh.write(header + "\n")
        for i in range(0, len(cleaned), args.width):
            out_fh.write("".join(cleaned[i:i + args.width]) + "\n")

    header = None
    seq_chars = []
    with open(args.in_fasta) as in_fh, open(args.out_fasta, "w") as out_fh:
        for line in in_fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush(out_fh, header, seq_chars)
                header = line
                seq_chars = []
            else:
                seq_chars.extend(line)
        flush(out_fh, header, seq_chars)

    sys.stderr.write(
        f"[{args.in_fasta}] {n_seqs} sequences, {n_seqs_changed} contained "
        f"non-standard characters, {n_chars_changed} character(s) replaced "
        f"with 'X', re-wrapped at {args.width} cols -> {args.out_fasta}\n"
    )
    if changed_chars:
        breakdown = ", ".join(f"{c!r}:{n}" for c, n in sorted(changed_chars.items()))
        sys.stderr.write(f"  replaced character counts: {breakdown}\n")
    if longest[1]:
        sys.stderr.write(f"  longest sequence: {longest[0]} aa, {longest[1]}\n")
    if n_seqs == 0:
        sys.stderr.write("WARNING: 0 sequences read -- check the input file.\n")


if __name__ == "__main__":
    main()
