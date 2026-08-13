#!/usr/bin/env python3
"""
Build a transcript_id -> gene_id map directly from a GTF file's attribute
column, instead of parsing species-specific conventions out of a
translation FASTA header (get_id_map.py's approach, which only works for
FlyBase-style `parent=FBgn...` headers).

Generic `key "value";` attribute parsing means this works on any
well-formed GTF, in particular:
  - FlyBase GTFs, e.g. gene_id "FBgn0031081"; transcript_id "FBtr0070000";
  - NCBI RefSeq/Gnomon GTFs, e.g. gene_id "LOC120284240";
    transcript_id "XM_039291137.1"; ...
so it drops in for D. simulans annotations sourced from NCBI (no FBgn IDs,
no FlyBase header conventions) just as well as FlyBase ones.

Pair with `gffread -y` protein extraction: gffread's default protein
FASTA header for a GTF-derived transcript is the transcript_id itself, so
this script's output lines up directly with those FASTA headers as
DIAMOND/BLAST query/subject IDs -- feed the result straight into
filter_rbh.py as an id_map.

Usage:
    python gtf_id_map.py annotation.gtf > id_map.tsv
"""
import re
import sys

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} annotation.gtf > id_map.tsv")

    gtf_path = sys.argv[1]
    transcript_to_gene = {}
    n_attr_rows = 0

    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            transcript_id = attrs.get("transcript_id")
            gene_id = attrs.get("gene_id")
            # Gene-level feature rows (e.g. FlyBase/NCBI "gene" rows) carry
            # gene_id but an empty transcript_id ("") -- skip those, we only
            # want rows that tie a transcript to its gene.
            if not transcript_id or not gene_id:
                continue
            n_attr_rows += 1
            transcript_to_gene.setdefault(transcript_id, gene_id)

    for transcript_id, gene_id in transcript_to_gene.items():
        print(f"{transcript_id}\t{gene_id}")

    sys.stderr.write(
        f"[{gtf_path}] {len(transcript_to_gene)} transcript_id -> gene_id "
        f"pairs (from {n_attr_rows} attribute rows carrying both fields)\n"
    )
    if not transcript_to_gene:
        sys.stderr.write(
            "WARNING: no transcript_id/gene_id pairs found -- check that "
            "this is a standard GTF with both attributes present on "
            "transcript/exon/CDS rows.\n"
        )


if __name__ == "__main__":
    main()
