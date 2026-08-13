#!/usr/bin/env python3
"""
Build a {protein_id, transcript_id} -> gene_id map directly from a GTF
file's attribute column, instead of parsing species-specific conventions
out of a translation FASTA header (get_id_map.py's approach, which only
works for FlyBase-style `parent=FBgn...` headers).

Generic `key "value";` attribute parsing means this works on any
well-formed GTF, in particular:
  - FlyBase GTFs, e.g. gene_id "FBgn0031081"; transcript_id "FBtr0070000";
    ... protein_id "FBpp0070000"; (on CDS rows)
  - NCBI RefSeq/Gnomon GTFs, e.g. gene_id "LOC120284240";
    transcript_id "XM_039291137.1"; ... protein_id "XP_039147071.1";
    (on CDS rows)
so it drops in for D. simulans annotations sourced from NCBI (no FBgn IDs,
no FlyBase header conventions) just as well as FlyBase ones.

Maps BOTH protein_id and transcript_id to gene_id, keyed off whichever ID
is present per row. This matters because `gffread -y` protein extraction
writes the FASTA header from `protein_id` when a GTF's CDS rows carry one
(true for both FlyBase and NCBI/Gnomon GTFs) and only falls back to
transcript_id otherwise -- if this script mapped transcript_id alone, the
protein_id-based headers gffread actually writes (FBpp### / XP_###...)
would never be found in the map, silently falling through filter_rbh.py's
`query_map.get(id, id)` fallback and leaking raw protein accessions into
the final ortholog table's gene_id columns instead of FBgn/LOC gene IDs.
Mapping both covers whichever ID gffread ends up using.

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
    id_to_gene = {}
    n_protein_rows = 0
    n_transcript_rows = 0

    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                continue
            attrs = dict(ATTR_RE.findall(fields[8]))
            gene_id = attrs.get("gene_id")
            if not gene_id:
                continue
            transcript_id = attrs.get("transcript_id")
            protein_id = attrs.get("protein_id")
            # Gene-level feature rows (e.g. FlyBase/NCBI "gene" rows) carry
            # gene_id but an empty transcript_id ("") -- skip those, we only
            # want rows that tie a transcript/protein to its gene.
            if transcript_id:
                n_transcript_rows += 1
                id_to_gene.setdefault(transcript_id, gene_id)
            if protein_id:
                n_protein_rows += 1
                id_to_gene.setdefault(protein_id, gene_id)

    for seq_id, gene_id in id_to_gene.items():
        print(f"{seq_id}\t{gene_id}")

    sys.stderr.write(
        f"[{gtf_path}] {len(id_to_gene)} total id -> gene_id pairs "
        f"({n_transcript_rows} transcript_id rows, {n_protein_rows} "
        f"protein_id rows seen)\n"
    )
    if not id_to_gene:
        sys.stderr.write(
            "WARNING: no transcript_id/protein_id/gene_id pairs found -- "
            "check that this is a standard GTF with these attributes "
            "present on transcript/exon/CDS rows.\n"
        )
    if n_protein_rows == 0:
        sys.stderr.write(
            "NOTE: no protein_id attribute found in this GTF -- gffread -y "
            "will fall back to transcript_id for FASTA headers here, which "
            "this map already covers.\n"
        )


if __name__ == "__main__":
    main()
