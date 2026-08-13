#!/usr/bin/env python3
"""Find rRNA genes in a genome annotation (GTF or GFF3) and print their gene
IDs, one per line.

Replaces hardcoding a fixed rRNA gene list per species/Wolbachia strain:
run this against whatever GTF/GFF3 a genome build actually uses and it
finds the rRNA-annotated genes in that file.

Works on:
  - FlyBase-style GTF (host species): feature type "rRNA" in column 3,
    gene id in the `gene_id` attribute (or `gene` after
    swap_gene_id_to_symbol.py has been run).
  - NCBI/Prokka-style bacterial GTF (Wolbachia): feature type "rRNA" in
    column 3, or `gbkey`/`gene_biotype` attribute == "rRNA", gene id in
    `locus_tag`, `gene`, or `gene_id`.
  - NCBI GFF3 (e.g. genomic.gff / genomic.fixed.gff): same feature-type/
    gbkey rules, but attributes are unquoted `key=value;key=value` instead
    of GTF's quoted `key "value";` -- both are parsed automatically.
    Some derived/converted GTFs (e.g. AGAT output) drop gbkey/gene_biotype
    and locus_tag entirely for non-coding features, so if a strain's own
    .gtf comes up empty, try running this against its source GFF3 instead.

Usage:
    find_rrna_genes.py annotation.gtf > rrna_genes.txt
    find_rrna_genes.py genomic.gff --id-attr locus_tag > rrna_genes.txt
"""
import argparse
import re
import sys

GTF_ATTR_RE  = re.compile(r'(\S+)\s+"([^"]*)"')
GFF3_ATTR_RE = re.compile(r'([^=;\s]+)=([^;]*)')


def parse_attrs(field: str) -> dict:
    """Parse a GTF attribute string ('key "value"; key "value";') or a
    GFF3 one ('key=value;key=value') -- tries GTF-style first since GFF3
    values could in principle contain a stray quoted substring, then falls
    back to GFF3-style if that finds nothing."""
    pairs = GTF_ATTR_RE.findall(field)
    if pairs:
        return dict(pairs)
    return dict(GFF3_ATTR_RE.findall(field))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf")
    ap.add_argument("--id-attr", default=None,
                     help="attribute to use as the gene id "
                          "(default: try locus_tag, then gene, then gene_id)")
    args = ap.parse_args()

    ids = []
    seen = set()

    with open(args.gtf) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature = fields[2]
            attrs = parse_attrs(fields[8])

            is_rrna = (
                feature == "rRNA"
                or attrs.get("gbkey", "").lower() == "rrna"
                or attrs.get("gene_biotype", "").lower() == "rrna"
            )
            if not is_rrna:
                continue

            if args.id_attr:
                gid = attrs.get(args.id_attr)
            else:
                gid = attrs.get("locus_tag") or attrs.get("gene") or attrs.get("gene_id")

            if gid and gid not in seen:
                seen.add(gid)
                ids.append(gid)

    for gid in ids:
        print(gid)

    print(f"Found {len(ids)} rRNA gene(s) in {args.gtf}", file=sys.stderr)


if __name__ == "__main__":
    main()
