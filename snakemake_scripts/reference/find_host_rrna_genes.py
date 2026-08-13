#!/usr/bin/env python3
"""Find host (Drosophila) rRNA genes in a raw FlyBase GTF and print their
identifiers, one per line.

find_rrna_genes.py identifies rRNA by feature type == "rRNA" or a
gbkey/gene_biotype attribute -- that works for NCBI bacterial GTFs but not
FlyBase GTFs, which don't reliably carry either. This instead matches on
FlyBase's own rRNA gene *symbol* naming convention: nuclear rRNA genes are
named "<N>SrRNA[-Psi]:CR#####" (2S, 5.8S, 18S, 28S -- e.g. FBgn0267496
"2SrRNA:CR45836"). Mitochondrial rRNA ("mt:lrRNA", "mt:srRNA") is excluded
by default since Wolbachia titer should be computed from nuclear, not
mitochondrial, rRNA (Wolbachia is known to alter host mtDNA copy number,
which would confound the titer).

Runs directly on the RAW FlyBase GTF (gene_id "FBgn..."; gene_symbol "...";
...) -- no need to run it through swap_gene_id_to_symbol.py first. By
default prints the gene_symbol, cleaned of the "Dsim\\"/"Dwil\\" species
prefix the same way swap_gene_id_to_symbol.py does (so it matches
adata.var_names when your transcriptome's t2g "gene" column is symbol-keyed).
Pass --id-attr gene_id to print raw FBgn IDs instead, if that's what your
var_names use.

Usage:
    find_host_rrna_genes.py dmel-all-r6.68.gtf > host_rrna_genes.txt
    find_host_rrna_genes.py dsim-all-r2.02.gtf --id-attr gene_id > host_rrna_genes.txt
    find_host_rrna_genes.py dmel-all-r6.68.gtf --include-mito > host_rrna_genes.txt
"""
import argparse
import re
import sys

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')
SPECIES_PREFIX_RE = re.compile(r'^D[a-z]+\\')
NUCLEAR_RRNA_RE = re.compile(r'^\d+(\.\d+)?SrRNA')
MITO_RRNA_RE = re.compile(r'^mt:.*rRNA', re.IGNORECASE)


def parse_attrs(field: str) -> dict:
    return dict(ATTR_RE.findall(field))


def clean_symbol(symbol: str) -> str:
    """Strip a FlyBase species prefix like 'Dsim\\' / 'Dwil\\' off a gene
    symbol, matching swap_gene_id_to_symbol.py's clean_symbol()."""
    return SPECIES_PREFIX_RE.sub('', symbol)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf")
    ap.add_argument("--id-attr", default="gene_symbol",
                     choices=["gene_symbol", "gene_id"],
                     help="attribute to print (default: gene_symbol, cleaned "
                          "of species prefix -- matches adata.var_names when "
                          "your transcriptome is symbol-keyed). Use gene_id "
                          "for raw FBgn IDs.")
    ap.add_argument("--include-mito", action="store_true",
                     help="also include mitochondrial rRNA (mt:lrRNA / "
                          "mt:srRNA); excluded by default since Wolbachia "
                          "can alter host mtDNA copy number, which would "
                          "confound the titer calculation.")
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

            attrs = parse_attrs(fields[8])
            symbol = attrs.get("gene_symbol")
            if not symbol:
                continue
            symbol = clean_symbol(symbol)

            is_nuclear_rrna = bool(NUCLEAR_RRNA_RE.match(symbol))
            is_mito_rrna = bool(MITO_RRNA_RE.match(symbol))
            if not (is_nuclear_rrna or (args.include_mito and is_mito_rrna)):
                continue

            if args.id_attr == "gene_id":
                gid = attrs.get("gene_id")
            else:
                gid = symbol

            if gid and gid not in seen:
                seen.add(gid)
                ids.append(gid)

    for gid in ids:
        print(gid)

    print(f"Found {len(ids)} host rRNA gene(s) in {args.gtf}", file=sys.stderr)


if __name__ == "__main__":
    main()
