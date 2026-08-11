#!/usr/bin/env python3
"""
Map each Drosophila simulans gene to its D. melanogaster ortholog using
FlyBase's own comparative annotation, rather than an external ortholog
database.

WHY THIS WORKS WITHOUT DIOPT/BLAST:
The current FlyBase D. simulans annotation (dsim r2.02, FB2017_04) was built
by projecting gene models from the D. melanogaster reference genome via
synteny. Where a Dsim gene has a confident 1:1 D. melanogaster counterpart,
FlyBase names it after that Dmel gene, with a species-prefixed symbol --
e.g. Dsim's ortholog of Dmel's CG1704 is annotated "Dsim\\CG1704", and the
ortholog of Adh is "Dsim\\Adh". Stripping that prefix and matching against
the D. melanogaster GTF's gene_symbol values recovers the ortholog call
FlyBase already made.

This is the same species-prefix-stripping logic swap_gene_id_to_symbol.py
uses when it rewrites gene_id -> symbol for the kallisto|bustools reference
build -- this script instead builds an explicit
dsim_gene_id -> dmel_gene_id ortholog table (FBgn to FBgn), independent of
that reference-build step.

(Note: as of FlyBase FB2026_02, the bulk "dmel_orthologs_in_drosophila_
species" precomputed file that used to make this a one-line download no
longer exists on the FlyBase FTP -- only dmel_paralogs and
dmel_human_orthologs_disease remain. Symbol-projection from the GTFs
themselves, done here, is the current best FlyBase-native substitute for
Dmel<->Dsim specifically.)

LIMITATIONS:
- Genes with no clear 1:1 Dmel counterpart (species-specific genes, or
  genes FlyBase hasn't confidently placed) won't have a gene_symbol that
  strips down to a real Dmel symbol, and are reported as "unmatched".
- A few Dmel symbols are shared by more than one Dsim gene (recent
  duplications / paralogs that both project to the same Dmel gene) --
  flagged as "ambiguous_dsim_paralogs" rather than silently picking one.
- This is a nomenclature/projection-based call, not a sequence-based one.
  For unmatched genes, or to independently verify matches, run
  run_reciprocal_blast_orthologs.sh in this same directory for a
  sequence-based (reciprocal-best-BLAST-hit) fallback.

Usage:
    find_dmel_dsim_orthologs.py --dsim-gtf dsim-all-r2.02.gtf \\
        --dmel-gtf dmel-all-r6.68.gtf \\
        -o dsim_to_dmel_orthologs.tsv \\
        [--unmatched-out dsim_unmatched_genes.txt]
"""
import argparse
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')
SPECIES_PREFIX_RE = re.compile(r'^D[a-z]+\\')


def parse_attrs(field: str) -> dict:
    return dict(re.findall(ATTR_RE, field))


def clean_symbol(symbol: str) -> str:
    """Strip a FlyBase species prefix like 'Dsim\\' / 'Dwil\\' off a gene symbol."""
    return SPECIES_PREFIX_RE.sub('', symbol)


def load_genes(gtf_path: str) -> dict:
    """Return {gene_id: gene_symbol} for a FlyBase GTF, one entry per gene_id.

    Prefers the 'gene' feature row for each gene_id when present (there's
    exactly one per gene in a well-formed FlyBase GTF); falls back to the
    first row carrying that gene_id otherwise.
    """
    symbol_by_id = {}
    from_gene_feature = set()
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            feature = fields[2]
            attrs = parse_attrs(fields[8])
            gene_id = attrs.get('gene_id')
            gene_symbol = attrs.get('gene_symbol')
            if not gene_id or not gene_symbol:
                continue
            if gene_id in from_gene_feature:
                continue
            if feature == 'gene':
                symbol_by_id[gene_id] = gene_symbol
                from_gene_feature.add(gene_id)
            elif gene_id not in symbol_by_id:
                symbol_by_id[gene_id] = gene_symbol
    return symbol_by_id


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--dsim-gtf', required=True, help='FlyBase D. simulans GTF (e.g. dsim-all-r2.02.gtf)')
    ap.add_argument('--dmel-gtf', required=True, help='FlyBase D. melanogaster GTF (e.g. dmel-all-r6.68.gtf)')
    ap.add_argument('-o', '--output', required=True, help='output TSV path')
    ap.add_argument('--unmatched-out', default=None,
                     help='optional: write unmatched dsim gene_ids here, one per line '
                          '(feed this into run_reciprocal_blast_orthologs.sh)')
    args = ap.parse_args()

    print(f"Loading D. simulans genes from {args.dsim_gtf}", file=sys.stderr)
    dsim_genes = load_genes(args.dsim_gtf)
    print(f"  {len(dsim_genes)} dsim genes", file=sys.stderr)

    print(f"Loading D. melanogaster genes from {args.dmel_gtf}", file=sys.stderr)
    dmel_genes = load_genes(args.dmel_gtf)
    print(f"  {len(dmel_genes)} dmel genes", file=sys.stderr)

    # dmel clean symbol -> [gene_id, ...] (list, in case >1 dmel gene shares a symbol)
    dmel_by_symbol = defaultdict(list)
    for gene_id, symbol in dmel_genes.items():
        dmel_by_symbol[clean_symbol(symbol)].append(gene_id)

    # dsim clean symbol -> [dsim_gene_id, ...], to flag >1 dsim gene projecting
    # to the same dmel symbol (paralogs / recent duplications)
    dsim_by_clean_symbol = defaultdict(list)
    for gene_id, symbol in dsim_genes.items():
        dsim_by_clean_symbol[clean_symbol(symbol)].append(gene_id)

    rows = []
    unmatched_ids = []
    counts = defaultdict(int)

    for dsim_id, dsim_symbol in dsim_genes.items():
        clean = clean_symbol(dsim_symbol)
        dmel_hits = dmel_by_symbol.get(clean, [])
        dsim_siblings = dsim_by_clean_symbol[clean]

        if not dmel_hits:
            status = 'unmatched'
            dmel_id = dmel_symbol_out = ''
            unmatched_ids.append(dsim_id)
        elif len(dmel_hits) > 1:
            status = 'ambiguous_dmel_symbol'
            dmel_id = ','.join(sorted(dmel_hits))
            dmel_symbol_out = clean
        elif len(dsim_siblings) > 1:
            status = 'ambiguous_dsim_paralogs'
            dmel_id = dmel_hits[0]
            dmel_symbol_out = dmel_genes[dmel_id]
        else:
            status = 'matched'
            dmel_id = dmel_hits[0]
            dmel_symbol_out = dmel_genes[dmel_id]

        counts[status] += 1
        rows.append((dsim_id, dsim_symbol, clean, dmel_id, dmel_symbol_out, status))

    rows.sort(key=lambda r: r[0])

    with open(args.output, 'w') as out:
        out.write('dsim_gene_id\tdsim_symbol\tdsim_symbol_clean\tdmel_gene_id\tdmel_symbol\tmatch_status\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')

    if args.unmatched_out:
        with open(args.unmatched_out, 'w') as out:
            for gid in sorted(unmatched_ids):
                out.write(gid + '\n')

    total = len(dsim_genes)
    print(f"\n[find_dmel_dsim_orthologs] wrote {args.output}", file=sys.stderr)
    print(f"  {total} dsim genes total", file=sys.stderr)
    for status in ('matched', 'ambiguous_dsim_paralogs', 'ambiguous_dmel_symbol', 'unmatched'):
        n = counts.get(status, 0)
        print(f"  {status}: {n} ({n/total:.1%})", file=sys.stderr)
    if args.unmatched_out:
        print(f"  unmatched gene_ids written to {args.unmatched_out}", file=sys.stderr)
        print("  -> for these, run run_reciprocal_blast_orthologs.sh for a "
              "sequence-based fallback", file=sys.stderr)


if __name__ == '__main__':
    main()
