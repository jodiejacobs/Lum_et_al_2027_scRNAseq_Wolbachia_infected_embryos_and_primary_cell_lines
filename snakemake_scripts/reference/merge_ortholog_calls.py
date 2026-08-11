#!/usr/bin/env python3
"""
Merge the symbol-projection ortholog table (find_dmel_dsim_orthologs.py)
with the sequence-based reciprocal-best-BLAST-hit table
(run_reciprocal_blast_orthologs.sh / reciprocal_best_hits.py) into one
final dsim_gene_id -> dmel_gene_id call per gene.

Most D. simulans genes in dsim-all-r2.02.gtf carry a "GD#####" systematic
annotation ID rather than a Dmel-projected symbol (GD is Dsim's equivalent
of Dmel's "CG" numbering -- FlyBase only gives a Dsim gene its Dmel
homolog's symbol when it has confidently identified that homolog). So the
symbol table alone leaves most genes unmatched; BLAST RBH is what actually
resolves the bulk of them. This script combines both, preferring the
symbol call where FlyBase already made one (it's curated, not just a
sequence-similarity guess) and falling back to the BLAST RBH call
everywhere else.

Precedence per dsim gene, using find_dmel_dsim_orthologs.py's match_status:
  - matched / ambiguous_dsim_paralogs
        -> keep the symbol call as-is (method = "symbol"). Both statuses
           carry a single confident dmel_gene_id; "ambiguous_dsim_paralogs"
           just means a sibling dsim gene also projects to it.
  - unmatched / ambiguous_dmel_symbol
        -> use the BLAST RBH call for this gene if one exists
           (method = "blast_rbh"). For ambiguous_dmel_symbol specifically,
           also flag whether the RBH pick matches one of the symbol-based
           candidate dmel genes (blast_confirms_symbol_candidate) or
           disagrees with all of them (blast_conflicts_symbol_candidates)
           -- worth a manual look either way.
  - no call from either method -> method = "none", dmel_gene_id left blank

Usage:
    merge_ortholog_calls.py --symbol-table dsim_to_dmel_orthologs.tsv \\
        --rbh-table dsim_dmel_rbh_orthologs.tsv \\
        --dmel-gtf dmel-all-r6.68.gtf \\
        -o dsim_to_dmel_orthologs_final.tsv
"""
import argparse
import csv
import re
import sys

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')

SYMBOL_CONFIDENT_STATUSES = {'matched', 'ambiguous_dsim_paralogs'}
SYMBOL_UNRESOLVED_STATUSES = {'unmatched', 'ambiguous_dmel_symbol'}


def parse_attrs(field: str) -> dict:
    return dict(re.findall(ATTR_RE, field))


def load_dmel_symbols(gtf_path: str) -> dict:
    """Return {gene_id: gene_symbol} from a FlyBase GTF (for labeling RBH calls)."""
    symbol_by_id = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            attrs = parse_attrs(fields[8])
            gene_id, symbol = attrs.get('gene_id'), attrs.get('gene_symbol')
            if gene_id and symbol and gene_id not in symbol_by_id:
                symbol_by_id[gene_id] = symbol
    return symbol_by_id


def load_symbol_table(path: str) -> dict:
    """Return {dsim_gene_id: row_dict} from find_dmel_dsim_orthologs.py's output."""
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows[row['dsim_gene_id']] = row
    return rows


def load_rbh_table(path: str) -> dict:
    """Return {dsim_gene_id: row_dict} from reciprocal_best_hits.py's output."""
    rows = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            rows[row['dsim_gene_id']] = row
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--symbol-table', required=True, help='output of find_dmel_dsim_orthologs.py')
    ap.add_argument('--rbh-table', required=True, help='output of reciprocal_best_hits.py')
    ap.add_argument('--dmel-gtf', required=True, help='FlyBase D. melanogaster GTF, to label RBH dmel symbols')
    ap.add_argument('-o', '--output', required=True)
    args = ap.parse_args()

    print(f"Loading symbol table from {args.symbol_table}", file=sys.stderr)
    symbol_rows = load_symbol_table(args.symbol_table)
    print(f"  {len(symbol_rows)} dsim genes", file=sys.stderr)

    print(f"Loading RBH table from {args.rbh_table}", file=sys.stderr)
    rbh_rows = load_rbh_table(args.rbh_table)
    print(f"  {len(rbh_rows)} RBH pairs", file=sys.stderr)

    print(f"Loading dmel gene symbols from {args.dmel_gtf}", file=sys.stderr)
    dmel_symbols = load_dmel_symbols(args.dmel_gtf)

    out_rows = []
    counts = {}

    for dsim_id, srow in symbol_rows.items():
        status = srow['match_status']
        rbh = rbh_rows.get(dsim_id)

        if status in SYMBOL_CONFIDENT_STATUSES:
            method = 'symbol'
            dmel_gene_id = srow['dmel_gene_id']
            dmel_symbol = srow['dmel_symbol']
            flag = status
            pident = evalue = bitscore = ''
        elif rbh:
            dmel_gene_id = rbh['dmel_gene_id']
            dmel_symbol = dmel_symbols.get(dmel_gene_id, '')
            pident, evalue, bitscore = rbh['pident'], rbh['evalue'], rbh['bitscore']
            method = 'blast_rbh'
            if status == 'ambiguous_dmel_symbol':
                candidates = set(srow['dmel_gene_id'].split(','))
                flag = ('blast_confirms_symbol_candidate' if dmel_gene_id in candidates
                        else 'blast_conflicts_symbol_candidates')
            else:
                flag = 'blast_only'
        else:
            method = 'none'
            dmel_gene_id = dmel_symbol = ''
            pident = evalue = bitscore = ''
            flag = status  # unmatched / ambiguous_dmel_symbol, and no RBH hit either

        counts[(method, flag)] = counts.get((method, flag), 0) + 1
        out_rows.append((dsim_id, srow['dsim_symbol'], dmel_gene_id, dmel_symbol,
                          method, flag, pident, evalue, bitscore))

    out_rows.sort(key=lambda r: r[0])

    with open(args.output, 'w') as out:
        out.write('dsim_gene_id\tdsim_symbol\tdmel_gene_id\tdmel_symbol\tmethod\tflag\tpident\tevalue\tbitscore\n')
        for row in out_rows:
            out.write('\t'.join(row) + '\n')

    total = len(out_rows)
    n_called = sum(1 for r in out_rows if r[2])
    print(f"\n[merge_ortholog_calls] wrote {args.output}", file=sys.stderr)
    print(f"  {total} dsim genes total, {n_called} with a dmel ortholog call ({n_called/total:.1%})", file=sys.stderr)
    for (method, flag), n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"  {method:10s} {flag:35s} {n}", file=sys.stderr)


if __name__ == '__main__':
    main()
