#!/usr/bin/env python3
"""
Rewrite a FlyBase GTF so gene_id holds the gene_symbol value instead of the
FBgn ID, stripping the species-prefix FlyBase puts on non-melanogaster
symbols (e.g. "Dsim\\CG1704" -> "CG1704", "Dwil\\Adh" -> "Adh") so orthologs
share the same base symbol across dmel/dsim/dwil. kb ref keys the t2g "gene"
column off gene_id, so running kb ref against the output of this script
produces a transcriptome/index/t2g keyed by symbol instead of FBgn.

Every GTF line already carries both gene_id and gene_symbol as redundant
attributes, so this is a straight line-by-line rewrite -- no need to build a
gene -> symbol table across lines first.

Gene symbols are not guaranteed unique the way FBgn IDs are (paralogs,
CG-number fallbacks, symbol reuse across releases). This prints a warning
listing any (stripped) symbol that maps back to more than one distinct
original gene_id -- those will collide in adata.var_names and need
`adata.var_names_make_unique()` downstream.

Usage: swap_gene_id_to_symbol.py <in.gtf> <out.gtf>
"""
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')
SPECIES_PREFIX_RE = re.compile(r'^D[a-z]+\\')


def parse_attributes(attr_str):
    """Return an ordered list of (key, value) pairs from a GTF attribute field."""
    return re.findall(ATTR_RE, attr_str)


def clean_symbol(symbol):
    """Strip a FlyBase species prefix like 'Dsim\\' / 'Dwil\\' off a gene symbol."""
    return SPECIES_PREFIX_RE.sub('', symbol)


def main(in_path, out_path):
    symbol_to_ids = defaultdict(set)
    n_lines = 0
    n_swapped = 0
    n_no_symbol = 0

    with open(in_path) as fin, open(out_path, 'w') as fout:
        for line in fin:
            if not line.strip() or line.startswith('#'):
                fout.write(line)
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                fout.write(line)
                continue

            n_lines += 1
            pairs = parse_attributes(fields[8])
            attr_dict = dict(pairs)

            gene_id = attr_dict.get('gene_id')
            gene_symbol = attr_dict.get('gene_symbol')

            if gene_symbol:
                symbol = clean_symbol(gene_symbol)
                if gene_id:
                    symbol_to_ids[symbol].add(gene_id)
                new_pairs = [
                    (k, symbol if k == 'gene_id' else v)
                    for k, v in pairs
                ]
                n_swapped += 1
            else:
                new_pairs = pairs
                n_no_symbol += 1

            fields[8] = ' '.join(f'{k} "{v}";' for k, v in new_pairs)
            fout.write('\t'.join(fields) + '\n')

    print(f"[swap_gene_id_to_symbol] {in_path} -> {out_path}")
    print(f"  {n_lines} feature lines, {n_swapped} with gene_id swapped to symbol, "
          f"{n_no_symbol} left unchanged (no gene_symbol attribute)")

    collisions = {s: ids for s, ids in symbol_to_ids.items() if len(ids) > 1}
    if collisions:
        print(f"  WARNING: {len(collisions)} symbol(s) map to more than one "
              f"original gene_id -- these will collide in adata.var_names:")
        for symbol, ids in sorted(collisions.items()):
            print(f"    {symbol}: {sorted(ids)}")
    else:
        print("  No symbol collisions found.")


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print("usage: swap_gene_id_to_symbol.py <in.gtf> <out.gtf>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
