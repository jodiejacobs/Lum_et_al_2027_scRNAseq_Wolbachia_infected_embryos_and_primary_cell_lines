#!/usr/bin/env python3
"""
Call reciprocal-best-BLAST-hit (RBH) gene orthologs from two BLASTP
outfmt-6 tables (dsim-query-vs-dmel-db and dmel-query-vs-dsim-db), aggregated
from transcript/protein hits up to the gene level.

Meant as the sequence-based counterpart to find_dmel_dsim_orthologs.py
(which calls orthologs from FlyBase's gene_symbol projection instead) --
run this on the genes that script left unmatched, or on everything to
cross-check the symbol-based calls. See run_reciprocal_blast_orthologs.sh,
which runs gffread + blastp and then calls this script.

BLASTP query/subject IDs are transcript IDs (gffread -y protein FASTA
headers), not gene IDs -- this script maps each hit back to its gene via
the two GTFs, then for each gene keeps its single best-bitscore hit gene
on the other side. A pair (dsim_gene, dmel_gene) is reported as an RBH
ortholog only when each is the other's best hit in both directions.

Usage:
    reciprocal_best_hits.py --forward dsim_vs_dmel.tsv --reverse dmel_vs_dsim.tsv \\
        --dsim-gtf dsim-all-r2.02.gtf --dmel-gtf dmel-all-r6.68.gtf \\
        -o dsim_dmel_rbh_orthologs.tsv [--restrict-ids dsim_unmatched_genes.txt]
"""
import argparse
import re
import sys
from collections import defaultdict

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')

BLAST_COLS = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
              'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']


def parse_attrs(field: str) -> dict:
    return dict(re.findall(ATTR_RE, field))


def load_transcript_to_gene(gtf_path: str) -> dict:
    """Return {transcript_id: gene_id} from a FlyBase GTF."""
    t2g = {}
    with open(gtf_path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if len(fields) != 9:
                continue
            attrs = parse_attrs(fields[8])
            tid, gid = attrs.get('transcript_id'), attrs.get('gene_id')
            if tid and gid:
                t2g[tid] = gid
    return t2g


def best_gene_hits(blast_tsv: str, query_t2g: dict, subject_t2g: dict) -> dict:
    """Collapse a BLASTP outfmt6 table to {query_gene: (best_subject_gene, pident, evalue, bitscore)},
    aggregating over all query/subject transcript pairs and keeping the
    highest-bitscore subject gene per query gene."""
    best = {}  # query_gene -> (bitscore, subject_gene, pident, evalue)
    n_lines = n_unmapped = 0
    with open(blast_tsv) as fh:
        for line in fh:
            if not line.strip():
                continue
            n_lines += 1
            fields = line.rstrip('\n').split('\t')
            row = dict(zip(BLAST_COLS, fields))
            qgene = query_t2g.get(row['qseqid'])
            sgene = subject_t2g.get(row['sseqid'])
            if not qgene or not sgene:
                n_unmapped += 1
                continue
            bitscore = float(row['bitscore'])
            cur = best.get(qgene)
            if cur is None or bitscore > cur[0]:
                best[qgene] = (bitscore, sgene, row['pident'], row['evalue'])
    if n_unmapped:
        print(f"  WARNING: {n_unmapped}/{n_lines} lines in {blast_tsv} had a "
              f"query/subject ID not found in the supplied GTF(s)", file=sys.stderr)
    return best


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--forward', required=True, help='dsim-query-vs-dmel-db BLASTP outfmt6 table')
    ap.add_argument('--reverse', required=True, help='dmel-query-vs-dsim-db BLASTP outfmt6 table')
    ap.add_argument('--dsim-gtf', required=True)
    ap.add_argument('--dmel-gtf', required=True)
    ap.add_argument('-o', '--output', required=True)
    ap.add_argument('--restrict-ids', default=None,
                     help='optional: only report RBH pairs for dsim gene_ids in this file '
                          '(one per line) -- e.g. find_dmel_dsim_orthologs.py --unmatched-out')
    args = ap.parse_args()

    dsim_t2g = load_transcript_to_gene(args.dsim_gtf)
    dmel_t2g = load_transcript_to_gene(args.dmel_gtf)
    print(f"dsim: {len(dsim_t2g)} transcripts, dmel: {len(dmel_t2g)} transcripts", file=sys.stderr)

    print(f"Collapsing {args.forward} (dsim -> dmel) to gene level", file=sys.stderr)
    forward_best = best_gene_hits(args.forward, dsim_t2g, dmel_t2g)
    print(f"  {len(forward_best)} dsim genes with a dmel hit", file=sys.stderr)

    print(f"Collapsing {args.reverse} (dmel -> dsim) to gene level", file=sys.stderr)
    reverse_best = best_gene_hits(args.reverse, dmel_t2g, dsim_t2g)
    print(f"  {len(reverse_best)} dmel genes with a dsim hit", file=sys.stderr)

    restrict = None
    if args.restrict_ids:
        with open(args.restrict_ids) as fh:
            restrict = {line.strip() for line in fh if line.strip()}
        print(f"Restricting output to {len(restrict)} dsim gene_ids from {args.restrict_ids}", file=sys.stderr)

    rows = []
    for dsim_gene, (bitscore, dmel_gene, pident, evalue) in forward_best.items():
        if restrict is not None and dsim_gene not in restrict:
            continue
        reverse_hit = reverse_best.get(dmel_gene)
        if reverse_hit and reverse_hit[1] == dsim_gene:
            rows.append((dsim_gene, dmel_gene, pident, evalue, f'{bitscore:.1f}'))

    rows.sort(key=lambda r: r[0])

    with open(args.output, 'w') as out:
        out.write('dsim_gene_id\tdmel_gene_id\tpident\tevalue\tbitscore\n')
        for row in rows:
            out.write('\t'.join(row) + '\n')

    print(f"\n[reciprocal_best_hits] wrote {args.output}: {len(rows)} RBH ortholog pairs", file=sys.stderr)


if __name__ == '__main__':
    main()
