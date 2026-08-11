#!/usr/bin/env python3
"""
Collapse two one-directional DIAMOND blastp best-hit tables into a
genome-wide, gene-level reciprocal-best-hit (RBH) ortholog table.

Inputs are DIAMOND results run with:
    --max-target-seqs 1 --outfmt 6 qseqid sseqid pident length evalue bitscore

so each is already restricted to one hit per query protein; this script
still re-derives the single best hit per query defensively, then maps
protein IDs -> gene IDs, then keeps only gene pairs that are each
other's best hit in both directions.

Usage:
    python filter_rbh.py \\
        dsim_vs_dmel.tsv dmel_vs_dsim.tsv \\
        dsim_id_map.tsv dmel_id_map.tsv \\
        > dmel_dsim_orthologs_rbh.tsv
"""
import csv
import sys


def load_id_map(path):
    mapping = {}
    with open(path) as fh:
        for line in fh:
            protein_id, gene_id = line.rstrip("\n").split("\t")
            mapping[protein_id] = gene_id
    return mapping


def load_best_hits_by_gene(path, query_map, subject_map):
    """Best hit per query gene, collapsed across isoforms, keyed by gene ID."""
    best = {}  # query_gene -> (subject_gene, bitscore, pident, evalue)
    with open(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        for row in reader:
            qseqid, sseqid, pident, length, evalue, bitscore = row
            bitscore = float(bitscore)
            q_gene = query_map.get(qseqid, qseqid)
            s_gene = subject_map.get(sseqid, sseqid)

            if q_gene not in best or bitscore > best[q_gene][1]:
                best[q_gene] = (s_gene, bitscore, pident, evalue)
    return best


def main():
    if len(sys.argv) != 5:
        sys.exit(
            f"Usage: {sys.argv[0]} A_vs_B.tsv B_vs_A.tsv A_id_map.tsv B_id_map.tsv"
        )

    a_vs_b_path, b_vs_a_path, a_map_path, b_map_path = sys.argv[1:5]

    a_map = load_id_map(a_map_path)  # e.g. dsim protein -> dsim gene
    b_map = load_id_map(b_map_path)  # e.g. dmel protein -> dmel gene

    # A queries against B database
    a_vs_b = load_best_hits_by_gene(a_vs_b_path, query_map=a_map, subject_map=b_map)
    # B queries against A database
    b_vs_a = load_best_hits_by_gene(b_vs_a_path, query_map=b_map, subject_map=a_map)

    writer = csv.writer(sys.stdout, delimiter="\t")
    writer.writerow(["geneA", "geneB", "pident", "evalue", "bitscore"])

    n_rbh = 0
    for gene_a, (gene_b, bitscore, pident, evalue) in a_vs_b.items():
        hit_back = b_vs_a.get(gene_b)
        if hit_back and hit_back[0] == gene_a:
            writer.writerow([gene_a, gene_b, pident, evalue, bitscore])
            n_rbh += 1

    sys.stderr.write(
        f"Genes in A: {len(a_vs_b)}, genes in B: {len(b_vs_a)}, "
        f"reciprocal best-hit gene pairs: {n_rbh}\n"
    )


if __name__ == "__main__":
    main()
