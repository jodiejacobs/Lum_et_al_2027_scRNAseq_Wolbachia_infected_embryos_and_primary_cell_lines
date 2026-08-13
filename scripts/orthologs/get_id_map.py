#!/usr/bin/env python3
"""
NOTE: only usable when you actually have a FlyBase translation FASTA for
both species (headers with a `parent=FBgn...` field). run_rbh_orthologs.sbatch
now builds proteins from genome FASTA + GTF via gffread instead (D. simulans
moved to an NCBI RefSeq/Gnomon annotation with no such FASTA or header
convention), and uses gtf_id_map.py for id mapping in that case. Keep this
script around for the pure-FlyBase-translation-FASTA case only.

Build a protein_id -> gene_id map from a FlyBase-style translation FASTA.

FlyBase protein FASTA headers look like:
    >FBpp0070001 type=protein; loc=2L:7529..8116,+; ID=FBpp0070001;
     name=CG12345-PA; parent=FBgn0031208,FBtr0070001; ...

We pull the FBgn ID out of the `parent=` field. If no FBgn is found
(e.g. you swapped in an NCBI protein FASTA instead), the protein ID
itself is used as a fallback gene ID and a warning is printed — in
that case do gene-level mapping from the GTF instead and re-run.

Usage:
    python get_id_map.py proteins.fasta > id_map.tsv
"""
import re
import sys

FBGN_RE = re.compile(r"parent=([^,;\s]*FBgn\d+)")


def main():
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} proteins.fasta > id_map.tsv")

    fasta_path = sys.argv[1]
    n_total = 0
    n_fallback = 0

    with open(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            n_total += 1
            header = line[1:].strip()
            protein_id = header.split()[0]

            m = FBGN_RE.search(header)
            if m:
                gene_id = m.group(1)
            else:
                gene_id = protein_id
                n_fallback += 1

            print(f"{protein_id}\t{gene_id}")

    sys.stderr.write(
        f"[{fasta_path}] {n_total} proteins, "
        f"{n_fallback} without a parsed FBgn (fell back to protein ID)\n"
    )
    if n_total and n_fallback == n_total:
        sys.stderr.write(
            "WARNING: no FBgn IDs found at all — this doesn't look like a "
            "FlyBase translation FASTA. Gene-level RBH results will actually "
            "be at the protein/transcript level. Map protein->gene from the "
            "matching GTF/GFF3 instead if you need true gene IDs.\n"
        )


if __name__ == "__main__":
    main()
