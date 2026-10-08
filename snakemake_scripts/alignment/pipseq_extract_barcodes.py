#!/usr/bin/env python
"""Rewrite PIPseq R1 into a fixed [barcode][UMI] layout for kallisto|bustools.

PIPseq R1 (confirmed on the T20 runs):
  [0-4 bp stagger][cb1 8][ATG][cb2 6][GAG][cb3 6][TCGAG][cb4 8][UMI 12][polyT]
Output R1:
  [cb1 cb2 cb3 cb4 = 28 bp][UMI 12]   ->  kb count -x 0,0,28:0,28,40:1,0,0

Only R1 is rewritten; kb reads the original R2 directly, so every R1 record is
kept to stay in sync with R2. Reads with no exact linker match get a random
28 bp barcode, which never reaches the whitelist and is discarded by
bustools correct. Writes a JSON log with matched counts and the stagger
distribution.
"""
import argparse
import json
import random
import re
import shutil
import subprocess
from collections import Counter

PAT = re.compile(rb"^[ACGTN]{0,4}?([ACGTN]{8})ATG([ACGTN]{6})GAG([ACGTN]{6})"
                 rb"TCGAG([ACGTN]{8})([ACGTN]{12})")

ap = argparse.ArgumentParser()
ap.add_argument("--r1", required=True)
ap.add_argument("--out-r1", required=True)
ap.add_argument("--log", required=True)
ap.add_argument("--threads", type=int, default=2)
a = ap.parse_args()

gz = shutil.which("pigz") or "gzip"
pthr = ["-p", str(a.threads)] if gz.endswith("pigz") else []
rng = random.Random(0)
NOQUAL = b"#" * 40

rp = subprocess.Popen([gz, "-dc", a.r1], stdout=subprocess.PIPE, bufsize=1 << 20)
wp = subprocess.Popen([gz, "-1", "-c", *pthr], stdin=subprocess.PIPE,
                      stdout=open(a.out_r1, "wb"), bufsize=1 << 20)
fin, out = rp.stdout, wp.stdin

total = matched = 0
stagger = Counter()
while True:
    h = fin.readline()
    if not h:
        break
    s, _, q = fin.readline(), fin.readline(), fin.readline()
    total += 1
    m = PAT.match(s)
    if m:
        matched += 1
        stagger[m.start(1)] += 1
        seq = b"".join(m.group(i) for i in range(1, 6))
        qual = b"".join(q[m.start(i):m.end(i)] for i in range(1, 6))
    else:
        seq = bytes(rng.choices(b"ACGT", k=40))
        qual = NOQUAL
    out.write(b"%s%s\n+\n%s\n" % (h, seq, qual))

out.close()
for p in (rp, wp):
    if p.wait() != 0:
        raise SystemExit(f"{p.args[0]} exited with {p.returncode}")

with open(a.log, "w") as fh:
    json.dump({"reads_total": total, "reads_matched": matched,
               "frac_matched": matched / total if total else 0,
               "stagger": {str(k): v for k, v in sorted(stagger.items())}},
              fh, indent=2)
