#!/usr/bin/env python
"""Rewrite PIPseq R1 into a fixed [barcode][UMI] layout for kallisto|bustools.

PIPseq R1 (confirmed on the T20 runs):
  [0-4 bp stagger][cb1 8][ATG][cb2 6][GAG][cb3 6][TCGAG][cb4 8][UMI 12][polyT]
Output R1:
  [cb1 cb2 cb3 cb4 = 28 bp][UMI 12]   ->  kb count -x 0,0,28:0,28,40:1,0,0

Read pairs whose R1 has no exact linker match are dropped; R2 is passed
through unchanged. Writes a JSON log with kept/dropped counts and the
stagger distribution.
"""
import argparse
import json
import re
import shutil
import subprocess
from collections import Counter

PAT = re.compile(rb"^[ACGTN]{0,4}?([ACGTN]{8})ATG([ACGTN]{6})GAG([ACGTN]{6})"
                 rb"TCGAG([ACGTN]{8})([ACGTN]{12})")

ap = argparse.ArgumentParser()
ap.add_argument("--r1", required=True)
ap.add_argument("--r2", required=True)
ap.add_argument("--out-r1", required=True)
ap.add_argument("--out-r2", required=True)
ap.add_argument("--log", required=True)
ap.add_argument("--threads", type=int, default=4)
a = ap.parse_args()

gz = shutil.which("pigz") or "gzip"
pthr = ["-p", str(a.threads)] if gz.endswith("pigz") else []


def reader(path):
    return subprocess.Popen([gz, "-dc", path], stdout=subprocess.PIPE,
                            bufsize=1 << 20)


def writer(path):
    return subprocess.Popen([gz, "-c", *pthr], stdin=subprocess.PIPE,
                            stdout=open(path, "wb"), bufsize=1 << 20)


r1p, r2p = reader(a.r1), reader(a.r2)
w1p, w2p = writer(a.out_r1), writer(a.out_r2)
f1, f2, o1, o2 = r1p.stdout, r2p.stdout, w1p.stdin, w2p.stdin

total = kept = 0
stagger = Counter()
while True:
    h1 = f1.readline()
    if not h1:
        break
    s1, p1, q1 = f1.readline(), f1.readline(), f1.readline()
    rec2 = b"".join((f2.readline(), f2.readline(), f2.readline(), f2.readline()))
    total += 1
    m = PAT.match(s1)
    if not m:
        continue
    kept += 1
    stagger[m.start(1)] += 1
    seq = b"".join(m.group(i) for i in range(1, 6))
    qual = b"".join(q1[m.start(i):m.end(i)] for i in range(1, 6))
    o1.write(b"%s%s\n+\n%s\n" % (h1, seq, qual))
    o2.write(rec2)

if f2.readline():
    raise SystemExit("R2 has more reads than R1; files are not paired")
for p in (w1p, w2p):
    p.stdin.close()
    p.wait()
for p in (r1p, r2p, w1p, w2p):
    if p.wait() != 0:
        raise SystemExit(f"{p.args[0]} exited with {p.returncode}")

with open(a.log, "w") as fh:
    json.dump({"reads_total": total, "reads_kept": kept,
               "frac_kept": kept / total if total else 0,
               "stagger": {str(k): v for k, v in sorted(stagger.items())}},
              fh, indent=2)
