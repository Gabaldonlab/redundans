#!/usr/bin/env python
# Skip self matches from MAF file

import sys

out = sys.stdout.write
err = sys.stderr.write

maf = []

for i, l in enumerate(sys.stdin, 1):
    #err(" %s  \n"%i)
    if l.startswith('#'):
        out(l)
        continue
    maf.append(l)
    if len(maf)==4:
        if maf[1].split()[:6]!=maf[2].split()[:6]:
            out("".join(maf))
        maf = []