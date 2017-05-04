#!/usr/bin/env python
# Skip self matches from MAF file

import sys

out, err = sys.stdout.write, sys.stderr.write

maf = []
minLength = 200
if len(sys.argv)>1 and sys.argv[1].isdigit():
    minLength = int(sys.argv[1])

for i, l in enumerate(sys.stdin, 1):
    #err(" %s  \n"%i)
    if l.startswith('#'):
        out(l)
        continue
    maf.append(l)
    if len(maf)==4:
        # report not self matches to larger reference and if query sequences >= minLength
        maf1, maf2 = maf[1].split()[:6], maf[2].split()[:6]
        if maf1!=maf2 and int(maf1[5])>=int(maf2[5]) and int(maf2[5])>=minLength:
            out("".join(maf))
        maf = []