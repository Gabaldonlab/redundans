#!/usr/bin/env python
# Produce symmetric tab for self-matches from lastal
# USAGE: f=$ref; time lastal -C2 -P 8 $ref $f | ~/src/redundans/bin/skip_selfmatches.py | last-split | maf-convert tab | ~/src/redundans/bin/last-split2symmetric.py | last-dotplot - $f.best.png

import sys
out = sys.stdout

for l in sys.stdin:
    if l.startswith('#'):
        out.write(l)
        continue
    ldata = l.split('\t')
    t, q = ldata[1:6], ldata[6:11]
    out.write("%s\t%s\t%s\t%s"%(ldata[0], "\t".join(t), "\t".join(q), "\t".join(ldata[11:])))
    out.write("%s\t%s\t%s\t%s"%(ldata[0], "\t".join(q), "\t".join(t), "\t".join(ldata[11:])))
