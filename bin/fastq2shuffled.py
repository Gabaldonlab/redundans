#!/usr/bin/env python

import os, sys, subprocess, gzip
if sys.version_info < (3,):
    from itertools import izip as zip

def fqparser(fn, stripNames=0, pair="", i=0):
    """Single process implementation of rawtrimmer.
    Open zcat subprocess and read from stdin."""
    if fn.endswith('.gz'):
        zcat = subprocess.Popen(['zcat', fn], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        handle = zcat.stdout
    else:
        handle = open(fn)
        
    fqlist = []
    for l in handle:
        if not l[:-1]:
            continue
        fqlist.append(l[:-1])
        if len(fqlist) == 4:
            i += 1
            name, seq, sep, quals = fqlist
            if stripNames: 
                yield "@%s%s\n%s\n%s\n%s\n"%(i, pair, seq, sep, quals)
            else:
                yield "%s%s\n%s\n%s\n%s\n"%(name, pair, seq, sep, quals)
            fqlist = []

def fastq2shuffled(fnames, out=sys.stdout, stripNames=1, limit=0):
    """Process FastQ files"""
    i = 0
    for fn1, fn2 in zip(fnames[::2], fnames[1::2]):
        for i, (fq1, fq2) in enumerate(zip(fqparser(fn1, stripNames, '/1', i), fqparser(fn2, stripNames, '/2', i)), i):
            if limit and i>limit:
                break
            out.write(fq1+fq2)
        
if __name__=="__main__":
    fastq2shuffled(sys.argv[1:])