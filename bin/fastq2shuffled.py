#!/usr/bin/env python3

#from msilib import sequence
import os, sys, subprocess, gzip
from unicodedata import name
#if sys.version_info < (3,):
    

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
            bname, bseq, bsep, bquals = fqlist

            #Convert byte stream to str
            name = bytes(bname).decode('utf-8')
            seq = bytes(bseq).decode('utf-8')
            sep = bytes(bsep).decode('utf-8')
            quals = bytes(bquals).decode('utf-8')
            if stripNames: 
                yield "@%s%s\n%s\n%s\n%s\n"%(i, pair, seq, sep, quals)
            else:
                yield "%s%s\n%s\n%s\n%s\n"%(name, pair, seq, sep, quals)
            fqlist = []

def fastq2shuffled(fnames, out, stripNames=1, limit=0):
    """Process FastQ files"""
    i = 0
    with open(out, 'w') as file:
        for fn1, fn2 in zip(fnames[0::2], fnames[1::2]):
            for i, (fq1, fq2) in enumerate(zip(fqparser(fn1, stripNames, '/1', i), fqparser(fn2, stripNames, '/2', i)), i):
                if limit and i>limit:
                    print("Breaked ", i, limit)
                    break
                file.write(fq1+fq2)
        
if __name__=="__main__":
    #Range of all arguments from 1 (0 is the script) till the last and the last, as it is the outfile path
    fastq2shuffled(sys.argv[1:-1], out=sys.argv[-1])