#!/usr/bin/env python
desc="""Split contigs on likely mis-assemblies
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Warsaw, 6/05/2017
"""

import math, os, sys, subprocess, subprocess
from datetime import datetime
from FastaIndex import FastaIndex
from fastq2sspace import _get_bwamem_proc, _get_snap_proc
from collections import Counter

def sam_proc(bam="600.bam", log=sys.stderr):

    args = ["samtools", "view", bam]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=log)
    return proc

def reads2window(faidx, iSize, wSize, mapQ):

    pref, pwi, c = "", 0, Counter()
    for l in sam_proc().stdout:
        if l.startswith('@'):
            continue
        flag, ref, pos, mapq, cigar, mref, mpos, isize = l.split('\t')[1:9] # q
        flag, pos, mapq, mpos, isize = list(map(int, (flag, pos, mapq, mpos, isize)))
        # skip low quality algs, short contigs and algs from contig ends 
        if mapq<mapQ or faidx.id2stats[ref][0]<3*iSize or pos<1.5*iSize or pos+1.5*iSize>faidx.id2stats[ref][0]:
            continue
        # get window id
        wi = pos / wSize    
        if pref!=ref or wi!=pwi:
            if c: yield pref, pwi, c
            pref, pwi, c = ref, wi, Counter()
        # count mates
        c[mref] += 1
    if pref!=ref or wi!=pwi:
        if c: yield pref, pwi, c
    
def fasta2qc(fasta, outfn="", fq1="", fq2="", iSize=600, wSize=100, minFrac=0.66, mapQ=10, minSize=200, verbose=1): 
    """Split contigs on likely mis-assemblies"""
    if not outfn:
        outfn = fasta+".split.fa"
        
    faidx = FastaIndex(fasta)

    contig2split = {}
    for ref, wi, c in reads2window(faidx, iSize, wSize, mapQ):
        mref, mcount = c.most_common(1)[0]
        if mref!="=" and 1.*mcount/sum(c.values()) > minFrac:
            if verbose: print(ref, wi, mref, c)
            if ref not in contig2split:
                contig2split[ref] = []
            contig2split[ref].append(wi)
            #break
    if verbose: print(contig2split)

    out = open(outfn, "w")
    for c in faidx:
        if c not in contig2split:
            out.write(faidx[c])
            continue
        i = pe = 0
        for i, wi in enumerate(contig2split[c], 1):
            e = wi*wSize
            if e-pe>minSize:
                out.write(faidx.get_fasta(contig=c, start=pe+1, stop=e, name="%s.%s"%(c,i)))
            pe = e+wSize
        e = faidx.id2stats[c][0]
        if e-pe>minSize:
            out.write(faidx.get_fasta(contig=c, start=pe+1, stop=e, name="%s.%s"%(c,i+1)))
    out.close()
        
def main():
    import argparse
    usage   = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("-f", "--fasta", required=True, help="genome fasta")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))
    
    fasta2qc(o.fasta)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!   \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
 