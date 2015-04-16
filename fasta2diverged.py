#!/usr/bin/env python
desc="""Return sequences diverged by given percentage compares to input.
LOH sizes are drawn from log-normal distribution. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 4/03/2014
"""

import os, random, sys
from datetime import datetime
from Bio import SeqIO, Seq
import numpy as np

aminos = 'ACDEFGHIKLMNPQRSTVWY'
nucleotides = 'ACGT'

def get_heterozygous(positions, divergence):
    """Return fraction of chromosome being heterozygous"""
    #count hetero SNPs within 100bp
    hetero = []
    pp = 0.0
    for p in sorted(positions):
        if   p - pp <=100:
            if not hetero or len(hetero[-1])==2:
                hetero.append([pp,])
            pp = p
        elif hetero and len(hetero[-1])==1:
            hetero[-1].append(pp)
        pp = p
    if hetero and len(hetero[-1])==1:
        hetero[-1].append(pp)
    return sum(e-s for s, e in hetero)

def seq2diverged(seq, divergence, loh, lohSizes, aminos, verbose):
    """Return diverged sequence"""
    seqlist = list(seq)
    #number of position to change
    k = int(round(divergence*len(seqlist)))
    positions = random.sample(xrange(len(seqlist)), k)
    #inclue LOHs
    if loh:
        lohs = []
        while get_heterozygous(positions, divergence)/len(seq) > 1-loh:
            #get LOH random start
            s = random.randint(0, len(seq)-1)
            #and random LOH length from negative binomial distribution
            lsize = 0
            #LOH has to be at least 2x larger than average distance between SNPs
            while lsize <= 2.0 / divergence:
                lsize = int(round(random.sample(lohSizes, 1)[0]))
            e = s + lsize
            lohs.append(lsize)
            #filter snp posiitons
            positions = filter(lambda p: p<s or p>e, positions)
        if verbose:
            hetero = 100.0 * get_heterozygous(positions, divergence)/len(seq)
            sys.stderr.write(" %.1f%s heterozygous with %s LOHs (%s - %sbp)\n" % \
                             (hetero, '%', len(lohs), min(lohs), max(lohs)))
    #change positions
    for i in positions:
        nb = seqlist[i]
        while nb == seqlist[i]:
            nb = random.choice(aminos)
        seqlist[i] = nb
    return "".join(seqlist)

def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0')   
    parser.add_argument("-i", "--input",     default=sys.stdin, type=file, 
                        help="fasta file(s)   [stdin]")
    parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream   [stdout]")
    parser.add_argument("-d", "--divergence",    default=0.01, type=float, 
                        help="divergence      [%(default)s]")
    parser.add_argument("--dna",             default=False, action='store_true',
                        help="DNA alphabet    [amino acids]")
    parser.add_argument("--loh",             default=0.00, type=float, 
                        help="level of LOH    [%(default)s]")
    parser.add_argument("--learn",    default=False, type=file, 
                        help="BED file to learn LOH sizes [%(default)s]")
    parser.add_argument("--power_a",         default=0.3, type=float, 
                        help="LOH power distribution a [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    
    #get LOH sizes distibution
    lohSizes = np.random.lognormal(0, 1.2, 1e5) * 1e3
    if o.verbose:
        sys.stderr.write("Processing chromosomes...\n")
    for i, r in enumerate(SeqIO.parse(o.input, 'fasta'), 1):
        if o.verbose:
            sys.stderr.write('%s %s %sbp %s\n'%(i, r.id, len(r), " "*20))
        if o.dna:
            alphabet = nucleotides
        else:
            alphabet = aminos
        seq = seq2diverged(r.seq, o.divergence, o.loh, lohSizes, alphabet, o.verbose)
        r.seq = Seq.Seq(seq)
        o.output.write(r.format('fasta'))

if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt  = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n" % dt)
