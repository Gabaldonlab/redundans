#!/usr/bin/env python
desc="""Generate paired-end / mate pairs libraries from long reads.
EXPERIMENTAL VERSION!

TBD:
- fix read orientation
- estimate which IS shall be generated
- need proprietary FastQ parser

More info at: https://github.com/lpryszcz/redundans
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 26/02/2016
"""

import os, sys, gzip
from datetime import datetime

def fastq_parser(fastq):
    """Parse fastq file"""
    pass

def read2mates(out, r, readlen=50, window=25, minLen=400, i=0):
    """Write mates from given read."""
    name = r.name
    r.description = ""
    for i, s in enumerate(xrange(0, len(r)-minLen-readlen, window), 1):
        r.id = "%s.r%s/1"%(name, i)
        # FR
        out[0].write(r[s:s+readlen].format('fastq'))
        # read2 reverse complement
        r2 = r[s+minLen:s+minLen+readlen].reverse_complement()
        r2.id = "%s.r%s/2"%(name, i)
        r2.description = ""
        out[1].write(r2.format('fastq'))
    return i

def fastq2mates(fastq, outbase, readlen, verbose, minLen=400):
    """Generate fake mate pair reads from long reads."""
    # get output fnames
    if not outbase:
        outbase = fastq
    fn1, fn2 = outbase+"_1.fq.gz", outbase+"_2.fq.gz"
    # skip if output file exists
    for fn in (fn1, fn2):
        if os.path.isfile(fn) or os.path.isfile(fn2):
            sys.stderr.write(" %s exists!\n"%fn)
            return
    # open output files
    out = (gzip.open(fn1, "w"), gzip.open(fn2, "w"))
    # get window size - this will be estimated in the future
    window = readlen / 2
    #if verbose:
    sys.stderr.write(" %s\n"%fastq)
    k = short = 0
    # open gzipped files
    if fastq.endswith('.gz'):
        fastq = gzip.open(fastq)
    for i, r in enumerate(fastq_parser(fastq), 1):
        if verbose and not i % 1e5:
            sys.stderr.write("  %s long: %s [%5.2f%s]  PE reads: %s  \r"%(i, i-short, 100.0*(i-short)/i, '%', k))
        # skip too short
        if len(r)<minLen+readlen:
            short += 1
            continue
        # generate reads
        k += read2mates(out, r, readlen=50, window=25, minLen=400)
    sys.stderr.write("  %s long: %s [%5.2f%s]  PE reads: %s  \n"%(i, i-short, 100.0*(i-short)/i, '%', k))

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.10a')   
    parser.add_argument("-i", "--fastq", nargs="+", 
                        help="FASTQ files with long reads")
    parser.add_argument("-o", "--outbase", default="", 
                        help="output base name [inputfn_1.fq.gz & inputfn_2.fq.gz]")
    parser.add_argument("-l", "--readlength", default=50, 
                        help="mates read length [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    for fn in o.fastq:
        fastq2mates(fn, o.outbase, o.readlength, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
