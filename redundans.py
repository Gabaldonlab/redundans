#!/usr/bin/env python
desc="""Heterozygous genome assembly pipeline.
It consists of three steps:
1. assembly reduction
2. scaffolding
3. gap closing

Note, FASTQ libraries need to be

PREREQUISITIES:
- genome assembly (fasta)
- BLAT
- BWA
- Biopython

To be done:
 - check if files exist
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow/Warsaw, 17/10/2014
"""

import os, resource, sys
from datetime import datetime
import numpy as np

from fasta2homozygous import fasta2homozygous
from fastq2sspace import fastq2sspace
from fastq2insert_size import fastq2insert_size

def get_timestamp():
    """Return formatted date-time string"""
    return "\n%s\n[%s] "%("#"*15, datetime.ctime(datetime.now()))

def symlink(file1, file2):
    """Create symbolic link taking care of real path."""
    os.symlink(os.path.join(os.path.realpath(os.path.curdir), file1), file2)

def run_scaffolding(outdir, scaffoldsFname, fastq, reducedFname, mapq, threads, \
                    joins, limit, iters, sspacebin, verbose, lib=""):
    """Execute scaffolding step"""
    # get libraries statistics using 2% of read limit
    libdata = fastq2insert_size(sys.stderr, fastq, reducedFname, mapq, threads, \
                                limit/50, verbose)
    # separate paired-end & mate pairs
    ## also separate 300 and 600 paired-ends
    
    # create symlink to reduced contigs
    tmpname = os.path.join(outdir,"sspace.%s.%s.final.scaffolds.fasta")
    symlink(reducedFname, tmpname%(0, iters))
    
    # run scaffolding using libraries with increasing insert size in multiple iterations
    for i, data in enumerate(libraries, 1):
        # init empty handles
        libnames, libFs, libRs, orientations, libIS, libISStDev = [], [], [], [], [], []
        for j in range(1, iters+1):
            out = "sspace.%s.%3s"%(libType, i)
            lib = ""
            # run fastq scaffolding
            fastq2sspace(out, open(reducedFname), lib, libnames, libFs, libRs, orientations, \
                         libIS, libISStDev, threads, mapq, upto, joins, \
                         sspacebin, verbose)
            
    # create symlink to final scaffolds
    symlink(tmpname%(i, iters), scaffoldsFname)
    
def redundants(libs, fasta, outdir="redundans", mapq=10, threads=1, \
               identity=0.8, overlap=0.75, minLength=200, \
               joins=5, limit=10e7, iters=3, spacebin, \
               reduction=1, scaffolding=1, gapclosing=1, \
               verbose=1, log=sys.stderr):
    """Launch redundans pipeline."""
    # redirect stderr
    sys.stderr = log
    
    # prepare outdir or quit if exists
    if os.path.isdir(outdir):
        log.write("Directory %s exists!\n"%outdir)
        #sys.exit(1)
    else:
        os.makedirs(outdir)

    # check if all files exists
    #_check_files(log, fasta, libs)
    
    # REDUCTION
    contigsFname = os.path.join(outdir, "contigs.fa")
    reducedFname = os.path.join(outdir, "contigs.reduced.fa")
    # link contigs & genome
    symlink(fasta, contigsFname)
    if reduction:
        if verbose:
            sys.stderr.write("%sReduction...\n"%timestamp())
            sys.stderr.write("#file name\tgenome size\tcontigs\theterozygous size\t[%]\theterozygous contigs\t[%]\tidentity [%]\tpossible joins\thomozygous size\t[%]\thomozygous contigs\t[%]\n")
        with open(reducedFname, "w") as out:
            fasta2homozygous(out, open(contigsFname), identity, overlap, minLength)
    else:
        symlink(contigsFname, reducedFname)

    # SCAFFOLDING
    scaffoldsFname = os.path.join(outdir, "scaffolds.fa")
    if scaffolding:
        if verbose:
            sys.stderr.write("%sScaffolding...\n"%timestamp())
        run_scaffolding()
    else:
        symlink(reducedFname, scaffoldsFname)
        
    # GAP CLOSING
    nogapsFname = os.path.join(outdir, "scaffolds.fa")
    if gapclosing:
        if verbose:
            sys.stderr.write("%sGap closing...\n"%timestamp())
        #run_gapclosing()
    else:
        symlink(scaffoldsFname, nogapsFname)

    # FASTA STATS
    
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-i", "--fastq", nargs="+", 
                        help="FASTQ PE/MP files")
    parser.add_argument("-f", "--fasta", 
                        help="assembly FASTA file")
    parser.add_argument("-o", "--outdir",  default="redundants", 
                        help="output directory")
    parser.add_argument("-q", "--mapq",    default=20, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("-t", "--threads", default=8, type=int, 
                        help="max threads to run [%(default)s]")
    parser.add_argument("--log",           default=None, type=argparse.FileType('w'), 
                        help="save log to file [stderr]")
    redu = parser.add_argument_group('Reduction options')
    redu.add_argument("--identity",        default=0.8, type=float,
                      help="min. identity [%(default)s]")
    redu.add_argument("--overlap",         default=0.75, type=float,
                      help="min. overlap  [%(default)s]")
    redu.add_argument("--minLength",       default=200, type=int, 
                      help="min. contig length [%(default)s]")
    ##missing redu options
    scaf = parser.add_argument_group('Scaffolding options')
    scaf.add_argument("-j", "--joins",  default=5, type=int, 
                      help="min k pairs to join contigs [%(default)s]")
    scaf.add_argument("-l", "--limit",  default=5000000, type=int, 
                      help="align at most l reads [%(default)s]")
    scaf.add_argument("-iters",         default=2, type=int, 
                      help="scaffolding iterations per library  [%(default)s]")
    scaf.add_argument("--sspacebin",    default="~/src/SSPACE/SSPACE_Standard_v3.0.pl", 
                       help="SSPACE path  [%(default)s]")
    gaps = parser.add_argument_group('Gap closing options')
    #gaps.add_argument("-l", "--limit",  default=7e, type=int, 
    #                  help="align l reads [%(default)s]")
    skip = parser.add_argument_group('Skip below steps (all performed by default)')
    skip.add_argument('--reduction',   action='store_false', default=True)   
    skip.add_argument('--scaffolding', action='store_false', default=True)   
    skip.add_argument('--gapclosing',  action='store_false', default=True)   
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    # initialise pipeline
    redundants(o.fastq, o.fasta, o.outdir, o.mapq, o.threads, \
               o.identity, o.overlap, o.minLength, \
               o.joins, o.limit, o.iters, o.spacebin, \
               o.reduction, o.scaffolding, o.gapclosing, \
               o.verbose, o.log)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
