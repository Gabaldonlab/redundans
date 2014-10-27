#!/usr/bin/env python
desc="""Heterozygous genome assembly pipeline.
It consists of three steps:
1. assembly reduction
2. scaffolding
3. gap closing

Note, FASTQ libraries need to be

To be done:
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow/Warsaw, 17/10/2014
"""

import os, resource, sys
try:
    import pysam
except:
    sys.stderr.write("Install pysam: `sudo easy_install -U pysam`!\n")
    #sys.exit(1)
from datetime import datetime
import numpy as np
#from scipy import stats, signal

class redundants(object):
    """Store BAM related data. Take BAM file name as input."""
    def __init__(self, *args):
        self._get_parameters(args)
        self._prepare_streams()

    def _get_params(self, libs, fasta, outdir, mapq, threads, \
                    ident, joins, limit, iters, \
                    reduction, scaffolding, gapclosing, \
                    verbose, log):
        """Set parameters"""
        #general
        self.libs   = libs
        self.fasta  = fasta
        self.outdir = outdir
        self.mapq   = mapq
        self.threads = threads
        #specific options
        self.ident  = ident
        self.joins  = joins
        self.limit  = limit
        self.iters  = iters
        #steps
        self.reduction   = reduction
        self.scaffolding = scaffolding
        self.gapclosing  = gapclosing
        
    def _prepare_streams(self):
        #prepare logging stream
        if log:
            self.log = log
        elif verbose:
            self.log = sys.stderr
        else:
            self.log     = None
        #prepare outdir
        if os.path.isdir(self.outdir):
            sys.stderr.write("Directory %s exists!\n")
            sys.exit(1)
        os.makedirs(self.outdir)
        
    def get_isize_stats(self, limit=1e5): 
        """Estimate insert size median, mean and stdev.
        Also count pair orientations and select main orientation.
        """
        if self.log:
            self.log.write("Estimating insert size stats...\n")
        isizes = []
        self.pairs = [0, 0, 0, 0]
        #read from stdin
        for alg in pysam.Samfile(self.bam):
            #take only reads with good alg quality and one read per pair
            if alg.mapq < self.mapq or alg.isize < 1:
                continue
            #store isize
            isizes.append(alg.isize)
            #store pair orientation
            self.pairs[self.alg2orientation(alg)] += 1
            #stop if limit reached
            if len(isizes) >= limit:
                break
        #get rid of right 5 percentile
        maxins = stats.scoreatpercentile(isizes, 100-self.q)
        minins = stats.scoreatpercentile(isizes, self.q)
        isizes = filter(lambda x: minins<x<maxins, isizes)
        #store
        self.isize_median = np.median(isizes)
        self.isize_mean   = np.mean(isizes)
        self.isize_stdev  = np.std(isizes)
                
    def parse(self, test=0):
        """Parse sam alignments and store info"""
        #parse algs
        if self.log:
            self.log.write("Parsing alignments...\n")
        pchrom = ""
        for i, alg in enumerate(self.sam, 1):
            if test and i > test:
                break
            #write log
            if self.log and not i % 1e5:
                info = " %s [%.1f%s]  reads for dels: %s dups: %s ins: %s invs: %s trans: %s [%s Mb]\r"
                self.log.write(info % (i, i*100.0/self.nalgs, '%', len(self.delReads), \
                                len(self.dupReads), len(self.insReads), len(self.invReads), len(self.traReads), \
                                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024))
            #skip unmapped and secondary alignments
            if alg.rname<0 or alg.is_secondary:
                continue
            #add read
            self.add_read(alg)
        if self.log:
            self.log.write(" %s alignments parsed. \n"%i)
        #dump all important info
        if not self.nodump and not os.path.isfile(self.bamdump):
            self.sv2bam()
        #get mean rlen
        if not self.rlen:
            self.rlen = np.mean([alg.rlen for alg in self.delReads])
            if self.log:
                self.log.write(" Mean read length: %.2f \n"%self.rlen)
        #call variants
        self.call_variants()        
            
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-i", "--fastq", nargs="+", 
                        help="FASTQ PE/MP files")
    parser.add_argument("-f", "--fasta", type=file, 
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
    redu.add_argument("--ident",        default=0.8, type=float,
                      help="align l reads [%(default)s]")
    ##missing redu options
    scaf = parser.add_argument_group('Scaffolding options')
    scaf.add_argument("-j", "--joins",  default=5, type=int, 
                      help="min k pairs to join contigs [%(default)s]")
    scaf.add_argument("-l", "--limit",  default=5000000, type=int, 
                      help="align l reads [%(default)s]")
    scaf.add_argument("-iters",         default=2, type=int, 
                      help="scaffolding iterations per library  [%(default)s]")
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
        
    #initialise pipeline
    redundants(o.fastq, o.fasta, o.outdir, o.mapq, o.threads, \
               o.ident, \
               o.joins, o.limit, o.iters, \
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
