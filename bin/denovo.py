#!/usr/bin/env python
desc="""De novo assembly module
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 20/11/2017
"""

import commands, os, sys, tempfile
from datetime import datetime
from subprocess import Popen, PIPE

def get_best_lib(fastq):
    """Select best library"""
    return fastq

def denovo(outdir, fastq, threads, verbose, log):
    """Select best libriaries and run de novo assembly using idba_ud"""
    # create missing outdir
    if not os.path.isdir(os.path.dirname(outdir)):
        os.makedirs(os.path.dirname(outdir))

    # select best library(ies)
    fastq = get_best_lib(fastq)
        
    # create named FIFO
    tmp = tempfile.mktemp()
    os.mkfifo(tmp)
    # run assembly
    #cmd = "platanus assemble -tmp %s -t %s -o %s -f <(zcat %s)" % (tmp, threads, prefix, " ".join(fastq))
    cmd1 = "idba_ud --num_threads %s -o %s -r %s"%(threads, outdir, tmp) 
    p = Popen(cmd1.split(), stdout=PIPE, stderr=PIPE)
    if verbose:
        log.write(" %s\n"%cmd1)
    # write FastA to fifo
    with open(tmp, 'wb') as pipe:
        fastq2fasta = Popen(["fastq2fasta.py", "-i"] + fastq, stdout=pipe, stderr=PIPE)
    # wait for process to finish
    fastq2fasta.wait()
    # rm fifo    
    os.unlink(tmp)
        
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", "--verbose",  action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.13c')
    parser.add_argument("-i", "--fastq", nargs="*", default=[], help="FASTQ PE / MP files")
    parser.add_argument("-o", "--outdir", default="denovo", help="output directory [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="max threads to run [%(default)s]")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), help="output log to [stderr]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
    o = parser.parse_args()
    if o.verbose:
        o.log.write("Options: %s\n"%str(o))

    denovo(o.outdir, o.fastq, o.threads, o.verbose, o.log)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
