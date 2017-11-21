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

def get_readlen_and_seqsize(fastq, limit=1000):
    """Return mean readlen and total seqsize in Mb"""
    handle = open(fastq)
    ratio = 1.
    # open gzip file as subprocess
    if handle.name.endswith('.gz'):
        zcat = Popen(['zcat', handle.name], bufsize=-1, stdout=PIPE, stderr=PIPE)
        handle = zcat.stdout
        ratio = 0.33 # compression rate 9201307/25077788=0.367 or 903148/3017787=0.299
    # calculate seq length vs total data length
    i = datalen = seqlen = 0
    for l in handle:
        i += 1
        if i/4>=limit:
            break
        datalen += len(l)
        l = l.strip()
        if not l: continue
        if i%4==2:
            seqlen += len(l)
    # stop handle
    handle.close()
    # mean read len    
    readlen = int(round(1.*seqlen/limit))
    # estimate size
    fsize = os.path.getsize(fastq) / ratio
    seqsize = round(seqlen * fsize / datalen / 1e6, 3)
    return readlen, seqsize

def get_best_lib(fastq, frac=0.66, verbose=0):
    """Select best library"""
    lib2data = {fn: get_readlen_and_seqsize(fn) for fn in fastq}#; print lib2data
    # select largest lib in 70-150 read len range
    _fastq = filter(lambda x: 70<=lib2data[x][0]<=150, fastq)
    if _fastq:
        fastq = _fastq
    # return libs having frac*max size
    sizes = [lib2data[fn][1] for fn in fastq]    
    maxsize = max(sizes)
    return filter(lambda x: lib2data[x][1] >= frac*maxsize, fastq), sum(sizes)

def denovo(outdir, fastq, threads, verbose, log, tmp='/tmp'):
    """Select best libriaries and run de novo assembly using idba_ud"""
    if verbose:    
        log.write(" %s libs: %s\n"%(len(fastq), ", ".join(fastq)))
        
    # create missing outdir
    if not os.path.isdir(os.path.dirname(outdir)):
        os.makedirs(os.path.dirname(outdir))

    # select best library(ies)
    fastq, seqsize = get_best_lib(fastq)
    if verbose:    
        log.write("  %s libs (~%.2f Mbases) selected for assembly: %s\n"%(len(fastq), seqsize, ", ".join(fastq)))
    # SPAdes <5 Gb sequence
    if seqsize < 5*1e3:
        cmd = "spades.py --only-assembler -t %s -o %s -s %s"%(threads, outdir, " -s ".join(fastq))
        if verbose:
            log.write(" %s\n"%cmd)
        p = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        outfn = os.path.join(outdir, "contigs.fasta")
    else:
        # platanus
        prefix = os.path.join(outdir, "out")
        # create named FIFO
        tmp = tempfile.mktemp()
        os.mkfifo(tmp)
        # run assembly
        #cmd = "idba_ud --mink 31 --maxk 101 --step 10 --num_threads %s -o %s -r %s"%(threads, outdir, tmp)
        cmd = "platanus assemble -t %s -o %s -f %s" % (threads, prefix, tmp)
        p = Popen(cmd.split(), stdout=PIPE, stderr=PIPE)
        if verbose:
            log.write(" %s\n"%cmd)
        # write FastA to fifo
        with open(tmp, 'wb') as pipe:
            fastq2fasta = Popen(["fastq2fasta.py", "-l 31", "-q 20", "-i"] + fastq, stdout=pipe, stderr=PIPE)
        # wait for process to finish
        fastq2fasta.wait()
        # rm fifo 
        os.unlink(tmp)
        outfn = prefix + "_contigs.fa"
    # read stdout to finish the process
    stdout = p.stdout.readlines()#; print stdout
    return outfn
        
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

    outfn = denovo(o.outdir, o.fastq, o.threads, o.verbose, o.log)
    sys.stderr.write("Final contigs in: %s\n"%outfn)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
