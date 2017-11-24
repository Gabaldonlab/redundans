#!/usr/bin/env python
desc="""De novo assembly module

TBA:
- catch combination of .fq & fq.gz
- specify IP / OP correctly
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 20/11/2017
"""

import commands, os, sys, tempfile
from datetime import datetime
from subprocess import Popen, PIPE
from fastq2insert_size import fastq2insert_size

# update environmental PATH with script dir
os.environ["PATH"] = "%s:%s"%(os.path.dirname(os.path.abspath(sys.argv[0])), os.environ["PATH"])

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

def get_named_fifo():
    """Return named FIFO"""
    tmpfn = tempfile.mktemp()
    os.mkfifo(tmpfn)
    return tmpfn
    
def run_assembly(prefix, fastq, threads, tmpdir, log, locallog):
    """Execute platanus assemble"""
    # create named FIFO
    tmp = get_named_fifo()
    # run assembly
    cmd = "platanus assemble -tmp %s -t %s -o %s -f %s" % (tmpdir, threads, prefix, tmp)
    p = Popen(cmd.split(), stdout=locallog, stderr=locallog)
    if log:
        log.write(" %s\n"%cmd)
    # write FastA to fifo
    with open(tmp, 'wb') as pipe:
        cmd = ["cat", ] + fastq
        if filter(lambda fn: fn.endswith('.gz'), fastq):
            cmd[0] = "zcat" 
        parser = Popen(cmd, stdout=pipe, stderr=locallog)
        parser.wait()
    # wait for process to finish & rm fifo
    p.wait()
    os.unlink(tmp)
    
def run_scaffolding(prefix, fastq, threads, tmpdir, log, locallog, limit=1.):
    """Execute platanus assemble"""
    tmp = get_named_fifo()
    cmd = "platanus scaffold -tmp %s -t %s -o %s -c %s_contig.fa -b %s_contigBubble.fa -ip1 %s" % (tmpdir, threads, prefix, prefix, prefix, tmp)
    p = Popen(cmd.split(), stdout=locallog, stderr=locallog)
    if log:
        log.write(" %s\n"%cmd)
    # write shuffled FastQ to fifo
    with open(tmp, 'wb') as pipe:
        parser = Popen(["fastq2shuffled.py", ] + fastq, stdout=pipe, stderr=locallog)
        parser.wait()
    # wait for process to finish & rm fifo
    p.wait()
    os.unlink(tmp)

def run_gapclosing(prefix, fastq, threads, tmpdir, log, locallog, limit=1.):
    """Execute platanus assemble"""
    tmp = get_named_fifo()
    cmd = "platanus gap_close -tmp %s -t %s -o %s -c %s_scaffold.fa -ip1 %s" % (tmpdir, threads, prefix, prefix, tmp)
    p = Popen(cmd.split(), stdout=locallog, stderr=locallog)
    if log:
        log.write(" %s\n"%cmd)
    # write shuffled FastQ to fifo
    with open(tmp, 'wb') as pipe:
        parser = Popen(["fastq2shuffled.py", ] + fastq, stdout=pipe, stderr=locallog)
        parser.wait()
    # wait for process to finish & rm fifo
    p.wait()
    os.unlink(tmp)
    
def denovo(outdir, fastq, threads, verbose, log, tmpdir='/tmp'):
    """Select best libriaries and run de novo assembly using idba_ud"""
    # create missing outdir
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    # select best library(ies)
    if verbose: 
        log.write(" %s libs: %s\n"%(len(fastq), ", ".join(fastq)))
    bestfastq, seqsize = get_best_lib(fastq)
    if verbose: 
        log.write("  %s libs (~%.2f Mbases) selected for assembly: %s\n"%(len(bestfastq), seqsize, ", ".join(bestfastq)))
    # platanus
    prefix = os.path.join(outdir, "out"); locallog = open(prefix+".log", "a")
    run_assembly(prefix, bestfastq, threads, tmpdir, log, locallog)
    outfn = prefix + "_contig.fa"
    
    # estimate insert size
    # fq1, fq2, readlen, ismedian, ismean, isstd, pairs, orientation
    libdata = fastq2insert_size(log, fastq, prefix+"_contig.fa", threads=threads)
    pelibs = filter(lambda x: x[-1] in ('FR', ) and 100<x[4]<1000, libdata) #'RF'
    if pelibs:
        pelibs = sorted(pelibs, key=lambda x: x[4])
        pefastq = []
        for x in pelibs: pefastq += x[:2]
        if verbose:
            log.write("  selected %s lib(s) for scaffolding & gap closing: %s\n"%(len(pelibs), ", ".join(pefastq)))
        run_scaffolding(prefix, pefastq, threads, tmpdir, log, locallog)
        run_gapclosing(prefix, pefastq, threads, tmpdir, log, locallog)
        outfn = prefix + "_gapClosed.fa"
    elif verbose:
        log.write("  No suitable libs for scaffolding & gap closing!\n")
    locallog.close()
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
    parser.add_argument("--tmp", default='/tmp', help="tmp directory [%(default)s]")
    
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
    o = parser.parse_args()
    if o.verbose:
        o.log.write("Options: %s\n"%str(o))

    outfn = denovo(o.outdir, o.fastq, o.threads, o.verbose, o.log, o.tmp)
    sys.stderr.write("Final contigs in: %s\n"%outfn)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
