#!/usr/bin/env python
desc="""Separate heterozygous contigs into two files based on identity to homozygous (one parental) reference
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Warsaw, 20/04/2017
"""

import gzip, math, os, sys, subprocess, tempfile
from datetime import datetime
from FastaIndex import FastaIndex

def run_last(ref, fasta, threads, verbose):
    """Start LAST with multi-threads"""
    if verbose:
        sys.stderr.write(" Running LAST...\n")
    # build db
    if not os.path.isfile(fasta+".suf"):
        os.system("lastdb %s %s" % (ref, ref))
    # run LAST
    fname = tempfile.mktemp()
    fifo = os.mkfifo(fname)
    args1 = ["lastal", "-P", str(threads), ref, fasta]
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=sys.stderr)
    proc2 = subprocess.Popen(["last-split",], stdin=proc1.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    proc3 = subprocess.Popen(["maf-convert", "tab"], stdin=proc2.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    proc4 = subprocess.Popen(["tee", "-a", fname], stdin=proc3.stdout, stdout=subprocess.PIPE, stderr=sys.stderr)
    #os.system("")
    proc4b = subprocess.Popen(["last-dotplot", fname, "%s.dotplot.png"%fasta])
    return proc4.stdout
    
def get_matches(ref, fasta, threads, verbose):
    """Return LASTal hits passing identity and overlap thresholds"""
    # execute last
    #last = 
    q2matches = {}
    for l in run_last(ref, fasta, threads, verbose): # gzip.open("pilon2.ABT301.fasta.tab.gz"): #
        if l.startswith('#'):
            continue
        # unpack
        (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
        (score, qstart, qalg, qsize, tstart, talg, tsize) = list(map(int, (score, qstart, qalg, qsize, tstart, talg, tsize)))
        if q not in q2matches:
            q2matches[q] = [0, 0, qsize]
        # count score and qalg
        q2matches[q][0] += score
        q2matches[q][1] += qalg
    return q2matches

def plot_scatter(fname, identities, algsizes, gcs):
    """Plot scatter with identity to reference and GC.

    Points are scaled proportionally to algsizes (https://stackoverflow.com/a/14860958/632242)
    """
    try:
        import numpy as np
        import matplotlib 
        matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
        import matplotlib.pyplot as plt
    except:
        sys.stderr.write("[WARNING] numpy or matplotlib missing! Cannot plot scatter plot\n")
        return

    s = [.1*(18/math.log(max(algsizes), 10))**math.log(x, 10) for x in algsizes]#; print s #(12/math.log(max(algsizes), 10))
        
    fig = plt.figure()
    axis = plt.scatter(identities, gcs, s=s, c=algsizes)
    #plt.legend(loc=2)
    plt.title("Identity vs GC")
    plt.xlabel("Identity to reference")
    plt.ylabel("GC [%]")
    #plt.ylim(38, 45)
    cbar = fig.colorbar(axis, ticks=[min(algsizes), np.mean(algsizes), max(algsizes)])
    cbar.ax.set_title("Contig size", size=12)
    #cbar.ax.set_yticklabels(["Low","Medium","High"])
    
    fig.savefig(fname+".GC.scatter.png", dpi=300)
    
def plot_histograms(fname, identities, algsizes):
    """Plot histogram for matches"""
    try:
        import numpy as np
        import matplotlib 
        matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend
        import matplotlib.pyplot as plt
    except:
        sys.stderr.write("[WARNING] numpy or matplotlib missing! Cannot plot histogram\n")
        return
        
    fig = plt.figure()
    # get bins
    bins = np.arange(.5, 1.01, 0.01)
    
    counts = [0]*len(bins)
    sizes = [0]*len(bins)
    for i, isize in zip(np.digitize(identities, bins, right=1), algsizes):
        counts[i] += 1
        sizes[i] += isize
    #print max(identities), bins
    bins -= 0.01    
    # plot no. of contigs at give identity
    plt.subplot(211)
    #plt.bar(bins*100, bestcounts, color="red", label="best", alpha=1.0)
    plt.bar(bins*100, counts, color="red", label="all", alpha=0.33)
    plt.xlim(50, 100)
    plt.legend(loc=2)
    plt.title("Identity between contigs")
    plt.ylabel("No. of contigs")

    # plot cumulative alignment size at give identity
    plt.subplot(212)
    #plt.bar(bins*100, np.array(bestsizes)/1e6, color="blue", label="best", alpha=1.0)
    plt.bar(bins*100, np.array(sizes)/1e6, color="blue", label="all", alpha=0.33)
    plt.xlim(50, 100)
    plt.legend(loc=2)
    plt.xlabel("Identity [%]")
    plt.ylabel("Cumulative alignment size [Mb]")
    fig.savefig(fname+".homo.hist.png", dpi=300)

def fasta2split(outs, ref, fasta, identity, overlap, minLength, threads, verbose, log=sys.stderr):
    """ """
    sizes = [0, 0]
    #create/load fasta index
    if verbose:
        log.write("Indexing fasta...\n")
    faidx = FastaIndex(fasta)
    
    # filter alignments & remove redundant
    if verbose:
        log.write("Parsing alignments...\n")
    # get matches
    q2matches = get_matches(ref, fasta, threads, verbose)

    identities, algsizes, qsizes, gcs = [], [], [], []
    for q, (score, algsize, qsize) in q2matches.items():
        _identity = 1.0 * (score+(algsize-score)/2) / algsize
        if 1.*algsize/qsize>=overlap:
            identities.append(_identity)
            algsizes.append(algsize)
            qsizes.append(qsize)
            (A, C, G, T) = faidx.id2stats[q][-4:]
            gc = 100.0*(G + C) / sum((A, C, G, T))
            gcs.append(gc) #stats[-4:] for stats in self.id2stats.itervalues()
        # split
        if _identity > identity and 1.*algsize/qsize>=overlap:
            i = 0
        else:
            i = 1
        sizes[i] += qsize #faidx.id2stats[q][0]
        outs[i].write(faidx[q])
    # check missing
    i = 1    
    for q in faidx:
        if q not in q2matches:
            sizes[i] += faidx.id2stats[q][0]
            outs[i].write(faidx[q])
    log.write("%s bp split into %s files: %s\n"%(sum(sizes), len(sizes), ', '.join(map(str, sizes))))
        
    # plot
    plot_scatter(fasta, identities, qsizes, gcs)
    plot_histograms(fasta, identities, algsizes)

def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.01d')   
    parser.add_argument("-v", "--verbose", action="store_true", help="verbose")    
    parser.add_argument("-f", "--fasta", type=file, help="contigs FastA")
    parser.add_argument("-r", "--ref", type=file, help="homozygous reference FastA")
    parser.add_argument("-i", "--identity", default=0.93, type=float, help="identity split [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, help="max threads to run [%(default)s]")
    parser.add_argument("--overlap",  default=0.33, type=float, help="min. overlap [%(default)s]")
    parser.add_argument("--minLength", default=200, type=int, help="min. contig length [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    outs = [open(o.fasta.name+".homo_%s.fa"%i, "w") for i in range(1, 3)]
    fasta2split(outs, o.ref.name, o.fasta.name, o.identity, o.overlap, o.minLength, o.threads, o.verbose)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
