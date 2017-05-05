#!/usr/bin/env python
desc="""Align genome onto itself (LAST) and keep only the longest
from heterozygous (redundant) contigs/scaffolds.

TO ADD:
- scaffold extension based on overlapping matches
- reporting of haplotypes
- recognise heterozygous contigs with translocations
- guess which identity cutoff will be the best
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 26/08/2014
"""

import gzip, os, sys, subprocess
from datetime import datetime
from FastaIndex import FastaIndex

import numpy as np

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
os.environ["PATH"] = "%s:%s"%(root, os.environ["PATH"])

def run_last(fasta, identity, threads, minLength=200, verbose=1):
    """Start LAST with multi-threads"""
    if verbose:
        sys.stderr.write(" Running LAST...\n")
    # build db
    ref = fasta
    if not os.path.isfile(ref+".suf"):
        os.system("lastdb -W 11 %s %s" % (ref, fasta))
    # run LAST
    args1 = [ "lastal", "-P", str(threads), "-f", "TAB", ref, fasta] 
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=sys.stderr)
    return proc1

def get_best_match(matches, q, qsize, identityTh, overlapTh):
    """Return best match for particular query"""
    if not matches:
        return
    # get best t mach - the one having max cumulative score
    t, (score, qalg) = sorted(matches.iteritems(), key=lambda x: x[1][0], reverse=1)[0]
    # get score, identity & overlap # LASTal is using +1/-1 for match/mismatch, while I need +1/0
    identity = 1.0 * (score+(sum(qalg)-score)/2) / sum(qalg)
    overlap  = 1.0 * sum(qalg>0) / qsize
    # filter by identity and overlap
    if identity >= identityTh and overlap >= overlapTh:
        return score, t, q, sum(qalg>0), identity, overlap
        
def fasta2hits(fasta, threads, identityTh, overlapTh, minLength, verbose):
    """Return LASTal hits passing identity and overlap thresholds"""
    # execute last
    last = run_last(fasta.name, identityTh, threads, minLength, verbose)
    pq, pqsize = '', 0
    matches = {}
    for l in last.stdout: 
        if l.startswith('#'): 
            continue
        # unpack
        (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
        (score, qstart, qalg, qsize, tstart, talg, tsize) = map(int, (score, qstart, qalg, qsize, tstart, talg, tsize))
        # skip reverse matches
        if tsize<minLength or qsize<minLength or t==q or tsize < qsize or tsize==qsize and t<q: 
            continue
        # report previous query
        if pq != q:
            if get_best_match(matches, pq, pqsize, identityTh, overlapTh): 
                yield get_best_match(matches, pq, pqsize, identityTh, overlapTh)
            # reset
            pq, pqsize = q, qsize
            matches = {}
        if t not in matches:
            matches[t] = [0, np.zeros(qsize, dtype='uint8')]
        # get qstart & qend
        if qstrand=="+":
            s = qstart
            e = s + qalg
        else:
            e = qsize - qstart
            s = qsize - qstart - qalg
        # allow partial query overlaps on target, but skip if more than 10% overlap with existing matches
        if sum(matches[t][1][s:e]>1) > 0.1*qalg:
            continue
        matches[t][0] += score
        matches[t][1][s:e] += 1
  
    # yield last bit
    if get_best_match(matches, pq, pqsize, identityTh, overlapTh): 
        yield get_best_match(matches, pq, pqsize, identityTh, overlapTh)
        
def fasta2skip(out, fasta, faidx, threads, identityTh, overlapTh, minLength, verbose):
    """Return dictionary with redundant contigs and their best alignments"""
    # get hits generator
    hits = fasta2hits(fasta, threads, identityTh, overlapTh, minLength, verbose)
    # iterate through hits
    identities, sizes = [], []
    contig2skip = {c: 0 for c in faidx} 
    for i, (score, t, q, algLen, identity, overlap) in enumerate(hits, 1):
        # skip contigs already marked as heterozygous
        if t not in contig2skip:
            sys.stderr.write(' [ERROR] `%s` (%s) not in contigs!\n'%(t, str(hits[i-1])))
            continue
        # catch missing contigs
        if q not in contig2skip:
            sys.stderr.write(' [ERROR] `%s` (%s) not in contigs!\n'%(q, str(hits[i-1])))
            continue
        # store first match or update best match
        if not contig2skip[q] or score > contig2skip[q][0]:
            contig2skip[q] = (score, t, algLen, identity, overlap)
        # store identity and alignment for plotting
        identities.append(identity)
        sizes.append(algLen)
    # plot histogram of identities
    plot_histograms(out.name, contig2skip, identities, sizes)
    return contig2skip

def plot_histograms(fname, contig2skip, identities, algsizes):
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
    contigs = contig2skip.keys()
    best = [contig2skip[c][3] for c in contigs if contig2skip[c]]
    bestalgsizes = [contig2skip[c][2] for c in contigs if contig2skip[c]]
    # get bins
    bins = np.arange(.5, 1.01, 0.01)
    
    # get counts
    bestcounts = [0]*len(bins)
    bestsizes = [0]*len(bins)
    for i, isize in zip(np.digitize(best, bins, right=1), bestalgsizes):
        bestcounts[i] += 1
        bestsizes[i] += isize
            
    counts = [0]*len(bins)
    sizes = [0]*len(bins)
    for i, isize in zip(np.digitize(identities, bins, right=1), algsizes):
        counts[i] += 1
        sizes[i] += isize
        
    bins -= 0.01
    # plot no. of contigs at give identity
    plt.subplot(211)
    plt.bar(bins*100, bestcounts, color="red", label="best", alpha=1.0)
    plt.bar(bins*100, counts, color="grey", label="all", alpha=0.33)
    plt.xlim(50, 100)
    plt.legend(loc=2)
    plt.title("Identity between contigs")
    plt.ylabel("No. of contigs")

    # plot cumulative alignment size at give identity
    plt.subplot(212)
    plt.bar(bins*100, np.array(bestsizes)/1e6, color="blue", label="best", alpha=1.0)
    plt.bar(bins*100, np.array(sizes)/1e6, color="grey", label="all", alpha=0.33)
    plt.xlim(50, 100)
    plt.legend(loc=2)
    plt.xlabel("Identity [%]")
    plt.ylabel("Cumulative alignment size [Mb]")
    fig.savefig(fname+".hist.png", dpi=300)
    
def fasta2homozygous(out, fasta, identity, overlap, minLength, threads=1, verbose=0, log=sys.stderr):
    """Parse alignments and report homozygous contigs.
    
    Return genomeSize, no. of contigs, removed contigs size & number
    and average identity between reduced contigs.
    """
    merged = []
    
    #create/load fasta index
    if verbose:
        log.write("Indexing fasta...\n")
    faidx = FastaIndex(fasta)
    genomeSize = faidx.genomeSize
    
    # filter alignments & remove redundant
    if verbose:
        log.write("Parsing alignments...\n")
    contig2skip = fasta2skip(out, fasta, faidx, threads, identity, overlap, minLength, verbose)
    
    #report homozygous fasta
    nsize, k, skipped, ssize, avgIdentity = save_homozygous(out, faidx, contig2skip, minLength, verbose)
    
    #summary    
    info = "%s\t%s\t%s\t%s\t%.2f\t%s\t%.2f\t%.3f\t%s\t%s\t%.2f\t%s\t%.2f\n"
    log.write(info%(fasta.name, genomeSize, len(faidx), ssize, 100.0*ssize/genomeSize, \
                    skipped, 100.0*skipped/len(faidx), avgIdentity, len(merged), \
                    nsize, 100.0*nsize/genomeSize, k, 100.0*k/len(faidx)))

    return genomeSize, len(faidx), ssize, skipped, avgIdentity

def save_homozygous(out, faidx, contig2skip, minLength, verbose):
    """Save homozygous contigs to out stream

    Here you could learn from distibution of identities,
    what really is the reasonable identity cut-off. 
    """
    k = skipped = ssize = nsize = identities = algLengths = 0
    # store skipped hetero contigs stats
    out2 = open(out.name+".hetero.tsv", "w")
    out2.write("#contig\tsize\ttarget\titentity\toverlap\n")
    # process contigs starting from the largest
    for i, c in enumerate(faidx, 1):
        # skip short
        if faidx.id2stats[c][0] < minLength:
            skipped += 1
            ssize   += faidx.id2stats[c][0]
        # skip hetero
        elif contig2skip[c]: 
            skipped += 1
            ssize   += faidx.id2stats[c][0]
            score, t, algLen, identity, overlap = contig2skip[c]
            out2.write("%s\t%s\t%s\t%s\t%s\n"%(c, faidx.id2stats[c][0], t, identity, overlap))
            # update identities and lengths
            identities += identity*algLen
            algLengths += algLen
        # save sequence
        else:
            out.write(faidx[c])
            # update counters
            k += 1
            nsize += faidx.id2stats[c][0]
    # close out2
    out2.close()
    # calculate average identity        
    avgIdentity = 0
    if algLengths:
        avgIdentity = 100.0 * identities / algLengths
    return nsize, k, skipped, ssize, avgIdentity
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.01d')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
    parser.add_argument("-i", "-f", "--fasta", nargs="+", type=file, help="FASTA file(s)")
    parser.add_argument("-t", "--threads", default=4, type=int, help="max threads to run [%(default)s]")
    parser.add_argument("--identity", default=0.51, type=float, help="min. identity [%(default)s]")
    parser.add_argument("--overlap", default=0.8, type=float, help="min. overlap [%(default)s]")
    parser.add_argument("--minLength",   default=200, type=int, help="min. contig length [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #process fasta
    sys.stderr.write("Homozygous assembly/ies will be written with input name + '.homozygous.fa.gz'\n")
    sys.stderr.write("#file name\tgenome size\tcontigs\theterozygous size\t[%]\theterozygous contigs\t[%]\tidentity [%]\tpossible joins\thomozygous size\t[%]\thomozygous contigs\t[%]\n")
    for fasta in o.fasta:
        out = gzip.open(fasta.name+".homozygous.fa.gz", "w")
        fasta2homozygous(out, fasta, o.identity, o.overlap, o.minLength, o.threads, o.verbose)
        out.close()

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
