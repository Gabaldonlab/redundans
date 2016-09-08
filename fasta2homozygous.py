#!/usr/bin/env python
desc="""Align genome onto itself (LAST) and keep only the longest
from heterozygous (redundant) contigs/scaffolds.

TO ADD:
- scaffold extension based on overlapping matches (overlapping already recognised)
- reporting of haplotypes
- recognise heterozygous contigs with translocations
- replace gzip with bgzip and indexing on the fly
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 26/08/2014
"""

import gzip, math, os, sys, subprocess
from datetime import datetime
from FastaIndex import FastaIndex

def run_last(fasta, identity, threads, verbose):
    """Start LAST with multi-threads. """
    if verbose:
        sys.stderr.write(" Running LAST...\n")
    # build db
    if not os.path.isfile(fasta+".suf"):
        os.system("lastdb %s %s" % (fasta, fasta))
    # run LAST
    args = ["lastal", "-T", "1", "-f", "TAB", "-P", str(threads), fasta, fasta]
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=sys.stderr)        
    return proc
    
def fasta2hits(fasta, threads, identityTh, overlapTh, joinOverlap, endTrimming, verbose):
    """Return valid hits. """
    overlapping = []
    hits = []
    added = set()
    # execute last
    last = run_last(fasta.name, identityTh, threads, verbose)
    for l in last.stdout: 
        if l.startswith('#'):
            continue
        # unpack
        (score, q, qstart, qalg, qstrand, qsize, t, tstart, talg, tstrand, tsize, blocks) = l.split()[:12]
        (score, qstart, qalg, qsize, tstart, talg, tsize) = map(int, (score, qstart, qalg, qsize, tstart, talg, tsize))
        # skip reverse matches
        if q==t or tsize<qsize or (t,q) in added: 
            continue
        added.add((q,t))
        #get score, identity & overlap # LASTal is using +1/-1 for match/mismatch, while I need +1/0
        identity = 1.0 * (score+(qalg-score)/2) / qalg
        overlap  = 1.0 * qalg / qsize
        #filter by identity and overlap
        if identity < identityTh or overlap < overlapTh:
            continue
        # store
        qend, tend = qstart + qalg, tstart + talg
        #data =(t, tsize, tstart, tend, q, qsize, qstart, qend, tstrand, identity, overlap, score)
        data = (score, t, q, qend-qstart, identity)
        hits.append(data)
        
    return hits, overlapping
    
def hits2skip(hits, faidx, verbose):
    """Return contigs to skip."""
    identities = algLengths = 0
    contig2skip = {}
    for c in faidx: 
        contig2skip[c] = 0
    # this may be slow! sorting by score
    for i, data in enumerate(sorted(hits, key=lambda x: x[0], reverse=1), 1):
    #for i, data in enumerate(hits, 1):
        #(t, tsize, tstart, tend, q, qsize, qstart, qend, strand, identity, overlap, score) = data
        (score, t, q, algLen, identity) = data
        if t not in contig2skip:
            sys.stderr.write(' [ERROR] `%s` (%s) not in contigs!\n'%(t, str(hits[i-1])))
            continue
        # skip alignments of contigs already removed
        if contig2skip[q]:
            # inform about matching already removed contig
            if verbose:
                info = " [WARNING]: Match to already removed conting: %s %s\n"
                sys.stderr.write(info%(q, str(hits[i-1])))
            continue
        # store
        contig2skip[q] += 1
        # update identities and lengths
        identities += identity*algLen
        algLengths += algLen
    # calculated identity
    identity = 0
    if algLengths:
        identity = 100.0*identities/algLengths
    return contig2skip, identity

'''def get_coverage(faidx, fasta, libraries, limit, verbose):
    """Align subset of reads and calculate per contig coverage"""
    # init c2cov make it python 2.6 compatible 
    c2cov = {} #c: 0 for c in faidx}
    covTh = 0
    
    return c2cov, covTh
'''
def fasta2homozygous(out, fasta, identity, overlap, minLength, \
                     libraries, limit, threads=1, joinOverlap=200, endTrimming=0,
                     verbose=0, log=sys.stderr):
    """Parse alignments and report homozygous contigs"""
    #create/load fasta index
    if verbose:
        log.write("Indexing fasta...\n")
    faidx = FastaIndex(fasta)
    genomeSize = faidx.genomeSize

    '''# depth-of-coverage info
    c2cov, covTh = None, None
    if libraries:
        c2cov, covTh = get_coverage(faidx, fasta.name, libraries, limit, verbose)
    '''
    if verbose:
        log.write("Parsing alignments...\n")
    #filter alignments
    hits, overlapping = fasta2hits(fasta, threads, identity, overlap, joinOverlap, endTrimming, verbose)

    #remove redundant
    ## maybe store info about removed also
    contig2skip, identity = hits2skip(hits, faidx, verbose)
    
    #report homozygous fasta
    nsize, k, skipped, ssize, merged = merge_fasta(out, faidx, contig2skip, \
                                                   overlapping, minLength, verbose)
    
    #summary    
    info = "%s\t%s\t%s\t%s\t%.2f\t%s\t%.2f\t%.3f\t%s\t%s\t%.2f\t%s\t%.2f\n"
    log.write(info%(fasta.name, genomeSize, len(faidx), ssize, 100.0*ssize/genomeSize, \
                    skipped, 100.0*skipped/len(faidx), identity, len(merged), \
                    nsize, 100.0*nsize/genomeSize, k, 100.0*k/len(faidx)))

    return genomeSize, len(faidx), ssize, skipped, identity

def get_name_abbrev(size, s, e):
    """Return s if s < size/2, otherwise return e."""
    if s + (e - s)/2.0 < size / 2.0:
       return "s"
    return "e"
    
def merge_fasta(out, faidx, contig2skip, overlapping, minLength, verbose):
    """Merged overlapping and report homozygous genome."""
    #merge
    joins = {}
    nsize = 0
    #ignore extensions with skipped contigs
    for data in filter(lambda x: not contig2skip[x[0]] and not contig2skip[x[4]], \
                       overlapping):
        (q, qsize, qstart, qend, t, tsize, tstart, tend, strand, identity, overlap, bitscore, evalue) = data
        qname = "%s%s"%(q, get_name_abbrev(qsize, qstart, qend))
        tname = "%s%s"%(t, get_name_abbrev(tsize, tstart, tend))
        #check if already present
        if qname in joins:
            if bitscore < joins[qname][1][-2]:
                continue
        if tname in joins:
            if bitscore < joins[tname][1][-2]:
                continue
        #rm joins
        if qname in joins:
            joins.pop(joins[qname][0])
        if tname in joins:
            joins.pop(joins[tname][0])
        #add
        joins[tname] = (qname, data)
        joins[qname] = (tname, data)
        
    #merging
    merged = {}
    
    #report not skipper, nor joined
    k = skipped = ssize = 0
    for i, c in enumerate(faidx, 1):
        # don't report skipped & merged
        if contig2skip[c] or faidx.id2stats[c][0]<minLength:
            skipped += 1
            ssize   += len(faidx[c])
            continue
        elif c in merged:
            continue
        out.write(faidx[c])
        k += 1
        nsize += faidx.id2stats[c][0] 
        
    #return nsize, k, skipped, ssize, merged
    return nsize, k, skipped, ssize, joins
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.01d')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "-f", "--fasta", nargs="+", type=file, 
                        help="FASTA file(s)")
    parser.add_argument("-t", "--threads", default=4, type=int, 
                        help="max threads to run [%(default)s]")
    #parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
    #                    help="output stream   [stdout]")
    parser.add_argument("--identity",    default=0.51, type=float, 
                        help="min. identity   [%(default)s]")
    parser.add_argument("--overlap",     default=0.66, type=float, 
                        help="min. overlap    [%(default)s]")
    parser.add_argument("--joinOverlap", default=200, type=int, 
                        help="min. end overlap to join two contigs [%(default)s]")
    parser.add_argument("--endTrimming", default=33, type=int, 
                        help="max. end trim on contig join [%(default)s]")
    parser.add_argument("--minLength",   default=200, type=int, 
                        help="min. contig length [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # allow depth-of-coverage
    libraries, limit = [], 0
        
    #process fasta
    sys.stderr.write("Homozygous assembly/ies will be written with input name + '.homozygous.fa.gz'\n")
    sys.stderr.write("#file name\tgenome size\tcontigs\theterozygous size\t[%]\theterozygous contigs\t[%]\tidentity [%]\tpossible joins\thomozygous size\t[%]\thomozygous contigs\t[%]\n")
    for fasta in o.fasta:
        out = gzip.open(fasta.name+".homozygous.fa.gz", "w")
        fasta2homozygous(out, fasta, o.identity, o.overlap, o.minLength, \
                         libraries, limit, o.threads, o.joinOverlap, o.endTrimming, o.verbose)
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
