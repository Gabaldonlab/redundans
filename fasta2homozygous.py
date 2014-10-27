#!/usr/bin/env python
desc="""Align genome onto itself (BLAT) and remove heterozygous (redundant) scaffolds.

TO ADD:
- scaffold extension based on overlapping matches (overlapping already recognised)
- reporting of haplotypes
- recognise heterozygous contigs with translocations
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 26/08/2014
"""

import gzip, math, os, sys
from datetime import datetime
from Bio import SeqIO

def blat(fasta, identity, verbose):
    """Start BLAT"""
    #prepare BLAT command
    identity = int(100*identity)
    args = ["-ooc=%s.11.ooc"%fasta, "-dots=1000", "-noHead", "-extendThroughN", \
            "-minScore=%s"%identity, "-minIdentity=%s"%identity]
    cmd = "blat %s %s %s %s.psl"%(" ".join(args), fasta, fasta, fasta)
    if not verbose:
        cmd += " > /dev/null"
    else:
        sys.stderr.write(cmd+'\n')
    #generate overepresented 11mers if not exists
    if not os.path.isfile(fasta+".11.ooc"):
        os.system(cmd.replace("-ooc=", "-makeOoc="))
    #run BLAT
    os.system(cmd)
    #sort and take into account only larger vs smaller
    cmd2 = "awk '$10!=$14 && $11>=$15' %s.psl | sort -k11nr,11 -k12n,12 -k13nr,13 | gzip > %s.psl.gz"%(fasta, fasta)
    if verbose:
        sys.stderr.write(cmd2+'\n')
    os.system(cmd2)
    #clean-up
    os.system("rm %s.psl %s.11.ooc"%(fasta, fasta))

def get_ranges(starts, sizes, offset=1):
    """Return str representation of alg ranges"""
    ranges = []
    for start, size in zip(starts.split(',')[:-1], sizes.split(',')[:-1]):
        start, size = int(start), int(size)
        start += offset
        end = start + size - offset
        coords = "%s-%s"%(start, end)
        ranges.append(coords)
    return " ".join(ranges)

def psl2hits(psl, identityTh, overlapTh, joinOverlap, endTrimming, \
             dblength=0, Lambda=0.318, K=0.13):
    """Return valid hits.
    PSL has to be sorted by q size: sort -k11nr,11 -k12n,12 -k13nr,13
    """
    overlapping = []
    hits = []
    added = set()
    for l in gzip.open(psl):
        if not l.strip() or not l.split()[0].isdigit():
            continue
        ##BLAT PSL without header
        (matches, mismatches, repm, Ns, Qgapc, Qgaps, Tgapc, Tgaps, strand, \
         q, qsize, qstart, qend, t, tsize, tstart, tend, blocks, bsizes, \
         qstarts, tstarts) = l.split('\t')
        #skip reverse matches - pairs tracking
        if q==t or (t,q) in added:
            continue
        added.add((q,t))
        #unpack batch
        matches, mismatches = int(matches), int(mismatches)
        Tgapc, Tgaps = int(Tgapc), int(Tgaps)
        qstart, qend = int(qstart), int(qend)
        qstart, qend, qsize = int(qstart), int(qend), int(qsize)
        tstart, tend, tsize = int(tstart), int(tend), int(tsize)
        #get score, identity & overlap
        score    = matches * 5 + mismatches * -3 + Tgapc * -4 + Tgaps * -1
        alglen   = int(tend) - int(tstart)
        identity = 1.0 * matches / alglen
        overlap  = 1.0 * alglen / tsize
        #bitscore & evalue
        bitscore = (Lambda*score-math.log(K))/math.log(2)
        pvalue = evalue = 0
        if dblength:
            pvalue = 2**-bitscore
            evalue = len(qseqs[0]) * dblength * pvalue
        data = (q, qsize, qstart, qend, t, tsize, tstart, tend, strand, \
                     identity, overlap, bitscore, evalue)
        #filter by identity and overlap
        if identity < identityTh:
            continue
        #capture overlapping #overlap
        elif overlap < overlapTh:
            if matches >= joinOverlap and qsize > tsize:
                if   strand == "+":
                    if qstart < endTrimming and tsize-tend < endTrimming \
                       or tstart < endTrimming and qsize-qend < endTrimming:
                        overlapping.append(data)
                elif strand == "-":
                    if qstart < endTrimming and tstart < endTrimming \
                       or tsize-tend < endTrimming and qsize-qend < endTrimming:
                        overlapping.append(data)
        #store
        else:
            hits.append(data)
        
    return hits, overlapping
    
def hits2skip(hits, faidx, verbose):
    """Return contigs to skip."""
    identities, algLengths = [], []
    contig2skip = {}
    for c in faidx.keys():
        contig2skip[c] = 0
    for i, data in enumerate(hits, 1):
        (q, qsize, qstart, qend, t, tsize, tstart, tend, strand, identity, overlap, bitscore, evalue) = data
        if q not in contig2skip:
            sys.stderr.write(' [ERROR] `%s` (%s) not in contigs!\n'%(q, str(hits[i-1])))
            continue
        #inform about matching already removed contig
        if verbose and contig2skip[q]:
            info = " [WARNING]: Match to already removed conting: %s %s\n"
            sys.stderr.write(info%(q, str(hits[i-1])))
        #store
        contig2skip[t] += 1
        #update identities and lengths
        algLen = tend-tstart
        identities.append(identity*algLen)
        algLengths.append(algLen)
    #calculated divergence
    if sum(algLengths):
        identity = 100.0*sum(identities)/sum(algLengths)
    else:
        identity = 0
    return contig2skip, identity

def fasta2homozygous(out, fasta, identity, overlap, joinOverlap, endTrimming, verbose):
    """Parse alignments and report homozygous contigs"""
    #create/load fasta index
    if verbose:
        sys.stderr.write("Indexing fasta...\n")
    faidx = SeqIO.index_db(fasta.name+".db3", fasta.name, "fasta")
    genomeSize = sum(len(faidx[c]) for c in faidx) 
    
    #run blat
    psl = fasta.name + ".psl.gz"
    if not os.path.isfile(psl):
        if verbose:
            sys.stderr.write("Running BLAT...\n")
        blat(fasta.name, identity, verbose)
    
    if verbose:
        sys.stderr.write("Parsing alignments...\n")
    #filter alignments
    hits, overlapping = psl2hits(psl, identity, overlap, joinOverlap, endTrimming)

    #remove redundant
    contig2skip, identity = hits2skip(hits, faidx, verbose)
    #print "\n".join("\t".join(map(str, x)) for x in overlapping[:100]); return
    
    #report homozygous fasta
    nsize, k, skipped, ssize, merged = merge_fasta(out, faidx, contig2skip, \
                                                   overlapping, verbose)
    
    #summary    
    info = "%s\t%s\t%s\t%s\t%.2f\t%s\t%.2f\t%.3f\t%s\t%s\t%.2f\t%s\t%.2f\n"
    sys.stderr.write(info%(fasta.name, genomeSize, len(faidx), ssize, 100.0*ssize/genomeSize, \
                           skipped, 100.0*skipped/len(faidx), identity, len(merged), \
                           nsize, 100.0*nsize/genomeSize, k, 100.0*k/len(faidx)))

def get_name_abbrev(size, s, e):
    """Return s if s < size/2, otherwise return e."""
    if s + (e - s)/2.0 < size / 2.0:
       return "s"
    return "e"
    
def merge_fasta(out, faidx, contig2skip, overlapping, verbose):
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

    """
    for k, v in sorted(joins.iteritems()):
        print k, v[0]
    print len(joins)
    sys.exit()
    #"""
        
    #merging
    merged = {}
    
    
    #report not skipper, nor joined
    k = skipped = ssize = 0
    for i, c in enumerate(faidx, 1):
        #don't report skipped & merged
        if contig2skip[c]:
            skipped += 1
            ssize   += len(faidx[c])
            continue
        elif c in merged:
            continue
        out.write(faidx[c].format('fasta'))
        k += 1
        nsize += len(faidx[c])
        
    #return nsize, k, skipped, ssize, merged
    return nsize, k, skipped, ssize, joins
        
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument('--version', action='version', version='1.0b')   
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")    
    parser.add_argument("-i", "-f", "--fasta", nargs="+", type=file, 
                        help="FASTA file(s)")
    #parser.add_argument("-o", "--output",    default=sys.stdout, type=argparse.FileType('w'), 
    #                    help="output stream   [stdout]")
    parser.add_argument("--identity",    default=0.8, type=float, 
                        help="min. identity   [%(default)s]")
    parser.add_argument("--overlap",     default=0.75, type=float, 
                        help="min. overlap    [%(default)s]")
    parser.add_argument("--joinOverlap", default=200, type=int, 
                        help="min. end overlap to join two contigs [%(default)s]")
    parser.add_argument("--endTrimming", default=33, type=int, 
                        help="max. end trim on contig join [%(default)s]")
    #parser.add_argument("-p", "--ploidy",    default=2, type=int, 
    #                    help="ploidy          [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #process fasta
    for fasta in o.fasta:
        out = gzip.open(fasta.name+".homozygous.fa.gz", "w")
        fasta2homozygous(out, fasta, o.identity, o.overlap, \
                         o.joinOverlap, o.endTrimming, o.verbose)
        out.close()

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
