#!/usr/bin/env python
desc="""Estimate insert size by aligning subset of reads onto contigs. 

PREREQUISITIES:
- BWA
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow, 10/04/2015
"""

import math, os, sys, commands, subprocess
from datetime import datetime
from FastaIndex import FastaIndex

def flag2orientation(flag):
    """Return orientation of read pair: 
    - FF: 0
    - FR: 1
    - RF: 2
    - RR: 4
    """
    # FR/RF # one F and one R
    if flag&16 != flag&32:
        # FR # first and not R OR second and R
        if flag&64 and not flag&16 or flag&128 and not flag&16:
            return 1
        # RF
        else:
            return 2
    # RR # first and R OR second and R
    if flag&64 and flag&16 or flag&128 and not flag&16:
        return 3
    # FF
    else:
        return 0
        
def get_bwa_subprocess(fq1, fq2, fasta, threads, verbose):
    """Return bwa subprocess"""
    # generate index if missing
    if not os.path.isfile(fasta+".bwt"):
        cmd = "bwa index %s"%fasta
        #if verbose: sys.stderr.write(" %s\n"%cmd)
        bwtmessage = commands.getoutput(cmd)
    # start BWA alignment stream
    ## -S skips mate rescue (faster)
    cmd = ["bwa", "mem", "-S", "-t %s"%threads, fasta, fq1, fq2]
    #if verbose: sys.stderr.write(" %s\n"%" ".join(cmd))
    bwa = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return bwa

def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values. 

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.

    @return - the percentile of the values

    From http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

def median(N, key=lambda x:x):
    """median is 50th percentile."""
    return percentile(N, 0.5, key)

def mean(data):
    """Return the sample arithmetic mean of data.
    http://stackoverflow.com/a/27758326/632242
    """
    n = len(data)
    if n < 1:
        raise ValueError('mean requires at least one data point')
    return sum(data)/float(n) 

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    if n < 2:
        raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5

def get_isize_stats(fq1, fq2, fasta, mapqTh=10, threads=1, limit=1e5, verbose=0, 
                    percent=5, stdfracTh=0.66, maxcfracTh=0.9): 
    """Return estimated insert size median, mean, stdev and
    pairing orientation counts (FF, FR, RF, RR). 
    Ignore bottom and up percentile for insert size statistics. 
    """
    orientations = ('FF', 'FR', 'RF', 'RR')
    # read dumped info
    if os.path.isfile(fq2+".is.txt") and os.path.getmtime(fq2) < os.path.getmtime(fq2+".is.txt"):
        ldata = open(fq2+".is.txt").readline().split("\t")
        if len(ldata) == 8:
            readlen, ismedian, ismean, isstd = map(float, ldata[:4])
            pairs = map(int, ldata[-4:])
            # select major orientation
            reads, orientation = sorted(zip(pairs, orientations), reverse=1)[0]
            # skip insert size estimation only if satisfactory previous estimate
            if isstd / ismean < stdfracTh and reads > maxcfracTh * sum(pairs): 
                return int(readlen), ismedian, ismean, isstd, pairs, orientation
    # run bwa
    bwa = get_bwa_subprocess(fq1, fq2, fasta, threads, verbose)
    # parse alignments
    isizes = [[], [], [], []] #0, 0, 0, 0]
    readlen = 0
    #read from stdin
    for i, sam in enumerate(bwa.stdout, 1):
        if sam.startswith("@"):
            continue
        if verbose and not i%1000:
            sys.stderr.write(' %s %s \r'%(i, sum(map(len, isizes))))
        # read sam entry
        rname, flag, chrom, pos, mapq, cigar, mchrom, mpos, isize, seq = sam.split('\t')[:10]
        flag, pos, mapq, mpos, isize = map(int, (flag, pos, mapq, mpos, isize))
        # take only reads with good alg quality and one read per pair
        # ignore not primary and supplementary alignments
        if mapq < mapqTh or isize < 1 or flag&256 or flag&2048: # or not flag&2:
            continue
        #store isize
        isizes[flag2orientation(flag)].append(isize)
        readlen += len(seq)
        #stop if limit reached
        if sum(map(len, isizes)) >= limit:
            break
    # terminate subprocess
    bwa.terminate()
    # catch cases with very few reads aligned
    pairs = map(len, isizes)
    readlen = int(round(1.*readlen/sum(pairs)))
    if sum(pairs) < 100:
        return 0, 0, 0, 0, [], ''
    # select major orientation - replace isizes by major isizes
    isizes, orientation = sorted(zip(isizes, orientations), key=lambda x: len(x[0]), reverse=1)[0]    
    # get frac of total reads
    maxcfrac = 1.0 * len(isizes) / sum(pairs)
    if maxcfrac < maxcfracTh:
        info = "[WARNING] Poor quality: Major orientation (%s) represent %s%s of pairs in %s - %s: %s\n"
        sys.stderr.write(info%(orientation, 100*maxcfrac, '%', fq1, fq2, str(pairs)))
    #get rid of 5 percentile from both sides
    isizes.sort()
    maxins = percentile(isizes, 0.01*(100-percent))
    minins = percentile(isizes, 0.01*percent)
    isizes = [x for x in isizes if x>minins and x<maxins]
    # get stats
    ismedian, ismean, isstd = median(isizes), mean(isizes), pstdev(isizes)
    # save info
    try:
        with open(fq2+".is.txt", "w") as out:
            out.write("%s\t%s\t%s\t%s\t%s\n"%(readlen, ismedian, ismean, isstd, "\t".join(map(str, pairs))))
    except:
        sys.stderr.write("[WARNING] Couldn't write library statistics to %s\n"%(fq2+".is.txt",))
    return readlen, ismedian, ismean, isstd, pairs, orientation

def prepare_genome(fasta, genomeFrac=0.05):
    """Prepare new fasta file that will contain
    genomefrac of the original reference in largest contigs.
    """
    faidx = FastaIndex(fasta)
    newfasta = "%s.%s"%(fasta, genomeFrac)
    with open(newfasta, "w") as out:
        for c in faidx.sort(genomeFrac=genomeFrac):
            out.write(faidx[c])
    return newfasta
    
def fastq2insert_size(out, fastq, fasta, mapq, threads, limit, verbose, log=sys.stderr, 
                      genomeFrac=0.05, stdfracTh=0.66, maxcfracTh=0.9):
    """Report insert size statistics and return all information."""
    # prepare genome
    fasta = prepare_genome(fasta, genomeFrac)
    # 
    header  = "Insert size statistics\t\t\t\tMates orientation stats\n"
    header += "FastQ files\tread length\tmedian\tmean\tstdev\tFF\tFR\tRF\tRR\n"
    if verbose:
        out.write(header)
    line = "%s %s\t%i\t%i\t%.2f\t%.2f\t%s\n"
    data = []
    for fq1, fq2 in zip(fastq[0::2], fastq[1::2]):
        # get IS stats
        isstats = get_isize_stats(fq1, fq2, fasta, mapq, threads, limit, verbose, stdfracTh, maxcfracTh)
        readlen, ismedian, ismean, isstd, pairs, orientation = isstats
        if not sum(pairs):
            log.write("[WARNING] No alignments for %s - %s!\n"%(fq1, fq2))
            continue
        # report
        if verbose:
            out.write(line%(fq1, fq2, readlen, ismedian, ismean, isstd, "\t".join(map(str, pairs))))
        # store data
        data.append((fq1, fq2, readlen, ismedian, ismean, isstd, pairs, orientation))
    return data
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='1.0a')   
    parser.add_argument("-i", "--fastq", nargs="+", 
                        help="FASTQ PE/MP files")
    parser.add_argument("-f", "--fasta", #type=file, 
                        help="reference assembly FASTA file")
    parser.add_argument("-o", "--output",  default=sys.stdout, type=argparse.FileType("w"), 
                        help="output stream [stdout]")
    parser.add_argument("-l", "--limit",  default=10000, type=int, 
                        help="align l reads [%(default)s]")
    parser.add_argument("-q", "--mapq",    default=10, type=int, 
                        help="min mapping quality for variants [%(default)s]")
    parser.add_argument("-t", "--threads", default=1, type=int, 
                        help="max threads to run [%(default)s]")
    parser.add_argument("-g", "--genomefrac",  default=0.05, type=float, 
                        help="use only this fraction of genome to speed up estimation [%(default)s]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
        
    fastq2insert_size(o.output, o.fastq, o.fasta, o.mapq, o.threads, \
                      o.limit, o.verbose, o.genomefrac)

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
