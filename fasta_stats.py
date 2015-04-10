#!/usr/bin/env python
desc="""Report FASTA statistics. Support gzipped files.

TBD:
Code need some polishing...
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow, 26/08/2014
"""

import os, sys
from Bio import SeqIO
import gzip, os, sys
from datetime import datetime
from Bio import SeqIO

def fasta_stats(f):
    """Report fasta statistics."""
    fn=f.name
    statsFn=fn+'.stats'
    #return content of .stats file if exists and younger than fasta
    if os.path.isfile(statsFn) and os.stat(fn).st_mtime < os.stat(statsFn).st_mtime:
        line=open(statsFn).readline()
        line=line.strip()
        line='%s\t%s\n' % (fn, '\t'.join( line.split('\t')[1:] ))
        return line
    #if not, generate that file
    lengths=[]; lengths1000=[]
    contigs=contigs1000=baseErrors=0
    #count bases frequencies
    bases={'A':0,'C':0,'G':0,'T':0,'N':0}
    for r in SeqIO.parse(f,'fasta' ):
        contigs+=1
        seq=str(r.seq)
        seq=seq.upper()
        lengths.append( len(seq) )
        if len(seq)>1000: 
            contigs1000+=1
            lengths1000.append( len(seq) )
        for base in seq:
            try:
                bases[base]+=1
            except:
                baseErrors+=1
    if not lengths:
        return fn+'\tError: No sequences!\n'
    #calculate GC
    if bases['A']+bases['T']:
        GC=(bases['G']+bases['C'])*100.0/(bases['A']+bases['C']+bases['G']+bases['T'])
    else:
        GC=0

    #N50 & N90
    size=sum(lengths) #sum(bases.itervalues())
    lengthSum=0
    n50=n90=0
    lengths.sort( reverse=True )
    for l in lengths:
        lengthSum += l
        if not n50 and lengthSum>=0.5*size:
            n50=l
        if lengthSum>=0.9*size:
            n90=l 
            break
            
    #print output
    line='%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\n' % ( fn,contigs,size,GC,contigs1000,sum(lengths1000),n50,n90,bases['N'],lengths[0] )
    try:
        out=open(statsFn,'wb'); out.write( line ); out.close()
    except IOError as e:
        sys.stderr.write("%s\n" % e)
    return line

def main():
    import argparse
    usage	 = "%(prog)s -i " #usage=usage, 
    parser	= argparse.ArgumentParser(description=desc, epilog=epilog, \
                                                                            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='1.1')	 
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")	
    parser.add_argument("-i", "--fasta", nargs="+", type=file, 
                        help="FASTA file(s)")
    parser.add_argument("-o", "--out",	 default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream	 [stdout]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    #header
    o.out.write('#fname\tcontigs\tbases\tGC [%]\tcontigs >1kb\tbases in contigs >1kb\tN50\tN90\tNs\tlongest\n')
    for f in o.fasta:
        if not os.path.isfile(f.name):
            continue
        if f.name.endswith('.gz'):
            f=gzip.open(f.name)
        o.out.write(fasta_stats(f))
	
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!		\n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)

