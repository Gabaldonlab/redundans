#!/usr/bin/env python
desc="""Report FASTA statistics. Support gzipped files.

Statistics are stored as .fai formatted file (http://www.htslib.org/doc/faidx.html),
with 4 extended columns, storing counts for A, C, G & T for each sequence. 
"""
epilog="""Author: l.p.pryszcz@gmail.com
Mizerow/Bratislava, 26/08/2014
"""

import gzip, os, sys
from FastaIndex import FastaIndex
from datetime import datetime
    
def fasta_stats(handle):
    """Report fasta statistics."""
    # load id2stats
    faidx = FastaIndex(handle)
    id2stats = faidx.id2stats
    # report stats
    contigs = len(id2stats)
    lengths = [stats[0] for stats in id2stats.itervalues()]
    size = sum(lengths)
    lengths1000 = [l for l in lengths if l>=1000]
    contigs1000 = len(lengths1000)
    A, C, G, T = map(sum, zip(*[stats[-4:] for stats in id2stats.itervalues()]))
    GC = 100.0*(G+C)/(A+T+G+C)
    nonACGT = size - A - C - G - T
    # N50 & N90
    n50 = n90 = lengthSum = 0
    lengths.sort( reverse=True )
    for l in lengths:
        lengthSum += l
        if not n50 and lengthSum >= 0.5*size:
            n50 = l
        if lengthSum >= 0.9*size:
            n90 = l 
            break
            
    _line = '%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\n'
    line = _line % (handle.name, contigs, size, GC, contigs1000, sum(lengths1000), n50, n90, nonACGT, lengths[0])
    return line
    
def main():
    import argparse
    usage	 = "%(prog)s -i " #usage=usage, 
    parser	= argparse.ArgumentParser(description=desc, epilog=epilog, \
                                          formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='1.2')	 
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
    header = '#fname\tcontigs\tbases\tGC [%]\tcontigs >1kb\tbases in contigs >1kb\tN50\tN90\tNs\tlongest\n'
    o.out.write(header)
    for f in o.fasta:
        if not os.path.isfile(f.name):
            continue
        if f.name.endswith('.gz'):
            f = gzip.open(f.name)
        o.out.write(fasta_stats(f))#, header))
	
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!		\n")
    except IOError as e:
        sys.stderr.write("I/O error({0}): {1}\n".format(e.errno, e.strerror))
    #[Errno 95] Operation not supported
    except OSError:
        sys.stderr.write("OS error({0}): {1}\n".format(e.errno, e.strerror))        
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)

