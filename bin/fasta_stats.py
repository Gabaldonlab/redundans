#!/usr/bin/env python
desc="""Report FASTA statistics. Support gzipped files.

Statistics are stored as .fai formatted file (http://www.htslib.org/doc/faidx.html),
with 4 extended columns, storing counts for A, C, G & T for each sequence. 
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Mizerow/Bratislava, 26/08/2014
"""

import gzip, os, sys
from FastaIndex import FastaIndex
from datetime import datetime
    
def main():
    import argparse
    usage	 = "%(prog)s -i " #usage=usage, 
    parser	= argparse.ArgumentParser(description=desc, epilog=epilog, \
                                          formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='1.2')	 
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")	
    parser.add_argument("-i", "--fasta", nargs="+", 
                        help="FASTA file(s)")
    parser.add_argument("-o", "--out",	 default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream	 [stdout]")

    o = parser.parse_args()
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)        
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))
    
    #header
    header = '#fname\tcontigs\tbases\tGC [%]\tcontigs >1kb\tbases in contigs >1kb\tN50\tN90\tNs\tlongest\n'
    o.out.write(header)
    for fname in o.fasta:
        if not os.path.isfile(fname):
            sys.stderr.write("[WARNING] No such file: %s\n"%fname)
            continue
        faidx = FastaIndex(fname)
        o.out.write(faidx.stats())
	
if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!		\n")
    except OSError:
        sys.stderr.write("OS error({0}): {1}\n".format(e.errno, e.strerror))        
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)

