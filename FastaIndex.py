#!/usr/bin/env python
desc="""FastA index (.fai) handler compatible with samtools faidx (http://www.htslib.org/doc/faidx.html)

TBA:
- get reverse complement if stop before start
- don't load entire sequence, just a slice; this is where Biopython is slow (https://www.biostars.org/p/58028/)

"""
epilog="""Author: l.p.pryszcz@gmail.com
Bratislava, 15/06/2016
"""

import os, sys
from datetime import datetime

class FastaIndex(object):
    """Facilitate Fasta index (.fai) operations compatible
    with samtools faidx (http://www.htslib.org/doc/faidx.html).
    """
    def __init__(self, handle, verbose=0):
        """ """
        self.verbose = verbose
        self.genomeSize = 0
        # guess handle
        if type(handle) is str and os.path.isfile(handle):
            self.handle = open(handle)
        elif type(handle) is file:
            if handle.name.endswith(('.gz','.bz')):
                raise Exception("Compressed files are currently not supported!")
            self.handle = handle
            
        self.fasta  = handle.name
        self.faidx  = self.fasta + ".fai"
        # create new index if no .fai or .fai younger than .fasta
        if not os.path.isfile(self.faidx) or \
           os.stat(self.fasta).st_mtime > os.stat(self.faidx).st_mtime: 
            self._generate_index()
        # otherwise load
        else:
            self._load_fai()
        # links
        self.get = self.get_fasta
            
    def _generate_index(self): 
        """Return fasta records"""
        if self.verbose:
            sys.stderr.write("Generating FastA index...\n")
        header = ""
        seq = []
        self.id2stats = {}
        with open(self.faidx, 'w') as out:
            for l in iter(self.handle.readline, ''): 
                if l.startswith(">"):
                    if header:
                        stats = self.get_stats(header, seq, offset)
                        seqid = self.get_id(header)
                        self.id2stats[seqid] = stats
                        out.write("%s\t%s\n"%(seqid, "\t".join(map(str, stats))))
                    header = l
                    offset = self.handle.tell() 
                    seq = []
                else:
                    seq.append(l)

            if header: 
                stats = self.get_stats(header, seq, offset)
                seqid = self.get_id(header)
                self.id2stats[seqid] = stats
                out.write("%s\t%s\n"%(seqid, "\t".join(map(str, stats))))

    def _load_fai(self):
        """Load stats from faidx file"""
        self.id2stats = {}
        for l in open(self.faidx):
            ldata = l[:-1].split('\t')
            if len(ldata)<9:
                return {}
            rid = ldata[0]
            stats = map(int, ldata[1:])
            self.id2stats[rid] = stats
            # update genomeSize
            self.genomeSize += stats[0]

    def __len__(self):
        """How many records are there?"""
        return len(self.id2stats)
            
    def __iter__(self):
        """Iterate over the keys."""
        for seqid in self.id2stats: #sorted(self.id2stats.keys(), key=lambda x: self.id2stats[x][1]):
            yield seqid
            
    def __getitem__(self, key, start=None, stop=None):
        """x.__getitem__(y) <==> x[y]"""
        if key not in self.id2stats:
            raise KeyError
        # get offset info
        size, offset, linebases, linebytes = self.id2stats[key][:4]
        # compute bytes to fetch
        linediff = linebytes - linebases
        bytesize = size / linebases * linebytes + size % linebases + linediff
        # load record
        self.handle.seek(offset)
        seq = self.handle.read(bytesize)
        seqid = key
        # get sequence slice
        if start and stop:
            # 1-base, inclusive end
            if start<1:
                start = 1
            seqid = "%s:%s-%s"%(key, start, stop)
            start -= 1
            seq = seq.replace('\n', '')[start:stop]
            seq = '\n'.join(seq[i:i+linebases] for i in range(0, len(seq), linebases))+'\n'
        record = ">%s\n%s"%(seqid, seq)
        return record

    def get_fasta(self, region="", contig="", start=None, stop=None):
        """Return FastA slice"""
        if region:
            if ':' in region:
                if '-' in region:
                    contig, startstop = region.split(':')
                    start, stop = map(int, startstop.split('-'))
            else:
                contig = region
        elif not contig:
            sys.stderr.write("Provide region or contig!\n")
            return            
        # get record
        record = self.__getitem__(contig, start, stop)
        return record

    def get_id(self, header):
        """Return seqid from header"""
        return header[1:].split()[0]

    def get_stats(self, header, seq, offset):
        """Return seq length, offset, linebases, linebyts and number of
        A, C, G and T in each sequence.
        Compatible with samtools faidx (http://www.htslib.org/doc/faidx.html).
        """
        errors = 0
        # get bases & bytes in line, ignoring last line
        if len(seq)>1:
            linebases = set(len(s.strip()) for s in seq[:-1])   
            linebytes = set(len(s) for s in seq[:-1])
            if len(linebases)>1:
                sys.stderr.write("[WARNING] Uneven line lengths in %s: %s\n"%(header, ",".join(map(str, linebases))))        
            linebases, linebytes = max(linebases), max(linebytes)
        else:
            linebases, linebytes = len(seq[0].strip()), len(seq[0])
        seq = "".join(s.strip() for s in seq)
        seqlen = len(seq)
        self.genomeSize += seqlen
        # count ACGT
        bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        for b in seq.upper():
            if b in bases:
                try:
                    bases[b] += 1
                except:
                    errors += 1
        return (seqlen, offset, linebases, linebytes, \
                bases['A'], bases['C'], bases['G'], bases['T'])    

def main():
    import argparse
    usage	 = "%(prog)s -i " #usage=usage, 
    parser	= argparse.ArgumentParser(description=desc, epilog=epilog, \
                                          formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='0.10a')	 
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")	
    parser.add_argument("-i", "--fasta", type=file, 
                        help="FASTA file(s)")
    parser.add_argument("-o", "--out",	 default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream	 [stdout]")
    parser.add_argument("-r", "--regions", nargs='*', default=[], 
                        help="contig or contig region to output [slices NOT IMPLEMENTED YET!]")

    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # init faidx
    faidx = FastaIndex(o.fasta, o.verbose)

    # report regions
    for region in o.regions:
        o.out.write(faidx.get_fasta(region))
	
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
    