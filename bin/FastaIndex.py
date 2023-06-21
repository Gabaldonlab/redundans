#!/usr/bin/env python3
desc="""FastA index (.fai) handler compatible with samtools faidx (http://www.htslib.org/doc/faidx.html).
.fai is extended with 4 columns storing counts for A, C, G & T for each sequence.
"""
epilog="""Author: l.p.pryszcz+git@gmail.com
Bratislava, 15/06/2016

Updated to Python3 by Diego Fuentes Palacios
Barcelona 07/07/2022
"""

import os, sys
from io import TextIOWrapper
from datetime import datetime
from functools import reduce

def symlink(file1, file2):
    """Create symbolic link taking care of real path."""
    if not os.path.isfile(file2):
        # check if need for absolute path
        file1abs = os.path.join(os.path.realpath(os.path.curdir), file1)
        if os.path.isfile(file1abs):
            os.symlink(file1abs, file2)
        # otherwise create symbolic link without full path
        else:
            os.symlink(file1, file2)

class FastaIndex(object):
    """Facilitate Fasta index (.fai) operations compatible
    with samtools faidx (http://www.htslib.org/doc/faidx.html).
    """
    def __init__(self, handle, verbose=0, log=sys.stderr):
        """ """
        ext = ".fai"
        self.verbose = verbose
        self.log = log.write
        self.genomeSize = 0
        self.whitespaces_in_headers = False
        # guess handle
        if type(handle) is str and os.path.isfile(handle): 
            handle = open(handle)
        if type(handle) is TextIOWrapper:
            if handle.name.endswith(('.gz','.bz')):
                raise Exception("Compressed files are currently not supported!")
            self.handle = handle
        else:
            sys.stderr.write("[ERROR] Couldn't guess handle for %s\n"%str(handle))
            sys.exit(1)
            
        self.fasta  = self.handle.name
        self.faidx  = self.fasta + ext
        # check if fasta is symlink
        if not os.path.isfile(self.faidx) and os.path.islink(self.fasta):
            _fasta = os.path.realpath(self.fasta)
            _faidx = _fasta+ext
            # symlink faidx if faidx exists and linked fasta is older than its faidx
            if os.path.isfile(_faidx) and os.stat(_fasta).st_mtime < os.stat(_faidx).st_mtime:
                symlink(_faidx, self.faidx)
        # create new index if no .fai, .fai loading failed or .fai younger than .fasta
        if not os.path.isfile(self.faidx) or not self._load_fai() or \
           os.stat(self.fasta).st_mtime > os.stat(self.faidx).st_mtime:
            self._generate_index()
        # links
        self.get = self.get_fasta
        # init storage
        self.base2rc= {"A": "T", "T": "A", "C": "G", "G": "C",
                       "a": "t", "t": "a", "c": "g", "g": "c",
                       "N": "N", "n": "n"}
        # Calculate total counts
        counts = []
        #Select first index, as it includes contig length
        base_stats = [self.id2stats[stats][0] for stats in self.id2stats]
        for count in base_stats[0:]:
            try:
                _, bp = count.split("_", 1)
            except:
                bp = count
            counts.append(bp)
        #Apply reduce to a lambda function
        self.totalcounts = reduce(lambda x, y: x+y, counts)
     
        #Store aswell all bases counted in a map function as well as Ns
        self.basecounts = list(map(sum, zip(*[stats[-4:] for stats in self.id2stats.values()])))
        
        #Ns
        self.Ns = int(self.genomeSize) - sum(self.basecounts)

    def __process_seqentry(self, out, header, seq, offset, pi):
        """Write stats to file and report any issues"""
        if header:
            # get seqid and sequence stats
            seqid = self.get_id(header)
            # catch empty headers
            if not seqid:
                self.log("[WARNING] No header at line: %s\n"%", ".join(map(str, (pi,seqid,header))))
                return
            stats = self.get_stats(header, seq, offset)
            # warn about empty sequences
            if not stats[0]:
                self.log("[WARNING] No sequence for: %s at line: %s\n"%(seqid, pi))
            # catch duplicates
            if seqid in self.id2stats:
                self.log("[WARNING] Duplicated sequence ID: %s at line: %s\n"%(seqid, pi))
            self.id2stats[seqid] = stats
            out.write("%s\t%s\n"%(seqid, "\t".join(map(str, stats))))
            
    def _generate_index(self): 
        """Return fasta records"""
        if self.verbose:
            self.log("Generating FastA index...\n")
        header, seq = "", []
        offset = pi = 0
        self.id2stats = {}
        with open(self.faidx, 'w') as out:
            for i, l in enumerate(iter(self.handle.readline, ''), 1): 
                if l.startswith(">"):
                    self.__process_seqentry(out, header, seq, offset, pi)
                    # mark that there is whitespace in headers
                    if len(l[:-1].split())>1:
                        self.whitespaces_in_headers = True
                    header = l
                    offset = self.handle.tell() 
                    seq = []
                    pi = i
                else:
                    seq.append(l)
            # process last entry
            self.__process_seqentry(out, header, seq, offset, pi)

    def _load_fai(self):
        """Load stats from faidx file.
        Return False if .fai is wrongly formatted.
        """
        self.id2stats = {}
        for l in open(self.faidx):
            ldata = l[:-1].split('\t')
            if len(ldata)<9:
                self.whitespaces_in_headers = False
                return 
            rid = ldata[0]
            stats = list(map(int, ldata[1:]))
            self.id2stats[rid] = stats
            # update genomeSize
            self.genomeSize += stats[0]
            if len(rid.split())>1:
                self.whitespaces_in_headers = True
        return True

    def __len__(self):
        """How many records are there?"""
        return len(self.id2stats)
            
    def __iter__(self):
        """Iterate over the keys."""
        for seqid in self.id2stats: 
            yield seqid

    def __getitem__(self, key, start=None, stop=None, name=None, seqonly=False):
        """x.__getitem__(y) <==> x[y]"""
        if key not in self.id2stats:
            raise KeyError
        # get offset info
        size, offset, linebases, linebytes = self.id2stats[key][:4]
        # compute bytes to fetch
        linediff = linebytes - linebases
        seqid = key
        # get sequence slice
        if start and stop:
            reverse_complement = 0
            if start<1:
                start = 1
            seqid = "%s:%s-%s"%(key, start, stop)
            if start>stop:
                reverse_complement = 1
                start, stop = stop, start
            if stop > size:
                stop = size
            # 1-base, inclusive end
            start -= 1
            # get bytesize and update offset
            offset += int(int(start / linebases) * linebytes) + start % linebases
            realsize = stop-start
            bytesize = int(int(realsize / linebases) * linebytes) + realsize % linebases
            # read sequence slice
            self.handle.seek(offset)
            seq = self.handle.read(bytesize).replace('\n', '')
            if reverse_complement:
                seq = self.get_reverse_complement(seq)
            if seqonly:
                return seq
            # format lines
            seq = '\n'.join(seq[i:i+linebases] for i in range(0, len(seq), linebases))+'\n'
        # load whole sequence record
        else:
            # get bytesize
            bytesize = int(int(size / linebases) * linebytes) + size % linebases
            ## add line diff for last line only for multiline fasta if last line is not complete
            if size / linebytes and size % linebases:
                bytesize += linediff 
            # read entire sequence
            self.handle.seek(offset)
            seq = self.handle.read(int(bytesize))
            if seqonly:
                return "".join(seq.split('\n'))
        # update name
        if not name:
            name = seqid
        record = ">%s\n%s"%(name, seq)
        return record

    def get_reverse_complement(self, seq):
        """Return reverse complement"""
        rc = []
        for seqsegment in seq.split():
            for b in seqsegment:
                if b in self.base2rc:
                    rc.append(self.base2rc[b])
                else:
                    rc.append(b)
        return "".join(reversed(rc))

    def get_sequence(self, contig, reverse=False):
        """Return sequence of given contig"""
        seq = self.__getitem__(contig, seqonly=True)
        if reverse:
            return self.get_reverse_complement(seq)
        return seq
        
    def get_fasta(self, region="", contig="", start=None, stop=None):
        """Return FastA slice"""
        if region:
            if ':' in region:
                #if '-' in region:
                try:
                    contig, startstop = region.split(':')
                    start, stop = map(int, startstop.split('-'))
                except Exception:
                    raise Exception("Malformed region definition: %s, while expected contig:start-stop"%region)
            else:
                contig = region
        elif not contig:
            self.log("Provide region or contig!\n")
            return            
        # get record
        record = self.__getitem__(contig, start, stop)
        return record

    def get_id(self, header):
        """Return seqid from header"""
        # catch empty headers
        if len(header.strip())<2:
            return
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
                self.log("[WARNING] Uneven line lengths in %s: %s\n"%(header, ",".join(map(str, linebases))))        
            linebases, linebytes = max(linebases), max(linebytes)
        elif len(seq)==1:
            linebases, linebytes = len(seq[0].strip()), len(seq[0])
        # handle empty sequences https://github.com/lpryszcz/redundans/issues/13
        else:
            linebases, linebytes = 60, 61 #len(seq[0].strip()), len(seq[0])
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

    def sort(self, reverse=1, minLength=0, genomeFrac=0):
        """Return list of contigs sorted by descending size (reverse=1).

        The list of returned contigs can be limited by: 
        - minLength  - return contigs longer than bases [0]
        - genomeFrac - return the longest contigs until genomeFrac is reached [all]
        """
        # get all contigs
        contigs = self.id2stats.keys()
        contigi = len(contigs)
        # filter by contig length
        if minLength:
            contigs = filter(lambda x: self.id2stats[x][0]>=minLength, self.id2stats)
        # sort by descending size
        sorted_contigs = sorted(contigs, key=lambda x: self.id2stats[x][0], reverse=reverse)
        # filter longest contigs by genome fraction
        if genomeFrac:
            totsize = 0
            for contigi, c in enumerate(sorted_contigs, 1):
                totsize += self.id2stats[c][0]
                if totsize >= genomeFrac*self.genomeSize:
                    break
        return sorted_contigs[:contigi]

    def get_N_and_L(self, genomeFrac=0.5, return_L=False, genomeSize=None):
        """Return N50 (and L50 if return_L) of given FastA.

        - genomeFrac - calculate N (and L) for this fraction of genome [0.5 for N50 & L50]
        - return NL  - return N50 (contig size) and L50 (number of contigs) [False]
        - genomeSize - if provided, it will use this size of the genome
                       instead of size of give assembly
        """
        if not genomeSize:
            genomeSize = self.genomeSize
        # parse contigs by descending size
        totsize = 0
        for i, x in enumerate(sorted(self.id2stats.values(), reverse=True), 1):
            size = x[0]
            totsize += size
            if totsize >= genomeFrac*genomeSize:
                break
        # return N & L
        if return_L:
            return size, i
        # return just N
        return size
        
    def N90(self):
        """Return N90"""
        return self.get_N_and_L(0.9)

    def L90(self):
        """Return N90"""
        return self.get_N_and_L(0.9, return_L=True)[1]

    def N50(self):
        """Return N90"""
        return self.get_N_and_L()

    def L50(self):
        """Return N90"""
        return self.get_N_and_L(return_L=True)[1]

    def GC(self):
        """Return GC and number of Ns"""
        # catch errors ie. empty files
        #if len(basecounts) != 4:
        #    return "%s\t[ERROR] Couldn't read file content\n"%handle.name
        #bases = set(self.basecounts)
        (A, C, G, T) = self.basecounts
        
        GC = 100.0*(int(G) + int(C)) / sum(self.basecounts)
        
        return GC

    def stats(self):
        """Return FastA statistics aka fasta_stats"""
        longest = max(stats[0] for stats in self.id2stats.values())
        lengths1000 = [x[0] for x in self.id2stats.values() if x[0]>=1000]
        contigs1000 = len(lengths1000)
        _line = '%s\t%s\t%s\t%.3f\t%s\t%s\t%s\t%s\t%s\t%s\n'
        line = _line % (self.fasta, len(self), self.genomeSize, self.GC(), contigs1000, sum(lengths1000),
                        self.N50(), self.N90(), self.Ns, longest)
        return line

def main():
    import argparse
    usage	 = "%(prog)s -i " #usage=usage, 
    parser	= argparse.ArgumentParser(description=desc, epilog=epilog, \
                                          formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--version', action='version', version='0.11c')	 
    parser.add_argument("-v", "--verbose", default=False, action="store_true",
                        help="verbose")	
    parser.add_argument("-i", "--fasta", type=str, 
                        help="FASTA file(s)")
    parser.add_argument("-o", "--out",	 default=sys.stdout, type=argparse.FileType('w'), 
                        help="output stream	 [stdout]")
    parser.add_argument("-r", "--regions", nargs='*', default=[], 
                        help="contig(s) or contig region(s) to output (returns reverse complement if end larger than start)")
    parser.add_argument("-N", default=0, type=int, 
                        help="calculate NXX and exit ie N50")
    parser.add_argument("-L", default=0, type=int, 
                        help="calculate LXX and exit ie L50")
    parser.add_argument("-S", "--stats", default=False, action="store_true",
                        help="return FastA stats aka fasta_stats")

    o = parser.parse_args()
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    if o.verbose:
        sys.stderr.write("Options: %s\n"%str(o))

    # init faidx
    faidx = FastaIndex(o.fasta, o.verbose)

    # report N & L
    if o.N:
        o.out.write("%s\n"%faidx.get_N_and_L(o.N/100.))
    if o.L:
        o.out.write("%s\n"%faidx.get_N_and_L(o.L/100., return_L=True)[1])
    # fasta_stats
    if o.stats:
        o.out.write(faidx.stats())

        print(faidx.stats())
        
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
