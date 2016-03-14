#/usr/bin/env python
"""


Note, memory usage can be dropped by half using str instead of tuple for links storage

from pympler import asizeof
l = g.links['NODE_1_length_29603_cov_57.5935_ID_1']['NODE_3_length_7388_cov_80.798_ID_5']
l2=['%s.%s.%s'%e for e in l]
1.0*asizeof.asizeof(l2)/asizeof.asizeof(l) # 0.45321229050279327

"""

import sys
import numpy as np

class SimpleGraph(object):
    """Graph class to represent scaffolds."""
    def __init__(self, contigs, sizes, mapq=10, limit=float('inf'), frac=1.5,
                 mingap=50, ratio=0.75, minlinks=5, printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "SimpleGraph"
        self.log = log
        self.printlimit = printlimit
        # prepare storage
        self.contigs = {c: s  for c, s in zip(contigs, sizes)}# if s>isize/frac}
        self.links   = {c: [0, 0] for c in self.contigs}
        self.ilinks  = 0
        # alignment options
        self.mapq  = mapq
        self.frac  = frac
        self.limit = limit
        # scaffolding
        self.ratio    = ratio
        self.minlinks = minlinks
        self.mingap  = mingap
        # read libraries storage
        self.libraries = []

    def shorter(self, v, i=4, sep="_"):
        """Return shortened contig name"""
        return sep.join(v.split(sep)[:i])
        
    def __str__(self):
        """Produce string representation of the graph"""
        out = '%s\n%s contigs: %s\n\n%s links:\n' % (self.name, len(self.contigs), self.contigs.keys()[:self.printlimit], self.ilinks)
        i = 0
        for v1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                if not self.links[v1][end1]:
                    continue
                v2, end2, links, gap = self.links[v1][end1]
                # skip if v2 longer than v1
                if self.contigs[v1] < self.contigs[v2]:
                    continue
                out += ' %s (%s) - %s (%s) with %s links; %s bp gap\n' % (self.shorter(v1), end1, self.shorter(v2), end2, links, gap)
                # print up to printlimit
                i += 1
                if i > self.printlimit:
                    break
        return out
        
    def _add_line(self, ref1, ref2, end1, end2, links, gap):
        """Add connection between contigs. """
        # store connection details
        if not self.links[ref1][end1]:
            self.links[ref1][end1] = (ref2, end2, links, gap)
            self.links[ref2][end2] = (ref1, end1, links, gap) 
            # update connection counter 
            self.ilinks += 1
        elif self.links[ref1][end1][0] != ref2:
            log.write("[WARNING] Overwritting existing connection %s with %s!\n"%(str(self.links[ref1][end1]), str((ref2, end2, links, gap))))

    def add_library(self, handle, name="lib1", isize=300, stdev=50, orientation="FR"):
        """Add sequencing library as ReadGraph. 

        Contigs are connected later. 
        """
        rg = ReadGraph(self.contigs, handle, name, isize, stdev, orientation, \
                       self.mapq, self.limit, self.frac, self.ratio, self.minlinks)
        print rg
        #self.libraries.append(rg)

        # populate graph
        for ref1, ref2, end1, end2, links, gap in rg.get_links():
            self._add_line(ref1, ref2, end1, end2, links, gap)
        
    def _simplify(self):
        """Simplify scaffold graph.

        For not it's doing nothing.
        """
        ## make sure there is not circles

    def save(self, out, genome, format='fasta'):
        """Report scaffolds 

        For not it's doing nothing.
        """
        # simplify graph
        self._simplify()
        # and report
        

class ReadGraph(SimpleGraph):
    """Graph class to represent paired alignments."""
    def __init__(self, contigs, handle, name, isize, stdev, orientation, \
                 mapq, limit, frac, ratio, minlinks, 
                 log=sys.stderr, printlimit=10):
        self.name = name
        self.log = log
        self.printlimit = printlimit
        # prepare storage
        self.contigs = contigs
        self.links   = {c: [{},{}] for c in self.contigs}
        self.ilinks  = 0
        # alignment options
        self.mapq  = mapq
        self.limit = limit
        # set distance function
        if   orientation in (1, "FR"):
            self.get_distance = self._get_distance_FR
        elif orientation in (2, "RF"):
            self.get_distance = self._get_distance_RF
        else:
            self.log.write("[WARNING] Provided orientation (%s) is not supported!\n"%orientation)
        # insert size
        self.isize = isize
        self.stdev = stdev
        self.orientation = orientation
        self.minlinks = minlinks
        self.ratio = ratio
        # max dist
        self.maxdist = self.isize + frac * self.stdev
        # load alignments
        self.load_from_SAM(handle)
                    
    def __str__(self):
        """Produce string representation of the graph"""
        # header
        out = '%s\n%s contigs: %s\n\n%s links' % (self.name, len(self.contigs), self.contigs.keys()[:self.printlimit], self.ilinks)
        i = 0
        for v1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                for v2 in sorted(self.links[v1][end1], key=lambda x: self.contigs[x], reverse=1):
                    for end2 in range(2):
                        # skip if v2 longer than v1
                        if self.contigs[v1] < self.contigs[v2]:
                            continue
                        # get positions of links
                        positions = self.links[v1][end1][v2][end2]
                        if not positions:
                            continue
                        # print up to printlimit
                        i += 1
                        if i > self.printlimit:
                            break
                        out += ' %s (%s) - %s (%s) with %s links: %s\n' % (self.shorter(v1), end1, self.shorter(v2), end2, len(positions), ", ".join(positions[:self.printlimit]))
        return out

    def _present(self, c):
        """Return True if v in present in the graph"""
        if c in self.links:
            return True
        if self.log:
            self.log.write("[WARNING] %s not in contigs!\n"%c)
            
    def _get_distance_FR(self, ref, pos, flag):
        """Return True if distance from the v end is smaller that frac * isize"""
        # FR - /1 and /2 doesn't matter, only F / R
        # reverse -> beginning of the contig
        if flag&16:
            return pos, 0
        # forward -> end of the contig
        else:
            return self.contigs[ref]-pos, 1

    def _get_distance_RF(self, ref, pos, flag):
        """Return distance and end from which the read is originating"""
        # RF - /1 and /2 doesn't matter, only F / R
        # forward -> beginning of the contig
        if not flag&16:
            return pos, 0
        # reverse -> end of the contig
        else:
            return self.contigs[ref]-pos, 1
                
    def add_line(self, ref1, ref2, pos1, pos2, flag1, flag2):
        """Add a line from v1 to v2"""
        if not self._present(ref1) or not self._present(ref2):
            return
        # check if distance is correct
        ## note, distance is corrected by pair orientation
        d1, end1 = self.get_distance(ref1, pos1, flag1)
        d2, end2 = self.get_distance(ref2, pos2, flag2)
        if d1 + d2 > self.maxdist:
            return
        # add connection
        if ref2 not in self.links[ref1][end1]:
            self.links[ref1][end1][ref2] = [[],[]]
        if ref1 not in self.links[ref2][end2]:
            self.links[ref2][end2][ref1] = [[],[]]
        # store connection details
        self.links[ref1][end1][ref2][end2].append('%s-%s'%(d1, d2))
        self.links[ref2][end2][ref1][end1].append('%s-%s'%(d2, d1))
        # update connection counter 
        self.ilinks += 1

    def sam_parser(self, handle):
        """Return tuple representing paired alignment from SAM format."""
        q1 = q2 = ""
        for l in handle:
            l = l.strip()
            if not l or l.startswith('@'):
                continue
            sam = l.split('\t')
            #first in pair
            if int(sam[1]) & 64:
                #skip multiple matches
                if sam[0] == q1:
                    continue
                q1, flag1, ref1, start1, mapq1, len1 = sam[0], int(sam[1]), sam[2], int(sam[3]), int(sam[4]), len(sam[9])
            else:
                #skip multiple matches
                if sam[0] == q2:
                    continue
                q2, flag2, ref2, start2, mapq2, len2 = sam[0], int(sam[1]), sam[2], int(sam[3]), int(sam[4]), len(sam[9])
            #report
            if q1 == q2:
                yield q1, flag1, ref1, start1, mapq1, len1, q2, flag2, ref2, start2, mapq2, len2

    def get_start_stop(self, start, length, flag):
        """Return start-end read boundaries.
        Return end-start if reverse aligned (flag & 16)."""
        if flag & 16:
            end    = start
            start += length
        else:
            end    = start + length
        return start, end
            
    #def load_SAM(handle, mapq=self.mapq, limit=self.limit):
    def load_from_SAM(self, handle):
        """Populate graph with links from SAM file.

        Note, SAM file has to be ordered by read name, F read always before R.
        Multiple alignments are allowed, but only the first (assuming best)
        alignment is taken into accound. 
        """
        # parse SAM
        i = j = k = pq1 = 0
        for q1, flag1, r1, s1, mapq1, len1, q2, flag2, r2, s2, mapq2, len2 in self.sam_parser(handle):
            i   += 1
            if self.limit and i > self.limit:
                break        
            #skip 0 quality pair
            if self.mapq:
                if mapq1 < self.mapq or mapq2 < self.mapq:
                    continue  
            j   += 1
            #skip self matches
            if r1==r2:
                #isizes.append(abs(s2-s1))
                continue
            k += 1
            # define start-stop ranges - correct distance by orientation and strand!!
            #start1,end1 = get_start_stop(s1,len1,flag1)
            #start2,end2 = get_start_stop(s2,len2,flag2)    
            #print output
            self.add_line(r1, r2, s1, s2, flag1, flag2)
        #return isizes
        
    def _get_major_link(self, links):
        """Return major link if any of the links full fill quality criteria.

        So far this is super simplistics!
        It works for PE, but for MP you need to allow for multi-joins,
        as many contigs will be shorter than i.e. 5kb insert... 
        """
        if not links:
            return
        # reorganise links into list
        links = [(c, e, pos) for c in links for e, pos in enumerate(links[c])]
        # sort starting from contig with most links
        best  = sorted(links, key=lambda x: x[-1], reverse=1)[0]
        # skip if not enough links or many links to more than one contigs
        sumi = sum(len(pos) for c, e, pos in links)#; print sumi, best
        if len(best[-1]) < self.minlinks or len(best[-1]) < self.ratio*sumi:
            return
        return best

    def _calculat_gap_size(self, positions):
        """Return estimated size of the gap"""
        # unload positions
        dists = [self.isize - sum(map(int, pos.split('-'))) for pos in positions]
        return np.median(dists)
        
    def get_links(self):
        """Generator of contig connections for given library."""
        # process starting from the longest
        for c in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end in range(2):
                # major connection
                major = self._get_major_link(self.links[c][end])
                if not major:
                    continue
                # check if reciprocal
                c2, end2, pos2 = major
                rmajor = self._get_major_link(self.links[c2][end2])
                if not rmajor:
                    continue
                c1, end1, pos1 = rmajor
                if c==c1 and end==end1:    
                    #print c, end, c2, end2, len(pos2)
                    # get distance
                    dist = self._calculat_gap_size(pos1)
                    # add connection
                    yield c1, c2, end1, end2, len(pos1), dist
    
if __name__=="__main__":
    """
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1
bwa index test/run1/contigs.reduced.fa
bwa mem -t4 test/run1/contigs.reduced.fa test/600_?.fq.gz > test/run1/600.sam
bwa mem -t4 test/run1/contigs.reduced.fa test/5000_?.fq.gz > test/run1/5000.sam

ipython

import pyScaf as ps

from Bio import SeqIO
sam = 'test/run1/600.sam'
sam2 = 'test/run1/5000.sam'
    
fasta = 'test/run1/contigs.reduced.fa'
contigs, sizes = [], []
for r in SeqIO.parse(fasta, 'fasta'):
    contigs.append(r.id)
    sizes.append(len(r))

reload(ps); g = ps.SimpleGraph(contigs, sizes, limit=19571);
g.load_from_SAM(open(sam), isize=600, stdev=100, orientation="FR"); print g
g.scaffold(); print g.scaffolds
# g.scaffolds.save()
g.load_from_SAM(open(sam), isize=5000, stdev=1000, orientation="FR"); print g 

    """
    from Bio import SeqIO
    sam = 'test/run1/600.sam'
    sam2 = 'test/run1/5000.sam'
    fasta = 'test/run1/contigs.reduced.fa'
    contigs, sizes = [], []
    for r in SeqIO.parse(fasta, 'fasta'):
        contigs.append(r.id)
        sizes.append(len(r))

    s = SimpleGraph(contigs, sizes, mapq=10, limit=19571);
    #s.add_library(open(sam), name=sam, isize=600, stdev=100, orientation="FR"); print s
    s.add_library(open(sam2), name=sam2, isize=5000, stdev=1000, orientation="FR"); print s


