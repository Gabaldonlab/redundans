#/usr/bin/env python
""" """

import sys

class SimpleGraph(object):
    """Graph class to represent scaffolds."""
    def __init__(self, contigs=[], printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.log = log
        self.printlimit = printlimit
        # prepare storage 
        self.contigs = set(contigs)
        self.links   = {c: {} for c in self.contigs}
        self.ilinks  = 0

    def add_line(self, ref1, ref2, links=0, rev1=0, rev2=0):
        # store connection details
        self.links[ref1][ref2] = (links, rev1, rev2) # it can be int int(""%)
        self.links[ref2][ref1] = (links, rev2, rev1) 
        # update connection counter 
        self.ilinks += 1
        
class Graph(object):
    """Graph class to represent scaffolds."""
    def __init__(self, contigs=[], sizes=[], frac=1.5, \
                 isize=300, stdev=50, orientation="FR", ratio=0.75, \
                 log=sys.stderr, printlimit=10):
        """Construct a graph with the given vertices & features"""
        self.log = log
        self.printlimit = printlimit
        # prepare storage 
        self.contigs = {c: s  for c, s in zip(contigs, sizes)}# if s>isize/frac}
        self.links   = {c: {} for c in self.contigs}
        self.ilinks  = 0
        # scaffolding options
        self.ratio = ratio
        # init orientation
        self._update_library_info(frac, isize, stdev, orientation)
                    
    def __str__(self):
        """Produce string representation of the graph"""
        # header
        out = '%s Vertices: %s\n%s Lines:\n' % (len(self.contigs), self.contigs.keys()[:self.printlimit], self.ilinks)
        i = 0
        for v1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for v2 in sorted(self.links[v1], key=lambda x: self.contigs[x], reverse=1):
                # skip if v2 longer than v1
                if self.contigs[v1] < self.contigs[v2]:
                    continue
                # get positions of links
                positions = self.links[v1][v2]
                # print up to printlimit
                i += 1
                if i > self.printlimit:
                    break
                out += '\t%s - %s with %s links: %s\n' % (v1, v2, len(positions), str(positions[:self.printlimit]))
        return out

    def _update_library_info(self, frac, isize, stdev, orientation):
        """Update orientatation dependent stuff"""
        # set _check function
        if   orientation in (1, "FR"):
            self.get_distance = self._get_distance_FR
        elif orientation in (2, "RF"):
            self.get_distance = self._get_distance_RF
        else:
            self.log.write("[WARNING] Provided orientation (%s) is not supported!\n"%orientation)
        # insert size
        self.frac  = frac
        self.isize = isize
        self.stdev = stdev
        self.orientation = orientation
        # max dist
        self.maxdist = self.isize + self.frac*self.stdev

    def _present(self, c):
        """Return True if v in present in the graph"""
        if c in self.links:
            return True
        if self.log:
            self.log.write("[WARNING] %s not in contigs!\n"%c)
            
    def _get_distance_FR(self, ref, pos, flag):
        """Return True if distance from the v end is smaller that frac * isize"""
        # FR - /1 and /2 doesn't matter, only F / R
        if self._present(ref):
            # reverse -> beginning of the contig
            if flag&16:
                return pos
            # forward -> end of the contig
            else:
                return self.contigs[ref]-pos

    def _get_distance_RF(self, ref, pos, flag):
        """Return True if distance from the v end is smaller that frac * isize"""
        # FR - /1 and /2 doesn't matter, only F / R
        if self._present(ref):
            # reverse -> beginning of the contig
            if not flag&16:
                return pos
            # forward -> end of the contig
            else:
                return self.contigs[ref]-pos
                
    def add_line(self, ref1, ref2, pos1, pos2, flag1, flag2):
        """Add a line from v1 to v2"""
        # check if distance is correct
        ## note, distance is corrected by pair orientation
        d1 = self.get_distance(ref1, pos1, flag1)
        d2 = self.get_distance(ref2, pos2, flag2)
        if not d1 or not d2 or d1 + d2 > self.maxdist:
            return
        # add connection
        if ref2 not in self.links[ref1]:
            self.links[ref1][ref2] = []
            self.links[ref2][ref1] = []
        # store connection details
        self.links[ref1][ref2].append((pos1, pos2))
        self.links[ref2][ref1].append((pos2, pos1))
        # update connection counter 
        self.ilinks += 1

    def _split_links(self, links):
        """Sepearte upsteam and downstream connections"""
        
    def scaffold(self):
        """Perform contig scaffolding"""
        # start simple graph
        s = SimpleGraph(self.contigs, log=self.log)
        # process starting from the longest
        added = set()
        for c1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            up, down = self._split_links(self.links(c1))
            # connect only contigs with no. of links above ratio
            #for v2 in sorted(self.links[v1], key=lambda x: self.contigs[x], reverse=1):
                
        
    def delete_line(self, v1, v2):
        """Delete line between v1 and v2 if any. 
        Return error if no line between them.
        """
        self._check((v1,v2))
        try:
            self.links[v1].remove(v2)
            self.links[v2].remove(v1)
        except KeyError:
            raise ValueError('No line between %s and %s' % (v1,v2))


def parse_sam(handle):
    """Return sam tuple"""
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
            q1,flag1,ref1,start1,mapq1,len1 = sam[0],int(sam[1]),sam[2],int(sam[3]),int(sam[4]),len(sam[9])
        else:
            #skip multiple matches
            if sam[0] == q2:
                continue
            q2,flag2,ref2,start2,mapq2,len2 = sam[0],int(sam[1]),sam[2],int(sam[3]),int(sam[4]),len(sam[9])
        #report
        if q1 == q2:
            yield q1,flag1,ref1,start1,mapq1,len1,q2,flag2,ref2,start2,mapq2,len2

def get_start_stop( start,length,flag ):
    """Return start-end read boundaries.
    Return end-start if reverse aligned (flag & 16)."""
    if flag & 16:
        end    = start
        start += length
    else:
        end    = start + length
    return start, end

def process_alignments(g, handle, orientation=1, mapqTh=0, upto=float('inf'), verbose=False, log=sys.stderr):
    """Convert SAM to SSPACE TAB file."""
    i = j = k = pq1 = 0
    isizes = []
    for q1, flag1, r1, s1, mapq1, len1, q2, flag2, r2, s2, mapq2, len2 in parse_sam(handle):
        i   += 1
        if upto and i>upto:
            break        
        #skip 0 quality pair
        if mapqTh:
            if mapq1<mapqTh or mapq2<mapqTh:
                continue  
        j   += 1
        #skip self matches
        if r1==r2:
            isizes.append(abs(s2-s1))
            continue
        k += 1
        #define start-stop ranges
        #start1,end1 = get_start_stop(s1,len1,flag1)
        #start2,end2 = get_start_stop(s2,len2,flag2)    
        #print output
        g.add_line(r1, r2, s1, s2, flag1, flag2)
    return isizes

    
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

reload(ps); g = ps.Graph(contigs, sizes, isize=600, orientation="FR");  isizes = ps.process_alignments(g, open(sam), upto=19571); print g
reload(ps); g = ps.Graph(contigs, sizes, isize=5000, orientation="FR"); isizes = ps.process_alignments(g, open(sam2), upto=19571); print g

    """
    from Bio import SeqIO
    sam = 'test/run1/600.sam'
    fasta = 'test/run1/contigs.reduced.fa'
    contigs, sizes = [], []
    for r in SeqIO.parse(fasta, 'fasta'):
        contigs.append(r.id)
        sizes.append(len(r))

    g = Graph(contigs, sizes, orientation="FR", isize=600)
    
    isizes = process_alignments(g, open(sam), upto=19571)
    print g

