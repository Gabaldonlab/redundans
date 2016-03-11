"""
Some graph functions. Usefull, to work with orthologs.
"""

import sys

class Graph(object):
    """Undirected Graph class."""
    def __init__(self, vertices=[], features=[], frac=0.95, \
                 isize=300, stdev=50, log=sys.stderr, printlimit=10):
        """Construct a graph with the given vertices & features"""
        self.printlimit = printlimit
        # populate vertices
        self._neighbours = {v: {} for v in vertices}
        self._features   = {v: f for v, f in zip(vertices, features)}
        # insert size
        self.log   = log
        self.frac  = frac
        self.isizes = []
        self.isize = isize
        self.stdev = stdev
        # links
        self.ilinks = 0
        self.links = self._neighbours
        self.sizes = self._features

    def __str__(self):
        """Produce string representation of the graph"""
        out = '%s Vertices: %s \n%s Lines:\n' % (len(self.sizes), self._neighbours.keys()[:self.printlimit], self.ilinks)
        i = 0
        for vertex, neighbours in self._neighbours.items():
            for neighbour, links in neighbours.iteritems():
                if vertex < neighbour:
                    i += 1
                    if i>self.printlimit:
                        break
                    out += '\t%s - %s [%s links]\n' % (vertex, neighbour, len(links))
        return out

    def _present(self, v):
        """Return True if v in present in the graph"""
        if v in self._features:
            return True
        if self.log:
            log.write("[WARNING] Vertex (%s) not in vertices\n"%(v))
            
    def _check_distance(self, v, p):
        """Return True if distance from the v end is smaller that frac * isize"""
        if self._present(v):
            if p < self.frac*self.isize or self.sizes[v]-p < self.frac*self.isize:
                return 
        return True
                
    def add_line(self, v1, v2, p1, p2):
        """Add a line from v1 to v2"""
        # check distance from the end
        if self._check_distance(v1, p1) or self._check_distance(v2, p2):
            return
        # add connection
        if v2 not in self._neighbours[v1]:
            self._neighbours[v1][v2] = []
            self._neighbours[v2][v1] = []
        self._neighbours[v1][v2].append((p1,p2))
        self._neighbours[v2][v1].append((p2,p1))
        self.ilinks += 1
       
    '''    
    def delete_line(self, v1, v2):
        """Delete line between v1 and v2 if any. 
        Return error if no line between them.
        """
        self._check((v1,v2))
        try:
            self._neighbours[v1].remove(v2)
            self._neighbours[v2].remove(v1)
        except KeyError:
            raise ValueError('No line between %s and %s' % (v1,v2))
    '''


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

def process_alignments(g, handle, mapqTh=0, upto=float('inf'), verbose=False, log=sys.stderr):
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
        g.add_line(r1, r2, s1, s2)
    return isizes

"""
from Bio import SeqIO
sam=open('test/run1/600.sam')
contigs, sizes = [], []
for r in SeqIO.parse('test/run1/contigs.reduced.fa', 'fasta'):
    contigs.append(r.id)
    sizes.append(len(r))

#g=Graph(contigs, sizes)
import pyScaf
reload(pyScaf); g=pyScaf.Graph(contigs, sizes); pyScaf.process_alignments(g, sam)
"""