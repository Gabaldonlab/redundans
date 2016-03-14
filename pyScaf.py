#/usr/bin/env python
"""Perform scaffolding of contigs using information from PE & MP libraries. 
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 12/03/2016
"""

import sys, resource
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class SimpleGraph(object):
    """Graph class to represent scaffolds."""
    def __init__(self, genome, mapq=10, limit=float('inf'), frac=1.5,
                 mingap=15, ratio=0.75, minlinks=5, printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "SimpleGraph"
        self.log = log
        self.printlimit = printlimit
        # load fasta into index
        self.sequences = SeqIO.index_db(genome+".db3", genome, 'fasta')
        self.seq = self.sequences
        # prepare storage
        self.contigs = {c: len(self.seq[c]) for c in self.seq}
        self.links   = {c: [{}, {}] for c in self.contigs}
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
        self.rlen = 100

    def shorter(self, v, i=4, sep="_"):
        """Return shortened contig name"""
        return sep.join(v.split(sep)[:i])
        
    def __str__(self):
        """Produce string representation of the graph"""
        out = '%s\n%s contigs: %s\n\n%s links:\n' % (self.name, len(self.contigs), self.contigs.keys()[:self.printlimit], self.ilinks)
        i = 0
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                '''if not self.links[ref1][end1]:
                    continue
                ref2, end2, links, gap = self.links[ref1][end1]'''
                for (ref2, end2), (links, gap) in self.links[ref1][end1].items():
                    # skip if v2 longer than v1
                    if self.contigs[ref1] < self.contigs[ref2]:
                        continue
                    out += ' %s (%s) - %s (%s) with %s links; %s bp gap\n' % (self.shorter(ref1), end1, self.shorter(ref2), end2, links, gap)
                    # print up to printlimit
                    i += 1
                    if i > self.printlimit:
                        break
        return out
        
    def _add_line(self, ref1, ref2, end1, end2, links, gap):
        """Add connection between contigs. """
        # store connection details
        #if not self.links[ref1][end1]:
        #    self.links[ref1][end1] = {}
        self.links[ref1][end1][(ref2, end2)] = (links, gap)
        # update connection counter 
        self.ilinks += 1
        '''
            self.links[ref1][end1] = (ref2, end2, links, gap)
            self.links[ref2][end2] = (ref1, end1, links, gap)
            
            # update connection counter 
            self.ilinks += 1
        elif self.log and self.links[ref1][end1][0] != ref2:
            self.log.write("[WARNING] Overwritting existing connection %s with %s!\n"%(str(self.links[ref1][end1]), str((ref2, end2, links, gap))))
        '''
            
    def info(self):
        """Add sequencing library as ReadGraph"""
        if self.log:
            self.log.write(" %s kb\n"%resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        
    def add_library(self, handle, name="lib1", isize=300, stdev=50, orientation="FR"):
        """Add sequencing library as ReadGraph. 

        Contigs are connected later. 
        """
        rg = ReadGraph(self.contigs, handle, name, isize, stdev, orientation, \
                       self.mapq, self.limit, self.frac, self.ratio, self.minlinks,
                       log=self.log, printlimit=self.printlimit)
        print rg
        #self.libraries.append(rg)

        # populate graph
        for ref1, ref2, end1, end2, links, gap in rg.get_links():
            self._add_line(ref1, ref2, end1, end2, links, gap)
        self.info()
        
    def _simplify(self):
        """Simplify scaffold graph by resolving circles"""
        # remove non-reciprocal connections
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                if not self.links[ref1][end1]:
                    continue
                for (ref2, end2), (links, gap) in self.links[ref1][end1].items(): 
                    if not (ref1, end1) in self.contigs[ref2][end2]:
                        # remove link and update counter
                        self.log.write("Removing connection %s -> %s\n"%(str(ref1, end1)(), str((ref2, end2))))
                        self.links[ref1][end1].pop((ref2, end2))
                        self.ilinks -= 1
                        
        ## make sure there are no circles

    def _populate_scaffold(self, links, pend, sid, scaffold, orientations, gaps, porientation):
        """Add links to scaffold representation"""
        ref, end, links, gap = links
        # skip if already added
        if ref in self.contig2scaffold:
            return scaffold, orientations, gaps, porientation
        # get orientation - get forward/reverse-complement signal by XOR
        ## if previous is 1 (end) & current is (0) start & orientation is 0 (forward) --> keep forward orientation
        orientation = (pend != end) != porientation
        # store at the end if previous contig was F and pend 1
        if porientation != pend:
            scaffold.append(ref)
            orientations.append(not orientation)
            gaps.append(gap)
        else:
            scaffold.insert(0, ref)
            orientations.insert(0, not orientation)
            gaps.insert(0, gap)
        # update contigs2scaffold info
        self.contig2scaffold[ref] = sid
        # populate further connections from another end
        links = self.links[ref][abs(end-1)]
        # skip if not links
        if not links:
            return scaffold, orientations, gaps, orientation
        # populate further connections from another end
        return self._populate_scaffold(links, end, sid, scaffold, orientations, gaps, orientation)

    def _get_scaffold(self, name, scaffold, orientations, gaps):
        """"Prepare fasta sequence for given scaffold"""
        # add empty gap at the end
        gaps.append(0)
        seqs = []
        for c, reverse, gap in zip(scaffold, orientations, gaps):
            if reverse:
                seq = self.seq[c].reverse_complement().seq
            else:
                seq = self.seq[c].seq
            # adjust gap size
            if gap and gap < self.mingap:
                strip = int(gap - self.mingap)
                seq = seq[:strip]
                gap = self.mingap
            seqs.append(str(seq)+"N"*gap)
        # generate seq record & report fasta
        r = SeqRecord(Seq("".join(seqs), IUPAC.ambiguous_dna), id=name)
        return r
        
    def save(self, out, format='fasta'):
        """Resolve & report scaffolds"""
        # simplify graph
        self._simplify()
        
        # build scaffolds
        self.scaffolds = []
        self.contig2scaffold = {}
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            # skip if already added
            if ref1 in self.contig2scaffold:
                continue
            # get scaffold id
            sid = len(self.scaffolds)
            self.contig2scaffold[ref1] = sid
            # store scaffold and its orientation (0-forward, 1-reverse-complement) and gaps
            ## consider blist instead of list!
            scaffold, orientations, gaps = [ref1], [0], []
            # populate scaffold with connections from both ends
            for end1 in range(2):
                links = self.links[ref1][end1]
                if not links:
                    continue
                scaffold, orientations, gaps, porientation = self._populate_scaffold(links, end1, sid, scaffold, orientations, gaps, 0)
            # store
            self.scaffolds.append((scaffold, orientations, gaps))
        # report
        #print self.contig2scaffold
        totsize = 0
        self.log.write("# name\tsize\tcontigs\n")
        for i, (scaffold, orientations, gaps) in enumerate(self.scaffolds, 1):
            # scaffold00001
            name = str("scaffold%5i"%i).replace(' ','0')
            # save scaffold
            r = self._get_scaffold(name, scaffold, orientations, gaps)
            out.write(r.format('fasta'))
            # report info
            totsize += len(r)
            self.log.write("%s\t%s\t%s\n"%(name, len(r), " ".join(scaffold)))
        # close output
        out.close()
        self.log.write("%s bp in %s scaffolds.\n"%(totsize, len(self.scaffolds)))
            
class ReadGraph(SimpleGraph):
    """Graph class to represent paired alignments."""
    def __init__(self, contigs, handle, name, isize, stdev, orientation, \
                 mapq, limit, frac, ratio, minlinks, 
                 log=sys.stderr, printlimit=10):
        self.name = name
        self.log = log
        self.printlimit = printlimit
        # prepare storage
        self.contigs = contigs # {c: s for c, s in contigs.items() if s>isize/frac}
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
        self.rlen = None
        # load alignments
        self.load_from_SAM(handle)
                    
    def __str__(self):
        """Produce string representation of the graph"""
        # header
        out = '%s\n%s contigs: %s\n\n%s links' % (self.name, len(self.contigs), self.contigs.keys()[:self.printlimit], self.ilinks)
        i = 0
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                for ref2 in sorted(self.links[ref1][end1], key=lambda x: self.contigs[x], reverse=1):
                    for end2 in range(2):
                        # skip if v2 longer than v1
                        if self.contigs[ref1] < self.contigs[ref2]:
                            continue
                        # get positions of links
                        positions = self.links[ref1][end1][ref2][end2]
                        if not positions:
                            continue
                        # print up to printlimit
                        i += 1
                        if i > self.printlimit:
                            break
                        out += ' %s (%s) - %s (%s) with %s links: %s\n' % (self.shorter(ref1), end1, self.shorter(ref2), end2, len(positions), str(positions[:self.printlimit]))
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
        self.links[ref1][end1][ref2][end2].append(float('%s.%s'%(d1, d2)))
        self.links[ref2][end2][ref1][end1].append(float('%s.%s'%(d2, d1)))
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
            # skip
            if self.limit and i > self.limit:
                break
            if not self.rlen:
                self.rlen = int(len1+len2)/2
            # skip alignments with low quality
            if self.mapq:
                if mapq1 < self.mapq or mapq2 < self.mapq:
                    continue  
            j   += 1
            #skip self matches
            if r1==r2:
                #isizes.append(abs(s2-s1))
                continue
            k += 1
            #print output
            self.add_line(r1, r2, s1, s2, flag1, flag2)
        #return isizes
        
    def _get_major_link(self, links, c1, end1):
        """Return major link if any of the links full fill quality criteria.

        So far this is too simplistic!
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

    def _cluster_links(self, links, maxdist=200):
        """Merge links from the same regions"""
        def get_median(positions):
            """return median value of first pos"""
            return np.median([int(str(pos).split('.')[0]) for pos in positions])
            
        alllinks = {}
        for (c, e, pos) in sorted(links, key=lambda x: len(x[-1]), reverse=1):
            # get median pos - maybe use mean & stdev?
            mpos = get_median(pos)
            if not alllinks:
                alllinks[mpos] = [(c, e, pos)]
                continue
            # sort by distance #filter(lambda x: abs(x-mpos)<maxdist, alllinks)
            closest = sorted(alllinks.keys(), key=lambda x: abs(x-mpos))[0]
            if abs(closest-mpos)<maxdist:
                alllinks[closest].append((c, e, pos))
            else:
                alllinks[mpos] = [(c, e, pos)]
        
        for mpos, links in alllinks.items():
            yield links
        
    def _select_links(self, links):
        """Return best links between set of contigs"""
        if not links:
            return []
        # reorganise links into list
        links = [(c, e, pos) for c in links for e, pos in enumerate(links[c])]
        # split links from different regions
        bests = []
        alllinks = self._cluster_links(links)
        for links in alllinks:
            best = links[0] #sorted(links, key=lambda x: len(x[-1]), reverse=1)[0]
            # skip if not enough links or many links to more than one contigs
            sumi = sum(len(pos) for c, e, pos in links)
            if len(best[-1]) < self.minlinks or len(best[-1]) < self.ratio*sumi:
                continue
            ### also check if the same connection in reciprocal bests - no need, will do it during scaffolding
            bests.append(best)
        # make sure best links are not overlapping too much
        
        return bests

    def _calculat_gap_size(self, positions):
        """Return estimated size of the gap"""
        # unload positions
        dists = [self.isize - sum(map(int, str(pos).split('.'))) for pos in positions]
        return np.median(dists)-self.rlen
        
    def _filter_links(self):
        """Filter links by removing contigs with too many connections ie. repeats"""
        # first identify contigs with too many matches
        
    def get_links(self):
        """Generator of contig connections for given library."""
        # memory profilling
        #from pympler import asizeof; print "Links: %s bases\n"%asizeof.asizeof(self.links)
        # filter links
        self._filter_links()
        # process starting from the longest 
        for c in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end in range(2):
                # best connections
                for c2, end2, pos in self._select_links(self.links[c][end]):
                    # get distance
                    dist = self._calculat_gap_size(pos)
                    # add connection
                    yield c, c2, end, end2, len(pos), dist
    
if __name__=="__main__":
    """
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1
bwa index test/run1/contigs.reduced.fa
bwa mem -t4 test/run1/contigs.reduced.fa test/600_?.fq.gz > test/run1/600.sam
bwa mem -t4 test/run1/contigs.reduced.fa test/5000_?.fq.gz > test/run1/5000.sam

    
ipython

import pyScaf as ps

sam = 'test/run1/600.sam'
sam2 = 'test/run1/5000.sam'
fasta = 'test/run1/contigs.reduced.fa'

reload(ps); s = ps.SimpleGraph(fasta, mapq=10, limit=19571, log=sys.stderr);
s.add_library(open(sam), name=sam, isize=600, stdev=100, orientation="FR"); print s
s.add_library(open(sam2), name=sam2, isize=5000, stdev=1000, orientation="FR"); print s

s.save(out=open(fasta+".scaffolds.fa", "w"))


bwa index test/run1/contigs.reduced.fa.scaffolds.fa
bwa mem -t4 test/run1/contigs.reduced.fa.scaffolds.fa test/5000_?.fq.gz > test/run1/5000.2.sam
    
sam2 = 'test/run1/5000.2.sam'
fasta = 'test/run1/contigs.reduced.fa.scaffolds.fa'

reload(ps); s = ps.SimpleGraph(fasta, mapq=10, limit=19571, log=sys.stderr);
s.add_library(open(sam2), name=sam2, isize=5000, stdev=1000, orientation="FR"); print s

    """
    sam = 'test/run1/600.sam'
    sam2 = 'test/run1/5000.sam'
    fasta = 'test/run1/contigs.reduced.fa'
    s = SimpleGraph(fasta, mapq=10, limit=19571);
    s.add_library(open(sam), name=sam, isize=600, stdev=50, orientation="FR"); print s
    #s.add_library(open(sam2), name=sam2, isize=5000, stdev=1000, orientation="FR"); print s
    s.save(out=open(fasta+".scaffolds.fa", "w"))


