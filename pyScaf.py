#!/usr/bin/env python
desc="""Perform scaffolding of contigs using information from (in this order):
- paired-end (PE) and/or mate-pair (MP) libraries (!!!NOT IMPLEMENTED YET!!!)
- synteny to reference genome

More info at: http://bit.ly/pyScaf
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Warsaw, 12/03/2016
"""

import os, sys
import resource, subprocess
from datetime import datetime
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

class Graph(object):
    """Graph class to represent scaffolds. It shouldn't be invoked directly,
    rather use its children: PairedGraph or SyntenyGraph.

    import pyScaf as ps
    
    #####
    # Scaffolding using PE/MP libraries
    reload(ps); s = ps.ReadGraph(fasta, mapq=10, limit=19571, log=sys.stderr);
    s.add_library(fastq, name="lib1", isize=600, stdev=100, orientation="FR"); print s
    s.add_library(fastq, name="lib2", isize=5000, stdev=1000, orientation="FR"); print s
    s.save(out=open(fasta+".scaffolds.fa", "w"))

    
    ######
    # Reference-based scaffolding
    reload(ps); s = ps.SyntenyGraph('test/run1/contigs.reduced.fa', 'test/ref.fa')
    s.save(out=open(fasta+".scaffolds.ref.fa", "w"))
    """
    def __init__(self, mingap=15, printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "Graph"
        self.log = log
        self.printlimit = printlimit
        self.ilinks  = 0

    def logger(self, mssg=""):
        """Logging function."""
        head = "\n%s"%("#"*50,)
        timestamp = "\n[%s]"% datetime.ctime(datetime.now())
        memusage  = "[%5i Mb] "%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, )
        if self.log:
            self.log.write("".join((head, timestamp, memusage, mssg)))
        
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
        
    def _init_storage(self, genome):
        """Load sequences from genome, their sizes and init links"""
        # load fasta into index
        self.sequences = SeqIO.index_db(genome+".db3", genome, 'fasta')
        self.seq = self.sequences
        # prepare storage
        self.contigs = {c: len(self.seq[c]) for c in self.seq}
        self.links   = {c: [{}, {}] for c in self.contigs}
        self.ilinks  = 0
        
    def _add_line(self, ref1, ref2, end1, end2, links, gap):
        """Add connection between contigs. """
        # store connection details
        self.links[ref1][end1][(ref2, end2)] = (links, gap)
        # update connection counter 
        self.ilinks += 1
                   
    def _get_seqrecord(self, name, scaffold, orientations, gaps):
        """"Return SeqRecord for given scaffold"""
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
        # generate scaffolds
        self._get_scaffolds()
        # report
        self.logger("Reporting scaffolds...\n")
        # log scaffold structure
        log = open(out.name+".tsv", "w")
        logline = "%s\t%s\t%s\t%s\t%s\t%s\n"
        log.write("# name\tsize\tno. of contigs\tordered contigs\tcontig orientations (0-forward; 1-reverse)\tgap sizes (negative gap size = adjacent contigs are overlapping)\n")
        totsize = 0
        for i, (scaffold, orientations, gaps) in enumerate(self.scaffolds, 1):
            # scaffold00001
            name = str("scaffold%5i"%i).replace(' ','0')
            # save scaffold
            r = self._get_seqrecord(name, scaffold, orientations, gaps)
            out.write(r.format('fasta'))
            # report info
            log.write(logline%(name, len(r), len(scaffold), " ".join(scaffold),
                               " ".join(map(str, (int(x) for x in orientations))), # x may be bool!
                               " ".join(map(str, (x for x in gaps)))))
            totsize += len(r)
        # close output & loge
        out.close()
        log.close()
        self.log.write(" %s bp in %s scaffolds. Details in %s\n"%(totsize, len(self.scaffolds), log.name))
        self.log.write("Scaffolds saved to: %s\n"%out.name)
        
class ReadGraph(Graph):
    """Graph class to represent scaffolds derived from NGS libraries.
    
    It stores connections from paired alignments in self.readlinks. 
    """
    def __init__(self, genome, mapq=10, load=0.2, threads=4, 
                 frac=1.5, mingap=15, ratio=0.75, minlinks=5,
                 log=sys.stderr, printlimit=10):
        self.name = "ReadGraph"
        self.log = log
        self.printlimit = printlimit
        self.threads = threads
        # load fasta into index
        self.ref = genome
        # prepare storage
        self._init_storage(genome)
        # alignment options
        self.mapq  = mapq
        self.limit = load * sum(self.contigs.itervalues())
        # scaffolding parameters
        self.ratio    = ratio
        self.minlinks = minlinks
        self.mingap   = mingap
        self.frac     = frac
        
    def _set_isize(self, isize=0, stdev=0, orientation=0):
        """Estimate isize stats if necessary"""
        # insert size
        self.isize = isize
        self.stdev = stdev
        self.orientation = orientation
        
    # read-specific part
    def show(self):
        """Produce string representation of the graph"""
        # header
        out = '%s\n%s contigs: %s\n\n%s links' % (self.name, len(self.contigs), self.contigs.keys()[:self.printlimit], self.ireadlinks)
        i = 0
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                for ref2 in sorted(self.readlinks[ref1][end1], key=lambda x: self.contigs[x], reverse=1):
                    for end2 in range(2):
                        # skip if v2 longer than v1
                        if self.contigs[ref1] < self.contigs[ref2]:
                            continue
                        # get positions of links
                        positions = self.readlinks[ref1][end1][ref2][end2]
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
        if c in self.readlinks:
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
                
    def add_readline(self, ref1, ref2, pos1, pos2, flag1, flag2):
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
        if ref2 not in self.readlinks[ref1][end1]:
            self.readlinks[ref1][end1][ref2] = [[],[]]
        if ref1 not in self.readlinks[ref2][end2]:
            self.readlinks[ref2][end2][ref1] = [[],[]]
        # store connection details
        self.readlinks[ref1][end1][ref2][end2].append(float('%s.%s'%(d1, d2)))
        self.readlinks[ref2][end2][ref1][end1].append(float('%s.%s'%(d2, d1)))
        # update connection counter 
        self.ireadlinks += 1

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
        
    def _bwamem(self, fastqs, bwalog=open("/dev/null", "w")):
        """Return pipe to stdout of BWA MEM subprocess"""
        # build db if index missing .sa .bwt .ann .amb .pac
        if not os.path.isfile(self.ref+".sa") or not os.path.isfile(self.ref+".bwt"):
            os.system("bwa index %s" % self.ref)
        # run BWA MEM
        args = ['bwa', 'mem', '-S', '-t', str(self.threads), self.ref] + fastqs
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=bwalog)
        return proc.stdout
            
    def load_from_SAM(self, handle):
        """Populate graph with readlinks from SAM file.

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
            self.add_readline(r1, r2, s1, s2, flag1, flag2)
        #return isizes
        
    def _get_major_link(self, links, c1, end1):
        """Return major readlink if any of the links full fill quality criteria.

        So far this is too simplistic!
        It works for PE, but for MP you need to allow for multi-joins, 
        as many contigs will be shorter than i.e. 5kb insert... 
        """
        if not links:
            return []
        # reorganise links into list
        links = [(c, e, pos) for c in links for e, pos in enumerate(links[c])]
        # sort starting from contig with most links
        best  = sorted(links, key=lambda x: x[-1], reverse=1)[0]
        # skip if not enough links or many links to more than one contigs
        sumi = sum(len(pos) for c, e, pos in links)#; print sumi, best
        if len(best[-1]) < self.minlinks or len(best[-1]) < self.ratio*sumi:
            return []
        return [best,]

    def _cluster_links(self, links, maxdist=200):
        """Merge readlinks from the same regions. Experimental!

        !!!Plenty of work needed here!!!
        """
        alllinks = {}
        for (c, e, pos) in sorted(links, key=lambda x: len(x[-1]), reverse=1):
            # get median pos - maybe use mean & stdev?
            mpos = np.median([int(str(p).split('.')[0]) for p in pos]) # get_median(pos)
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
        """Return best readlinks between set of contigs"""
        if not links:
            return []
        # reorganise links into list
        links = [(c, e, pos) for c in links for e, pos in enumerate(links[c]) if pos]
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
        """Filter readlinks by removing contigs with too many connections ie. repeats"""
        # first identify contigs with too many matches
        
    def get_links(self):
        """Generate connections contig connections for given library."""
        # memory profilling
        #from pympler import asizeof; print "Links: %s bases\n"%asizeof.asizeof(self.links)
        # filter links
        self._filter_links()
        # process starting from the longest 
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                # best connections
                for ref2, end2, pos in self._select_links(self.readlinks[ref1][end1]): # self._get_major_link(self.readlinks[ref1][end1], ref1, end1):
                    # get distance
                    gap = self._calculat_gap_size(pos)
                    # add connection
                    self._add_line(ref1, ref2, end1, end2, len(pos), gap)

    def add_library(self, fastqs, isize=0, stdev=0, orientation=0):
        """Add connections from sequencing library to the graph"""
        self.logger("Adding library...\n")
        # for read info
        self.readlinks  = {c: [{},{}] for c in self.contigs}
        self.ireadlinks = 0
        # set isize
        self._set_isize(isize, stdev, orientation)
        # set distance function
        if   orientation in (1, "FR"):
            self.get_distance = self._get_distance_FR
        elif orientation in (2, "RF"):
            self.get_distance = self._get_distance_RF
        else:
            self.log.write("[WARNING] Provided orientation (%s) is not supported!\n"%orientation)
        # set parameters
        self.maxdist = self.isize + self.frac * self.stdev
        self.rlen = None
        # run alignments; parse alignments & get links
        #self.log.write(" Library stats...\n")
        handle = self._bwamem(fastqs)
        self.load_from_SAM(handle)
        self.get_links()
        
    # scaffolding part
    def _simplify(self):
        """Simplify scaffold graph by resolving circles"""
        # remove non-reciprocal connections
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                if not self.links[ref1][end1]:
                    continue
                for (ref2, end2), (links, gap) in self.links[ref1][end1].items(): 
                    if not (ref1, end1) in self.links[ref2][end2]:
                        # remove link and update counter
                        #self.log.write("Removing connection %s -> %s\n"%(str((ref1, end1)), str((ref2, end2))))
                        self.links[ref1][end1].pop((ref2, end2))
                        self.ilinks -= 1
                        
        ## make sure there are no circles

    def _populate_scaffold(self, links, pend, sid, scaffold, orientations, gaps, porientation):
        """Add links to scaffold representation. Experimental!

        !!!Plenty of work needed here!!!
        """
        # self.links[ref1][end1][(ref2, end2)] = (links, gap)
        #ref, end, links, gap = links
        ## there may be many connections now, but only one is processed so far!!!
        #print links
        for (ref, end), (links, gap) in links.iteritems():
            break
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

    def _get_scaffolds(self):
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
                    
class SyntenyGraph(Graph):
    """Graph class to represent scaffolds derived from synteny information"""
    def __init__(self, genome, reference, identity=0.51, overlap=0.66, threads=4, 
                 mingap=15, maxgap=0, printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "ReferenceGraph"
        self.log = log
        self.printlimit = printlimit
        # vars
        self.genome = genome
        # don't load reference genome - maybe we can avoid that
        self.reference = reference
        self.ref = self.reference
        # prepare storage
        self._init_storage(genome)
        # alignment options
        self.identity = identity
        self.overlap  = overlap
        self.threads  = threads
        # scaffolding options
        self.mingap  = mingap
        self._set_maxgap(maxgap)
        
    def _set_maxgap(self, maxgap=0, frac=0.01, min_maxgap=10000):
        """Set maxgap to 0.01 of assembly size, 0.01 of assembly size"""
        # set to 0.01 of assembly size
        if not maxgap:
            maxgap = int(round(frac * sum(self.contigs.itervalues())))
        # check if big enough
        if maxgap < min_maxgap:
            maxgap = min_maxgap
        # set variable
        self.maxgap = maxgap
        #self.log.write(" maxgap cut-off of %s bp\n"%self.maxgap)

    def _lastal(self):
        """Start LAST"""
        # build db
        if not os.path.isfile(self.ref+".suf"):
            os.system("lastdb %s %s" % (self.ref, self.ref))
        # run LAST
        args = ["lastal", "-T", "1", "-f", "TAB", "-P", str(self.threads), self.ref, self.genome]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=sys.stderr)        
        return proc.stdout
        
    def _get_hits(self):
        """Resolve & report scaffolds"""
        ## consider splitting into two functions
        ## to facilitate more input formats
        t2hits = {}
        t2size = {}
        q2hits = {}
        for l in self._lastal():
            if l.startswith('#'):
                continue
            # unpack
            (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
            (score, qstart, qalg, qsize, tstart, talg, tsize) = map(int, (score, qstart, qalg, qsize, tstart, talg, tsize))
            #get score, identity & overlap
            identity = 1.0 * score / qalg
            overlap  = 1.0 * qalg / qsize
            #filter by identity and overlap
            if identity < self.identity or overlap < self.overlap:
                continue
            # keep only best match to ref for each q
            if q in q2hits:
                pt, pscore = q2hits[q]
                # skip if worse score
                if score < q2hits[q][1]:
                    continue
                # update if better score
                q2hits[q] = (t, score)
                # remove previous q hit
                # normally should be at the end of matches to pt
                t2hits[pt].pop(-1)
            else:
                q2hits[q] = (t, score)
            # store
            if t not in t2hits:
                t2hits[t] = []
                t2size[t] = tsize            
            # For - strand matches, coordinates in the reverse complement of the 2nd sequence are used.
            strand = 0 # forward
            if qstrand == "-":
                strand = 1 # reverse
            qend, tend = qstart + qalg, tstart + talg
            t2hits[t].append((tstart, tend, q, qstart, qend, strand))
        
        # remove q that overlap too much on t
        # sort by r pos
        for t in t2hits:
            hits = t2hits[t]
            # remove hits overlapping too much
            t2hits[t] = []
            for tstart, tend, q, qstart, qend, strand in sorted(hits):
                # overlap with previous above threshold
                if t2hits[t] and t2hits[t][-1][1]-tstart > self.overlap*self.contigs[q]:
                    phit = t2hits[t][-1]
                    # do nothing if previous hit is better
                    if tend - tstart <= phit[1]-phit[0]:
                        continue
                    # remove previous match
                    t2hits[t].pop(-1)
                # add match only if first,
                # no overlap with previous or better than previous
                t2hits[t].append((tstart, tend, q, qstart, qend, strand))
        return t2hits, t2size

    def _estimate_gap(self, data, pdata):
        """Return estimated gap size"""
        # current start - previous end
        # this needs to be corrected by qstart and qend !!
        gap = data[0] - pdata[1]
        return gap
        
    def _get_scaffolds(self):
        """Resolve & report scaffolds"""
        self.logger("Aligning contigs on reference...\n")
        # get best ref-match to each contig
        t2hits, t2size = self._get_hits()

        # store scaffold structure
        ## [(contigs, orientations, gaps), ]
        self.scaffolds = []
        added = set()
        for t in sorted(t2size, key=lambda x: t2size[x], reverse=1):
            # skip one-to-one matches
            if len(t2hits[t])<2:
                continue
            # add empty scaffold
            scaffold, orientations, gaps = [], [], []
            #print t, t2size[t]
            for data in t2hits[t]:
                gap = 0
                tstart, tend, q, qstart, qend, strand = data
                # calculate gap only if scaffold has at least one element
                if scaffold:
                    gap = self._estimate_gap(data, pdata)
                    if gap > self.maxgap:
                        # add to scaffolds
                        if len(scaffold)>1:
                            self.scaffolds.append([scaffold, orientations, gaps])
                            added.update(scaffold)
                        # reset storage
                        scaffold, orientations, gaps = [], [], []
                    else:    
                        gaps.append(gap)
                # store contig & orientation
                scaffold.append(q)
                orientations.append(strand)
                # keep track of previous data
                pdata = data
                #print " %s %s-%s %s:%s-%s %s %s bp"%(len(self.scaffolds)+1, tstart, tend, self.shorter(q), qstart, qend, strand, gap)
                
            # add to scaffolds
            if len(scaffold)>1:
                self.scaffolds.append([scaffold, orientations, gaps])
                added.update(scaffold)
                
        # add missing
        for c in filter(lambda x: x not in added, self.contigs):
            self.scaffolds.append([(c,),(0,),[]])
            
def main():
    import argparse
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-f", "--fasta", required=1, 
                        help="assembly FASTA file")
    #parser.add_argument("-o", "--outdir",  default="redundans", 
    #                    help="output directory [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, 
                        help="max no. of threads to run [%(default)s]")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), 
                        help="output log to [stderr]")
    
    refo = parser.add_argument_group('Reference-based scaffolding options')
    refo.add_argument("-r", "--ref", "--reference", default='', 
                      help="reference FASTA file")
    refo.add_argument("--identity",        default=0.51, type=float,
                      help="min. identity [%(default)s]")
    refo.add_argument("--overlap",         default=0.66, type=float,
                      help="min. overlap  [%(default)s]")
    refo.add_argument("-g", "--maxgap",   default=0, type=int,
                      help="max. distance between adjacent contigs [0.01 * assembly_size]")
    
    scaf = parser.add_argument_group('Scaffolding options')
    scaf.add_argument("-i", "--fastq", nargs="+",
                      help="FASTQ PE/MP files")
    scaf.add_argument("-j", "--joins",  default=5, type=int, 
                      help="min pairs to join contigs [%(default)s]")
    scaf.add_argument("-a", "--linkratio", default=0.7, type=float,
                       help="max link ratio between two best contig pairs [%(default)s]")    
    scaf.add_argument("-l", "--load",  default=0.2, type=float, 
                      help="align subset of reads [%(default)s]")
    scaf.add_argument("-q", "--mapq",    default=10, type=int, 
                      help="min mapping quality [%(default)s]")
    # standard
    parser.add_argument("-v", dest="verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.10a')   
    o = parser.parse_args()
    if o.verbose:
        o.log.write("Options: %s\n"%str(o))

    # check logic
    if not o.ref and not o.fastq:
        sys.stderr.write("Provide FastQ files or reference genome (or both)!")
        sys.exit(1)
        
    # check if input files exists
    fnames = [o.fasta, o.ref]
    if o.fastq:
        fnames += o.fastq
    for fn in fnames: 
        if fn and not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)
            
    fasta = o.fasta
    log = o.log

    # perform PE/MP scaffolding if NGS provided
    if o.fastq:
        # NOT IMPLEMENTED
        # get library statistics
        #sys.stderr.write("PE/MP-based scaffolding is not implemented yet!\n"); sys.exit(1)
        # init
        s = ReadGraph(fasta, mapq=o.mapq, load=o.load, threads=o.threads, log=log)
        # add library
        s.add_library(o.fastq, isize=600, stdev=100, orientation="FR"); s.show()
        # add another library
        
        # save output
        s.save(out=open(fasta+".scaffolds.fa", "w"))
        
        # update fasta at the end
        fasta = fasta
    
    # perform referece-based scaffolding only if ref provided
    if o.ref:
        # init
        s = SyntenyGraph(fasta, o.ref, identity=o.identity, overlap=o.overlap, \
                         maxgap=o.maxgap, threads=o.threads, log=log)
        # save output
        s.save(out=open(fasta+".scaffolds.ref.fa", "w"))

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
    