#!/usr/bin/env python
desc="""Perform scaffolding of contigs using information from (in this order):

- long reads that are going to be assembled through miniasm and then scaffolded through synteny
- synteny to reference genome
More info at: https://github.com/lpryszcz/pyScaf
"""
epilog="""Author:
l.p.pryszcz@gmail.com
Warsaw, 12/03/2016

Updated to Python3 and repurposed by Diego Fuentes Palacios
Barcelona 08/18/2022
"""

import math, os, sys
from re import T
import time
import subprocess, resource, subprocess
from datetime import datetime
from FastaIndex import FastaIndex


def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values. 
    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @return - the percentile of the values
    From http://code.activestate.com/recipes/511478-finding-the-percentile-of-the-values/
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1

def median(N, key=lambda x:x):
    """median is 50th percentile."""
    return percentile(N, 0.5, key)

def mean(data):
    """Return the sample arithmetic mean of data.
    http://stackoverflow.com/a/27758326/632242
    """
    n = len(data)
    #if n < 1:
    #    raise ValueError('mean requires at least one data point')
    return sum(data)/float(n) 

def _ss(data):
    """Return sum of square deviations of sequence data."""
    c = mean(data)
    ss = sum((x-c)**2 for x in data)
    return ss

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    #if n < 2:
    #    raise ValueError('variance requires at least two data points')
    ss = _ss(data)
    pvar = ss/n # the population variance
    return pvar**0.5

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

    def logger(self, mssg="", decorate=1):
        """Logging function."""
        head = "\n%s"%("#"*50,)
        timestamp = "\n[%s]"% datetime.ctime(datetime.now())
        memusage  = "[%5i Mb] "%(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, )
        if self.log:
            if decorate:
                self.log.write("".join((head, timestamp, memusage, mssg)))
            else:
                self.log.write(mssg)
        
    def shorter(self, v, i=4, sep="_"):
        """Return shortened contig name"""
        return sep.join(v.split(sep)[:i])
        
    def __str__(self):
        """Produce string representation of the graph"""
        out = '%s\n%s contigs: %s\n\n%s links:\n' % (self.name, len(self.contigs), list(self.contigs.keys())[:self.printlimit], self.ilinks)
        i = 0
        for ref1 in sorted(self.contigs, key=lambda x: self.contigs[x], reverse=1):
            for end1 in range(2):
                '''if not self.links[ref1][end1]:
                    continue
                ref2, end2, links, gap = self.links[ref1][end1]'''
                for (ref2, end2), (links, gap) in list(self.links[ref1][end1].items()):
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
        self.sequences = FastaIndex(genome)
        self.seq = self.sequences
        # prepare storage
        self.contigs = {c: self.seq.id2stats[c][0] for c in self.seq}
        self.links   = {c: [{}, {}] for c in self.contigs}
        self.ilinks  = 0
        
    def _add_line(self, ref1, ref2, end1, end2, links, gap):
        """Add connection between contigs. """
        # store connection details
        self.links[ref1][end1][(ref2, end2)] = (links, gap)
        # update connection counter 
        self.ilinks += 1
        
    def _get_seqrecord(self, name, scaffold, orientations, gaps):
        """"Return name & seq for given scaffold"""
        # add empty gap at the end
        gaps.append(0)
        seqs = []
        for c, reverse, gap in zip(scaffold, orientations, gaps):
            seq = self.seq.get_sequence(c, reverse)
            # adjust gap size
            if gap and gap < self.mingap:
                strip = int(gap - self.mingap)
                seq = seq[:strip]
                gap = self.mingap
            seqs.append(str(seq)+"N"*gap)
        return name, "".join(seqs)

    def _format_fasta(self, name, seq, linebases=60):
        """Return formatted FastA entry"""
        seq = '\n'.join(seq[i:i+linebases] for i in range(0, len(seq), linebases))
        return ">%s\n%s\n"%(name, seq)
        
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
        added = set()
        for i, (scaffold, orientations, gaps) in enumerate(self.scaffolds, 1):
            # scaffold00001
            name = str("scaffold%5i"%i).replace(' ','0')
            # save scaffold
            name, seq = self._get_seqrecord(name, scaffold, orientations, gaps)
            out.write(self._format_fasta(name, seq))
            # report info
            log.write(logline%(name, len(seq), len(scaffold), " ".join(scaffold),
                               " ".join(map(str, (int(x) for x in orientations))), # x may be bool!
                               " ".join(map(str, (x for x in gaps)))))
            totsize += len(seq)
            # store info about added
            for c in scaffold:
                added.add(c)
        # add contigs that were not scaffolded
        for contig in [x for x in self.seq if x not in added]:
            name = "unscaffolded.%s"%contig
            seq = self.seq.get_sequence(contig)
            out.write(self._format_fasta(name, seq))
        # close output & loge
        out.close()
        log.close()
        self.logger(" %s bp in %s scaffolds. Details in %s\n"%(totsize, len(self.scaffolds), log.name), 0)
        self.logger("Scaffolds saved to: %s\n"%out.name, 0)

    def save_dotplot(self, query, readstdin=False): # open('/dev/null','w')): #
        """Produce query to reference dotplot"""
        ### TBD: generate dotplot only if less than several hundreds of seqs in query
        outfn = "%s.%s"%(query, self.dotplot)
        self.logger("Saving dotplots to: %s\n"%outfn) 
        args = ["last-dotplot", "-", outfn]
        if readstdin:
            proc = subprocess.Popen(args, stdin=subprocess.PIPE, stderr=self.log)
        else:
            proc = subprocess.Popen(args, stdin=self._lastal(query), stderr=self.log)
        return proc

    def _minimap2(self, queries=[], index="4G"):
        """Minimap2 for long read alignment using presets"""

        #Provide necessary arguments:
        if not queries:
            queries = self.fastq
        # convert filename to list of filenames
        if type(queries) is str:
            queries = [queries, ]
        args0 = ["cat", ] + queries
        if queries[0].endswith('.gz'):
            args0[0] = "zcat"
        proc0 = subprocess.Popen(args0, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        args1 = ["minimap2", "-x", str(self.preset), "-t", str(self.threads),  "-I", index, "--cs=long", self.ref, "-"]
        proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stdin=proc0.stdout, stderr=subprocess.DEVNULL)
        args2 = ["k8-Linux", "bin/minimap2/misc/paftools.js", "view", "-f", "maf", "-"]
        proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stdin=proc1.stdout, stderr=subprocess.DEVNULL)
        #Added maf converter from LAST to keep same format
        args3 = ["maf-convert", "tab", "-"]
        proc3 = subprocess.Popen(args3, stdout=subprocess.PIPE, stdin=proc2.stdout, stderr=subprocess.DEVNULL)
        return proc3.stdout

    def _lastal(self, queries=[]):
        """Start LAST in local mode and with FastQ input (-Q 1)."""
        # build db
        if not os.path.isfile(self.ref+".suf"):
            os.system("lastdb %s %s" % (self.ref, self.ref))
        # decide on input
        if not queries:
            queries = self.fastq
        # convert filename to list of filenames
        if type(queries) is str:
            queries = [queries, ]
        # run LAST aligner, split and maf-convert in pipe
        seqformati = 1
        args0 = ["cat", ] + queries
        if queries[0].endswith('.gz'):
            args0[0] = "zcat"
            seqformati += 1
        # deduce sequence format
        seqformat = queries[0].split(".")[-seqformati].lower()
        if seqformat in ("fasta", "fa"):
            seqformatcode = "0" # FASTA
        elif seqformat in ("fastq", "fq"):
            seqformatcode = "1" # FastQ
        else:
            self.logger("[WARNING] Unsupported sequence format `%s` in %s\n"%(seqformat, queries[0]), 0)
            sys.exit(1)
        # combine processes
        proc0 = subprocess.Popen(args0, stdout=subprocess.PIPE, stderr=sys.stderr)
        args1 = ["lastal", "-Q", seqformatcode, "-P", str(self.threads), self.ref, "-"]
        proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stdin=proc0.stdout, stderr=sys.stderr)
        args2 = ["last-split", "-"]
        proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stdin=proc1.stdout, stderr=sys.stderr)
        args3 = ["maf-convert", "tab", "-"]
        proc3 = subprocess.Popen(args3, stdout=subprocess.PIPE, stdin=proc2.stdout, stderr=sys.stderr)
        #print " ".join(args1), queries
        return proc3.stdout

    def _lastal_global(self, query=''):
        """Start LAST in overlap mode. Slightly faster than local mode,
        but much less sensitive."""
        # build db
        if not os.path.isfile(self.ref+".suf"):
            os.system("lastdb %s %s" % (self.ref, self.ref))
        # select query
        if not query:
            query = self.genome
        # run LAST
        args = ["lastal", "-T", "1", "-f", "TAB", "-P", str(self.threads), self.ref, query]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=sys.stderr)
        return proc.stdout

    ###
    # SCAFFOLDING PART
    ###

    #Commented the recursive version of the function
    #def _populate_scaffold(self, links, pend, sid, scaffold, orientations, gaps, porientation):
    #    """Add links to scaffold representation. Experimental!
    #    !!!Plenty of work needed here!!!
    #    """
    #    # self.links[ref1][end1][(ref2, end2)] = (links, gap)
    #    #ref, end, links, gap = links
    #    ## there may be many connections now, but only one is processed so far!!!
    #    #print links
    #    links_sorted = sorted(iter(links.items()), key=lambda x: x[1][0], reverse=1)
    #    if len(links_sorted)>1:
    #        self.logger(" multi connections: %s %s\n"%(scaffold, links), 0)
    #    (ref, end), (links, gap) = links_sorted[0]
    #    #for (ref, end), (links, gap) in links.iteritems(): break
    #    # skip if already added
    #    if ref in self.contig2scaffold:
    #        return scaffold, orientations, gaps, porientation
    #    # get orientation - get forward/reverse-complement signal by XOR
    #    ## if previous is 1 (end) & current is (0) start & orientation is 0 (forward) --> keep forward orientation
    #    orientation = (pend != end) != porientation
    #    # store at the end if previous contig was F and pend 1
    #    if porientation != pend:
    #        scaffold.append(ref)
    #        orientations.append(not orientation)
    #        gaps.append(gap)
    #    else:
    #        scaffold.insert(0, ref)
    #        orientations.insert(0, not orientation)
    #        gaps.insert(0, gap)
    #    # update contigs2scaffold info
    #    self.contig2scaffold[ref] = sid
    #    # populate further connections from another end
    #    links = self.links[ref][abs(end-1)]
    #    # skip if not links
    #    if not links:
    #        return scaffold, orientations, gaps, orientation
    #    # populate further connections from another end
    #    return self._populate_scaffold(links, end, sid, scaffold, orientations, gaps, orientation)
    
    def _populate_scaffold(self, links, pend, sid, scaffold, orientations, gaps, porientation):
        """Add links to scaffold representation.
        
        Iterative version
        """

        links_sorted = sorted(iter(links.items()), key=lambda x: x[1][0], reverse=1)

        while len(links_sorted) > 0:
            if len(links_sorted) > 1:
                self.logger(" multi connections: %s %s\n"%(scaffold, links), 0)

            (ref, end), (links, gap) = links_sorted[0]

            # skip if already added
            if ref in self.contig2scaffold:
                return scaffold, orientations, gaps, porientation

            # get orientation - get forward/reverse-complement signal by XOR
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

            links_sorted = sorted(iter(links.items()), key=lambda x: x[1][0], reverse=1)

        return scaffold, orientations, gaps, orientation

           
class SyntenyGraph(Graph):
    """Graph class to represent scaffolds derived from synteny information"""
    def __init__(self, genome, reference, identity=0.51, overlap=0.66, preset="asm10", index="4G", useminimap2=0, norearrangements=0, 
                 threads=4, mingap=15, maxgap=0, nofilter_overlaps=0, 
                 dotplot="png", printlimit=10, log=sys.stderr, uselongreads=0, preset_long="ava-ont", fastq=None):
        """Construct a graph with the given vertices & features"""
        self.name = "ReferenceGraph"
        self.log = log
        self.printlimit = printlimit
        self.threads  = threads
        # vars
        self.genome = genome
        # don't load reference genome - maybe we can avoid that
        if uselongreads and fastq:
            self.fastq = fastq
            try:
                self.reference = _miniasm(os.path.join(os.path.dirname(self.genome), "miniasm.assm.fa"), self.fastq, self.threads, preset_long)
            except ValueError as e:
                sys.stderr.write("\n[ERROR] something went wrong with the provided data, resulting in error : %s"%e)
        else:
            self.reference = reference
        self.ref = self.reference
        # prepare storage
        self._init_storage(genome)
        # alignment options
        self.identity = identity
        self.overlap  = overlap
        self.dotplot  = dotplot
        self.useminimap2 = useminimap2
        self.index=index
        self.preset = preset
        # 0-local alignment; 
        if norearrangements:
            # 1-global/overlap - simpler and faster
            self._get_hits = self._get_hits_global
            self._lastal   = self._lastal_global # needed by save_dotplot
        # scaffolding options
        self.mingap  = mingap
        self._set_maxgap(maxgap)
        self.nofilter_overlaps = nofilter_overlaps
        
    def _set_maxgap(self, maxgap=0, frac=0.01, min_maxgap=10000):
        """Set maxgap to 0.01 of assembly size, 0.01 of assembly size"""
        # set to 0.01 of assembly size
        if not maxgap:
            maxgap = int(round(frac * sum(self.contigs.values())))
        # check if big enough
        if maxgap < min_maxgap:
            maxgap = min_maxgap
        # set variable
        self.maxgap = maxgap
        self.logger(" maxgap cut-off of %s bp\n"%self.maxgap, 0)

    def _clean_overlaps_on_reference(self, _t2hits):
        """Remove hits that overlap on reference too much"""
        t2hits = {}
        for t, hits in _t2hits.items():
            t2hits[t] = []
            # remove hits overlapping too much # sort by r pos
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
        return t2hits

    def _get_hits_global(self):
        """Resolve & report scaffolds"""
        dotplot = None
        if self.dotplot:
            dotplot = self.save_dotplot(self.genome, readstdin=True)
        ## consider splitting into two functions
        ## to facilitate more input formats
        t2hits = {}
        t2size = {}
        q2hits = {}
        if self.useminimap2:
            handle = self._minimap2(index=self.index)
        else:
            handle = self._lastal_global()
        for le in handle:
            l = le.decode("utf-8")
            if dotplot:
                try:
                    dotplot.stdin.write(l)
                except:
                    self.logger("[WARNING] dotplot generation failed!\n", 0)
                    dotplot = None
            if l.startswith('#'):
                continue
            # unpack
            (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
            (score, qstart, qalg, qsize, tstart, talg, tsize) = list(map(int, (score, qstart, qalg, qsize, tstart, talg, tsize)))
            #get score, identity & overlap
            identity = 1.0 * score / qalg
            overlap  = 1.0 * qalg / qsize
            #filter by identity and overlap.
            if identity < self.identity or overlap < self.overlap:
                continue
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
                # reverse
                strand = 1 
                qstart = qsize - qstart - qalg
            qend, tend = qstart + qalg, tstart + talg
            t2hits[t].append((tstart, tend, q, qstart, qend, strand))

        # remove q that overlap too much on t
        t2hits = self._clean_overlaps_on_reference(t2hits)
        return t2hits, t2size

    def _get_best_global_match(self, hits):
        """Return best, longest match for given q-t pair"""
        newHits = [[hits[0]]]
        for hit in hits[1:]:
            # break synteny if too large gap
            #if hit[2]=='scaffold22|size195699': print hit, hit[0]-newHits[-1][-1][0]
            if hit[0]-newHits[-1][-1][0] > self.maxgap:
                newHits.append([])
            newHits[-1].append(hit)
        # sort by the longest consecutive alg
        newHits = sorted(newHits, key=lambda x: sum(y[1] for y in x), reverse=1)
        return newHits[0]
        
    def _calculate_global_algs(self, t2hits):
        """Return simplified, global alignments"""
        t2hits2 = {}
        for t in t2hits:
            if t not in t2hits2:
                t2hits2[t] = []
            for q in t2hits[t]:
                # sort by r pos
                hits = self._get_best_global_match(sorted(t2hits[t][q]))
                #print(hits)
                #get score, identity & overlap
                score = sum(x[-1] for x in hits)
                qalg  = sum(x[4] for x in hits)
                identity = 1.0 * score / qalg # this local identity in hits ignoring gaps
                overlap  = 1.0 * qalg / self.contigs[q]

                #filter by identity and overlap. 
                if identity < self.identity or overlap < self.overlap:
                    continue

                ## this needs work and bulletproofing!!!
                tstart = hits[0][0]
                tend   = hits[-1][0] + hits[-1][1]
                qstart = hits[0][3]
                qend   = hits[-1][3] + hits[-1][4]
                # and report rearrangements
                
                # get strand correctly - by majority voting
                strand = int(round(1.0 * sum(x[4]*x[5] for x in hits) / qalg))
                t2hits2[t].append((tstart, tend, q, qstart, qend, strand))
                
        return t2hits2
        
    def _get_hits(self):
        """Resolve & report scaffolds"""
        dotplot = None
        if self.dotplot:
            dotplot = self.save_dotplot(self.genome, readstdin=True)
        ## consider splitting into two functions
        ## to facilitate more input formats
        t2hits, t2size = {}, {}
        index=self.index
        q2hits = {}
        if self.useminimap2:
            handle = self._minimap2(self.genome, index=index)
        else:
            handle = self._lastal(self.genome)
        for le in handle:
            # strip leading spaces
            l = le.decode("utf-8")
            l = l.lstrip()
            if dotplot:
                try:
                    dotplot.stdin.write(l)
                except:
                    self.logger("[WARNING] dotplot generation failed!\n", 0)
                    dotplot = None
            if l.startswith('#'):
                continue
            # unpack
            (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
            (score, qstart, qalg, qsize, tstart, talg, tsize) = list(map(int, (score, qstart, qalg, qsize, tstart, talg, tsize)))
            # prepare storage
            if t not in t2hits:
                t2hits[t] = {q: []}
                t2size[t] = tsize
            elif q not in t2hits[t]:
                t2hits[t][q] = []
            if q not in q2hits:
                q2hits[q] = []
            # For - strand matches, coordinates in the reverse complement of the 2nd sequence are used.
            strand = 0 # forward
            if qstrand == "-":
                # reverse -> adjust start
                strand = 1 
                qstart = qsize - qstart - qalg
            t2hits[t][q].append((tstart, talg, q, qstart, qalg, strand, score))
            q2hits[q].append((qstart, qalg, strand, t, tstart, talg))
        
        #tstart, tend, q, qstart, qend, strand

        # get simplified global alignments

        t2hits = self._calculate_global_algs(t2hits)
        #for t, hits in t2hits.iteritems(): print t, hits

        # clean overlaps on reference
        if not self.nofilter_overlaps:
            t2hits = self._clean_overlaps_on_reference(t2hits)

        #print "after clean-up"
        #for t, hits in t2hits.iteritems(): print t, hits
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
                
            # add to scaffolds
            if len(scaffold)>1:
                self.scaffolds.append([scaffold, orientations, gaps])
                added.update(scaffold)
                
        # add missing
        for c in [x for x in self.contigs if x not in added]:
            self.scaffolds.append([(c,),(0,),[]])

def _check_executable(cmd):
    """Check if executable exists."""
    p = subprocess.Popen("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return "".join(p.stdout.readlines())

def _check_dependencies(dependencies):
    """Return error if wrong software version"""
    warning = 0
    info = "[WARNING] Old version of %s: %s. Update to version %s+!\n"
    for cmd, version in dependencies.items():
        out = _check_executable(cmd)
        if "not found" in out:
            warning = 1
            sys.stderr.write("[ERROR] %s\n"%out)
        elif version:
            out = subprocess.getoutput("%s --version"%cmd)
            curver = out.split()[-1]
            if not curver.isdigit():
                warning = 1
                sys.stderr.write("[WARNING] Problem checking %s version: %s\n"%(cmd, out))
            elif int(curver)<version:
                warning = 1
                sys.stderr.write(info%(cmd, curver, version))
                
    message = "Make sure you have installed all dependencies from https://github.com/lpryszcz/pyScaf#dependencies !"
    if warning:
        sys.stderr.write("\n%s\n\n"%message)
        sys.exit(1)

#TO DO#
#Here we define the function to generate a new reference based on long reads

def overlap(longreads, threads, preset = "ava-ont"):

    if type(longreads) is not str:
                longreads = str(longreads[0])
    
    args1 = ["minimap2", "-x", preset, "-t", str(threads), longreads, longreads, "-"]
    #print(args1)
    sys.stderr.write(str(args1))
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=sys.stderr)
    #proc1.communicate()

    return proc1.stdout

def _miniasm(outref, longreads, threads, preset = "ava-ont"):
    """Minimap2 for long read alignment using presets"""

    if type(longreads) is not str:
        longreads = str(longreads[0])

    #Here we define locations to generate intermediate files
    basedir = os.path.dirname(outref)

    sam = os.path.join(basedir, "minimap.sam")
    args1 = "minimap2 -x %s -t %s %s %s > %s"%(preset, str(threads), longreads, longreads, sam)
    proc1 = subprocess.call(args1, shell=True, stderr=subprocess.DEVNULL)
    time.sleep(10)
    

    gfa= os.path.join(basedir, "miniasm.gfa")
    sys.stderr.write("\nRunning miniasm to generate a primary assembly: %s"%outref)
    args2 = "miniasm -e2 -n2 -f %s %s > %s"%(longreads, sam, gfa)
    proc2 = subprocess.call(args2, shell=True, stderr=subprocess.DEVNULL)

    with open(outref, "wb") as file:
        args3 =  "awk \'/^S/{print \">\"$2\"\\n\"$3}\' %s | fold > %s"%(gfa, outref)
        proc3 = subprocess.call(args3, shell=True, stdout=file)

    return outref

"""def get_proc_status(pid):
    ""Get the status of the process which has the specified process id.""

    proc_status = None
    try:
        proc_status = psutil.Process(pid).status()
    except psutil.NoSuchProcess as no_proc_exc:
        print(no_proc_exc)
    except psutil.ZombieProcess as zombie_proc_exc:  
        # For Python 3.0+ in Linux (and MacOS?).
        print(zombie_proc_exc)
    return proc_status"""

class LongReadGraph(Graph):
    """Graph class to represent scaffolds derived from synteny after using long reads to assemble a reference"""
    def __init__(self, genome, fastq, identity=0.51, overlap=0.66, preset="map-ont", index="4G", useminimap2=0, norearrangements=0, 
                 threads=4, dotplot="png", mingap=15, maxgap=0, printlimit=10, log=sys.stderr):
        """Construct a graph with the given vertices & features"""
        self.name = "ReferenceGraph"
        self.log = log
        self.printlimit = printlimit
        # vars
        self.genome = genome
        self.ref = self.genome
        self.fastq = fastq
        # prepare storage
        self._init_storage(genome)
        # alignment options

        self.preset = preset
        self.useminimap2 = useminimap2
        self.index=index
        self.identity = identity
        self.overlap  = overlap
        self.threads  = threads
        self.dotplot  = dotplot
        # scaffolding options
        self.mingap  = mingap
        self._set_maxgap(maxgap)
        self.maxoverhang = 0.1
        # store long links
        self.longlinks = {c: [{}, {}] for c in self.contigs}
        self.ilonglinks = 0
        
    def _set_maxgap(self, maxgap=0, frac=0.01, min_maxgap=10000):
        """Set maxgap to 0.01 of assembly size, 0.01 of assembly size"""
        # set to 0.01 of assembly size
        if not maxgap:
            maxgap = int(round(frac * sum(self.contigs.values())))
        # check if big enough
        if maxgap < min_maxgap:
            maxgap = min_maxgap
        # set variable
        self.maxgap = maxgap
        self.logger(" maxgap cut-off of %s bp\n"%self.maxgap, 0)
        
    def _get_hits(self):
        """Resolve & report scaffolds"""
        # maybe instead of last-split, get longest, non-overlapping matches here
        q2hits, q2size = {}, {}

        #Add minimap processing
        if self.useminimap2:
            handle = self._minimap2
        else:
            handle = self._lastal
        #for le in self._lastal(): # open('contigs.reduced.fa.tab7'): #
        for le in handle():
            # strip leading spaces
            l = le.decode("utf-8")
            l = l.lstrip()
            
            if l.startswith('#') or l.startswith('@'):
                continue

            # unpack
            (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
            (score, qstart, qalg, qsize, tstart, talg, tsize) = list(map(int, (score, qstart, qalg, qsize, tstart, talg, tsize)))
            if q not in q2hits:
                q2hits[q] = []
                q2size[q]  = qsize
            # For - strand matches, coordinates in the reverse complement of the 2nd sequence are used.
            strand = 0 # forward
            if qstrand == "-":
                # reverse -> adjust start
                strand = 1 
                qstart = qsize - qstart - qalg
            q2hits[q].append((qstart, qalg, strand, t, tstart, talg))

        return q2hits, q2size
        
    def _get_best_global_match(self, hits):
        """Return best, longest match for given q-t pair"""
        newHits = [[hits[0]]]
        for hit in hits[1:]:
            # break synteny if too large gap
            #if hit[2]=='scaffold22|size195699': print hit, hit[0]-newHits[-1][-1][0]
            if hit[0]-newHits[-1][-1][0] > self.maxgap:
                newHits.append([])
            newHits[-1].append(hit)
        # sort by the longest consecutive alg
        newHits = sorted(newHits, key=lambda x: sum(y[1] for y in x), reverse=1)
        return newHits[0]
        
    def _add_longread_line(self, ref1, ref2, end1, end2, gap):
        """Add connection between contigs from long reads."""
        if (ref2, end2) not in self.longlinks[ref1][end1]:
            self.longlinks[ref1][end1][(ref2, end2)] = []
            self.longlinks[ref2][end2][(ref1, end1)] = []
        # store connection details
        self.longlinks[ref1][end1][(ref2, end2)].append(gap)
        self.longlinks[ref2][end2][(ref1, end1)].append(gap)
        # update connection counter 
        self.ilonglinks += 1

    def _contained_hits(self, hits):
        """Return True if cointaned hits ie A, B, A"""
        added = set([hits[0][3]])
        for i, (qstart, qalg, strand, t, tstart, talg) in enumerate(hits[1:]):
            if t!=hits[i][3] and t in added:
                return True
            added.add(t)
        
    def _hits2longlinks(self, q2hits, q2size, score=0):
        """Filter alignments and populate links.
        Skip: 
        - long reads aligning to only one contig
        - check read overlap
        - mixed alignments ie c1, c2, c1, c2
        - clearly wrong alignment ie c1s - c2s
        """
        for q in list(q2hits.keys()):
            qsize, hits = q2size[q], q2hits[q]
            
            # check if more than 2 contigs aligned
            targets = set(t for qstart, qalg, strand, t, tstart, talg in hits)
            if len(targets)<2:
                continue
                
            # check if enough overlap
            aligned = sum(qalg for qstart, qalg, strand, t, tstart, talg in hits)
            if aligned < self.overlap*qsize:
                continue
            
            hits.sort()

            # check if contained hit (ie A, B, A)
            if self._contained_hits(hits):
                self.logger("contained hits" + "\n".join("\t".join(map(str, (qsize, qstart, qalg, strand, "_".join(t.split("_")[:2]), tstart, talg, self.contigs[t]))) for qstart, qalg, strand, t, tstart, talg in hits) + "\n", 0)
                continue
            
            # combine hits for the same pair
            t2hits = {}
            for qstart, qalg, strand, t, tstart, talg in hits:
                if t not in t2hits:
                    t2hits[t] = []
                t2hits[t].append((tstart, talg, q, qstart, qalg, strand, score))

            uhits = []
            for t in t2hits:
                # get best target hit
                hits = self._get_best_global_match(sorted(t2hits[t]))
                
                # and report rearrangements
                
                # get strand correctly - by majority voting
                talg   = sum(x[1] for x in hits)
                strand = int(round(1.0 * sum(x[1]*x[5] for x in hits) / talg))
                #if strand>1: self.logger("[ERROR] Wrong strand from majority voting: %s %s\n"%(strand, str(hits)), 0)

                if strand:
                    qstart = hits[-1][3]
                    qend   = hits[0][3] + hits[0][4]
                else:
                    qstart = hits[0][3]
                    qend   = hits[-1][3] + hits[-1][4]
                # get global start & end
                tstart = hits[0][0]
                tend   = hits[-1][0] + hits[-1][1]
                uhits.append((q, qstart, qend, strand, t, tstart, tend))

            uhits = sorted(uhits, key=lambda x: x[1])
            self.logger("\t".join(map(str, uhits[0]))+"\n", 0)
            for i, (q, qstart, qend, strand, t, tstart, tend) in enumerate(uhits[1:]):
                dist = qstart - uhits[i][2]
                c1, c2 = t, uhits[i][4]
                # get contig orientation
                end1, end2 = 0, 1
                pos1 = tstart
                if strand:
                    end1 = 1
                    pos1 = self.contigs[c1] - tend
                pos2 = self.contigs[c2] - uhits[i][6]
                if uhits[i][3]: 
                    end2 = 0
                    pos2 = uhits[i][5]
                # calculate gap
                gap = dist - pos1 + pos2
                overhang = gap - dist
                self.logger("\t".join(map(str, (q, qstart, qend, strand, t, tstart, tend, gap)))+"\n", 0)
                self.logger(" %s:%s %s -> %s:%s %s  %s  %s bp\n"%(c1, pos1, end1, c2, pos2, end2, dist, gap), 0)
                # skip if too big overhang on contig edge
                if gap > self.maxgap or overhang > self.maxoverhang*(self.contigs[c1]+self.contigs[c2]):
                    self.logger(" too big contig overhang (%s) or gap (%s)!\n\n"%(overhang, gap), 0)
                    continue
                self._add_longread_line(c1, c2, end1, end2, gap)
            self.logger("\n", 0)
                
        # get links
        self._get_links()
        self.logger("%s %s\n"%(self.ilonglinks, self.ilinks), 0)
                
    def _get_links(self):
        """Combine longlinks into links"""
        for ref1 in self.longlinks:
            for end1, data in enumerate(self.longlinks[ref1]):
                for (ref2, end2), gaps in data.items():
                    if ref1 > ref2:
                        continue
                    links = len(gaps)
                    gap = int(round(mean(gaps)))
                    gapstd = pstdev(gaps)
                    if gap and gapstd>100 and gapstd / gap > 0.25:
                        self.logger("highly variable gap size at %s %s -> %s %s: %s +- %.f %s\n"%(ref1, end1, ref2, end2, gap, gapstd, str(gaps)), 0)
                    self._add_line(ref1, ref2, end1, end2, links, gap)
                    self._add_line(ref2, ref1, end2, end1, links, gap)
                
    def _get_scaffolds(self):
        """Resolve & report scaffolds"""
        self.logger("Aligning long reads on contigs...\n")
        # get best ref-match to each contig
        q2hits, q2size = self._get_hits()

        # get simplified global alignments
        #self.longlinks = {c: [{}, {}] for c in self.contigs}
        self._hits2longlinks(q2hits, q2size)
        
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
            