#!/usr/bin/env python
desc="""Filter QSEQ/FastQ reads. Store output as FastQ.
Reads are clipped at first undetermined base (. in qseq or N in fastq)
and at first base having qual below -q.
Reads (and their pairs if -p) not passing filtering are discarded.
Orphaned reads may be store optionally (-u).
"""
epilog="""Author:
Leszek Pryszcz/l.p.pryszcz@gmail.com
Modified by Salvador Capella-Gutierrez/salcagu@gmail.com

Barcelona, 01/01/2014
"""

"""
Fixes:
-0.23:
--Stable version after different new options

-0.22:
--Added Phix removal possibility

-0.21:
--Added Encoding quality check
--Added output file prefixes

-0.2
--preformance:
---multiprocessing tested, but very little improvement
---zcat subprocess speeds up by 30% (1M of paired 100bp reads processed in 40s in neptune)
--tested biopython, but ~4x slower than raw fq parser

-0.1
--wrong sep in _clipSeq solved

-0.11
--output always PHRED+33 quals (Sanger, CASAVA1.8+)
--include reads with '.' bases -> 'N'
"""

import argparse, locale, subprocess, gzip, os, sys
from datetime import datetime
#locale.setlocale(locale.LC_ALL, 'en_US.utf8')

def checkQualityEncoding(inFile, number_reads, qual64offset, qseq):

    iFile =  gzip.open(inFile, "rb") if inFile.endswith(".gz") else \
      open(inFile, "rU")

    ## We analyze the first 1000 reads to determine the format
    counter, qualities = number_reads, []
    while counter > 0:
        read = (qseqparser(iFile) if qseq else fqparser(iFile)).next()
        if read == None:
            continue
        name, seq, quals = read
        ## We only consider quality scores between [33, 126] to compute whether
        ## input quality encoding parameter have been properly selected or not
        qualities += [ord(q) for q in quals if ord(q) < 127 and ord(q) > 32]
        counter -= 1
    iFile.close()

    minQ, maxQ = min(qualities), max(qualities)

    ## for PHRED+33, values will be 33-126 and for PHRED+64: 64/59-126.
    if qual64offset and minQ < 59:
        msg = ("\nERROR: Check quality encoding. Selected [PHRED+64]. MinQ: %d "
            + "- MaxQ: %d\n") % (minQ, maxQ)
        return False, msg

    if not qual64offset and minQ > 58 and (maxQ - 68) < minQ:
        msg = ("\nERROR: Check quality encoding. Selected [PHRED+33]. MinQ: %d "
            + "- MaxQ: %d\n") % (minQ, maxQ)
        return False, msg

    return True, ""

def qseqparser(handle, limit=0):
    """Parse QSEQ fromat and yield name, sequence and qualities."""

    for l in handle:
        if limit and i>limit:
            break
        ## Example SOLEXA read
        ## SOLEXA 90403 4 1 23 1566 0 1 ACCGCTCTCGTGCTCGTCGCTGCGTTGAGGCTTGCG `aaaaa```aZa^`]a``a``a]a^`a\Y^`^^]V` 1
        qseq_element = l[:-1].split('\t')

        if len(qseq_element) != 11 or qseq_element[-1] != '1':
            yield
            continue
        ## Format input read following FASTQ standard
        name = '@%s:%s:%s:%s:%s#%s/%s' % (qseq_element[0], qseq_element[2], \
            qseq_element[3], qseq_element[4], qseq_element[5], qseq_element[6],\
            qseq_element[7])
        seq, quals = qseq_element[8].replace(".", "N"), qseq_element[9]

        yield name, seq, quals

def fqparser(handle, limit=0):
    """Parse FASTQ format and yield name, sequence and qualities."""
    fqlist = []
    for i, l in enumerate(handle, 1):
        if limit and i>4*limit:
            break
        fqlist.append(l[:-1])
        if len(fqlist) != 4:
            continue
        #unload & reset fqlist
        name, seq, sep, quals = fqlist
        fqlist = []
        yield name, seq, quals

def _clipSeq(seq, quals, sep='.'):
    """Clip sequence at first sep base (. or N). Clip quals accordingly."""
    if sep in seq:
        pos = seq.index(sep)
        seq, quals = seq[:pos], quals[:pos]
    return seq, quals

def rawtrimmer(infile, minlen, maxlen, limit, minqual, \
               qual64offset, qseq, stripHeaders, outformat, \
               pi, pair="", phixReads=[], logFile=sys.stderr):
    """Single process implementation of rawtrimmer.
    Open zcat subprocess and read from stdin."""
    handle = infile
    if infile.name.endswith('.gz'):
        zcat = subprocess.Popen(['zcat', infile.name], bufsize=-1, \
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        handle = zcat.stdout

    ## Get parser
    parser = qseqparser(handle, limit) if qseq else fqparser(handle, limit)

    ## Process entries
    phixs = 0
    for seqi, read in enumerate(parser, pi+1):

        ## Discard any unvalid read
        if not read:
            yield
            continue
        name, seq, quals = read

        ## Clip seq & quals @ N ( unknown base )
        seq, quals = _clipSeq(seq, quals, 'N')
        #check if correct length
        if not seq or len(seq) < minlen:
            yield
            continue

        ## Return PHRED+33 quals (Sanger encoding)
        if qual64offset:
            quals=''.join([chr(ord(q)-31) for q in quals])
        #cut sequence & quals @ quality
        if minqual:
            for pos, qual in enumerate(quals):
                if ord(qual)-33 < minqual:
                    seq = seq[:pos]
                    break
            #check if correct length
            if not seq or len(seq) < minlen:
                yield
                continue
            #clip seq and qual
            quals = quals[:len(seq)]

        ## Check wether current read is not part of PhiX Sequence
        if phixReads and seq[:minqual] in phixReads and phixSeq.find(seq) != -1:
            phixs += 1
            logFile.write(("%s Detected a read [%d] mapping to PhiX sequence"
                + "\r") % (locale.format("%d", i, grouping=True), phixs))
            logFile.flush()
            yield
            continue

        # hard-trim read
        if maxlen:
            seq, quals = seq[:maxlen], quals[:maxlen]
            
        ## Define fastQ line
        if stripHeaders:
            name = "@%s%s" % (seqi, pair)
        if outformat == "fasta":
            fastq = '>%s\n%s\n' % (name[1:], seq)
        else:
            fastq = '%s\n%s\n+\n%s\n' % (name, seq, quals)
        yield fastq

def filter_paired(fpair, outfiles, minlen, maxlen, limit, minqual, \
                  qual64offset=0, qseq=0, stripHeaders=1, outformat='fastq', \
                  pi=0, logFile=0):
    """Filter paired reads."""
    inF, inR = fpair
    outF, outR, outCombined, outUnpaired = outfiles

    ## Define parsers rawtrimmer fqtrimmer
    fqparser1 = rawtrimmer(inF, minlen, maxlen, limit, minqual, qual64offset, qseq, \
                           stripHeaders, outformat, pi)
    fqparser2 = rawtrimmer(inR, minlen, maxlen, limit, minqual, qual64offset, qseq, \
                           stripHeaders, outformat, pi)

    ## Process
    i = both = fori = revi = filtered = 0
    for i, rec1 in enumerate(fqparser1, pi+1):
        # will crash if len(fq1) > len(fq2)
        rec2 = fqparser2.next()
        if rec1 and rec2:
            #store paired output
            if outF:
                outF.write(rec1)
                outR.write(rec2)
            #store combined output
            if outCombined:
                outCombined.write(rec1+rec2)
            both += 1
        elif outUnpaired and rec1:
            ## Store F read if R didn't pass filtering and orphans requested
            fori+=1
            outUnpaired.write(rec1)
        elif outUnpaired and rec2:
            ## Store R read if F didn't pass filtering and orphans requested
            revi+=1
            outUnpaired.write(rec2)
        else:
            #count skipped reads
            filtered += 1
        #print stats
        if logFile and not i % 10e3:
            info = "%9s processed [" % (locale.format("%d", i, grouping=True))
            info += '%6.2f%s ok] Both passed: %s Orphans F/R: %s/%s      \r'\
                % ((i-filtered)*100.0/i, '%', both, fori, revi)
            logFile.write(info)
            logFile.flush()

    if logFile:
        info = "%9s processed [" % (locale.format("%d", i, grouping=True))
        info += '%6.2f%s ok] Both passed: %s Orphans F/R: %s/%s      \r' \
                % ((i-filtered)*100.0/i, '%', both, fori, revi)
        logFile.write(info)
        logFile.flush()

    return i, filtered, fori+revi

def process_paired(inputs, qseq, outdir, outprefix, unpaired, minlen, maxlen, limit, minqual, \
                   noSeparate, combined, qual64offset, replace, \
                   stripHeaders, fasta, verbose):
    """Process paired libraries."""

    ## Define output fnames
    fnend = outformat = 'fasta' if fasta else 'fastq'
    prefix = ("%sq%s_l%s") % (outprefix, minqual, minlen)

    outfnF     = os.path.join(outdir, '%s.1.%s'        % (prefix, fnend))
    outfnR     = os.path.join(outdir, '%s.2.%s'        % (prefix, fnend))
    unpairedfn = os.path.join(outdir, '%s.unpaired.%s' % (prefix, fnend))
    combinedfn = os.path.join(outdir, '%s.combined.%s' % (prefix, fnend))

    ## Check if outfiles exists
    if not replace:
        if os.path.isfile(outfnF) or os.path.isfile(outfnR) or \
            os.path.isfile(unpairedfn) or os.path.isfile(combinedfn):
            logFile.write("At least one of the output files is present. Remove "
                "them or run with --replace parameter. Exiting!\n")
            logFile.flush()
            exit(-3)

    #open files for writting
    outF = outR = outCombined = outUnpaired = False
    if not noSeparate:
        outF = open(outfnF, 'w')
        outR = open(outfnR, 'w')
    #open out file for unpaired reads
    if unpaired:
        outUnpaired = open(unpairedfn, 'w')
    #open out file for combined FastQ
    if combined:
        outCombined = open(combinedfn, 'w')
    outfiles = (outF, outR, outCombined, outUnpaired)

    #process all input files
    fpair = []
    i = pi = filtered = single = 0
    for fn in inputs:
        fpair.append(fn)
        if len(fpair) != 2:
            continue

        ## Process QSEQ files: GERALD->FASTA
        i, pfiltered, psingle = filter_paired(fpair, outfiles, minlen, maxlen, limit, minqual,\
            qual64offset, qseq, stripHeaders, outformat, pi)
        ## Print info
        if verbose:
            logFile.write('[%s]  %s  %s  %s  %s\n' % \
                (datetime.ctime(datetime.now()), fpair[0].name, fpair[1].name,\
                i-pi, pfiltered))
            logFile.flush()
        #update read counts
        pi        = i
        filtered += pfiltered
        single   += psingle
        #reset fnames
        fpair = []

    ## Close outfiles
    for outfile in outfiles:
        if outfile:
            outfile.close()
    ## Print info
    ratio = (i-filtered)*(100.0/i)
    logFile.write('Processed pairs: %s. Filtered: %s. Reads ' % (i, filtered))
    logFile.write('pairs included: %s [%.2f%c]. ' % (i-filtered, ratio, '%'))
    logFile.write('Orphans: %s [%.2f%c]\n' % (single, single*(100.0/i), '%'))
    logFile.flush()

def filter_single(infile, out, minlen, maxlen, limit, minqual, qual64offset, qseq, \
                  stripHeaders, outformat, pi):
    """Filter single reads."""
    #define parser
    fqparser = fqtrimmer(infile, minlen, maxlen, limit, minqual, qual64offset, qseq, stripHeaders, outformat, pi)
    #process output of subprocesses
    both = filtered = 0
    for i, rec in enumerate(fqparser, pi+1):
        #write pair if F & R passed filtering
        if rec:
            out.write(rec)
            both += 1
        #nothing if both didn't pass filtering
        else:
            filtered += 1
        #print stats
        if not i % 10e3:
            info = "%9s processed [%6.2f%s ok]      \r" % (locale.format("%d", \
                i, grouping=True), (i-filtered)*100.0/i, '%')
            logFile.write(info)
            logFile.flush()
    info = "%9s processed [%6.2f%s ok]       \n" % (locale.format("%d", i, \
        grouping=True), (i-filtered)*100.0/i, '%')
    logFile.write(info)
    logFile.flush()

    return i, filtered

def process_single(inputs, qseq, outdir, outprefix, minlen, maxlen, limit, minqual, \
                   qual64offset, replace, stripHeaders, fasta, verbose):
    """Process single end libraries."""

    ## Define output fnames
    fnend = outformat = 'fasta' if fasta else 'fastq'
    prefix = ("%sq%s_l%s") % (outprefix, minqual, minlen)
    outfn     = os.path.join(outdir, '%s.%s'        % (prefix, fnend))

    #check if outfiles exists
    if not replace:
        if os.path.isfile(outfn):
            logFile.write("File exists: %s. Remove them or run with --replace "
                + "parameter. Exiting!\n"%outfn)
            logFile.flush()
            exit(-3)
    #process input files
    i = pi = filtered = 0
    out = open(outfn, 'w')
    for fn in inputs:
        ## Process QSEQ files: GERALD->FASTA
        i, pfiltered = filter_single(fn, out, minlen, maxlen, limit, minqual, qual64offset, \
            qseq, stripHeaders, outformat, pi)
        ## Print info
        if verbose:
            logFile.write('[%s]   %s  %s  %s\n' % \
                (datetime.ctime(datetime.now()), fn, i-pi, pfiltered))
            logFile.flush()
        #update read counts
        pi        = i
        filtered += pfiltered
    #close outfile
    out.close()
    #print info
    logFile.write('Processed: %s. Filtered: %s. Reads included: %s [%.2f%s].\n'\
        % (i, filtered, i-filtered, (i-filtered)*100.0/i, '%'))
    logFile.flush()

def main():

    usage  = "%(prog)s -p -u -q10 -l31 -i sample_read1.fastq.gz sample_read2.fastq.gz -o sample [options]"
    parser = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("--version", action="version", version='%(prog)s 0.23')
    parser.add_argument("-i", "--inputs", nargs="+", type=str,
                        help="input file(s)")
    parser.add_argument("-o", "--outdir", default='outdir',
                        help="define where to store output files")
    parser.add_argument("-g", "--qseq", action="store_true",  default=False,
                        help="input is QSEQ, not FastQ")
    parser.add_argument("-l", "--minlen", default=31, type=int,
                        help="min read length (shorter after quality trimming are removed) [%(default)s]" )
    parser.add_argument("-m", "--maxlen", default=0, type=int,
                        help="max read length [entire read, unless quality trimming]" )
    parser.add_argument("--limit", default=0, type=int,
                        help="process subset of reads [all reads]" )
    parser.add_argument("-q", "--minqual", default=None, type=int,
                        help="read is clipped @ first base having PHRED quality lower than [%(default)s]" )
    parser.add_argument("-t", dest="qual64offset", default=False, action='store_true',
                        help="use PHRED+64 (illumina/solexa) quality encoding [Sanger PHRED+33 encoding]")
    parser.add_argument("-p", "--paired", default=False, action="store_true",
                        help="paired-end reads (qXX_1.fastq & qXX_2.fastq)")
    parser.add_argument("-u", "--unpaired", default=False, action="store_true",
                        help="store orphaned reads > PREFIXqXX_mYY.unpaired.fastq")
    parser.add_argument("-r", "--replace", default=False, action="store_true",
                        help="overwrite output files")
    parser.add_argument("-b", "--noSeparate", default=False, action="store_true",
                        help="don't store separate fastQ for F & R reads > PREFIXqXX_mYY.1.fastq PREFIXqXX_mYY.2.fastq" )
    parser.add_argument("-c", "--combined", default=False, action="store_true",
                        help="store combined fastQ for paired reads > PREFIXqXX_mYY.combined.fastq [%(default)s]" )
    parser.add_argument("-H", "--stripHeaders", default=False, action="store_true",
                        help="replace headers by int [%(default)s]" )
    parser.add_argument("--fasta", default=False, action="store_true",
                        help="report fasta, not FastQ" )

    parser.add_argument("--log", default="", type=str,
                        help="Dump log into a file rather than to the stdout/stderr")
    parser.add_argument("--prefix", default="", type=str,
                        help="Add a prefix to output files")
    parser.add_argument("--phix_seq", default=None, type=str,
                        help="Input the PHIX Sequence to remove any suspicious read")
    parser.add_argument("--ignore_quality", default=True, action="store_false",
                        help="Ignore Encoding Quality Check" )

    o = parser.parse_args()

    global logFile
    ## Define the output stream
    logFile = open(o.log, "w") if o.log != "" else sys.stderr

    if o.verbose:
        logFile.write("Options: %s\n" % str(o))
        logFile.flush()

    ## Check if input files exist
    for inFile in o.inputs:
        if not os.path.isfile(inFile):
            logFile.write(("\nERROR: Check input file '%s'\n") % (inFile))
            logFile.flush()
            exit(-1)

    ## Create output directory if not present already
    if not os.path.isdir(o.outdir):
        os.makedirs(o.outdir)

    ## If define, read the PHIX Sequence files and fragment it into chunks of
    ## min_length size
    global phixReads, phixSeq
    phixReads, phixSeq = set(), None
    '''
    if o.phix_seq and os.path.isfile(o.phix_seq):
        phixSeq = str((SeqIO.read(open(o.phix_seq, "rU"), "fasta")).seq)
        end = len(phixSeq) - o.minlen + 1
        phixReads = set([phixSeq[pos:pos+o.minlen] for pos in range(0, end)])'''

    ## Check if quality encoding correct
    ## Quality offset checking need to be done!
    if o.ignore_quality:
        state, msg = checkQualityEncoding(o.inputs[0], 1000, o.qual64offset, \
            o.qseq)
        if not state:
            logFile.write(msg)
            logFile.flush()
            exit(-2)

    ## Check if prefix is endend with '.' or '_'. Otherwise, add '.' at the end
    if o.prefix and not o.prefix[-1] in ['.', '_']:
        o.prefix += '.'

    ## We pass file descriptors to the main functions
    inputs = [open(inFile, "r") for inFile in o.inputs]

    ## Gerald2fastq
    if o.paired:
        process_paired(inputs, o.qseq, o.outdir, o.prefix, o.unpaired, \
            o.minlen, o.maxlen, o.limit, o.minqual, o.noSeparate, o.combined, o.qual64offset, \
            o.replace, o.stripHeaders, o.fasta, o.verbose)
    else:
        process_single(inputs, o.qseq, o.outdir, o.prefix, o.minlen, \
            o.minqual, o.maxlen, o.limit, o.qual64offset, o.replace, o.stripHeaders, o.fasta, \
            o.verbose)

if __name__=='__main__':
    t0=datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!                \n")
    dt=datetime.now()-t0
    logFile.write("## Time elapsed: %s\n" % dt)
    logFile.flush()
