#!/usr/bin/env python
desc="""Align pairs/mates onto contigs and run SSPACE scaffolder.
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

19/06/2012 Dublin/Warsaw

Updated to Python3 by Diego Fuentes Palacios
Barcelona 08/18/2022
"""

import os, subprocess, sys, tempfile
from datetime import datetime
def _unload_sam(sam):
    return sam[0], int(sam[1]), sam[2], int(sam[3]), int(sam[4]), sam[9], sam[10]

def parse_sam(handle):
    """Return tuple representing entries from SAM."""
    q1 = q2 = ""
    for le in handle:
        l= le.decode('utf-8')
        l = l.strip()
        if not l or l.startswith('@'):
            continue
        sam = l.split('\t')
        ## first in pair

        #However, ingore if the number of fields does not reach minimum
        if len(sam) < 11:
            continue

        if int(sam[1]) & 64:
            # skip multiple matches
            if sam[0] == q1:
                continue
            q1, flag1, ref1, start1, mapq1, seq1, qual1 = _unload_sam(sam)
        else:
            # skip multiple matches
            if sam[0] == q2:
                continue
            q2, flag2, ref2, start2, mapq2, seq2, qual2 = _unload_sam(sam)
        # report
        if q1 == q2:
            yield q1, flag1, ref1, start1, mapq1, seq1, qual1, q2, flag2, ref2, start2, mapq2, seq2, qual2

def get_start_stop(start, length, flag):
    """Return start-end read boundaries.
    Return end-start if reverse aligned (flag & 16)."""
    if flag & 16:
        end = start
        start += length
    else:
        end = start + length
    return start, end

base2rc = {"A": "T", "T": "A", "C": "G", "G": "C", "a": "t", "t": "a", "c": "g", "g": "c", "N": "N", "n": "n"}
def sam2fastq(name, seq, qual, flag):
    """Return FastQ"""
    if flag & 16:
        seq = "".join(reversed([base2rc[s] for s in seq]))
        qual = "".join(reversed(qual))
    return "@%s\n%s\n+\n%s\n"%(name, seq, qual)
    
def sam2sspace_tab(inhandle, outhandle, mapqTh=0, upto=float('inf'), verbose=False, log=sys.stderr, ref="", cores=4):
    """Convert SAM to SSPACE TAB file."""
    i = j = k = pq1 = 0
    if ref:
        _tmpfile = tempfile.NamedTemporaryFile(delete=False)
    else:
        _tmpfile = ""
    for q1, flag1, ref1, start1, mapq1, seq1, qual1, q2, flag2, ref2, start2, mapq2, seq2, qual2 in parse_sam(inhandle):
        i   += 1
        if upto and i>upto:
            break        
        #skip 0 quality pair
        if mapqTh:
            if mapq1 < mapqTh or mapq2 < mapqTh:
                if _tmpfile:
                    _tmpfile.write((sam2fastq("%s/1"%i, seq1, qual1, flag1).encode('utf-8')+sam2fastq("%s/2"%i, seq2, qual2, flag2).encode('utf-8')))
                continue
        if q1!=q2:
            log.write("[Warning] Queries have different names: %s vs %s\n" % (q1, q2))
            continue
        j   += 1
        #skip self matches
        if ref1==ref2:
            continue
        k += 1
        #define start-stop ranges
        start1, end1 = get_start_stop(start1, len(seq1), flag1)
        start2, end2 = get_start_stop(start2, len(seq2), flag2)
        #print output
        outhandle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ref1, start1, end1, ref2, start2, end2))
    if i:
        info = "   %s pairs. %s passed filtering [%.2f%s]. %s in different contigs [%.2f%s].\n" % (i, j, j*100.0/i, '%', k, k*100.0/i, '%')
    else:
        info = "   No pairs were aligned!\n"
    log.write(info)
    # run last
    if _tmpfile:
        _tmpfile.close()
        lastproc = _get_last_proc(_tmpfile.name, ref, cores)
        last_tab2sspace_tab(lastproc.stdout, outhandle, mapqTh, upto, verbose, log)
        lastproc.kill()
        os.unlink(_tmpfile.name)
    
def _get_bwamem_proc(fn1, fn2, ref, cores, verbose, log=sys.stderr):
    """Return bwamem subprocess.
    bufsize: 0 no buffer; 1 buffer one line; -1 set system default.
    """
    # create genome index
    idxfn = ref + ".pac"
    if not os.path.isfile(idxfn):
        cmd = "bwa index %s" % (ref,)
        if verbose:
            log.write(" Creating index...\n  %s\n" % cmd)
        subprocess.Popen(cmd.split(), stdout=log, stderr=log).wait()
    # skip mate rescue
    bwaArgs = ['bwa', 'mem', '-S', '-t', str(cores), ref, fn1, fn2]
    if verbose:
        log.write( "  %s\n" % " ".join(bwaArgs))
    #select ids
    bwaProc = subprocess.Popen(bwaArgs, stdout=subprocess.PIPE, stderr=log)
    return bwaProc
    
def _get_snap_proc(fn1, fn2, ref, cores, verbose, log=sys.stderr):
    """Return snap-aligner subprocess.
    bufsize: 0 no buffer; 1 buffer one line; -1 set system default.
    """
    # create genome index
    idxfn = ref + ".snap"
    idxcmd = "snap-aligner index %s %s -bSpace" % (ref, idxfn)
    if not os.path.isdir(idxfn):
        if verbose:
            log.write(" Creating index...\n  %s\n" % idxcmd)
        subprocess.Popen(idxcmd.split(), stdout=log, stderr=log).wait()
    # -d maxEditDist should be set based on readlen and expected divergence
    args = ['snap-aligner', 'paired', idxfn, fn1, fn2, '-d', '30', '--b', '-t', str(cores), '-o', '-sam', '-']
    if verbose:
        log.write( "  %s\n" % " ".join(args))
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=log)
    return proc

def _get_last_proc(fqfname, ref, cores, verbose=0, log=sys.stderr):
    """Run last process"""
    # create genome index
    if not os.path.isfile(ref+".suf"):
        cmd = "lastdb -uNEAR -W 11 %s %s" % (ref, ref)
        subprocess.Popen(cmd.split(), stdout=log, stderr=log).wait()
    # skip mate rescue 
    args1 = ['cat', fqfname]
    args2 = ['lastal', '-T1', '-Q1', '-fTAB', '-P%s'%cores, ref, fqfname] 
    if verbose:
        log.write("  %s\n"%(" ".join(args2),))
    #select ids
    #proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=log) # stdin=proc1.stdout, 
    proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stderr=log) #, preexec_fn=os.setpgrp)
    return proc2

def _last2pairs(handle):
    """Yield pairs from LASTal"""
    pq = ""
    hits = [[]]
    for le in handle:
        l= le.decode('utf-8')
        if l.startswith('#'): 
            continue
        # unpack
        (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
        score, tstart, talg = list(map(int, (score, tstart, talg)))
        # report previous query
        if pq != q:
            if list(map(len, hits)) == [1, 1] and hits[0][0][-1][:-1] == hits[1][0][-1][:-1]:
                yield hits[0][0][1:4], hits[1][0][1:4]
                hits = [[]]
            elif len(hits)==2:
                hits = hits[1:]
            if hits[-1]:
                hits.append([])
        # store current hit
        s, e = tstart, tstart+talg
        if qstrand=="-":
            s, e = e, s
        data = (score, t, s, e, q)
        if not hits[-1]: # or score>hits[-1][0][0]:
            hits[-1] = [data]
        elif score==hits[-1][0][0]:
            hits[-1].append(data)
        pq = q
    # yield last bit
    if list(map(len, hits)) == [1, 1] and hits[0][0][-1][:-1] == hits[1][0][-1][:-1]:
        yield hits[0][0][1:4], hits[1][0][1:4]
    
def last_tab2sspace_tab(handle, out, mapqTh, upto, verbose, log, proc=""):
    """Generate TAB based on LASTal alignments"""
    i = k = 0
    for i, ((ref1, start1, end1), (ref2, start2, end2)) in enumerate(_last2pairs(handle), 1):
        if upto and i>upto:
            break
        if ref1==ref2:
            continue
        k += 1
        out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ref1, start1, end1, ref2, start2, end2))
    if i:
        info = "    %s pairs. %s in different contigs [%.2f%s].\n" % (i, k, k*100.0/i, '%')
    else:
        info = "    No pairs were aligned!\n"
    log.write(info)
    
def get_tab_files(outdir, reffile, libNames, fReadsFnames, rReadsFnames, inserts, iBounds, libreadlen, \
                  cores, mapqTh, upto, verbose, usebwa=False, log=sys.stderr):
    """Prepare genome index, align all libs and save TAB file"""
    ref = reffile.name
    tabFnames = []
    #if max(libreadlen)<=500 and min(libreadlen)>40:
    if usebwa:
        _get_aligner_proc = _get_bwamem_proc
        ref2 = '' # disable lastal
    else:
        _get_aligner_proc = _get_snap_proc
        ref2 = ref
    # process all libs
    for libName, f1, f2, iSize, iFrac in zip(libNames, fReadsFnames, rReadsFnames, inserts, iBounds):
        if verbose:
            log.write("[%s] [lib] %s\n" % (datetime.ctime(datetime.now()), libName))
        # define tab output
        outfn = "%s.%s.tab" % (outdir, libName)
        # skip if file exists
        if os.path.isfile(outfn):
            log.write("  File exists: %s\n" % outfn)
            tabFnames.append(outfn)
            continue
        # run alignment for all libs
        out = open(outfn, "w")
        bwalog = open(outfn+".log", "w")
        proc = _get_aligner_proc(f1.name, f2.name, ref, cores, verbose, bwalog)
        # parse botwie output
        sam2sspace_tab(proc.stdout, out, mapqTh, upto, verbose, log, ref2, cores) #
        # close file & terminate subprocess
        out.close()
        tabFnames.append(outfn)
        proc.kill()
    return tabFnames
    
def get_libs(outdir, libFn, libNames, tabFnames, inserts, iBounds, orientations, libreadlen, verbose, log=sys.stderr):
    """Save lib fname and return it's path"""
    # load libs from file
    lines = []
    if libFn:
        if verbose:
            log.write(" Reading libs from %s\n" % libFn)
        lines = open(libFn).readlines()
    # add TAB libs
    tabline = "%s\tTAB\t%s\t%s\t%s\t%s\n"
    for libname, tabfn, isize, isfrac, orient, rlen in zip(libNames, tabFnames, inserts, iBounds, orientations, libreadlen):
        lines.append(tabline%(libname, os.path.basename(tabfn), isize, isfrac, orient))

    outfn = "%s.libs.txt" % outdir 
    if verbose:
        log.write( " Updated libs saved to: %s\n" % outfn )
    with open(outfn, "w") as out:
        out.write("".join(lines))
    return outfn

def fastq2sspace(out, fasta, lib, libnames, libFs, libRs, orientations,  \
                 libIS, libISStDev, libreadlen, cores, mapq, upto, minlinks, linkratio, \
                 sspacebin, verbose, usebwa=False, log=sys.stderr):
    """Map reads onto contigs, prepare library file and execute SSPACE2"""
    # get dir variables
    curdir = os.path.abspath(os.path.curdir)
    outfn  = os.path.basename(out)
    outdir = os.path.dirname(out)
    # generate outdirs if out contain dir and dir not exists
    if outdir:
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    
    # get tab files
    if verbose:
        log.write("[%s] Generating TAB file(s) for %s library/ies...\n" % (datetime.ctime(datetime.now()),len(libnames)) )
    tabFnames = get_tab_files(out, fasta, libnames, libFs, libRs, libIS, libISStDev, libreadlen, \
                              cores, mapq, upto, verbose, usebwa, log)
    
    # generate lib file
    if  verbose:
        log.write("[%s] Generating library file(s)...\n" % datetime.ctime(datetime.now()))
    libFn = get_libs(out, lib, libnames, tabFnames, libIS, libISStDev, orientations, libreadlen, verbose, log)

    # run sspace
    ## change dir to outdir - sspace output intermediate files always in curdir
    os.chdir(outdir)
    CMD = "perl %s -l %s -a %s -k %s -s %s -b %s > %s.sspace.log"
    cmd = CMD%(sspacebin, os.path.basename(libFn), minlinks, linkratio, os.path.basename(fasta.name), outfn, outfn)
    if verbose:
        log.write(" %s\n"%cmd)
    os.system(cmd)
    ## change to basal dir
    os.chdir(curdir)
    
def main():
    import argparse
    usage   = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser(usage=usage, description=desc, epilog=epilog)

    parser.add_argument("-v", "--verbose", default=False, action="store_true")
    parser.add_argument("-f", "--fasta", required=True, type=file, help="genome fasta")
    parser.add_argument("-k", "--minlinks", default=5, type=int, help="min number of links [%(default)s]")    
    parser.add_argument("-a", "--linkratio", default=0.7, type=float, help="max link ratio between two best contig pairs [%(default)s]") 
    parser.add_argument("-l", "--lib", default="", help="lib file [No libs]")    
    parser.add_argument("-o", "--out", default="sspace_out", help="output basename [%(default)s]")
    parser.add_argument("-n", "--libnames", nargs="+", help="libraries names [%(default)s]")
    parser.add_argument("-1", "--libFs", nargs="+", type=file, help="libs forward reads [%(default)s]")
    parser.add_argument("-2", "--libRs", nargs="+", type=file, help="libs reverse reads [%(default)s]")
    parser.add_argument("-i", "--libIS", nargs="+", type=int, help="libs insert sizes [%(default)s]")
    parser.add_argument("-s", "--libISStDev", nargs="+", type=float, help="libs IS StDev [%(default)s]")    
    parser.add_argument("-t", "--orientations", nargs="+", help="libs orientations [%(default)s]")    
    parser.add_argument("-c", "--cores", default=2, type=int, help="no. of cpus [%(default)s]")
    parser.add_argument("-q", "--mapq", default=10, type=int, help="min map quality [%(default)s]")
    parser.add_argument("-u", "--upto", default=0, type=int, help="process up to pairs [all]")
    parser.add_argument("-b", "--usebwa", action='store_true', help="use bwa mem for alignment [use snap-aligner]")
    parser.add_argument("--sspacebin", default="~/src/SSPACE/SSPACE_Standard_v3.0.pl", help="SSPACE perl script [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write("Options: %s\n" % str(o))
    
    if len(o.libnames)*6 != len(o.libnames)+len(o.libFs)+len(o.libRs)+len(o.libIS)+len(o.libISStDev)+len(o.orientations):
        parser.error("Wrong number of arguments!")

    fastq2sspace(o.out, o.fasta, o.lib, o.libnames, o.libFs, o.libRs, o.orientations, \
                 o.libIS, o.libISStDev, o.cores, o.mapq, o.upto, o.minlinks, o.linkratio, \
                 o.sspacebin, o.usebwa, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!   \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    