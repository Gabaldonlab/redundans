#!/usr/bin/env python
desc="""Align pairs/mates onto contigs and run SSPACE scaffolder.
Example:
fastq2sspace.py -v -f contigs.fa -n pe300 pe600 pe5000 -1 ../../_archives/PL429.{3,6,50}00_read1*.fastq.gz -2 ../../_archives/PL429.{3,6,50}00_read2*.fastq.gz -i 300 600 5000 -s 0.15 0.25 0.5 -t FR FR RF -u 5000000
"""
epilog="""Author:
l.p.pryszcz@gmail.com

19/06/2012 Dublin
"""

import commands, os, subprocess, sys
from datetime import datetime

def _unload_sam(sam):
    return sam[0], int(sam[1]), sam[2], int(sam[3]), int(sam[4]), len(sam[9])

def parse_sam(handle):
    """Return tuple representing entries from SAM."""
    q1 = q2 = ""
    for l in handle:
        l = l.strip()
        if not l or l.startswith('@'):
            continue
        sam = l.split('\t')
        # first in pair
        if int(sam[1]) & 64:
            # skip multiple matches
            if sam[0] == q1:
                continue
            q1, flag1, ref1, start1, mapq1, len1 = _unload_sam(sam)
        else:
            # skip multiple matches
            if sam[0] == q2:
                continue
            q2, flag2, ref2, start2, mapq2, len2 = _unload_sam(sam)
        # report
        if q1 == q2:
            yield q1, flag1, ref1, start1, mapq1, len1, q2, flag2, ref2, start2, mapq2, len2

def get_start_stop(start, length, flag):
    """Return start-end read boundaries.
    Return end-start if reverse aligned (flag & 16)."""
    if flag & 16:
        end    = start
        start += length
    else:
        end    = start + length
    return start, end

def sam2sspace_tab(inhandle, outhandle, mapqTh=0, upto=float('inf'), verbose=False, log=sys.stderr):
    """Convert SAM to SSPACE TAB file."""
    i = j = k = pq1 = 0
    for q1, flag1, ref1, start1, mapq1, len1, q2, flag2, ref2, start2, mapq2, len2 in parse_sam(inhandle):
        i   += 1
        if upto and i>upto:
            break        
        #skip 0 quality pair
        if mapqTh:
            if mapq1 < mapqTh or mapq2 < mapqTh:
                continue  
        if q1!=q2:
            log.write("Warning: Queries have different names: %s vs %s\n" % (q1, q2))
            continue
        j   += 1
        #skip self matches
        if ref1==ref2:
            continue
        k += 1
        #define start-stop ranges
        start1, end1 = get_start_stop(start1, len1, flag1)
        start2, end2 = get_start_stop(start2, len2, flag2)
        #print output
        outhandle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (ref1, start1, end1, ref2, start2, end2))
    if i:
        info = "   %s pairs. %s passed filtering [%.2f%s]. %s in different contigs [%.2f%s].\n" % (i, j, j*100.0/i, '%', k, k*100.0/i, '%')
    else:
        info = "   No pairs were aligned!\n"
    log.write(info)
    
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
        bwtmessage = commands.getoutput(cmd)
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
    idxcmd = "snap-aligner index %s %s" % (ref, idxfn)
    if not os.path.isdir(idxfn):
        if verbose:
            log.write(" Creating index...\n  %s\n" % idxcmd)
        idxmessage = commands.getoutput(idxcmd)
        log.write(idxmessage)
    # skip mate rescue
    args = ['snap-aligner', 'paired', idxfn, fn1, fn2, '--b', '-t', str(cores), '-o', '-sam', '-']
    if verbose:
        log.write( "  %s\n" % " ".join(args))
    #select ids
    proc = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=log)
    return proc
    
def get_tab_files(outdir, reffile, libNames, fReadsFnames, rReadsFnames, inserts, iBounds, libreadlen, \
                  cores, mapqTh, upto, verbose, log=sys.stderr):
    """Prepare genome index, align all libs and save TAB file"""
    ref = reffile.name
    tabFnames = []
    _get_aligner_proc = _get_bwamem_proc
    if max(libreadlen)<=300 and min(libreadlen)>40:
        _get_aligner_proc = _get_snap_proc
    # process all libs
    for libName, f1, f2, iSize, iFrac in zip(libNames, fReadsFnames, rReadsFnames, inserts, iBounds):
        if verbose:
            log.write( "[%s] [lib] %s\n" % (datetime.ctime(datetime.now()), libName))
        # define tab output
        outfn = "%s.%s.tab" % (outdir, libName)
        # skip if file exists
        if os.path.isfile(outfn):
            log.write("  File exists: %s\n" % outfn)
            tabFnames.append(outfn)
            continue
        out = open(outfn, "w")
        # define max insert size allowed
        maxins = (1.0+iFrac) * iSize
        # run alignment for all libs        
        bwalog = open(outfn+".log", "w")
        proc = _get_aligner_proc(f1.name, f2.name, ref, cores, verbose, bwalog)
        # parse botwie output
        sam2sspace_tab(proc.stdout, out, mapqTh, upto, verbose, log)
        # close file
        out.close()
        tabFnames.append(outfn)
        # terminate subprocess
        proc.terminate()
    return tabFnames
    
def get_libs(outdir, libFn, libNames, tabFnames, inserts, iBounds, orientations, verbose, log=sys.stderr):
    """Save lib fname and return it's path"""
    lines = []
    # load libs from file
    if libFn:
        if verbose:
            log.write(" Reading libs from %s\n" % libFn)
        lines = open(libFn).readlines()
    # add TAB libs
    tabline = "%s\tTAB\t%s\t%s\t%s\t%s\n"
    for libname, tabfn, isize, isfrac, orient in zip(libNames, tabFnames, inserts, \
                                                     iBounds, orientations):
        lines.append(tabline%(libname, os.path.basename(tabfn), \
                              isize, isfrac, orient))

    outfn = "%s.libs.txt" % outdir #os.path.join( outdir,"libs.txt" )
    if verbose:
        log.write( " Updated libs saved to: %s\n" % outfn )
    out   = open( outfn,"w" ); out.write( "".join(lines) ); out.close()
    return outfn

def fastq2sspace(out, fasta, lib, libnames, libFs, libRs, orientations,  \
                 libIS, libISStDev, libreadlen, cores, mapq, upto, minlinks, linkratio, \
                 sspacebin, verbose, log=sys.stderr):
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
    tabFnames = get_tab_files(out, fasta, libnames, libFs, libRs, libIS, libISStDev, libreadlen, cores, mapq, upto, verbose, log)

    # generate lib file
    if  verbose:
        log.write("[%s] Generating libraries file...\n" % datetime.ctime(datetime.now()) )
    libFn = get_libs(out, lib, libnames, tabFnames, libIS, libISStDev, orientations, verbose, log)

    # run sspace
    ## change dir to outdir - sspace output intermediate files always in curdir
    os.chdir(outdir)
    CMD = "perl %s -l %s -a %s -k %s -s %s -b %s > %s.sspace.log"
    cmd = CMD%(sspacebin, os.path.basename(libFn), minlinks, linkratio, \
               os.path.basename(fasta.name), outfn, outfn)
    if verbose:
        log.write(" %s\n"%cmd)
    os.system(cmd)
    ## change to basal dir
    os.chdir(curdir)
    
def main():
    import argparse
    usage   = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose", default=False, action="store_true")
    parser.add_argument("-f", dest="fasta",      required=True, type=file,
                       help="genome fasta        [mandatory]")
    parser.add_argument("-k", dest="minlinks",   default=5, type=int,
                       help="min number of links [%(default)s]")    
    parser.add_argument("-a", "--linkratio",     default=0.7, type=float,
                       help="max link ratio between two best contig pairs [%(default)s]")    
    parser.add_argument("-l", dest="lib",        default="",
                       help="lib file            [No libs]")    
    parser.add_argument("-o", dest="out",        default="sspace_out",
                       help="output basename     [%(default)s]")
    parser.add_argument("-n", dest="libnames",   nargs="+",
                       help="libraries names     [%(default)s]")
    parser.add_argument("-1", dest="libFs",      nargs="+", type=file,
                       help="libs forward reads  [%(default)s]")
    parser.add_argument("-2", dest="libRs",      nargs="+", type=file,
                       help="libs reverse reads  [%(default)s]")
    parser.add_argument("-i", dest="libIS",      nargs="+", type=int,
                       help="libs insert sizes   [%(default)s]")
    parser.add_argument("-s", dest="libISStDev", nargs="+", type=float,
                       help="libs IS StDev       [%(default)s]")    
    parser.add_argument("-t", dest="orientations", nargs="+", #type=float,
                       help="libs orientations   [%(default)s]")    
    parser.add_argument("-c", dest="cores",      default=2, type=int,
                       help="no. of cpus         [%(default)s]")
    parser.add_argument("-q", dest="mapq",       default=10, type=int,
                       help="min map quality     [%(default)s]")
    parser.add_argument("-u", dest="upto",       default=0,  type=int,
                       help="process up to pairs [all]")
    parser.add_argument("--sspacebin", default="~/src/SSPACE/SSPACE_Standard_v3.0.pl", 
                       help="SSPACE2 perl script [%(default)s]")
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )
    
    if len(o.libnames)*6 != len(o.libnames)+len(o.libFs)+len(o.libRs)+len(o.libIS)+len(o.libISStDev)+len(o.orientations):
        parser.error("Wrong number of arguments!")

    fastq2sspace(o.out, o.fasta, o.lib, o.libnames, o.libFs, o.libRs, o.orientations, \
                 o.libIS, o.libISStDev, o.cores, o.mapq, o.upto, o.minlinks, o.linkratio, \
                 o.sspacebin, o.verbose)
    
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!   \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    