#!/usr/bin/env python
desc="""Heterozygous genome assembly pipeline. It consists of three steps:
reduction, scaffolding and gap closing.

More info at: http://bit.ly/Redundans

TBA:
- add exception if lastdb or lastal doesn't finish successfully
"""
epilog="""Author:
l.p.pryszcz@gmail.com

Mizerow/Warsaw/Bratislava/Barcelona, 17/10/2014
"""

import commands, os, resource, sys
import glob, subprocess, time
from datetime import datetime

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["bin", "bin/bwa", "bin/snap", "bin/last/build", "bin/last/scripts", "bin/last/src", "bin/pyScaf"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])

from fasta2homozygous import fasta2homozygous
from fastq2sspace import fastq2sspace
from fastq2insert_size import fastq2insert_size
from filterReads import filter_paired
#from fasta_stats import fasta_stats
from FastaIndex import FastaIndex, symlink
from pyScaf import LongReadGraph, SyntenyGraph

def timestamp():
    """Return formatted date-time string"""
    return "\n%s\n[%s] "%("#"*50, datetime.ctime(datetime.now()))

def get_libraries(fastq, fasta, mapq, threads, verbose, log=sys.stderr, limit=0,
                  libraries=[], stdfracTh=0.66, maxcfracTh=0.9, genomeFrac=0.05):
    """Return libraries"""
    # skip if all libs OKish
    ## max stdfrac cannot be larger than stdfracTh in any of the libraries
    if libraries and not filter(lambda x: x>stdfracTh, (max(lib[5]) for lib in libraries)):
        return libraries
    
    # otherwise process all reads
    if not limit or limit<10e5:
        limit = 10e5
    
    # get libraries statistics using 1% of mapped read limit
    libdata = fastq2insert_size(log, fastq, fasta, mapq, threads, limit/100, verbose, log, \
                                genomeFrac, stdfracTh, maxcfracTh)
    # separate paired-end & mate pairs
    ## also separate 300 and 600 paired-ends
    libraries = []
    # add libraries strating from lowest insert size
    for data in sorted(libdata, key=lambda x: x[3]):
        fq1, fq2, readlen, ismedian, ismean, isstd, pairs, orientation = data
        # add new library set if 
        if not libraries or ismean > 1.5*libraries[-1][4][0]:
            # libnames, libFs, libRs, orientations, libIS, libISStDev
            libraries.append([[], [], [], [], [], [], []])
            i = 1
        # add libname & fastq files
        libraries[-1][0].append("lib%s"%i)
        libraries[-1][1].append(open(fq1))
        libraries[-1][2].append(open(fq2))
        libraries[-1][6].append(readlen)
        # orientation
        #orientation = get_orientation(pairs, fq1, fq2, log, maxcfracTh)
        libraries[-1][3].append(orientation)
        # insert size information
        libraries[-1][4].append(int(ismean))
        stdfrac = isstd / ismean
        # capture large stdev
        if stdfrac > stdfracTh:
            log.write("[WARNING] Highly variable insert size (%.f +- %.2f) in %s - %s!\n"%(ismean, isstd, fq1, fq2))
        # SSSPACE accepts stdfrac 0-1.0
        if stdfrac > 1:
            stdfrac = 1.0
        libraries[-1][5].append(stdfrac)
        # update counter
        i += 1
    return libraries

def get_read_limit(fasta, readLimit, verbose, log=sys.stderr):
    """Return read limit and libraries."""
    # limit no. of reads to align as fraction of genome size
    limit = 0
    if readLimit:
        faidx = FastaIndex(fasta)
        fastaSize = faidx.genomeSize
        limit = int(readLimit * fastaSize)
        if verbose:
            log.write(" Aligning %s mates per library...\n"%limit)
    return limit
    
def run_scaffolding(outdir, scaffoldsFname, fastq, libraries, reducedFname, mapq, threads, \
                    joins, linkratio, limit, iters, sspacebin, gapclosing, verbose, log, \
                    identity, overlap, minLength, resume, lib=""):
    """Execute scaffolding step using libraries with increasing insert size
    in multiple iterations.
    """
    pout = reducedFname
    i = 0
    while i < len(libraries):
        libnames, libFs, libRs, orients, libIS, libISStDev, libreadlen = libraries[i]
        i += 1
        for j in range(1, iters+1):
            out = os.path.join(outdir, "_sspace.%s.%s"%(i, j))
            # resume if files don't exist
            if resume>1 or _corrupted_file(out+".fa"):
                resume += 1
                if verbose:
                    log.write(" iteration %s.%s of %s.%s ...\n"%(i, j, len(libraries), iters))
                lib = ""
                # run fastq scaffolding
                fastq2sspace(out, open(pout), lib, libnames, libFs, libRs, orients, \
                             libIS, libISStDev, libreadlen, threads, mapq, limit, linkratio, joins, \
                             sspacebin, verbose=0, log=log)
            # store out info
            pout = out+".fa"
            # link output ie out/_sspace.1.1/_sspace.1.1.scaffolds.fasta --> out/_sspace.1.1.scaffolds.fasta
            targetout = os.path.join(os.path.basename(out), os.path.basename(out+".final.scaffolds.fasta"))
            symlink(targetout, pout)
            # if number of gaps larger than 1%, run gap closer & reduction
            stats     = FastaIndex(pout).stats()
            fastaSize = int(stats.split('\t')[2])
            gapSize   = int(stats.split('\t')[-2])
            if gapclosing and 1.0 * gapSize / fastaSize > 0.01:
                nogapsFname = ".".join(pout.split(".")[:-1]) + ".filled.fa"
                if resume>1 or _corrupted_file(nogapsFname):
                    resume += 1
                    # close gaps
                    if verbose:
                        log.write("  closing gaps ...\n")
                    basename    = "_sspace.%s.%s._gapcloser"%(i, j)
                    run_gapclosing(outdir, mapq, [libraries[i-1],], nogapsFname, pout, \
                                   threads, limit, iters=1, resume=resume, verbose=0, log=log, basename=basename)
                pout = nogapsFname
        # update library insert size estimation, especially for mate-pairs
        libraries = get_libraries(fastq, pout, mapq, threads, verbose=0,log=log,
                                  libraries=libraries)
    # create symlink to final scaffolds or pout
    symlink(os.path.basename(pout), scaffoldsFname)
    symlink(os.path.basename(pout+".fai"), scaffoldsFname+".fai")
    return libraries, resume

def filter_reads(outdir, fq1, fq2, minlen, maxlen, limit, minqual):
    """Filter FastQ files and return output fnames."""
    fastq = (fq1, fq2)
    # generate output files
    fn1 = os.path.join(outdir, "_reads.%s"%(os.path.basename(fq1.name.rstrip('.gz'))))
    fn2 = os.path.join(outdir, "_reads.%s"%(os.path.basename(fq2.name.rstrip('.gz'))))
    # skip if fq files already generated
    if os.path.isfile(fn1) and os.path.isfile(fn2):
        # notify about empty trimmed libs
        if not os.path.getsize(fn1) or not os.path.getsize(fn2):
            return fn1, fn2, 0
        return fn1, fn2, 1
    # open output files
    out1 = open(fn1, "w")
    out2 = open(fn2, "w")
    outfiles = (out1, out2, 0, 0)
    # run filtering
    i, filtered, orphans = filter_paired(fastq, outfiles, minlen, maxlen, limit, minqual)
    out1.close()
    out2.close()
    return fn1, fn2, i-filtered
    
def prepare_gapcloser(outdir, mapq, configFn, libFs, libRs, orientations,
                      libIS, libISStDev, minlen, maxlen, limit, verbose, log): 
    """Return SOAPdenovo2 config file needed by GapCloser."""
    lines  = "[LIB]\navg_ins=%s\nreverse_seq=%s\nasm_flags=3\nrank=%s\npair_num_cutoff=5\nmap_len=35\nq1=%s\nq2=%s\n"
    config = ["max_rd_len=%s"%maxlen]
    for i, (fq1, fq2, orient, iSize, iFrac) in enumerate(zip(libFs, libRs, orientations, libIS, libISStDev), 1):
        # consider skipping mate-pairs is libIS>1kb
        # skip orientations other than FR RF
        if orient == "FR":
            reverse_seq = 0
        elif orient == "RF":
            reverse_seq = 1
        else:
            if verbose:
                info = "  Skipping unsupported library %s in: %s - %s!\n"
                log.write(info%(orient, fq1, fq2))
            continue
        # filter reads
        fn1, fn2, passed = filter_reads(outdir, fq1, fq2, minlen, maxlen, limit, mapq)
        #store config info only if some reads
        if passed:
            config.append(lines%(iSize, reverse_seq, i, fn1, fn2))
        
    # store config only if some libs passed filtering
    if len(config)>1:
        with open(configFn, "w") as out:
            out.write("\n".join(config))
        return True
    
def run_gapclosing(outdir, mapq, libraries, nogapsFname, scaffoldsFname, \
                   threads, limit, iters, resume, verbose, log, basename="_gapcloser", \
                   overlap=25, minReadLen=40):
    """Execute gapclosing step."""
    pout = scaffoldsFname
    
    for i, (libnames, libFs, libRs, orientations, libIS, libISStDev, libreadlen) in enumerate(libraries, 1):
        # prepare config file and filter reads
        configFn = os.path.join(outdir, "%s.%s.conf"%(basename, i))
        # skip if not suitable libraries
        maxReadLen = max(libreadlen)
        if not prepare_gapcloser(outdir, mapq, configFn, libFs, libRs, orientations, libIS, libISStDev, \
                                 minReadLen, maxReadLen, limit, verbose, log):
            continue
        # run iterations
        for j in range(1, iters+1):
            out = os.path.join(outdir, "%s.%s.%s.fa"%(basename, i, j))
            # skip if file exists
            if resume>1 or _corrupted_file(out):
                resume += 1
                # run GapCloser
                cmd = ["GapCloser", "-t %s"%threads, "-p %s"%overlap, "-l %s"%maxReadLen, \
                       "-a", pout, "-b", configFn, "-o", out]
                if verbose:
                    log.write(" iteration %s.%s ...\n"%(i,j))
                # run GapCloser and save stdout/err to log file
                with open(out+".log", "w") as gapcloselog:
                    GapCloser = subprocess.Popen(cmd, stdout=gapcloselog, stderr=gapcloselog)
                    GapCloser.wait()
            # store out info
            pout = out
    # create symlink to final scaffolds or pout
    symlink(os.path.basename(pout), nogapsFname)
    symlink(os.path.basename(pout+".fai"), nogapsFname+".fai")
    return resume

def _corrupted_file(fname):
    """Return True if output file doesn't exists or is corrupted."""
    if not os.path.isfile(fname)        or not os.path.getsize(fname) or \
       not os.path.isfile(fname+".fai") or not os.path.getsize(fname+".fai"):
        return True

def prepare_contigs(fasta, contigsFname, minLength=200):
    """Sort contigs starting from the longest and remove too short"""
    with open(contigsFname, "w") as out:
        # init fasta index
        faidx = FastaIndex(fasta)
        # filter out sequences shorter than minLength
        for i, c in enumerate(faidx.sort(minLength=minLength), 1):
            if i%1e5 == 1:
                sys.stderr.write(' %s   \r'%i)
            #seq = faidx.__getitem__(c, name=str(i))
            out.write(faidx[c])
        sys.stderr.write(' %s sequences stored.\n'%i)
        
def redundans(fastq, longreads, fasta, reference, outdir, mapq, 
              threads, resume, identity, overlap, minLength, \
              joins, linkratio, readLimit, iters, sspacebin, \
              reduction=1, scaffolding=1, gapclosing=1, cleaning=1, \
              norearrangements=0, verbose=1, log=sys.stderr):
    """Launch redundans pipeline."""
    fastas = [fasta, ]
    # update fasta list
    lastOutFn = os.path.join(outdir, "contigs.fa")
    fastas.append(lastOutFn)
    # check resume
    orgresume = resume
    if resume:
        log.write("%sResuming previous run from %s...\n"%(timestamp(), outdir))
        if not os.path.isdir(outdir):
            log.write("No such directory: %s!\n"%outdir)
            sys.exit(1)
    # prepare outdir or quit if exists
    elif os.path.isdir(outdir):
        log.write("Directory %s exists!\n"%outdir)
        sys.exit(1)
    else:
        os.makedirs(outdir)
    
    # REDUCTION
    # prepare contigs
    if verbose:
        log.write("%sPreparing contigs...\n"%timestamp())
    prepare_contigs(fasta, lastOutFn, minLength)
    # reduce
    outfn = os.path.join(outdir, "contigs.reduced.fa")
    if reduction and _corrupted_file(outfn):
        resume += 1
        if verbose:
            log.write("%sReduction...\n"%timestamp())
            log.write("#file name\tgenome size\tcontigs\theterozygous size\t[%]\theterozygous contigs\t[%]\tidentity [%]\tpossible joins\thomozygous size\t[%]\thomozygous contigs\t[%]\n")
        # reduce
        with open(outfn, "w") as out:
            info = fasta2homozygous(out, open(lastOutFn), identity, overlap, \
                                    minLength, threads, verbose=0, log=log)
        # update fasta list
        lastOutFn = outfn
        fastas.append(lastOutFn)

    # get read limit & libraries
    if fastq:
        if verbose:
            log.write("%sEstimating parameters of libraries...\n"%timestamp())
        limit     = get_read_limit(lastOutFn, readLimit, verbose, log)
        libraries = get_libraries(fastq, lastOutFn, mapq, threads, verbose, log)
    
    # SCAFFOLDING
    outfn = os.path.join(outdir, "scaffolds.fa")
    if fastq and scaffolding: 
        if verbose:
            log.write("%sScaffolding...\n"%timestamp())
        libraries, resume = run_scaffolding(outdir, outfn, fastq, libraries, lastOutFn, mapq, threads, joins, \
                                            linkratio, limit, iters, sspacebin, gapclosing, verbose, log, \
                                            identity, overlap, minLength, resume)
        # update fasta list
        fastas += filter(lambda x: "_gapcloser" not in x, sorted(glob.glob(os.path.join(outdir, "_sspace.*.fa"))))
        lastOutFn = outfn
        fastas.append(lastOutFn)

    # SCAFFOLDING WITH LONG READS
    outfn = os.path.join(outdir, "scaffolds.longreads.fa")
    if longreads and _corrupted_file(outfn):
        # here maybe sort reads by increasing median read length
        resume += 1
        if verbose:
            log.write("%sScaffolding with long reads...\n"%timestamp())
        poutfn = lastOutFn
        for i, fname in enumerate(longreads, 1):
            if verbose:
                log.write(" iteration %s out of %s ...\n"%(i, len(longreads)))
            s = LongReadGraph(lastOutFn, fname, identity, overlap, maxgap=0, threads=threads, \
                              dotplot="", norearrangements=norearrangements, log=0)
            # save output
            _outfn = os.path.join(outdir, "scaffolds.longreads.%s.fa"%i)
            with open(_outfn, "w") as out:
                s.save(out)
            # store fname
            fastas.append(_outfn)
            poutfn = _outfn
        # symlink last iteration
        symlink(poutfn, outfn)
        # update fasta list
        lastOutFn = outfn
        fastas.append(lastOutFn)

    # REFERENCE-BASED SCAFFOLDING
    outfn = os.path.join(outdir, "scaffolds.ref.fa")
    if reference and _corrupted_file(outfn):
        resume += 1
        if verbose:
            log.write("%sScaffolding based on reference...\n"%timestamp())        
        s = SyntenyGraph(lastOutFn, reference, identity, overlap, maxgap=0, threads=threads, \
                         dotplot="", norearrangements=norearrangements, log=0)
        # save output
        with open(outfn, "w") as out:
            s.save(out)
        # update fasta list
        lastOutFn = outfn
        fastas.append(lastOutFn)
        
    # GAP CLOSING
    outfn = os.path.join(outdir, "scaffolds.filled.fa")
    if fastq and gapclosing: 
        if verbose: 
            log.write("%sGap closing...\n"%timestamp())
        resume = run_gapclosing(outdir, mapq, libraries, outfn, lastOutFn, threads, \
                                limit, iters, resume, verbose, log)
        # update fasta list
        fastas += sorted(glob.glob(os.path.join(outdir, "_gap*.fa")))
        lastOutFn = outfn
        fastas.append(lastOutFn)
    
    # FASTA STATS
    if verbose:
        log.write("%sReporting statistics...\n"%timestamp())
    # report stats
    log.write('#fname\tcontigs\tbases\tGC [%]\tcontigs >1kb\tbases in contigs >1kb\tN50\tN90\tNs\tlongest\n')
    for fn in fastas:
        log.write(FastaIndex(fn).stats())
    
    # Clean-up
    if cleaning:
        if verbose:
            log.write("%sCleaning-up...\n"%timestamp())
        for root, dirs, fnames in os.walk(outdir):
            for i, fn in enumerate(filter(lambda x: not x.endswith(('.fa', '.fasta', '.fai', '.tsv', '.png')), fnames), 1):
                os.unlink(os.path.join(root, fn))
            # rmdir of snap index
            if root.endswith('.snap') and i==len(fnames):
                os.rmdir(root)

    if orgresume:
        log.write("%sResume report: %s step(s) have been recalculated.\n"%(timestamp(), resume-1))
    
def _check_executable(cmd):
    """Check if executable exists."""
    p = subprocess.Popen("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return "".join(p.stdout.readlines())

def _check_dependencies(dependencies):
    """Return error if wrong software version"""
    warning = 0
    info = "[WARNING] Old version of %s: %s. Update to version %s+!\n"
    for cmd, version in dependencies.iteritems():
        out = _check_executable(cmd)
        if "not found" in out:
            warning = 1
            sys.stderr.write("[ERROR] %s\n"%out)
        elif version:
            out = commands.getoutput("%s --version"%cmd)
            curver = out.split()[-1]
            if not curver.isdigit():
                warning = 1
                sys.stderr.write("[WARNING] Problem checking %s version: %s\n"%(cmd, out))
            if int(curver)<version:
                warning = 1
                sys.stderr.write(info%(cmd, curver, version))
                
    message = "Make sure you have installed all dependencies from https://github.com/lpryszcz/redundans#manual-installation !"
    if warning:
        sys.stderr.write("\n%s\n\n"%message)
        sys.exit(1)
    
def main():
    import argparse
    usage   = "%(prog)s -v" #usage=usage, 
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, \
                                      formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", "--verbose",  default=False, action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.13a')
    
    parser.add_argument("-i", "--fastq", nargs="*", default=[], 
                        help="FASTQ PE/MP files")
    parser.add_argument("-f", "--fasta", required=1, 
                        help="FASTA file with contigs")
    parser.add_argument("-o", "--outdir",  default="redundans", 
                        help="output directory [%(default)s]")
    parser.add_argument("-t", "--threads", default=4, type=int, 
                        help="max threads to run [%(default)s]")
    parser.add_argument("--resume",  default=False, action="store_true",
                        help="resume previous run")
    parser.add_argument("--log",           default=sys.stderr, type=argparse.FileType('w'), 
                        help="output log to [stderr]")
    
    redu = parser.add_argument_group('Reduction options')
    redu.add_argument("--identity",        default=0.51, type=float,
                      help="min. identity [%(default)s]")
    redu.add_argument("--overlap",         default=0.66, type=float,
                      help="min. overlap  [%(default)s]")
    redu.add_argument("--minLength",       default=200, type=int, 
                      help="min. contig length [%(default)s]")
    
    scaf = parser.add_argument_group('Scaffolding options')
    scaf.add_argument("-j", "--joins",  default=5, type=int, 
                      help="min pairs to join contigs [%(default)s]")
    scaf.add_argument("-a", "--linkratio", default=0.7, type=float,
                       help="max link ratio between two best contig pairs [%(default)s]")    
    scaf.add_argument("--limit",  default=0.2, type=float, 
                      help="align subset of reads [%(default)s]")
    scaf.add_argument("-q", "--mapq",    default=10, type=int, 
                      help="min mapping quality [%(default)s]")
    scaf.add_argument("--iters",         default=2, type=int, 
                      help="scaffolding iterations per library [%(default)s]")
    scaf.add_argument("-l", "--longreads", nargs="*", default=[], 
                      help="FastQ/FastA files with long reads")
    scaf.add_argument("-r", "--reference", default='', 
                      help="reference FastA file")
    scaf.add_argument("--norearrangements", default=False, action='store_true', 
                      help="high identity mode (rearrangements not allowed)")
    
    gaps = parser.add_argument_group('Gap closing options')
    
    skip = parser.add_argument_group('Skip below steps (all performed by default)')
    skip.add_argument('--noreduction',   action='store_false', default=True)   
    skip.add_argument('--noscaffolding', action='store_false', default=True)   
    skip.add_argument('--nogapclosing',  action='store_false', default=True)   
    skip.add_argument('--nocleaning',    action='store_false', default=True)   
    
    o = parser.parse_args()
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    if o.verbose:
        o.log.write("Options: %s\n"%str(o))

    # check if input files exists
    for fn in [o.fasta,] + o.fastq + o.longreads: 
        if not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)

    # patch sspacebin
    sspacebin = os.path.join(root, "bin/SSPACE/SSPACE_Standard_v3.0.pl")

    # check if all executables exists & in correct versions
    dependencies = {'lastal': 700, 'lastdb': 700, 'bwa': 0, sspacebin: 0, 'GapCloser': 0, 'snap-aligner': 0}
    _check_dependencies(dependencies)
    
    # initialise pipeline
    redundans(o.fastq, o.longreads, o.fasta, o.reference, o.outdir, o.mapq, \
              o.threads, o.resume, o.identity, o.overlap, o.minLength,  \
              o.joins, o.linkratio, o.limit, o.iters, sspacebin, \
              o.noreduction, o.noscaffolding, o.nogapclosing, o.nocleaning, \
              o.norearrangements, o.verbose, o.log)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    #except IOError as e:
    #    log.write("I/O error({0}): {1}\n{2}\n".format(e.errno, e.strerror, str(e)))
    #[Errno 95] Operation not supported ie symlinks over samba or in NFS shares
    #except OSError as e:
    #    log.write("%s\n"%str(e))
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
