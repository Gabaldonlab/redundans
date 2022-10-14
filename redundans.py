#!/usr/bin/env python3
desc="""Heterozygous genome assembly pipeline. It consists of three steps:
reduction, scaffolding and gap closing.

More info at: http://bit.ly/Redundans

TBA:
- reduction
 - split contigs on likely problems (fasta2qc.py)
 - add exception if lastdb or lastal doesn't finish successfully
- pyScaf short reads
- fastq2sspace.py - maybe prefilter reads that are not worth to be aligned?
 - create sparse index for snap?
 - filter last hits by e-value?
"""
epilog="""Author:
l.p.pryszcz+git@gmail.com

Mizerow/Warsaw/Bratislava/Barcelona, 17/10/2014

Updated to Python3 by Diego Fuentes Palacios
Barcelona 08/18/2022
"""

import os, resource, sys, re
from ssl import PROTOCOL_TLSv1_2
import glob, subprocess, time
from datetime import datetime
from io import TextIOWrapper, StringIO
from traceback import print_list

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["bin/", "bin/bwa/", "bin/snap/", "bin/pyScaf/", "bin/last/build/",
    "bin/last/scripts/", "bin/last/src/", "bin/minimap2/", "bin/merqury"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])

# make sure using Python 3
assert sys.version_info >= (3,), "Only Python 3 is supported!"

from fasta2homozygous import fasta2homozygous
from fastq2sspace import fastq2sspace
from fastq2insert_size import fastq2insert_size
from filterReads import filter_paired
from FastaIndex import FastaIndex, symlink
from pyScaf import LongReadGraph, SyntenyGraph
from denovo import denovo

def timestamp():
    """Return formatted date-time string"""
    return "\n%s\n[%s] "%("#"*50, datetime.ctime(datetime.now()))

def get_libraries(fastq, fasta, mapq=10, threads=4, verbose=1, log=sys.stderr, limit=0,
                  libraries=[], stdfracTh=0.66, maxcfracTh=0.9, genomeFrac=0.05, usebwa=0):
    """Return libraries"""
    # skip if all libs OKish
    ## max stdfrac cannot be larger than stdfracTh in any of the libraries
    if libraries and not [x for x in (max(lib[5]) for lib in libraries) if x>stdfracTh]:
        return libraries
    
    # otherwise process all reads
    if not limit or limit<10e5:
        limit = 10e5
    
    # get libraries statistics using 1% of mapped read limit
    libdata = fastq2insert_size(log, fastq, fasta, mapq, threads, limit/100, genomeFrac, stdfracTh, maxcfracTh, usebwa=usebwa)
    # separate paired-end & mate pairs
    ## also separate 300 and 600 paired-ends
    libraries = []
    # add libraries strating from lowest insert size
    for data in sorted(libdata, key=lambda x: x[3]):
        fq1, fq2, readlen, ismedian, ismean, isstd, pairs, orientation = data
        # add new library set if 
        if not libraries or ismean > 1.5*libraries[-1][4][0]:
            # libnames, readlen, libFs, libRs, orientations, libIS, libISStDev
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
                    joins, linkratio, limit, iters, sspacebin, gapclosing, verbose, usebwa, log, \
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
                    log.write(" iteration %s.%s: %s"%(i, j, FastaIndex(pout).stats()))
                lib = ""
                # run fastq scaffolding
                fastq2sspace(out, open(pout), lib, libnames, libFs, libRs, orients, \
                             libIS, libISStDev, libreadlen, threads, mapq, limit, linkratio, joins, \
                             sspacebin, verbose=0, usebwa=usebwa, log=log)
            # store out info
            pout = out+".fa"
            # link output ie out/_sspace.1.1/_sspace.1.1.scaffolds.fasta --> out/_sspace.1.1.scaffolds.fasta
            targetout = os.path.join(os.path.basename(out), os.path.basename(out+".final.scaffolds.fasta"))
            symlink(targetout, pout)
        # if number of gaps larger than 1%, run gap closer & reduction
        stats     = FastaIndex(pout).stats()
        fastaSize = int(stats.split('\t')[2])
        gapSize   = int(stats.split('\t')[-2])
        if i<len(libraries) and gapclosing and 1.0 * gapSize / fastaSize > 0.02:
            nogapsFname = ".".join(pout.split(".")[:-1]) + ".filled.fa"
            if resume>1 or _corrupted_file(nogapsFname):
                resume += 1
                # close gaps & reduce
                if verbose:
                    log.write("  closing gaps ...\n")
                basename    = "_sspace.%s.%s._gapcloser"%(i, j)
                run_gapclosing(outdir, [libraries[i-1],], nogapsFname, pout, threads, limit, \
                               iters=1, resume=resume, verbose=0, log=log, basename=basename)
            reducedFname = ".".join(pout.split(".")[:-1]) + ".reduced.fa"
            if resume>1 or _corrupted_file(reducedFname):
                if verbose:
                    log.write("  reducing ...\n")
                with open(reducedFname, "w") as out:
                    info = fasta2homozygous(out, open(nogapsFname), identity, overlap, minLength, threads, verbose=0, log=log)
            pout = reducedFname            
        # update library insert size estimation, especially for mate-pairs
        libraries = get_libraries(fastq, pout, mapq, threads, verbose=0,log=log, libraries=libraries, usebwa=usebwa)
    # create symlink to final scaffolds or pout
    symlink(os.path.basename(pout), scaffoldsFname)
    symlink(os.path.basename(pout+".fai"), scaffoldsFname+".fai")
    return libraries, resume

def filter_reads(outdir, fq1, fq2, minlen, maxlen, limit, minqual=10):
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
    
def prepare_gapcloser(outdir, configFn, libFs, libRs, orientations,
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
        fn1, fn2, passed = filter_reads(outdir, fq1, fq2, minlen, maxlen, limit)
        #store config info only if some reads
        if passed:
            config.append(lines%(iSize, reverse_seq, i, fn1, fn2))
        
    # store config only if some libs passed filtering
    if len(config)>1:
        with open(configFn, "w") as out:
            out.write("\n".join(config))
        return True
    
def run_gapclosing(outdir, libraries, nogapsFname, scaffoldsFname,  threads, limit, iters, 
                   resume, verbose, log, basename="_gapcloser", overlap=25, minReadLen=40):
    """Execute gapclosing step."""
    pout = scaffoldsFname
    stop = 1
    for i, (libnames, libFs, libRs, orientations, libIS, libISStDev, libreadlen) in enumerate(libraries, 1):
        # prepare config file and filter reads
        configFn = os.path.join(outdir, "%s.%s.conf"%(basename, i))
        # skip if not suitable libraries
        maxReadLen = max(libreadlen)
        if not prepare_gapcloser(outdir, configFn, libFs, libRs, orientations, libIS, libISStDev, \
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
                    log.write(" iteration %s.%s: %s"%(i, j, FastaIndex(pout).stats()))
                # run GapCloser and save stdout/err to log file
                with open(out+".log", "w") as gapcloselog:
                    GapCloser = subprocess.Popen(cmd, stdout=gapcloselog, stderr=gapcloselog)
                    GapCloser.wait()
            # store out info
            pout = out
            # skip if number of gaps smaller than 0.1%
            stats     = FastaIndex(pout).stats()
            fastaSize = int(stats.split('\t')[2])
            gapSize   = int(stats.split('\t')[-2])
            if 1.0 * gapSize / fastaSize < 0.001:
                stop = 1
                break
        if stop:
            break
    # create symlink to final scaffolds or pout
    symlink(os.path.basename(pout), nogapsFname)
    symlink(os.path.basename(pout+".fai"), nogapsFname+".fai")
    return resume

def _corrupted_file(fname):
    """Return True if output file doesn't exists or is corrupted."""
    if not os.path.isfile(fname) or not os.path.getsize(fname) or \
       not os.path.isfile(fname+".fai") or not os.path.getsize(fname+".fai"):
        return True

def _check_fasta(lastOutFn, minSize=1000, log=sys.stderr):
    """Exit if empty FastA file"""
    faidx = FastaIndex(open(lastOutFn))
    if not len(faidx) or faidx.genomeSize<minSize:
        log.write("[ERROR] Empty FastA file encountered: %s !\n"%lastOutFn)
        sys.exit(1)

def _build_meryldb(outpath, fastq, threads, mem, kmer=21):
    """Build meryldb using Illumina reads"""
    

    outname = os.path.join(outpath, "complete.meryl")

    for fq1, fq2 in zip(fastq[0::2], fastq[1::2]):
        # meryl k=$k count output read$i.meryl read$i.fastq.gz
        k = "k=%s"%kmer

        #First copy fastq files to outpath to generate the meryl db there. Those files are going to be deleted later
        
        args0 = ["cp", "-t", outpath, fq1, fq2]
        proc0 = subprocess.Popen(args0, stderr=sys.stderr)
        proc0.communicate()

        basefq1=os.path.basename(fq1)
        basefq2=os.path.basename(fq2)

        fq1 = os.path.join(outpath, basefq1)
        fq2 = os.path.join(outpath, basefq2)
        meryldb1 = os.path.join(outpath, "%s.meryl"%basefq1)
        meryldb2 = os.path.join(outpath, "%s.meryl"%basefq2)
        t = "threads=%s"%threads
        m = "memory=%s"%mem

        #Generate a meryl db for each reads file
        args1 = ["meryl", k, "count", "output", t, m, meryldb1, fq1]
        proc1 = subprocess.Popen(args1, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        args2 = ["meryl", k, "count", "output", t, m, meryldb2, fq2]
        proc2 = subprocess.Popen(args2, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

        #Communicate to enforce process has ended to avoid memory/thread issues
        proc1.communicate()
        proc2.communicate()


    path = os.path.join(outpath, "*.meryl")
    
    #Tried with code below, unknown error
    #args4 = ["meryl", "union-sum", "output", t, m, outname, path]
    args4 = "meryl union-sum output "+t+" "+m+" "+outname+" "+path

    proc4 = subprocess.Popen(args4, stderr=sys.stderr, shell=True)

    proc4.communicate()

    return outname

def merqury_statistics(outdir, meryldb, outfile, threads, mem, kmer=21, verbose=0):

    """This function recreates merqury functionality in python instead of a bash wrapper to apply subprocesses"""

    #####Define base variables
    mqout = os.path.join(outdir, "merqury_results")
    t = "threads=%s"%threads
    m = "memory=%s"%mem
    k = "k=%s"%kmer
    fname = os.path.basename(outfile)+".meryl"
    assembly_db = os.path.join(mqout, fname)

    ##Generate results directory
    args_pre = ["mkdir", mqout]
    proc_pre = subprocess.Popen(args_pre)

    #Glob all read meryl + complete.meryl dbs and create a symbolic link
    pattern = os.path.dirname(os.path.abspath(meryldb))+"/*.meryl"
    list_meryl = glob.glob(pattern) 
    
    for readdb in list_meryl:
        args_tmp = ["ln", "-s", os.path.abspath(readdb), mqout]
        proc_tmp = subprocess.Popen(args_tmp)

    #Recreating merqury:

    #Step 1: find filter kmer size and filter complete.meryl

    args0= "meryl histogram %s %s %s > %s"%(t, m, meryldb, os.path.join(mqout, "filt.hist"))
    proc0 = subprocess.Popen(args0, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc0.communicate()

    #java -jar -Xmx1g $MERQURY/eval/kmerHistToPloidyDepth.ja
    
    args1 = "java -jar -Xmx1g bin/merqury/eval/kmerHistToPloidyDepth.jar %s > %s"%(os.path.join(mqout, "filt.hist"), os.path.join(mqout, ".result"))
    proc1 = subprocess.Popen(args1, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc1.communicate()

    #arg2 = ["java", "-jar", "-Xmx1g","bin/merqury/eval/kmerHistToPloidyDepth.jar"]
    #proc2 = subprocess.Popen(arg2, stdout=subprocess.PIPE, stdin=proc1.stdout, stderr=sys.stderr)
    #proc2.communicate()

    args2 = "sed -n 2p %s | awk '{print $NF}'"%os.path.join(mqout, ".result")
    proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    kmer_filt = proc2.stdout.read().decode("utf-8").strip()
    proc2.communicate()
    if verbose:
        sys.stderr.write("[INFO] Filtering meryldb with kmers below %s kmer size\n"%kmer_filt)

    args3 = "meryl greater-than %s output %s %s"%(kmer_filt, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt), meryldb)
    args3 = ["meryl", "greater-than", kmer_filt, "output", t, m, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt), meryldb]
    proc3 = subprocess.Popen(args3, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    proc3.communicate()
    #meryl greater-than $filt output $db.gt$filt.meryl $db.meryl

    #Step 2:
    
    if verbose:
        sys.stderr.write("=== Generate spectra-cn plots per assemblies and get QV, k-mer completeness ===\n")

    args4 = ["meryl", k, "count", "output", t, m, assembly_db, outfile]
    proc4 = subprocess.Popen(args4, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    proc4.communicate()


    #Generate the spectra hist
    #with open(os.path.join(mqout, "spectra-cn.hist"), "a+") as file:
        #file.write("Copies\tkmer_multiplicity\tCount")

    #First the ones specific to reads
    name = "asm.k%s.0.%s"%(kmer, os.path.basename(assembly_db))
    path = os.path.join(mqout, name)

    cmd = "echo \"Copies\tkmer_multiplicity\tCount\" > %s"%os.path.join(mqout, "spectra-cn.hist")
    subprocess.Popen(cmd, shell=True)

    args5 = "meryl difference output %s %s %s %s %s"%(t, m, path, meryldb, assembly_db)
    proc5 = subprocess.Popen(args5, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc5.communicate()

    args6 =  "meryl histogram %s %s %s | awk '{print \"read-only\t\"$0}' >> %s"%(t, m, path, os.path.join(mqout, "spectra-cn.hist"))
    proc6 = subprocess.Popen(args6, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc6.communicate()

    #group kmers based on copy number from 1 to 4
    if verbose:
        sys.stderr.write("Group kmers based on copy number from 1 to 4 and larger than 4.\n")
    for i in range(1, 5, 1):
        name_i = "asm.k%s.%s.%s"%(kmer, i, os.path.basename(assembly_db))
        path_i = os.path.join(mqout, name_i)
        tree_cmd = "[equal-to %s %s ]"%(i, assembly_db)
        #args7=["meryl", "intersect", t, m, "output", path_i, meryldb, tree_cmd]
        args7= "meryl intersect output %s %s %s %s %s"%(t, m, path_i, meryldb, tree_cmd)
        proc7=subprocess.Popen(args7, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
        proc7.communicate()
        args8="meryl histogram %s %s %s | awk -v cn=%s \'{print cn\"\t\"$0}\' >> %s"%(t, m, path_i, i, os.path.join(mqout, "spectra-cn.hist"))
        proc8 = subprocess.Popen(args8, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
        proc8.communicate()
        args9=["rm", "-r", path_i]
        proc9 = subprocess.Popen(args9)

    #For copy number >4
    name_cn = "asm.k%s.gt4.%s"%(kmer, os.path.basename(assembly_db))
    path_cn = os.path.join(mqout, name_cn)
    tree_cmd = "[greater-than %s %s ]"%(i, assembly_db)

    args10= "meryl intersect output %s %s %s %s %s"%(t, m, path_cn, meryldb, tree_cmd)
    proc10=subprocess.Popen(args10, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc10.communicate()

    args11="meryl histogram %s %s %s | awk -v cn=\">4\" \'{print cn\"\t\"$0}\' >> %s"%(t, m, path_cn, os.path.join(mqout, "spectra-cn.hist"))
    proc11 = subprocess.Popen(args11, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc11.communicate()

    args12=["rm", "-r", path_cn]
    proc12 = subprocess.Popen(args12)
    
    # Copy numbers in k-mers found only in asm")

    args13 = "meryl difference output %s %s %s %s %s"%(t, m, path, assembly_db, meryldb)
    proc13 = subprocess.Popen(args13, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc13.communicate()

    #Extract the ones that are are distinct and multi

    args14 = "meryl statistics %s %s %s | head -n4 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc14 = subprocess.Popen(args14, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    PRESENT = proc14.stdout.read().decode("utf-8").strip()
    proc14.communicate()

    args15 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc15 = subprocess.Popen(args15, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    DISTINCT = proc15.stdout.read().decode("utf-8").strip()
    proc15.communicate()
    MULTI = str(int(PRESENT)-int(DISTINCT))

    #Store results in a complementary hist file
    cmd = "echo \"1\t0\t%s\" > %s"%(DISTINCT, os.path.join(mqout, "only.hist"))
    subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    cmd = "echo \"2\t0\t%s\" >> %s"%(MULTI, os.path.join(mqout, "only.hist"))
    subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


    #Generate the QV statistics
    if verbose:
        sys.stderr.write("Generate the QV statistics here: %s\n"%os.path.join(mqout, (os.path.basename(outfile)+".qv")))

    args16="meryl statistics %s %s %s| head -n4 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc16 = subprocess.Popen(args16, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    ASM_ONLY=proc16.stdout.read().decode("utf-8").strip()
    proc16.communicate()

    args17="meryl statistics %s %s %s | head -n4 | tail -n1 | awk '{print $2}'"%(t, m, assembly_db)
    proc17 = subprocess.Popen(args17, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    TOTAL=proc17.stdout.read().decode("utf-8").strip()
    proc17.communicate()

    args18= "echo \"%s %s\" | awk -v k=%s '{print (1-(1-$1/$2)^(1/k))}'"%(ASM_ONLY, TOTAL, kmer)
    proc18 = subprocess.Popen(args18, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    ERROR=proc18.stdout.read().decode("utf-8").strip()

    args19="echo \"%s %s\" | awk -v k=%s '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'"%(ASM_ONLY, TOTAL, kmer)
    proc19 = subprocess.Popen(args19, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    QV=proc19.stdout.read().decode("utf-8").strip()

    #Store it in QV file
    qv_cmd="echo \"%s\t%s\t%s\t%s\t%s\" >> %s"%(os.path.basename(outfile), ASM_ONLY, TOTAL, QV, ERROR, os.path.join(mqout, (os.path.basename(outfile)+".qv")))
    subprocess.Popen(qv_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    # Per seq QV statistics"

    args20="meryl-lookup -existence -sequence %s -mers %s/ | awk -v k=%s '{print $1\"\t\"$4\"\t\"$2\"\t\"(-10*log(1-(1-$4/$2)^(1/k))/log(10))\"\t\"(1-(1-$4/$2)^(1/k))}' > %s"%(outfile, path, kmer, os.path.join(mqout, (os.path.basename(outfile)+".perseq.qv")))
    proc20 = subprocess.Popen(args20, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc20.communicate()


    # k-mer completeness (recoveray rate) with solid k-mers for $asm with > $filt counts"

    if verbose:
        sys.stderr.write("Generate the k-mer completeness (recovery rate) with solid k-mers over %s.\n"%kmer)

    args21 = "meryl intersect output %s %s %s %s %s"%(t, m, os.path.join(mqout, "assembly_filtered_k%s.meryl"%kmer_filt), assembly_db, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt))
    proc21 = subprocess.Popen(args21, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc21.communicate()

    args22 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt))
    proc22 = subprocess.Popen(args22, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    TOTAL=proc22.stdout.read().decode("utf-8").strip()

    args23 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, os.path.join(mqout, "assembly_filtered_k%s.meryl"%kmer_filt))
    proc23 = subprocess.Popen(args23, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    ASM=proc23.stdout.read().decode("utf-8").strip()

    try:
        completeness_cmd = "echo \"%s\tall\t%s\t%s\" | awk '{print $0\"\t\"((100*$3)/$4)}' >> %s"%(os.path.basename(outfile), ASM, TOTAL, os.path.join(mqout, os.path.basename(outfile)+".completeness.stats"))
        subprocess.Popen(completeness_cmd, shell=True)
    except ValueError as e:
        sys.stderr.write("[WARNING] Couldn't generate kmer completeness plot")
        pass


    cmd_rm=["rm", "-r", os.path.join(mqout, "assembly_filtered_k%s.meryl"%kmer_filt)]
    subprocess.Popen(cmd_rm)

    #Generate the assembly histogram for plotting

    args24 = "meryl intersect output %s %s %s %s %s"%(t, m, os.path.join(mqout, "assembly.k%s.meryl"%kmer_filt), meryldb, assembly_db)
    proc24 = subprocess.Popen(args24, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc24.communicate()
    
    histname=os.path.join(mqout, os.path.basename(outfile)+".spectra-asm.hist")
    args25= "echo \"Assembly\tkmer_multiplicity\tCount\" > %s"%(histname)
    subprocess.Popen(args25, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    args26 = "meryl histogram %s %s %s | awk '{print \"read-only\t\"$0}' >> %s"%(t, m, path, histname)
    proc26 = subprocess.Popen(args26, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc26.communicate()

    args27 = "meryl histogram %s %s %s | awk -v hap=\"%s\" '{print hap\"\t\"$0}' >> %s"%(t, m, os.path.join(mqout, "assembly.k%s.meryl"%kmer_filt), os.path.basename(outfile), histname)
    proc27 = subprocess.Popen(args27, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, shell=True)
    proc27.communicate()

    # Get asm only for spectra-asm"

    args28 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc28 = subprocess.Popen(args28, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    ASM_ONLY=proc28.stdout.read().decode("utf-8").strip()

    hist_only=os.path.join(mqout, os.path.basename(outfile)+".dist_only.hist")
    args25= "echo \"%s\t0\t%s\" >  %s"%(os.path.basename(outfile), ASM_ONLY, hist_only)
    subprocess.Popen(args25, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)

    if verbose:
        sys.stderr.write("[INFO] Plotting the results in spectra plots here %s\n"%mqout)

    rcommand = "Rscript bin/merqury/plot/plot_spectra_cn.R -f %s -o %s -z %s"%(os.path.join(mqout, "spectra-cn.hist"), os.path.join(mqout, "results.spectra-cn"), os.path.join(mqout, "only.hist"))
    rproc = subprocess.Popen(rcommand, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)

    rcommand2 = "Rscript bin/merqury/plot/plot_spectra_cn.R -f %s -o %s -z %s"%(histname, os.path.join(mqout, "results.spectra-asm"), hist_only)
    rproc2 = subprocess.Popen(rcommand2, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)

    rproc.communicate()
    rproc2.communicate()

    return True


#TODO Add a check for R version, if fails, notify.
        
def redundans(fastq, longreads, fasta, reference, outdir, mapq, 
              threads, mem, resume, identity, overlap, minLength, \
              joins, linkratio, readLimit, iters, sspacebin, refpreset, \
              reduction=1, scaffolding=1, gapclosing=1, usemerqury=1, kmer=21, cleaning=1, \
              norearrangements=0, verbose=1, usebwa=0, useminimap2=0, log=sys.stderr, tmp="/tmp"):
    """Launch redundans pipeline."""
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

    #return True

    # DE NOVO CONTIGS
    lastOutFn = os.path.join(outdir, "contigs.fa") 
    if not fasta and _corrupted_file(lastOutFn):
        resume += 1
        if verbose:
            log.write("%sDe novo assembly...\n"%timestamp())        
        fasta = denovo(os.path.join(outdir, "denovo"), fastq, threads, mem, verbose, log, tmp)
    elif not fasta:
        fasta = lastOutFn
        
    # REDUCTION
    fastas = [fasta, ]; _check_fasta(fasta)
    symlink(fasta, lastOutFn)
    fastas.append(lastOutFn); _check_fasta(lastOutFn)
    # update fasta list
    outfn = os.path.join(outdir, "contigs.reduced.fa")
    if reduction and _corrupted_file(outfn):
        resume += 1
        if verbose:
            log.write("%sReduction...\n"%timestamp())
            log.write("#file name\tgenome size\tcontigs\theterozygous size\t[%]\theterozygous contigs\t[%]\tidentity [%]\tpossible joins\thomozygous size\t[%]\thomozygous contigs\t[%]\n")
        with open(outfn, "w") as out:
            info = fasta2homozygous(out, open(fastas[-1]), identity, overlap, minLength, threads, verbose=0, log=log)
        # update fasta list
        lastOutFn = outfn
        fastas.append(lastOutFn); _check_fasta(lastOutFn)


    # get read limit & libraries
    if fastq:
        if verbose:
            log.write("%sEstimating parameters of libraries...\n"%timestamp())
        limit     = get_read_limit(lastOutFn, readLimit, verbose, log)
        libraries = get_libraries(fastq, lastOutFn, mapq, threads, verbose, log, usebwa=usebwa)
    
    # SCAFFOLDING
    outfn = os.path.join(outdir, "scaffolds.fa")
    if fastq and scaffolding: 
        if verbose:
            log.write("%sScaffolding...\n"%timestamp())
        libraries, resume = run_scaffolding(outdir, outfn, fastq, libraries, lastOutFn, mapq, threads, joins, \
                                            linkratio, limit, iters, sspacebin, gapclosing, verbose, usebwa, log, \
                                            identity, overlap, minLength, resume)
        # update fasta list
        fastas += [x for x in sorted(glob.glob(os.path.join(outdir, "_sspace.*.fa"))) if "_gapcloser" not in x]
        lastOutFn = outfn
        fastas.append(lastOutFn); _check_fasta(lastOutFn)


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
                log.write(" iteration %s...\n"%i)
            
            #Add a check to change preset based on filename if provided a list of files:
            
            if re.search("ont", fname, flags=re.IGNORECASE) or re.search("nanopore", fname, flags=re.IGNORECASE) or re.search("oxford", fname, flags=re.IGNORECASE):
                preset = "map-ont"
            elif re.search("pb", fname, flags=re.IGNORECASE) or re.search("pacbio", fname, flags=re.IGNORECASE) or re.search("smrt", fname, flags=re.IGNORECASE):
                preset = "map-pb"
            elif re.search("hifi", fname, flags=re.IGNORECASE) or re.search("hi_fi", fname, flags=re.IGNORECASE) or re.search("hi-fi", fname, flags=re.IGNORECASE):
                preset = "map-hifi"
            else:
                #Added this to default to ONT
                preset = "map-ont"

            if useminimap2 and verbose:
                log.write("Using minimap2 preset %s for file %s...\n"%(preset, fname))
            elif verbose:
                log.write("Using LAST for file %s...\n"%fname)

            s = LongReadGraph(lastOutFn, fname, identity, overlap, preset=preset, useminimap2=useminimap2, maxgap=0, threads=threads, \
                        dotplot="", norearrangements=norearrangements, log=0)
            #else:
            #    s = LongReadGraph(lastOutFn, fname, identity, overlap, maxgap=0, threads=threads, \
            #                  dotplot="", norearrangements=norearrangements, log=0)
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
        fastas.append(lastOutFn); _check_fasta(lastOutFn)

    # REFERENCE-BASED SCAFFOLDING
    outfn = os.path.join(outdir, "scaffolds.ref.fa")
    if reference and _corrupted_file(outfn):
        resume += 1
        if verbose:
            log.write("%sScaffolding based on reference...\n"%timestamp())        
        s = SyntenyGraph(lastOutFn, reference, identity=0.51, overlap=0.66, maxgap=0, threads=threads, \
                         dotplot="", norearrangements=norearrangements, useminimap2=useminimap2, preset=refpreset, log=0)
        # save output
        with open(outfn, "w") as out:
            s.save(out)
        # update fasta list
        lastOutFn = outfn
        fastas.append(lastOutFn); _check_fasta(lastOutFn)

    
        
    # GAP CLOSING
    outfn = os.path.join(outdir, "scaffolds.filled.fa")
    if fastq and gapclosing: 
        if verbose: 
            log.write("%sGap closing...\n"%timestamp())
        resume = run_gapclosing(outdir, libraries, outfn, lastOutFn, threads, limit, iters, resume, verbose, log)
        # update fasta list
        fastas += sorted(glob.glob(os.path.join(outdir, "_gap*.fa")))
        lastOutFn = outfn
        fastas.append(lastOutFn); _check_fasta(lastOutFn)

    # FINAL REDUCTION
    outfn = os.path.join(outdir, "scaffolds.reduced.fa")
    if reduction and _corrupted_file(outfn):
        resume += 1
        if verbose:
            log.write("%sFinal reduction...\n"%timestamp())
            log.write("#file name\tgenome size\tcontigs\theterozygous size\t[%]\theterozygous contigs\t[%]\tidentity [%]\tpossible joins\thomozygous size\t[%]\thomozygous contigs\t[%]\n")
        # reduce
        with open(outfn, "w") as out:
            info = fasta2homozygous(out, open(lastOutFn), identity, overlap,  minLength, threads, verbose=0, log=log)
        # update fasta list
        lastOutFn = outfn
        fastas.append(lastOutFn); _check_fasta(lastOutFn)

    # MERQURY ANALYSIS

    if usemerqury:

        # MERYL DB
        if verbose:
            log.write("\n%sGenerating a complete Meryl database for your samples...\n"%(timestamp()))
        meryldb = _build_meryldb(outdir, fastq, threads, mem, kmer)

        if verbose:
            log.write("%sGenerating Merqury statistics...\n"%timestamp())
        merqury_statistics(outdir, meryldb, outfn, threads, mem, kmer, verbose)
        
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
            #If you want to check meryl data, add '.merylData', '.merylIndex', 'merylIndex'
            endings = ('.fa', '.fasta', '.fai', '.tsv', '.qv', '.stats', '.png', '.log', '.hist')
            for i, fn in enumerate([x for x in fnames if not x.endswith(endings)], 1):
                os.unlink(os.path.join(root, fn))
            # rmdir of snap index
            if root.endswith('.snap') and i==len(fnames):
                os.rmdir(root)

    if orgresume:
        log.write("%sResume report: %s step(s) have been recalculated.\n"%(timestamp(), resume-1))
    
def _check_executable(cmd):
    """Check if executable exists."""
    p = subprocess.Popen("type " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout = p.stdout.readlines()[0].decode("utf-8")
    return "".join(stdout)

def _check_dependencies(dependencies):
    """Return error if wrong software version"""
    warning = 0
    # check dependencies
    info = "[WARNING] Old version of %s: %s. Update to version %s+!\n"
    for cmd, version in list(dependencies.items()):
        #print("Checking %s version %s"%(cmd, version))
        out = _check_executable(cmd)
        if "not found" in out:
            warning = 1
            sys.stderr.write("[ERROR] %s\n"%out)
        elif version:
            p = subprocess.Popen([cmd, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            if cmd == "R":
                out = p.stdout.readline().decode("utf-8")
                _, _, ver, _ = out.split(" ", 3)
                curver = "".join(ver.split(".", 2))
                #For clarity in the warning message
                out = ver
            else:
                out = "".join(p.stdout.readlines()[0].decode("utf-8"))
                curver = out.split()[-1]
            if not curver.isdigit():
                warning = 1
                sys.stderr.write("[WARNING] Problem checking %s version: %s\n"%(cmd, out))
            elif int(curver)<version:
                warning = 1
                sys.stderr.write(info%(cmd, curver, version))
                
    message = "Make sure you have installed all dependencies from https://github.com/lpryszcz/redundans#manual-installation !"
    if warning:
        sys.stderr.write("\n%s\n\n"%message)
        sys.exit(1)
    
def main():
    import argparse
    usage   = "%(prog)s -v"
    parser  = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
  
    parser.add_argument("-v", "--verbose",  action="store_true", help="verbose")    
    parser.add_argument('--version', action='version', version='0.14a')
    parser.add_argument("-i", "--fastq", nargs="*", default=[], help="FASTQ PE / MP files")
    parser.add_argument("-f", "--fasta", default="", help="FASTA file with contigs / scaffolds")
    parser.add_argument("-o", "--outdir", default="redundans", help="output directory [%(default)s]")
    parser.add_argument("-t", "--threads", default=16, type=int, help="max threads to run [%(default)s]")
    parser.add_argument("--resume",  default=False, action="store_true", help="resume previous run")
    parser.add_argument("--log", default=sys.stderr, type=argparse.FileType('w'), help="output log to [stderr]")
    parser.add_argument('--nocleaning', action='store_false', help="keep intermediate files")
    
    denovo = parser.add_argument_group('De novo assembly options')
    denovo.add_argument("-m", "--mem", default=2, type=int, help="max memory to allocate (in GB) [%(default)s]")
    denovo.add_argument("--tmp", default='/tmp', help="tmp directory [%(default)s]")
    
    redu = parser.add_argument_group('Reduction options')
    redu.add_argument("--identity", default=0.51, type=float, help="min. identity [%(default)s]")
    redu.add_argument("--overlap", default=0.80, type=float, help="min. overlap [%(default)s]")
    redu.add_argument("--minLength", default=200, type=int, help="min. contig length [%(default)s]")
    redu.add_argument('--noreduction', action='store_false', help="Skip reduction")
    
    scaf = parser.add_argument_group('Short-read scaffolding options')
    scaf.add_argument("-j", "--joins", default=5, type=int, help="min pairs to join contigs [%(default)s]")
    scaf.add_argument("-a", "--linkratio", default=0.7, type=float,
                      help="max link ratio between two best contig pairs [%(default)s]")    
    scaf.add_argument("--limit", default=0.2, type=float, help="align subset of reads [%(default)s]")
    scaf.add_argument("-q", "--mapq", default=10, type=int, help="min mapping quality [%(default)s]")
    scaf.add_argument("--iters", default=2, type=int, help="iterations per library [%(default)s]")
    scaf.add_argument('--noscaffolding', action='store_false', help="Skip short-read scaffolding")
    scaf.add_argument("-b", "--usebwa", action='store_true', help="use bwa mem for alignment [use snap-aligner]")
     
    longscaf = parser.add_argument_group('Long-read scaffolding options')
    longscaf.add_argument("-l", "--longreads", nargs="*", default=[], help="FastQ/FastA files with long reads. By default LAST")
    longscaf.add_argument("--useminimap2", action='store_true', help="Use Minimap2 for aligning long reads. Preset usage dependant on file name convention (case insensitive): ont, nanopore, pb, pacbio, hifi, hi_fi, hi-fi. ie: s324_nanopore.fq.gz.")
    
    refscaf = parser.add_argument_group('Reference-based scaffolding options')
    refscaf.add_argument("-r", "--reference", default='', help="reference FastA file")
    refscaf.add_argument("--norearrangements", default=False, action='store_true', 
                         help="high identity mode (rearrangements not allowed)")
    refscaf.add_argument("-p", "--preset", default='asm10', help="Preset option for minimap2 Reference-based scaffolding. Possible options: asm5 (5%. sequence divergence), asm10 (10%. sequence divergence) and asm20(20%. sequence divergence). Default [%(default)s]")
    
    gaps = parser.add_argument_group('Gap closing options')
    gaps.add_argument('--nogapclosing',  action='store_false', default=True)

    stats = parser.add_argument_group('Meryl and Merqury options')
    stats.add_argument('--nomerqury',  action='store_false', default=True, help="Skip meryldb and merqury assembly stats.")
    stats.add_argument('-k', "--kmer",  default=21, type=int, help="K-mer size for meryl [%(default)s]")    
        
    # print help if no parameters
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
        
    o = parser.parse_args()
    if o.verbose:
        o.log.write("Options: %s\n"%str(o))

    # need contigs or at least PE or MP libs to compute those
    if not o.fasta and not o.fastq:
        sys.stderr.write("Provide contigs and/or paired-end/mate pairs libraries\n")
        sys.exit(1)
        
    # check if input files exists [o.fasta,] +
    for fn in o.fastq + o.longreads: 
        if not os.path.isfile(fn):
            sys.stderr.write("No such file: %s\n"%fn)
            sys.exit(1)

    # patch sspacebin
    sspacebin = os.path.join(root, "bin/SSPACE/SSPACE_Standard_v3.0.pl")

    # check if all executables exists & in correct versions
    #If using merqury, check for R version, else do not bother
    if o.nomerqury:
        dependencies = {'lastal': 800, 'lastdb': 800, 'GapCloser': 0, 'paste': 0, 'tr': 0, 'zcat': 0, 'platanus': 0, 'R' : 360}
    else:
        dependencies = {'lastal': 800, 'lastdb': 800, 'GapCloser': 0, 'paste': 0, 'tr': 0, 'zcat': 0, 'platanus': 0}
    _check_dependencies(dependencies)
    
    # initialise pipeline
    redundans(o.fastq, o.longreads, o.fasta, o.reference, o.outdir, o.mapq, \
              o.threads, o.mem, o.resume, o.identity, o.overlap, o.minLength,  \
              o.joins, o.linkratio, o.limit, o.iters, sspacebin, \
              o.preset, o.noreduction, o.noscaffolding, o.nogapclosing, o.nomerqury, o.kmer, o.nocleaning, \
              o.norearrangements, o.verbose, o.usebwa, o.useminimap2, o.log, o.tmp)

if __name__=='__main__': 
    t0 = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!      \n")
    dt = datetime.now()-t0
    sys.stderr.write("#Time elapsed: %s\n"%dt)
