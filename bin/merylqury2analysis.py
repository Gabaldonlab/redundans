#!/usr/bin/env python3

"""
Module in charge of the meryl db creation and the merqury analysis. Readapted from the original Merqury bash wrapper
Author:
Diego Fuentes Palacios
Barcelona 08/18/2022
"""

import os, resource, sys, re
import glob, subprocess, time
from datetime import datetime
from io import TextIOWrapper, StringIO
from traceback import print_list

# update sys.path & environmental PATH
root = os.path.dirname(os.path.abspath(sys.argv[0]))
src = ["bin/", "bin/merqury"]
paths = [os.path.join(root, p) for p in src]
sys.path = paths + sys.path
os.environ["PATH"] = "%s:%s"%(':'.join(paths), os.environ["PATH"])


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
    
    args4 = "meryl union-sum output "+t+" "+m+" "+outname+" "+path

    proc4 = subprocess.Popen(args4, stderr=subprocess.DEVNULL, shell=True)

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
    proc2 = subprocess.Popen(args2, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    kmer_filt = proc2.stdout.read().decode("utf-8").strip()
    proc2.communicate()
    if verbose:
        sys.stderr.write("[INFO] Filtering meryldb with kmers below %s kmer size\n"%kmer_filt)

    args3 = "meryl greater-than %s output %s %s"%(kmer_filt, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt), meryldb)
    args3 = ["meryl", "greater-than", kmer_filt, "output", t, m, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt), meryldb]
    proc3 = subprocess.Popen(args3, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    proc3.communicate()
    #meryl greater-than $filt output $db.gt$filt.meryl $db.meryl

    #Step 2:
    
    if verbose:
        sys.stderr.write("=== Generate spectra-cn plots per assemblies and get QV, k-mer completeness ===\n")

    args4 = ["meryl", k, "count", "output", t, m, assembly_db, outfile]
    proc4 = subprocess.Popen(args4, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
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
    proc5 = subprocess.Popen(args5, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc5.communicate()

    args6 =  "meryl histogram %s %s %s | awk '{print \"read-only\t\"$0}' >> %s"%(t, m, path, os.path.join(mqout, "spectra-cn.hist"))
    proc6 = subprocess.Popen(args6, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
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
        proc7=subprocess.Popen(args7, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        proc7.communicate()
        args8="meryl histogram %s %s %s | awk -v cn=%s \'{print cn\"\t\"$0}\' >> %s"%(t, m, path_i, i, os.path.join(mqout, "spectra-cn.hist"))
        proc8 = subprocess.Popen(args8, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
        proc8.communicate()
        args9=["rm", "-r", path_i]
        proc9 = subprocess.Popen(args9)

    #For copy number >4
    name_cn = "asm.k%s.gt4.%s"%(kmer, os.path.basename(assembly_db))
    path_cn = os.path.join(mqout, name_cn)
    tree_cmd = "[greater-than %s %s ]"%(i, assembly_db)

    args10= "meryl intersect output %s %s %s %s %s"%(t, m, path_cn, meryldb, tree_cmd)
    proc10=subprocess.Popen(args10, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc10.communicate()

    args11="meryl histogram %s %s %s | awk -v cn=\">4\" \'{print cn\"\t\"$0}\' >> %s"%(t, m, path_cn, os.path.join(mqout, "spectra-cn.hist"))
    proc11 = subprocess.Popen(args11, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc11.communicate()

    args12=["rm", "-r", path_cn]
    proc12 = subprocess.Popen(args12)
    
    # Copy numbers in k-mers found only in asm")

    args13 = "meryl difference output %s %s %s %s %s"%(t, m, path, assembly_db, meryldb)
    proc13 = subprocess.Popen(args13, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc13.communicate()

    #Extract the ones that are are distinct and multi

    args14 = "meryl statistics %s %s %s | head -n4 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc14 = subprocess.Popen(args14, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    PRESENT = proc14.stdout.read().decode("utf-8").strip()
    proc14.communicate()

    args15 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc15 = subprocess.Popen(args15, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    DISTINCT = proc15.stdout.read().decode("utf-8").strip()
    proc15.communicate()
    MULTI = str(int(PRESENT)-int(DISTINCT))

    #Store results in a complementary hist file
    cmd = "echo \"1\t0\t%s\" > %s"%(DISTINCT, os.path.join(mqout, "only.hist"))
    subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    cmd = "echo \"2\t0\t%s\" >> %s"%(MULTI, os.path.join(mqout, "only.hist"))
    subprocess.Popen(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


    #Generate the QV statistics
    if verbose:
        sys.stderr.write("Generate the QV statistics here: %s\n"%os.path.join(mqout, (os.path.basename(outfile)+".qv")))

    args16="meryl statistics %s %s %s| head -n4 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc16 = subprocess.Popen(args16, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    ASM_ONLY=proc16.stdout.read().decode("utf-8").strip()
    proc16.communicate()

    args17="meryl statistics %s %s %s | head -n4 | tail -n1 | awk '{print $2}'"%(t, m, assembly_db)
    proc17 = subprocess.Popen(args17, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    TOTAL=proc17.stdout.read().decode("utf-8").strip()
    proc17.communicate()

    args18= "echo \"%s %s\" | awk -v k=%s '{print (1-(1-$1/$2)^(1/k))}'"%(ASM_ONLY, TOTAL, kmer)
    proc18 = subprocess.Popen(args18, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    ERROR=proc18.stdout.read().decode("utf-8").strip()

    args19="echo \"%s %s\" | awk -v k=%s '{print (-10*log(1-(1-$1/$2)^(1/k))/log(10))}'"%(ASM_ONLY, TOTAL, kmer)
    proc19 = subprocess.Popen(args19, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    QV=proc19.stdout.read().decode("utf-8").strip()

    #Store it in QV file
    qv_cmd="echo \"%s\t%s\t%s\t%s\t%s\" >> %s"%(os.path.basename(outfile), ASM_ONLY, TOTAL, QV, ERROR, os.path.join(mqout, (os.path.basename(outfile)+".qv")))
    subprocess.Popen(qv_cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    # Per seq QV statistics"

    args20="meryl-lookup -existence -sequence %s -mers %s/ | awk -v k=%s '{print $1\"\t\"$4\"\t\"$2\"\t\"(-10*log(1-(1-$4/$2)^(1/k))/log(10))\"\t\"(1-(1-$4/$2)^(1/k))}' > %s"%(outfile, path, kmer, os.path.join(mqout, (os.path.basename(outfile)+".perseq.qv")))
    proc20 = subprocess.Popen(args20, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc20.communicate()


    # k-mer completeness (recoveray rate) with solid k-mers for $asm with > $filt counts"

    if verbose:
        sys.stderr.write("Generate the k-mer completeness (recovery rate) with solid k-mers over %s.\n"%kmer)

    args21 = "meryl intersect output %s %s %s %s %s"%(t, m, os.path.join(mqout, "assembly_filtered_k%s.meryl"%kmer_filt), assembly_db, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt))
    proc21 = subprocess.Popen(args21, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc21.communicate()

    args22 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, os.path.join(mqout, "filtered_k%s.complete.meryl"%kmer_filt))
    proc22 = subprocess.Popen(args22, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    TOTAL=proc22.stdout.read().decode("utf-8").strip()

    args23 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, os.path.join(mqout, "assembly_filtered_k%s.meryl"%kmer_filt))
    proc23 = subprocess.Popen(args23, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
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
    proc24 = subprocess.Popen(args24, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc24.communicate()
    
    histname=os.path.join(mqout, os.path.basename(outfile)+".spectra-asm.hist")
    args25= "echo \"Assembly\tkmer_multiplicity\tCount\" > %s"%(histname)
    subprocess.Popen(args25, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    args26 = "meryl histogram %s %s %s | awk '{print \"read-only\t\"$0}' >> %s"%(t, m, path, histname)
    proc26 = subprocess.Popen(args26, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc26.communicate()

    args27 = "meryl histogram %s %s %s | awk -v hap=\"%s\" '{print hap\"\t\"$0}' >> %s"%(t, m, os.path.join(mqout, "assembly.k%s.meryl"%kmer_filt), os.path.basename(outfile), histname)
    proc27 = subprocess.Popen(args27, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, shell=True)
    proc27.communicate()

    # Get asm only for spectra-asm"

    args28 = "meryl statistics %s %s %s | head -n3 | tail -n1 | awk '{print $2}'"%(t, m, path)
    proc28 = subprocess.Popen(args28, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)
    ASM_ONLY=proc28.stdout.read().decode("utf-8").strip()

    hist_only=os.path.join(mqout, os.path.basename(outfile)+".dist_only.hist")
    args25= "echo \"%s\t0\t%s\" >  %s"%(os.path.basename(outfile), ASM_ONLY, hist_only)
    subprocess.Popen(args25, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    if verbose:
        sys.stderr.write("[INFO] Plotting the results in spectra plots here %s\n"%mqout)

    rcommand = "Rscript bin/merqury/plot/plot_spectra_cn.R -f %s -o %s -z %s"%(os.path.join(mqout, "spectra-cn.hist"), os.path.join(mqout, "results.spectra-cn"), os.path.join(mqout, "only.hist"))
    rproc = subprocess.Popen(rcommand, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)

    rcommand2 = "Rscript bin/merqury/plot/plot_spectra_cn.R -f %s -o %s -z %s"%(histname, os.path.join(mqout, "results.spectra-asm"), hist_only)
    rproc2 = subprocess.Popen(rcommand2, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, shell=True)

    rproc.communicate()
    rproc2.communicate()

    return True