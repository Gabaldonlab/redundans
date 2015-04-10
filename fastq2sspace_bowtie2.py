#!/usr/bin/env python
desc="""Align pairs/mates onto contigs and run SSPACE scaffolder.
Example:
fastq2sspace.py -v -f contigs.fa -n pe300 pe600 pe5000 -1 ../../_archives/PL429.{3,6,50}00_read1*.fastq.gz -2 ../../_archives/PL429.{3,6,50}00_read2*.fastq.gz -i 300 600 5000 -s 0.15 0.25 0.5 -t FR FR RF -u 5000000
"""
epilog="""Author:
l.p.pryszcz@gmail.com

19/06/2012 Dublin
"""

import argparse, commands, os, subprocess, sys
from datetime import datetime

def parse_sam( handle ):
    """Return sam tuple"""
    for l in handle:
        l = l.strip()
        if not l or l.startswith('@'):
            continue
        yield l.split('\t')

def get_start_stop( start,length,flag ):
    """Return start-end read boundaries.
    Return end-start if reverse aligned (flag & 16)."""
    if flag & 16:
        end    = start
        start += length
    else:
        end    = start + length
    return start,end

def sam2sspace_tab( inhandle,outhandle,mapqTh=0,verbose=False ):
    """Convert SAM to SSPACE TAB file."""
    i = j = k = pq1 = 0
    sam = parse_sam( inhandle )
    while 1:
        try:
            #read pair sam
            sam1 = sam.next()
            sam2 = sam.next()
            #get variables
            q1,flag1,ref1,start1,mapq1,len1 = sam1[0],int(sam1[1]),sam1[2],int(sam1[3]),int(sam1[4]),len(sam1[9])
            q2,flag2,ref2,start2,mapq2,len2 = sam2[0],int(sam2[1]),sam2[2],int(sam2[3]),int(sam2[4]),len(sam2[9])
            i   += 1
            #gem uses entire fasta header as seq name
            ref1 = ref1.split()[0]
            ref2 = ref2.split()[0]
            #skip 0 quality pair
            if mapqTh:
                if mapq1<mapqTh or mapq2<mapqTh:
                    continue  
            if q1!=q2:
                sys.stderr.write("Warning: Queries has different names: %s vs %s\n" % (q1,q2) )
                continue
            j   += 1
            #skip self matches
            if ref1==ref2:
                continue
            k += 1
            #define start-stop ranges
            start1,end1 = get_start_stop( start1,len1,flag1 )
            start2,end2 = get_start_stop( start2,len2,flag2 )    
            #print output
            outhandle.write( "%s\t%s\t%s\t%s\t%s\t%s\n" % ( ref1,start1,end1,ref2,start2,end2 ) )
        except StopIteration:
            break
    sys.stderr.write( "   %s pairs. %s passed filtering [%.2f%s]. %s in different contigs [%.2f%s].\n" % (i,j,j*100.0/i,'%',k,k*100.0/i,'%') )

def _get_bowtie2_proc( fn1,fn2,ref,maxins,cores,upto,verbose,bufsize=-1):
    """Return bowtie subprocess.
    bufsize: 0 no buffer; 1 buffer one line; -1 set system default.
    """
    fformat = "-q"
    if fn1.endswith(('.fa','.fasta','.fa.gz','.fasta.gz')):
        fformat = "-f"
    bwtArgs = ['bowtie2','--quiet','--very-fast-local','-p',str(cores),'-x',ref, fformat,"-1", fn1, "-2", fn2,"--maxins",str(maxins) ]
    if upto:
        bwtArgs += [ "--upto",str(upto) ]
    if verbose:
        sys.stderr.write( "  %s\n" % " ".join(bwtArgs) )
    #select ids
    bwtProc = subprocess.Popen( bwtArgs,bufsize=bufsize,stdout=subprocess.PIPE )
    return bwtProc
    
def _get_gem_proc( fn1,fn2,ref,maxins,upto,cores,verbose,bufsize=-1):
    """Return bowtie subprocess.
    bufsize: 0 no buffer; 1 buffer one line; -1 set system default.
    """
    preArgs = [ 'fastq2shuffledFastQ.py', fn1, fn2 ]
    if upto:
        preArgs += [ '-u',str(upto) ]    
    gemArgs = [ 'gem-mapper','-I',ref+'.gem','--unique-mapping','--threads',str(cores),'-q','offset-33','--max-insert-size',str(maxins)]#,'2>','/dev/null' ]
    gem2samArgs = [ 'gem-2-sam','-I',ref,'--expect-paired-end-reads','-q','offset-33']#,'2>','/dev/null']
    sam2uniqArgs = [ 'samgem2unique.py', ]
    if verbose:
        sys.stderr.write( "%s | %s | %s | %s\n" % (" ".join(preArgs)," ".join(gemArgs)," ".join(gem2samArgs)," ".join(sam2uniqArgs)) )
    #select ids
    preProc = subprocess.Popen( preArgs,bufsize=bufsize,stdout=subprocess.PIPE )
    gemProc = subprocess.Popen( gemArgs,bufsize=bufsize,stdout=subprocess.PIPE,stdin=preProc.stdout )
    gem2samProc = subprocess.Popen( gem2samArgs,bufsize=bufsize,stdout=subprocess.PIPE,stdin=gemProc.stdout )
    sam2uniqProc = subprocess.Popen( sam2uniqArgs,bufsize=bufsize,stdout=subprocess.PIPE,stdin=gem2samProc.stdout )
    return sam2uniqProc

def get_tab_files( outdir,reffile,libNames,fReadsFnames,rReadsFnames,inserts,iBounds,cores,mapqTh,upto,verbose ):
    """Prepare genome index, align all libs and save TAB file"""
    #create genome index
    ref = reffile.name
    #'''
    idxfn = ref + ".1.bt2"
    if not os.path.isfile( idxfn ):
        cmd = "bowtie2-build %s %s" % (ref,ref)
        if verbose:
            sys.stderr.write( " Creating index...\n  %s\n" % cmd )
        bwtmessage = commands.getoutput( cmd )
    '''
    idxfn = ref + ".gem"
    if not os.path.isfile( idxfn ):
        cmd = "gem-indexer -i %s -o %s" % (ref,ref)
        if verbose:
            sys.stderr.write( " Creating index...\n  %s\n" % cmd )
        bwtmessage = commands.getoutput( cmd )#'''
        
    tabFnames = []
    #process all libs
    for libName,f1,f2,iSize,iFrac in zip( libNames,fReadsFnames,rReadsFnames,inserts,iBounds ):
        if verbose:
            sys.stderr.write( "[%s] [lib] %s\n" % (datetime.ctime(datetime.now()),libName) )
        #define tab output
        outfn = "%s.%s.tab" % ( outdir,libName ) 
        #skip if file exists
        if os.path.isfile( outfn ):
            sys.stderr.write( "  File exists: %s\n" % outfn )
            tabFnames.append( outfn )
            continue
        out   = open( outfn,"w" )
        #define max insert size allowed
        maxins = ( 1.0+iFrac ) * iSize
        #run bowtie2 for all libs        
        proc = _get_bowtie2_proc( f1.name,f2.name,ref,maxins,cores,upto,verbose )
        #proc = _get_gem_proc( f1.name,f2.name,ref,maxins,upto,cores,verbose )
        #parse botwie output
        sam2sspace_tab( proc.stdout,out,mapqTh )
        #close file
        out.close()
        tabFnames.append( outfn )

    return tabFnames
    
def get_libs( outdir,libFn,libNames,tabFnames,inserts,iBounds,orientations,verbose ):
    """Save lib fname and return it's path"""
    lines = []
    #load libs from file
    if libFn:
        if verbose:
            sys.stderr.write( " Reading libs from %s\n" % libFn )
        lines = open(libFn).readlines()
    #add TAB libs
    for t in zip( libNames,tabFnames,inserts,iBounds,orientations ):
        lines.append( "%s\tTAB\t%s\t%s\t%s\t%s\n" % t )

    outfn = "%s.libs.txt" % outdir #os.path.join( outdir,"libs.txt" )
    if verbose:
        sys.stderr.write( " Updated libs saved to: %s\n" % outfn )
    out   = open( outfn,"w" ); out.write( "".join(lines) ); out.close()
    return outfn
    
def main():

    usage   = "%(prog)s [options]" 
    parser  = argparse.ArgumentParser( usage=usage,description=desc,epilog=epilog )

    parser.add_argument("-v", dest="verbose", default=False, action="store_true")
    parser.add_argument("-f", dest="fasta",      required=True, type=file,
                       help="genome fasta        [mandatory]")
    parser.add_argument("-k", dest="minlinks",   default=5, type=int,
                       help="min number of links [%(default)s]")    
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
    
    o = parser.parse_args()
    if o.verbose:
        sys.stderr.write( "Options: %s\n" % str(o) )
    
    if len(o.libnames)*6 != len(o.libnames)+len(o.libFs)+len(o.libRs)+len(o.libIS)+len(o.libISStDev)+len(o.orientations):
        parser.error("Wrong number of arguments!")

    #generate outdirs if out contain dir and dir not exists
    if os.path.dirname(o.out):
        if not os.path.isdir( os.path.dirname(o.out) ):
            os.makedirs( os.path.dirname(o.out) )
    
    #get tab files
    if o.verbose:
        sys.stderr.write("[%s] Generating TAB file(s) for %s library/ies...\n" % (datetime.ctime(datetime.now()),len(o.libnames)) )
    tabFnames = get_tab_files( o.out,o.fasta,o.libnames,o.libFs,o.libRs,o.libIS,o.libISStDev,o.cores,o.mapq,o.upto,o.verbose )

    #generate lib file
    if o.verbose:
        sys.stderr.write("[%s] Generating libraries file...\n" % datetime.ctime(datetime.now()) )
    libFn = get_libs( o.out,o.lib,o.libnames,tabFnames,o.libIS,o.libISStDev,o.orientations,o.verbose )

    #print sspace cmd
    cmd = "perl /users/tg/lpryszcz/src/SSPACE-BASIC-2.0_linux-x86_64/SSPACE_Basic_v2.0.pl -l %s -a 0.7 -k %s -s %s -b %s > %s.sspace.log" % ( libFn,o.minlinks,o.fasta.name,o.out,o.out ); print cmd
    os.system( cmd )
          
if __name__=='__main__': 
    t0  = datetime.now()
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nCtrl-C pressed!   \n")
    dt  = datetime.now()-t0
    sys.stderr.write( "#Time elapsed: %s\n" % dt )
    