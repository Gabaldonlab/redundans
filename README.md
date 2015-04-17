### Table of Contents
**[Redundans](#redundans)**  
 **[Prerequisites](#prerequisites)**  
 **[Running the pipeline](#running-the-pipeline)**  
  **[Parameters](#parameters)**  
  **[Test run](#test-run)**  
 **[FAQ](#faq)**  
 **[Citation](#citation)**  

# Redundans

Redundans pipeline assists **an assembly of heterozygous genomes**.  
Program takes as input **assembled contigs**, **paired-end and/or mate pairs 
sequencing libraries** and returns **scaffolded homozygous genome assembly**, 
that should be **less fragmented** and with total **size smaller** than the input contigs. 
In addition, Redundans will automatically **close the gaps** resulting from genome assembly or scaffolding. 

The pipeline consists of three steps/modules: 
- **redundancy reduction**: detection and selectively removal of redundant contigs from an initial *de novo* assembly
- **scaffolding**: joining of genome fragments using paired-end and/or mate-pairs reads
- **gap closing**

Redundans is: 
- **fast** & **lightweight**, multi-core support and memory-optimised, 
so it can be run even on then laptop for small-to-medium size genomes
- **flexible** toward many sequencing technologies (Illumina, 454 or Sanger) and library types (paired-end, mate pairs, fosmids)
- **modular**: every step can be ommited or replaced by another tools

For more information have a look at the [poster](https://github.com/lpryszcz/redundans/blob/master/docs/poster.pdf).

## Prerequisites
- Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
- [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
- [BWA](http://bio-bwa.sourceforge.net/)
- [SSPACE3](http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
- [Gap2Seq](http://www.cs.helsinki.fi/u/lmsalmel/Gap2Seq/)

## Running the pipeline
Redundans input consists of **assembled contigs** (FastA) and **paired-end and/or mate pairs reads** (FastQ). Gzipped FastQ files are also accepted. 
Redundans will return **homozygous genome assembly** in `scaffolds.filled.fa` (FastA).  
In addition, the program reports statistics for every pipeline step, including number of contigs that were removed, GC content, N50, N90 and size of gap regions.   

### Parameters
For the user convinience, Redundans is equipped with a wrapper that **automatically estimates run parameters** and executes all steps/modules. 
The only parameters required at the runtime are **assembled contigs** (FastA) and **paired-end and/or mate pairs reads** (FastQ).  
Nevertheless, most of the pipeline parameters can be adjusted manually (default values are given in square brackets []):  
- Genral options:
```
  -h, --help            show this help message and exit
  -v                    verbose
  --version             show program's version number and exit
  -i FASTQ [FASTQ ...], --fastq FASTQ [FASTQ ...]
                        FASTQ PE/MP files
  -f FASTA, --fasta FASTA
                        assembly FASTA file
  -o OUTDIR, --outdir OUTDIR
                        output directory [redundans]
  -t THREADS, --threads THREADS
                        max threads to run [2]
  --log LOG             output log to [stderr]
```
- Reduction options:
```
  --identity IDENTITY   min. identity [0.8]
  --overlap OVERLAP     min. overlap  [0.75]
  --minLength MINLENGTH
                        min. contig length [200]
```
- Scaffolding options:
```
  -j JOINS, --joins JOINS
                        min k pairs to join contigs [5]
  -l LIMIT, --limit LIMIT
                        align subset of reads [0.2]
  -q MAPQ, --mapq MAPQ  min mapping quality [10]
  -iters ITERS          scaffolding iterations per library  [2]
  --sspacebin SSPACEBIN
                        SSPACE path  [~/src/SSPACE/SSPACE_Standard_v3.0.pl]
```

Redundans is **extremely flexible**. All steps of the pipeline can be ommited using: `--noreduction`, `--noscaffolding` and/or `--nogapclosing` parameters. 

### Test run
To run the test example, execute: 
```bash
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1 
```

Note, the **order of libraries (`-i/--input`) is not important**, as long as `read1` and `read2` from each library are given one after another 
i.e. `-i 600_1.fq.gz 600_2.fq.gz 5000_1.fq.gz 5000_2.fq.gz` would be interpreted the same as `-i 5000_1.fq.gz 5000_2.fq.gz 600_1.fq.gz 600_2.fq.gz`.

For more details have a look at [test/ directory](/lpryszcz/redundans/tree/master/test). 

## FAQ
- SSPACE fails with an error `Can't locate getopts.pl in @INC`.  
This is due to missing getops in recent Perl. Just do:
```bash
sudo cpan
install Perl4::CoreLibs
```

## Citation
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes 
