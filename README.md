### Table of Contents
- **[Redundans](#redundans)**  
  - **[Prerequisites](#prerequisites)**  
  - **[Running the pipeline](#running-the-pipeline)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[FAQ](#faq)**  
  - **[Citation](#citation)**  

# Redundans

Redundans pipeline assists **an assembly of heterozygous genomes**.  
Program takes as input **assembled contigs**, **paired-end and/or mate pairs 
sequencing libraries** and returns **scaffolded homozygous genome assembly**, 
that should be **less fragmented** and with total **size smaller** than the input contigs. 
In addition, Redundans will automatically **close the gaps** resulting from genome assembly or scaffolding [more details](/test#redundans-pipeline). 

The pipeline consists of three steps/modules: 
- **redundancy reduction**: detection and selectively removal of redundant contigs from an initial *de novo* assembly
- **scaffolding**: joining of genome fragments using paired-end and/or mate-pairs reads
- **gap closing**

Redundans is: 
- **fast** & **lightweight**, multi-core support and memory-optimised, 
so it can be run even on then laptop for small-to-medium size genomes
- **flexible** toward many sequencing technologies (Illumina, 454 or Sanger) and library types (paired-end, mate pairs, fosmids)
- **modular**: every step can be ommited or replaced by another tools

For more information have a look at the [poster](/docs/poster.pdf) or [manuscript](/docs/manuscript.pdf).

![Flowchart](/docs/redundans_flowchart.png)

## Prerequisites
- Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
- [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3) & [LAST](http://last.cbrc.jp/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [SSPACE3](http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
- [GapCloser](http://sourceforge.net/projects/soapdenovo2/files/GapCloser/)

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
                        max threads to run [4]
  --log LOG             output log to [stderr]
```
- Reduction options:
```
  --identity IDENTITY   min. identity [0.51]
  --overlap OVERLAP     min. overlap  [0.66]
  --minLength MINLENGTH
                        min. contig length [200]
```
- Scaffolding options:
```
  -j JOINS, --joins JOINS
                        min k pairs to join contigs [5]
  -l LIMIT, --limit LIMIT
                        align subset of reads [0.2]; this means 0.2*genome size reads will be aligned; so for 100Mb genome, redundans will process 20M reads per library
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

For more details have a look in [test directory](/test). 

## FAQ
### SSPACE fails with an error `Can't locate getopts.pl in @INC`.  
This is due to missing getops in recent Perl. Just do:
```bash
sudo cpan
install Perl4::CoreLibs
```

### Reduction step takes a lot of time.   
Reduction step execute all-vs-all similarity search on your contigs. This may take some time if your assembly is heavily fragmented (>100k contigs).
In order to speed-up the analysis, you can use threads (i.e. `-t 8` for 8 cores). Make sure, you are running Python 2.7, as threading in reduction step is disabled in Python 2.6 and earlier.

### Estimation of my library statistics is incorrect. Can I specify these values manually?   
This can happen for highly fragmented assemblies or poor quality libraries. You can specify library statistics manually. To do so, look for *.is.txt file specific for your library i.e. for `-i 5000_1.fq.gz 5000_2.fq.gz` you will have to enter requested values into `5000_2.fq.gz.is.txt`. Make sure you specify some large number of mates for requested orientation (ie 100,000). For example if you want mate-pairs with RF orientation and 5kb insert size +/- 1.5kb, enter into respective *.is.txt file:
```bash
5000.0        5000.0  1500.0  0       0      100000   0
```

### No alignments error.  
If you see warning messages like the ones below while running the test set: 
```bash
[WARNING] No alignments for test/5000_1.fq.gz - test/5000_2.fq.gz!
[WARNING] No alignments for test/600_1.fq.gz - test/600_2.fq.gz!
```

The problem is because older version of BWA lack MEM algorithm (check it by executing `bwa mem`). If it gives you an error, download the latest [BWA](http://bio-bwa.sourceforge.net/).  

### Redundans fails with `OSError` or `maf-convert: not found`.  
Make sure you are using the latest version of [LAST](http://last.cbrc.jp/) aligner and that all dependencies are accessible through your PATH environmental variable. 

```bash
OSError: [Errno 13] Permission denied

# or
OSError: [Errno 2] No such file or directory

# or
[ERROR] maf-convert: not found
```

### Why does Redundans use two similarity search algorithms, [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3) & [LAST](http://last.cbrc.jp/)?   
BLAT is lightweight & very fast, but lack sensitivity for more diverged sequences. If you specify `--identity` below 0.85, the pipeline will use LAST, that is ~4x slower, but more sensitive than BLAT.
Our simulations shows LAST is capable of correctly reducing heterozygous assemblies with up to 45% divergence between haplotypes.   
To limit speed difference between these two algorithms, LAST **runs in multiple threads**, so using `-t 4` you shouldn't see any difference in runtime between runs for `--identity 0.9` or `--identity 0.5`. Note, this only works in Python 2.7! 

### How is multiple redundancy handled #8 ? 
Redundans removes all contigs, but the longest one, that fullfill identity & overlap critaria during reduction step. 

## Citation
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes. 
