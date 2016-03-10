### Table of Contents
- **[Redundans](#redundans)**  
  - **[Prerequisites](#prerequisites)**  
    - **[UNIX installer](#unix-installer)**  
    - **[Docker image](#docker-image)** 
    - **[Manual installation](#manual-installation)**  
  - **[Running the pipeline](#running-the-pipeline)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[Support](#support)**
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
so it can be run even on the laptop for small-to-medium size genomes
- **flexible** toward many sequencing technologies (Illumina, 454 or Sanger) and library types (paired-end, mate pairs, fosmids)
- **modular**: every step can be ommited or replaced by another tools

For more information have a look at the [poster](/docs/poster.pdf) or [manuscript](/docs/manuscript.pdf).

![Flowchart](/docs/redundans_flowchart.png)

## Prerequisites

### UNIX installer
UNIX installer will automatically fetch, compile and configure Redundans together with all dependencies. It should work on most UNIX systems, but was only tested on some platforms.
It will install all dependencies from the scratch, ignoring versions already installed. 
**This is EXPERIMENTAL version, so you may want to create new user for installation process, to avoid data loss!**   
```bash
# sudo adduser test && su test
bash <(curl -Ls http://bit.ly/redundans_installer)
```

### Docker image
First, you  need to install [docker](https://www.docker.com/): `wget -qO- https://get.docker.com/ | sh`  
Then, you can run the test example by executing: 
```bash
# process the data inside the image - all data will be lost at the end
docker run -it -w /root/src/redundans lpryszcz/redundans ./redundans.py -v -i test/{600,5000}_{1,2}.fq.gz -f test/contigs.fa -o test/run1

# if you wish to process local files, you need to mount the volume with -v
## make sure you are in redundans repo directory (containing test/ directory)
docker run -v `pwd`/test:/test:rw -it lpryszcz/redundans /root/src/redundans/redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1
```
Docker images are very handy, but they have certain limitation. 
The most annoying for me is the **lack of autocompletion**, unless you specify the path in host and container in the exactly same manner as in the example above.
In addition, the volume needs to be mounted every time, leading to a bit complex commands. 

### Manual installation
Alternatively, you can download and configure all dependencies manually: 
- Python 2.7+ & dependencies `sudo pip install -U biopython numpy scipy`
 - biopython requires [sqlite3](https://www.sqlite.org/)
- [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
- [LAST](http://last.cbrc.jp/)
- [BWA](http://bio-bwa.sourceforge.net/)
- [SSPACE3](http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
 - SSPACE require Perl; for perl5+ you will need to copy [getopts.pl](http://cpansearch.perl.org/src/GBARR/perl5.005_03/lib/getopts.pl) into SSPACE/dotlib
 - by default Redundans looks for SSPACE in ~/src/SSPACE directory (`--sspacebin` parameter)
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

## Support 
If you have any issues or doubts check [documentation](/docs) and [FAQ (Frequently Asked Questions)](/doc#faq).  
You may want also to sign to [our forum](https://groups.google.com/d/forum/redundans).


## Citation
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes. NAR
