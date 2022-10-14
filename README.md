### Table of Contents
- **[Redundans](#redundans)**  
  - **[Prerequisites](#prerequisites)**  
    - **[UNIX installer](#unix-installer)**  
    - **[Docker image](#docker-image)** 
  - **[Running the pipeline](#running-the-pipeline)**  
    - **[Parameters](#parameters)**  
    - **[Test run](#test-run)**  
  - **[Support](#support)**
  - **[Citation](#citation)**  

# Redundans
  
Redundans pipeline assists **an assembly of heterozygous genomes**.  
Program takes [as input](#parameters) **assembled contigs**, **sequencing libraries** and/or **reference sequence** and returns **scaffolded homozygous genome assembly**. Final assembly should be **less fragmented** and with total **size smaller** than the input contigs. In addition, Redundans will automatically **close the gaps** resulting from genome assembly or scaffolding. 

<img align="right" src="/docs/redundans_flowchart.png">

The pipeline consists of several steps (modules):  
1. **de novo contig assembly** (optional if no contigs are given)
2. **redundancy reduction**: detection and selective removal of redundant contigs from an initial *de novo* assembly 
3. **scaffolding**: joining of genome fragments using paired-end reads, mate-pairs, long reads and/or reference chromosomes 
4. **gap closing**: filling the gaps after scaffolding using paired-end and/or mate-pair reads 

Redundans is: 
- **fast** & **lightweight**, multi-core support and memory-optimised, 
so it can be run even on the laptop for small-to-medium size genomes
- **flexible** toward many sequencing technologies (Illumina, 454, Sanger, PacBio & Nanopore) and library types (paired-end, mate pairs, fosmids, long reads)
- **modular**: every step can be omitted or replaced by other tools
- **reliable**: it has been already used to improve genome assemblies varying in size (several Mb to several Gb) and complexity (fungal, animal & plants)

For more information have a look at the [documentation](/docs), [poster](/docs/poster.pdf), [publication](http://nar.oxfordjournals.org/content/44/12/e113), [test dataset](/test) or [manual](http://bit.ly/redundans_manual). 

## Prerequisites
Redundans uses several programs (all provided within this repository): 
- [Platanus](http://platanus.bio.titech.ac.jp/?page_id=14) 
- [LAST](http://last.cbrc.jp/) v800+
- [BWA](http://bio-bwa.sourceforge.net/) v0.7.12+
- [Minimap2](https://github.com/lh3/minimap2)
- [SNAP aligner](https://github.com/amplab/snap)
- [SSPACE3](http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
- [GapCloser](http://sourceforge.net/projects/soapdenovo2/files/GapCloser/)
- [pyScaf](https://github.com/lpryszcz/pyScaf)
- [FastaIndex](https://github.com/lpryszcz/FastaIndex)
- [Meryl](https://github.com/marbl/meryl)
- [Merqury](https://github.com/marbl/merqury)
- [k8](https://github.com/attractivechaos/k8/) v0.2.4+

On most Linux distros, the installation should be as easy as:
```
git clone --recursive https://github.com/lpryszcz/redundans.git
cd redundans && bin/.compile.sh
```

If it fails, make sure you have below dependencies installed: 
- Python >= 3
- Perl [SSPACE3]
- make, gcc & g++ [BWA, Minimap2 & LAST] ie. `sudo apt-get install make gcc g++`
- R >= 3.6, ggplot2, scales, argparser for Merqury ie. `sudo apt-get install r-base r-cran-ggplot2 r-cran-scales r-cran-argparse`
- [zlib including zlib.h headers](http://zlib.net/) [BWA] ie. `sudo apt-get install zlib1g-dev`
- optionally for plotting `numpy` and `matplotlib` ie. `sudo -H pip install -U matplotlib numpy`

For user convenience, we provide [UNIX installer](#unix-installer) and [Docker image](#docker-image), that can be used instead of manually installation.  

### UNIX installer
UNIX installer will automatically fetch, compile and configure Redundans together with all dependencies.
It should work on all modern Linux systems, given Python 2.7, commonly used programmes (ie. wget, curl, git, perl, gcc, g++, ldconfig) and libraries (zlib including zlib.h) are installed. 
```bash
source <(curl -Ls http://bit.ly/redundans_installer)
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

[![](https://images.microbadger.com/badges/version/lpryszcz/redundans.svg)](https://hub.docker.com/r/lpryszcz/redundans/)
![](https://images.microbadger.com/badges/image/lpryszcz/redundans.svg)

## Running the pipeline
Redundans input consists of any combination of:
- **assembled contigs** (FastA)
- **paired-end and/or mate pairs reads** (FastQ*)
- **long reads** (FastQ/FastA*) - both PacBio and Nanopore are supported
- and/or **reference chromosomes/contigs** (FastA). 
* gzipped files are also accepted.

Redundans will return **homozygous genome assembly** in `scaffolds.filled.fa` (FastA).  
In addition, the program reports [statistics for every pipeline step](/test#summary-statistics), including number of contigs that were removed, GC content, N50, N90 and size of gap regions. 

### Parameters
For the user convenience, Redundans is equipped with a wrapper that **automatically estimates run parameters** and executes all steps/modules.
You should specify some sequencing libraries (FastA/FastQ) or reference sequence (FastA) in order to perform scaffolding. 
If you don't specify `-f` **contigs** (FastA), Redundans will assemble contigs *de novo*, but you'll have to provide **paired-end and/or mate pairs reads** (FastQ). 
Most of the pipeline parameters can be adjusted manually (default values are given in square brackets []):  
**HINT**: If you run fails, you may try to resume it, by adding `--resume` parameter. 
- General options:
```
  -h, --help            show this help message and exit
  -v, --verbose         verbose
  --version             show program's version number and exit
  -i FASTQ, --fastq FASTQ
                        FASTQ PE / MP files
  -f FASTA, --fasta FASTA
                        FASTA file with contigs / scaffolds
  -o OUTDIR, --outdir OUTDIR
                        output directory [redundans]
  -t THREADS, --threads THREADS
                        no. of threads to run [4]
  --resume              resume previous run
  --log LOG             output log to [stderr]
  --nocleaning
```
- Reduction options:
```
  --identity IDENTITY   min. identity [0.51]
  --overlap OVERLAP     min. overlap  [0.80]
  --minLength MINLENGTH
                        min. contig length [200]
  --noreduction         Skip reduction
```
- Short-read scaffolding options:
```
  -j JOINS, --joins JOINS
                        min pairs to join contigs [5]
  -a LINKRATIO, --linkratio LINKRATIO
                        max link ratio between two best contig pairs [0.7]
  --limit LIMIT         align subset of reads [0.2]
  -q MAPQ, --mapq MAPQ  min mapping quality [10]
  --iters ITERS         iterations per library [2]
  --noscaffolding       Skip short-read scaffolding
  -b, --usebwa          use bwa mem for alignment [use snap-aligner]
```
- Long-read scaffolding options:
```
  -l LONGREADS, --longreads LONGREADS
                        FastQ/FastA files with long reads
  --useminimap2         Use Minimap2 for aligning long reads. Preset usage dependant on file name convention (case insensitive): ont, nanopore, pb, pacbio, hifi, hi_fi, hi-fi. ie: s324_nanopore.fq.gz.
  --identity IDENTITY   min. identity [0.51]
  --overlap OVERLAP     min. overlap  [0.80]
```
- Reference-based scaffolding options:
```
  -r REFERENCE, --reference REFERENCE
                        reference FastA file
  --norearrangements    high identity mode (rearrangements not allowed)
  --identity IDENTITY   min. identity [0.51]
  --overlap OVERLAP     min. overlap  [0.80]
```
- Gap closing options:
```
  --iters ITERS         iterations per library [2]
  --nogapclosing                        
```
- Meryl and Merqury options:
```
  --nomerqury           Skip meryldb and merqury assembly stats.
  -k KMER, --kmer KMER  K-mer size for meryl [21]
```

Redundans is **extremely flexible**. All steps of the pipeline can be ommited using: `--noreduction`, `--noscaffolding`, `--nogapclosing` and/or `--nomerqury` parameters. 

### Test run
To run the test example, execute: 
```bash
./redundans.py -v -i test/*_?.fq.gz -f test/contigs.fa -o test/run1

# if your run failed for any reason, you can try to resume it
rm test/run1/_sspace.2.1.filled.fa
./redundans.py -v -i test/*_?.fq.gz -f test/contigs.fa -o test/run1 --resume

# if you have no contigs assembled, just run without `-f`
./redundans.py -v -i test/*_?.fq.gz -o test/run.denovo
```

Note, the **order of libraries (`-i/--input`) is not important**, as long as `read1` and `read2` from each library are given one after another 
i.e. `-i 600_1.fq.gz 600_2.fq.gz 5000_1.fq.gz 5000_2.fq.gz` would be interpreted the same as `-i 5000_1.fq.gz 5000_2.fq.gz 600_1.fq.gz 600_2.fq.gz`.

You can play with **any combination of inputs** ie. paired-end, mate pairs, long reads and / or reference-based scaffolding, for example:
```bash
# reduction, scaffolding with paired-end, mate pairs and long reads, and gap closing with paired-end and mate pairs
./redundans.py -v -i test/*_?.fq.gz -l test/pacbio.fq.gz test/nanopore.fa.gz -f test/contigs.fa -o test/run_short_long

# scaffolding and gap closing with paired-end and mate pairs (no reduction)
./redundans.py -v -i test/*_?.fq.gz -f test/contigs.fa -o test/run_short-scaffolding-closing --noreduction

# reduction, reference-based scaffolding and gap closing with paired-end reads (--noscaffolding disables only short-read scaffolding)
./redundans.py -v -i test/600_?.fq.gz -r test/ref.fa -f test/contigs.fa -o test/run_ref_pe-closing --noscaffolding
```

For more details have a look in [test directory](/test). 

## Support 
If you have any issues or doubts check [documentation](/docs) and [FAQ (Frequently Asked Questions)](/docs#faq). 
You may want also to sign to [our forum](https://groups.google.com/d/forum/redundans).

## Citation
Leszek P. Pryszcz and Toni Gabaldón (2016) Redundans: an assembly pipeline for highly heterozygous genomes. NAR. [doi: 10.1093/nar/gkw294](http://nar.oxfordjournals.org/content/44/12/e113)
