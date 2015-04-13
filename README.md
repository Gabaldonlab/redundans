# Redundans

Redundans pipeline **assists assembly of heterozygous genomes**. Redundans takes as input **assembled contigs**, **paired-end and/or mate pairs sequencing libraries**. Redundans will return **homozygous genome assembly** (FastA), that should be **less fragmented** and with total **size smaller** than the input contigs. In addition, Redundans will automatically **close the gaps** resulting from genome assembly or scaffolding. 

It consists of three steps/modules: 

- redundancy reduction: detection and selectively removal of redundant contigs from an initial **de novo** assembly
- scaffolding: joining of genome fragments using paired-end and/or mate-pairs reads
- gap closing

Redundans is: 

- **modular**: every step can be ommited or replaced with another tools i.e. if you wish to skip reduction, execute the program with `--noreduction` parameter,     
- **flexible** toward many sequencing technology i.e. Illumina, 454 or Sanger. 

For more information have a look at the [poster](https://github.com/lpryszcz/redundans/blob/master/docs/poster.pdf).

# Prerequisites

- Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
- [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
- [BWA](http://bio-bwa.sourceforge.net/)
- [SSPACE3](http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
- [Gap2Seq](http://www.cs.helsinki.fi/u/lmsalmel/Gap2Seq/)

# Running
Redundans input consists of: 
- FastA-formatted **assembled contigs**
- FastQ-formatted **paired-end and/or mate pairs reads** - gzipped reads are supported ie .fastq.gz or .fq.gz 
Redundans will return FastA-formatted **homozygous genome assembly** as `scaffolds.filled.fa`. In addition, the program by default reports statistics of every step and iteration of its process.  

## Parameters
Mandatory parameters are marked **in bold**: 
- Genral options:
  -h, --help            show this help message and exit
  -v                    verbose
  --version             show program's version number and exit
  **-i FASTQ [FASTQ ...], --fastq FASTQ [FASTQ ...]**
                        FASTQ PE/MP files
  **-f FASTA, --fasta FASTA**
                        assembly FASTA file
  -o OUTDIR, --outdir OUTDIR
                        output directory
  -q MAPQ, --mapq MAPQ  min mapping quality for variants [10]
  -t THREADS, --threads THREADS
                        max threads to run [2]
  --log LOG             output log to [stderr]

- Reduction options:
  --identity IDENTITY   min. identity [0.8]
  --overlap OVERLAP     min. overlap  [0.75]
  --minLength MINLENGTH
                        min. contig length [200]

- Scaffolding options:
  -j JOINS, --joins JOINS
                        min k pairs to join contigs [5]
  -l LIMIT, --limit LIMIT
                        align at most l reads [5000000]
  -iters ITERS          scaffolding iterations per library  [2]
  --sspacebin SSPACEBIN
                        SSPACE path  [~/src/SSPACE/SSPACE_Standard_v3.0.pl]

You can skip some pipeline steps (all performed by default):
  --noreduction
  --noscaffolding
  --nogapclosing

## Test set run
In the folder [./test](https://github.com/lpryszcz/redundans/tree/master/test) you can find test dataset with 100kb genomic region and three simulated Illumina libraries: 
- paired-end with 300bp insert (300_?.fastq.gz), 
- paired-end with 600bp insert (600_?.fastq.gz),  
- mate pairs with 5kb insert (5000_?.fastq.gz). 
To run the test example, just execute: 
```bash
cd test
mkdir redundans 
cd redundans
../../redundans.py -v -i ../**.fastq.gz -f ../contigs.fasta
```

Note, the **order of libraries is not important**, as long as `_read1` and `_read2` from each library are given one after another i.e. `-i 600_1.fastq.gz 600_2.fastq.gz 300_1.fastq.gz 300_2.fastq.gz` would be interpreted the same as `-i 300_1.fastq.gz 300_2.fastq.gz 600_1.fastq.gz 600_2.fastq.gz`. 

# FAQ

- SSPACE fails with an error `Can't locate getopts.pl in @INC`.  
This is due to missing getops in recent Perl. Just do:
```bash
sudo cpan
install Perl4::CoreLibs
```

# Citation
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes 
