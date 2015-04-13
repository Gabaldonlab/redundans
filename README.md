================
 Redundans
================


Summary
================

Redundans pipeline assist assembly of heterozygous genomes. 
It consists of three steps/modules: 

1. redundancy reduction: detection and selectively removal of redundant contigs from an initial *de novo* assembly
2. scaffolding 
3. gap closing

While just genome assembly is enough to proceed with reduction step, paired-end and/or mate-pair libraries are required for scaffolding and gap closing. 

Redundans is: 

* modular - every step can be ommited or replaced with another tools. 
* flexible toward any sequencing technology (i.e. Illumina, 454 or Sanger). 

Every step of Redundans pipeline may be skipped i.e. if you wish to skip reduction, execute the program with `--noreduction` parameter.  

For more information have a look at the [poster](https://github.com/lpryszcz/redundans/blob/master/docs/poster.pdf). 

Prerequisites
================

* Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
* [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
* [BWA](http://bio-bwa.sourceforge.net/)
* [SSPACE3](http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
* [Gap2Seq](http://www.cs.helsinki.fi/u/lmsalmel/Gap2Seq/)

Test
================
In the folder [./test](https://github.com/lpryszcz/redundans/tree/master/test) you can find test dataset with 100kb genomic region and three simulated Illumina libraries: 
* paired-end with 300bp insert (300_?.fastq.gz), 
* paired-end with 600bp insert (600_?.fastq.gz),  
* mate pairs with 5kb insert (5000_?.fastq.gz). 
To run the test example, just execute: 

```bash
cd test
mkdir redundans 
cd redundans
../../redundans.py -v -i ../*.fastq.gz -f ../contigs.fasta
```

Note, the *order of libraries is not important*, as long as `_read1` and `_read2` from each library are given one after another i.e. `-i 600_1.fastq.gz 600_2.fastq.gz 300_1.fastq.gz 300_2.fastq.gz` would be interpreted the same as `-i 300_1.fastq.gz 300_2.fastq.gz 600_1.fastq.gz 600_2.fastq.gz`. 

FAQ
================

* SSPACE fails with an error `Can't locate getopts.pl in @INC`.  
This is due to missing getops in recent Perl. Just do:

```bash
sudo cpan
install Perl4::CoreLibs
```

Citation
================
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes 


