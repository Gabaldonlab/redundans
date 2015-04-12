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

Prerequisites
================

* Python 2.7+ & Biopython 1.6+ `sudo easy_install -U biopython`
* [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3)
* [BWA](http://bio-bwa.sourceforge.net/)
* [SSPACE3](http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
* [Gap2Seq](http://www.cs.helsinki.fi/u/lmsalmel/Gap2Seq/)

Test
================
In the folder ./test you can find test dataset with 100kb genomic region and three libraries: pe300bp, pe600bp and mp5000bp. 
To run the test example, just execute: 

```bash
cd test
mkdir redundans 
cd redundans
../../redundans.py -v -i ../*.fastq.gz -f ../contigs.fasta
```

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


