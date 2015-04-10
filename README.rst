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
* samtools (https://github.com/samtools/htslib)
* BWA (http://bio-bwa.sourceforge.net/)
* SSPACE2 (http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE)
* some *de novo* genome assembler ie. SOAPdenovo2 or SPAdes
* Python 2.7+
* Biopython 1.6+ ```sudo easy_install -U biopython```
* pysam ```sudo easy_install -U pysam```

Test
================

```
cd test
mkdir redundans
cd redundans
~/src/redundans/redundans.py -v -i ../*.fastq.gz -f ../contigs.fasta
```


Citation
================
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes 


