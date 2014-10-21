================
 Redundans
================


Summary
================

Many genomes display high levels of heterozygosity (i.e. presence of different alleles at the same loci in homologous chromosomes), being those of hybrid organisms an extreme such case. The assembly of highly heterozygous genomes from short sequencing reads is a challenging task because it is difficult to accurately recover the different haplotypes. When confronted with highly heterozygous genomes, the standard assembly process tends to collapse homozygous regions and reports heterozygous regions in alternative contigs. The boundaries between homozygous and heterozygous regions result in multiple paths that are hard to resolve, which leads to highly fragmented assemblies with a total size larger than expected. This, in turn, causes numerous problems in downstream analyses i.e. fragmented gene models, wrong gene copy number, broken synteny. 

Redundans pipeline assist assembly of heterozygous genomes. 
It consists of three steps/modules:
1. redundancy reduction: detection and selectively removal of redundant contigs from an initial *de novo* assembly
2. scaffolding 
3. gap closing

While just genome assembly is enough to proceed with reduction step, paired-end and/or mate-pair libraries are required for scaffolding and gap closing. 

Redundans is: 
* **modular**, thus every step can be ommited or replaced with another tools. 
* **flexible** toward any sequencing technology (i.e. Illumina, 454 or Sanger). 

Prerequisites
================
* Python 2.7+
* Biopython 1.6+
* BWA
* pysam
* SSPACE2
* some de novo genome assembler ie SOAPdenovo2 or SPAdes


Citation
================
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes 


