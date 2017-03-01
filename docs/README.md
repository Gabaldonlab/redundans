### Table of Contents
- **[Redundans docs](#redundans-docs)**  
  - **[FAQ](#faq)**
  - **[FAQ - INSTALL.sh](#faq---installsh)**
  
# Redundans docs

## FAQ
### Estimation of my library statistics is incorrect. Can I specify these values manually?   
This can happen for highly fragmented assemblies or poor quality libraries. You can specify library statistics manually.
To do so, look for *.is.txt file specific for your library i.e. for `-i 5000_1.fq.gz 5000_2.fq.gz` you will have to enter requested values into `5000_2.fq.gz.is.txt`.
.is.txt file is tab-delimited, with following columns:
- read length
- median insert size
- mean insert size
- standard deviation of insert size
- and number of reads with FF, FR, RF & RR orientation. 
Make sure you specify some large number of mates for requested orientation (ie 100,000).
For example if you want mate-pairs with RF orientation and 5kb insert size +/- 1.5kb, enter into respective *.is.txt file:
```bash
60  5000.0        5000.0  1500.0  0       0      100000   0
```

### How is multiple redundancy handled? 
Redundans removes all contigs, but the longest one, that fullfill identity & overlap critaria during reduction step. For more info see [issue #8](https://github.com/lpryszcz/redundans/issues/8).

### Why there are two github repositories for Redundans?
https://github.com/Gabaldonlab/redundans is the official repository for Redundans, but we keep https://github.com/lpryszcz/redundans for back-compatibility, as some of the very first users of Redundans use it.
**Both Redundans repositories contain the same code and are regularly updated.**

### How to cite Redundans?
Leszek P. Pryszcz and Toni Gabaldón (2016) Redundans: an assembly pipeline for highly heterozygous genomes. NAR. [doi: 10.1093/nar/gkw294](http://nar.oxfordjournals.org/content/44/12/e113)

## FAQ - INSTALL.sh


