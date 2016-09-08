
####0.12d
- added all dependencies to github
 - INSTALL.sh downloads & compiles everything
 - all necessary paths are defined by `redundans.py`
- README.md updated
- subprocess is closed when finished to lower memory
- contigs FastA is sorted by descending contig size & contigs below `--minLength` are removed (this speeds up reduction greatly)
- reduction step (`fasta2homozygous.py`) was polished & optimised
  - speed-optimised: avoided LASTal results sorting by sorting contigs FastA
  - memory-optimised ie. generator instead of list (thanks to sorted contigs FastA)
  - contigs FastA file has to be ordered by descending contig size
/- changed default parameters:
/  - `--identity 0.75`
/  - `--overlap 0.75`
- prints how many iterations in total ie. `iteration 1.1 of 2.2 ...`
- libraries are estimated after reduction (avoiding double estimations for crappy libs)
TBD
- check if something need recompilation and try to recompile
- instead of processing 100M reads for 4 libs with similar IS, split it by 25M each

####0.12c
- added `--resume` option
- added warning if unable to fetch dependencies versions
- reduction step `fasta2homozygous.py`
 - produce more accurate identity between recognised heterozygous contigs
 - polished code

####0.12-beta
- LASTal version checked on runtime
- added FastaIndex.py
  - generate stats into .fai file - `samtools faidx` compatible
- simplified dependencies
  - Biopython, scipy, numpy & SQLite no longer needed


####0.12-alpha
- improved reduction step performance (fasta2homozygous.py)
  - no sorting - greatly improves performance on large and fragmented genomes
    - removed `-S / --sortopt` parameter
  - no BLAT - unable to multi-thread BLAT, thus relying on LAST completely
    - global alignment `-T 1`, instead of local - speed-up and less permissive (more accurate) reduction
  - enabled [multi-threading in LAST from v693](http://last.cbrc.jp/last/index.cgi/rev/4174fdbdb9a1)
  - no temp files
    - all is processed through pipes
    - removed unnecessary maf-convert, gzip and file parsing
- code polished

####0.11-beta
- [UNIX installer](https://github.com/lpryszcz/redundans#unix-installer) and [docker image](https://github.com/lpryszcz/redundans#docker-image)
- added new parameters
  - `--log`
  - `-S / --sortopt`
  - `-a / --linkratio`
- code polished
- solved minor issues

####0.11b
- two similarity search algorithms: BLAT for `--identity 0.85`+ and LAST for `--identity` < 0.85
- corrected error in fastq2sspace.py, so now libraries are merged based on mean insert size, not median

####0.11a
- iterative insert size estimation refining for mate pairs
- fasta2diverged.py deprecated
- `-o/--output` now contains also sspace intermediate files
- reduction (fasta2homozygous.py) uses LAST instead of BLAT
  - LAST multi-threaded (only in Python 2.7+)
- cleaning-up output directory from intermediate files

####0.10b
- `-l/--limit` now takes number of reads limit as fraction of homozygous genome size. To process all reads set -l to 0.
- Gap2Seq is default gap closing software
- BLAT is default reduction software
- BWA MEM is default scaffolding software
