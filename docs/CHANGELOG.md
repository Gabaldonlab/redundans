####0.12 (planned...)
- depth-of-coverage is used in reduction step
- calculate insert size from global distribution maximun, not mean
- write propietary FastA indexing module
  - get rid of SeqIO.index_db
  - incorporate fasta stats

####0.11d (work in progress...)
- improved reduction step performance (fasta2homozygous.py)
  - no sorting - greatly improves performance on large and fragmented genomes
    - removed `-S / --sortopt` parameter
  - no BLAT - unable to multi-thread BLAT, thus relying on LAST completely
  - enabled [multi-threading in LAST from v693](http://last.cbrc.jp/last/index.cgi/rev/4174fdbdb9a1)
  - no temp files
    - all is processed through pipes
    - removed unnecessary maf-convert, gzip and file parsing
- code polished

####0.11c
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
