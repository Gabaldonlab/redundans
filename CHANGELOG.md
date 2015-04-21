####0.11c
	-depth-of-coverage is used in reduction step
	-calculate insert size from global distribution maximun, not mean

####0.11b
    - two similarity search algorithms: BLAT for `--identity 0.85`+ and LAST for `--identity` < 0.85

####0.11a
	- iterative insert size estimation refining for mate pairs
	- fasta2diverged.py deprecated
    - `-o/--output` now contains also sspace intermediate files
    - reduction (fasta2homozygous.py) uses LAST instead of BLAT & multi-threaded

####0.10b
	- `-l/--limit` now takes number of reads limit as fraction of homozygous genome size. To process all reads set -l to 0.
    - Gap2Seq is default gap closing software
    - BLAT is default reduction software
    - BWA MEM is default scaffolding software
