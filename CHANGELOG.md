####0.11a
	-depth-of-coverage is used in reduction step
	- iterative insert size estimation refining for mate pairs
	-calculate insert size from global distribution maximun, not mean
	- fasta2diverged.py deprecated
    - `-o/--output` now contains also sspace intermediate files
    - fasta2homozygous.py uses LAST instead of BLAT & multi-threaded

####0.10b
	- `-l/--limit` now takes number of reads limit as fraction of homozygous genome size. To process all reads set -l to 0. 
