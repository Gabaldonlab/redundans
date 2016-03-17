### Table of Contents
- **[Redundans docs](#redundans-docs)**  
  - **[FAQ](#faq)**
  - **[FAQ - INSTALL.sh](#faq---installsh)**
  
- **[pyScaf docs](#pyscaf-docs)**  
  - **[NGS-based scaffolding](#ngs-based-scaffolding)**  
  - **[Reference-based scaffolding](#reference-based-scaffolding)**  

# Redundans docs

## FAQ
### SSPACE fails with an error `Can't locate getopts.pl in @INC`.  
This is due to missing getops in recent Perl. Just do:
```bash
sudo cpan
install Perl4::CoreLibs
```

### Reduction step takes a lot of time.   
Reduction step execute all-vs-all similarity search on your contigs. This may take some time if your assembly is heavily fragmented (>100k contigs).
In order to speed-up the analysis, you can use threads (i.e. `-t 8` for 8 cores). Make sure, you are running Python 2.7, as threading in reduction step is disabled in Python 2.6 and earlier.

### Estimation of my library statistics is incorrect. Can I specify these values manually?   
This can happen for highly fragmented assemblies or poor quality libraries. You can specify library statistics manually. To do so, look for *.is.txt file specific for your library i.e. for `-i 5000_1.fq.gz 5000_2.fq.gz` you will have to enter requested values into `5000_2.fq.gz.is.txt`. Make sure you specify some large number of mates for requested orientation (ie 100,000). For example if you want mate-pairs with RF orientation and 5kb insert size +/- 1.5kb, enter into respective *.is.txt file:
```bash
5000.0        5000.0  1500.0  0       0      100000   0
```

### No alignments error.  
If you see warning messages like the ones below while running the test set: 
```bash
[WARNING] No alignments for test/5000_1.fq.gz - test/5000_2.fq.gz!
[WARNING] No alignments for test/600_1.fq.gz - test/600_2.fq.gz!
```

The problem is because older version of BWA lack MEM algorithm (check it by executing `bwa mem`). If it gives you an error, download the latest [BWA](http://bio-bwa.sourceforge.net/).  

### Redundans fails with `OSError` or `maf-convert: not found`.  
Make sure you are using the latest version of [LAST](http://last.cbrc.jp/) aligner and that all dependencies are accessible through your PATH environmental variable. 

```bash
OSError: [Errno 13] Permission denied

# or
OSError: [Errno 2] No such file or directory

# or
[ERROR] maf-convert: not found
```

### Redundans fails with `IOError`.  
Make sure the `--sspacebin` point to the path containing SSPACE3 binaries ie. /home/user/bin/SSPACE-STANDARD-3.0_linux-x86_64/SSPACE_Standard_v3.0.pl

```bash
IOError: [Errno 2] No such file or directory: '.../test/run11/_sspace.1.1.fa'
```

### Why does Redundans use two similarity search algorithms, [BLAT](https://genome.ucsc.edu/FAQ/FAQblat.html#blat3) & [LAST](http://last.cbrc.jp/)?   
Prior to version v0.11d, Redundans was using BLAT or LAST for reduction step, in newer version only multi-threaded LAST is used. 
~~BLAT is lightweight & very fast, but lack sensitivity for more diverged sequences. If you specify `--identity` below 0.85, the pipeline will use LAST, that is ~4x slower, but more sensitive than BLAT.
Our simulations shows LAST is capable of correctly reducing heterozygous assemblies with up to 45% divergence between haplotypes.   
To limit speed difference between these two algorithms, LAST **runs in multiple threads**, so using `-t 4` you shouldn't see any difference in runtime between runs for `--identity 0.9` or `--identity 0.5`.~~

### How is multiple redundancy handled? 
Redundans removes all contigs, but the longest one, that fullfill identity & overlap critaria during reduction step. For more info see [issue #8](https://github.com/lpryszcz/redundans/issues/8).

### Why there are two github repositories for Redundans?
https://github.com/Gabaldonlab/redundans is the official repository for Redundans, but we keep https://github.com/lpryszcz/redundans for back-compatibility, as some of the very first users of Redundans use it.
**Both Redundans repositories contain the same code and are regularly updated.**

### How to cite Redundans?
Leszek P. Pryszcz and Toni Gabaldón (Submitted) Redundans: an assembly pipeline for highly heterozygous genomes. NAR

## FAQ - INSTALL.sh
### Installation succeeded, but redundans fails with `ImportError: No module named Bio`
Make sure you opened new terminal window after installation finished. 

### Installation succeeded, but redundans fails with `Bio.MissingPythonDependencyError: Requires sqlite3, which is included Python 2.5+`
Most likely you didn't install libsqlite3-dev before running installer. Try this:
```bash
# install missing library
sudo apt-get install sqlite3 libsqlite3-dev
# uninstall all
rm -rI ~/.pythonbrew ~/src/{*SSPACE,bwa,GapCloser,last,redundans}*
cp ~/.bashrc_bak ~/.bashrc
# open new terminal and relaunch installer
bash <(curl -Ls http://bit.ly/redundans_installer)
```

# pyScaf docs
pyScaf orders contigs from genome assemblies utilising several types of information:
- paired-end (PE) and/or mate-pair libraries ([NGS-based mode](#ngs-based-scaffolding))
- synteny to the genome of some related species ([reference-based mode](#reference-based-scaffolding))

## NGS-based scaffolding
This is under development... Stay tuned. 

## Reference-based scaffolding
In reference-based mode, pyScaf uses synteny to the genome of closely related species in order to order contigs and estimate distances between adjacent contigs.

Contigs are aligned globally (end-to-end) onto reference chromosomes, ignoring:
- matches not satisfying cut-offs (`--identity` and `--overlap`)
- suboptimal matches (only best match of each query to reference is kept) 
- and removing overlapping matches on reference. 

In preliminary tests, pyScaf performed superbly on simulated heterozygous genomes based on *C. parapsilosis* (13 Mb; CANPA) and *A. thaliana* (119 Mb; ARATH) chromosomes, reconstructing correctly all chromosomes always for CANPA and nearly always for ARATH ([Figures in dropbox](https://www.dropbox.com/sh/bb7lwggo40xrwtc/AAAZ7pByVQQQ-WhUXZVeJaZVa/pyScaf?dl=0), [CANPA table](https://docs.google.com/spreadsheets/d/1InBExy-qKDLj-upd8tlPItVSKc4mLepZjZxB31ii9OY/edit#gid=2036953672), [ARATH table](https://docs.google.com/spreadsheets/d/1InBExy-qKDLj-upd8tlPItVSKc4mLepZjZxB31ii9OY/edit#gid=1920757821)).  
Runs took ~0.5 min for CANPA on `4 CPUs` and ~2 min for ARATH on `16 CPUs`. 

**Important remarks:**
- Reduce your assembly before (fasta2homozygous.py) as any redundancy will likely break the synteny.
- pyScaf works better with contigs than scaffolds, as scaffolds are often affected by mis-assemblies (no *de novo assembler* / scaffolder is perfect...), which breaks synteny. 
- pyScaf works very well if divergence between reference genome and assembled contigs is below 20% at nucleotide level. 
- pyScaf does not allow any rearrangements (ie. translocations). I don't know also how the tool will behave if there are big gaps in contigs... 
- Consider closing gaps after scaffolding. 


```bash
./pyScaf.py -f test/run1/contigs.reduced.fa -r test/ref.fa
```

