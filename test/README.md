# Test dataset
This directory contains data from simulated hybrid genome with 5% divergence between parentals and 40% of loss of heterozygosity (LOH).

To run the test example, execute: 
```bash
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1 
```

**Table of contents**  
- [Generation of test set](#Generation-of-test-set)
  - [Simulations](#simulations)
    - [Heterozygous genome](#heterozygous-genome)
    - [Short reads](#short-reads)
  - [*De novo* genome assembly](#de-novo-genome-assembly)
  - [Redundans pipeline](#redundans-pipeline)
    - [Run statistics](#run-statistics)
      - [Parameters estimation](#parameters-estimation)
      - [Reduction](#reduction)
      - [Scaffolding](#scaffolding)
      - [Gap closing](#gap-closing)
      - [Summary statistics](#summary-statistics)
- [Accuracy estimation](#accuracy-estimation)

## Generation of test set
Below, you can find detailed description for regeneration of test dataset with any input genome. 

### Simulations
Test dataset is based on 100 Kb region from [*C. parapsilosis* CDC317](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2834264/?report=classic) [HE605202:100000-200000]. 
```bash
samtools faidx CANPA.fa HE605202:100000-200000 > ref.fa
```

#### Heterozygous genome
Heterozygous regions, 5% divergence and 40% LOH, have been simulated with [fasta2diverged.py](https://github.com/lpryszcz/bin). 
```bash
fasta2diverged.py -v -i ref.fa --loh 0.40 --dna -d 0.05 > ref.d05.l40.fa
```

#### Short reads
Two shotgun sequencing libraries including typical Illumina-related errors were simulated using [GemSIM](http://sourceforge.net/projects/gemsim/): 
- 100 bp paired-end reads with 600 bp insert [200X coverage]:
 - 600_1.fq.gz
 - 600_2.fq.gz
- 50 bp mate-pair reads with 5 kb insert [20X coverage]:
 - 5000_1.fq.gz
 - 5000_2.fq.gz

```bash
# generate SNPs
~/src/GemSIM_v1.6/GemHaps.py -r ref.fa -g '.50,0 .50,6'

# 2x 50,000 reads
## 100 bp paired-end with 600 bp insert +/- 50 bp
~/src/GemSIM_v1.6/GemReads.py -r ref.fa -n 50000 -g ref.txt -l 100 -m /home/lpryszcz/src/GemSIM_v1.6/models/ill100v5_p.gzip -q 33 -o 600 -u 600 -s 50 -p
~/src/GemSIM_v1.6/GemReads.py -r ref.d05.l40.fa -n 50000 -g ref.txt -l 100 -m /home/lpryszcz/src/GemSIM_v1.6/models/ill100v5_p.gzip -q 33 -o 600.d05.l40 -u 600 -s 50 -p

# 2x 10,000 reads
## 50 bp mate-pairs with 5 kb insert +/- 500 bp
~/src/GemSIM_v1.6/GemReads.py -r ref.fa -n 10000 -g ref.txt -l 50 -m /home/lpryszcz/src/GemSIM_v1.6/models/ill100v5_p.gzip -q 33 -o 5000 -u 5000 -s 500 -p
~/src/GemSIM_v1.6/GemReads.py -r ref.d05.l40.fa -n 10000 -g ref.txt -l 50 -m /home/lpryszcz/src/GemSIM_v1.6/models/ill100v5_p.gzip -q 33 -o 5000.d05.l40 -u 5000 -s 500 -p

# combine both libraries in random order
## 600_
paste <(zcat 600_1.fastq.gz 600.d05.l40_1.fastq.gz) <(zcat 600_2.fastq.gz 600.d05.l40_2.fastq.gz) | paste - - - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 | "gzip > random.600.d05.l40_1.fq.gz"; print $2,$4,$6,$8 | "gzip > random.600.d05.l40_2.fq.gz"}'
## 5000_
paste <(zcat 5000_1.fastq.gz 600.d05.l40_1.fastq.gz) <(zcat 5000_2.fastq.gz 5000.d05.l40_2.fastq.gz) | paste - - - - | shuf | awk -F'\t' '{OFS="\n"; print $1,$3,$5,$7 | "gzip > random.5000.d05.l40_1.fq.gz"; print $2,$4,$6,$8 | "gzip > random.5000.d05.l40_2.fq.gz"}'

# symlinks with handier names
for f in *.fq.gz; do ln -s $f `echo $f | sed 's/random.//;s/.d05.l40//'`; done

# mate-pairs contaminated with paired-end reads
#for i in 1 2; do zcat <(zcat test/5000_$i.fq.gz test/5000_$i.fq.gz | fastx2reverse_complement.py fastq | gzip) test/600_$i.fq.gz > test/mixed.$i.fq.gz; done
```

### *De novo* genome assembly
Simulated reads were assembled with SPAdes v.3.5.0. 
```bash
dipspades.py --only-assembler -t 4 -o 600 -1 600_1.fq.gz -2 600_2.fq.gz
```

### Redundans pipeline
Finally, heterozygous genome assembly pipeline can be applied to simulated heterozygous contigs and paired/mate-pairs reads.

```bash
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1 
```

#### Run statistics
##### Parameters estimation
At the beginning, Redundans estimates number of parameters:
- number of reads that it's going to align (based on `-l/--limit` parameter) 
```bash
 Aligning 32779 mates per library...
```
- library statistics:
 - insert size: median, mean, stdev
 - mates orientation: FF, FR, RF, RR

###### Table 1: Library statistics

FastQ files | read length | median | mean | stdev | FF | FR | RF | RR
----- | -----: | -----: | -----: | -----: | -----: | -----: | -----: | -----: 
test/5000_1.fq.gz test/5000_2.fq.gz | 50 | 4998 | 4990.20 | 721.47 | 0 | 4,674 | 0 | 0
test/600_1.fq.gz test/600_2.fq.gz | 100 | 599 | 598.63 | 47.68 | 0 | 10,000 | 0 | 0

Above, you can see, there are two libraries:
- test/5000_1.fq.gz & test/5000_2.fq.gz:
  - FR orientation
  - mate-pairs with 5000 bp insert size +/- 600 bp
- test/600_1.fq.gz & test/600_2.fq.gz:
  - FR orientation
  - paired-end with 600 bp insert size +/- 39 bp

Based on these statistics, the program will perform scaffolding and gap closing starting from libraries with the smallest insert size.  
Also at this stage, the user is notified about libraries with poor statistics i.e. large stdev or not consisted orientation. For these cases, the statistics are updated before every scaffolding iteration. This is done, because for libraries with large insert size, the estimation of insert size may be highly affected by the fragmentation of the assembly. 

##### Reduction
In the first step, Redundans will recognise and selectively remove heterozygous contigs from initial assembly.
The details about all removed heterozygous contigs are stored in tab-delimited file `contigs.reduced.fa.hetero.tsv`.
Additionally, the histogram with identities of reduced contigs are saved as [`contigs.reduced.fa.hist.png`](/docs/contigs.reduced.fa.hist.png). 

At the end of this step, the user will be notified about:
- initial assembly size and fragmentation
- what portion of the the initial assembly have been recognised as heterozygous
- identity between heterozygous regions
- size and fragmentation of reduced assembly

###### Table 2: Reduction statistics

file name | genome size | contigs | heterozygous size | [%] | heterozygous contigs | [%] | identity [%] | possible joins | homozygous size | [%] | homozygous contigs | [%]
----- | -----: | -----: | -----: | ----- | -----: | ----- | -----: | ----- | -----: | ----- | -----: | ----- 
test/run6/contigs.fa | 163897 | 245 | 76363 | 46.59 | 222 | 90.61 | 94.837 | 0 | 97319 | 59.38 | 23 | 9.39

Above, you can see that:
- initial SPAdes assembly was fragmented (245 contigs) and larger (163.8 Kb) than reference sequence (100 Kb)
- Redundans:
 - removed 220 heterozygous contigs, resulting in reduced assembly of 97.8 Kb (very similar to reference) in 28 contigs
 - estimated 94.837% identity between heterozygous regions, which is very similar to the divergence that was used for simulations (5%)

##### Scaffolding
Redundans repeat scaffolding on each library several times (`--iters`). For each scaffolding iteration, the program reports number of:
- processed pairs
- pairs passing alignement quality
- pairs that aligned in different contigs

```bash
 iteration 1.1 ...
   19572 pairs. 18505 passed filtering [94.55%]. 1993 in different contigs [10.18%].
 iteration 1.2 ...
   19572 pairs. 18569 passed filtering [94.88%]. 340 in different contigs [1.74%].
 iteration 2.1 ...
   19572 pairs. 17009 passed filtering [86.90%]. 755 in different contigs [3.86%].
  closing gaps ...
 iteration 2.2 ...
   19572 pairs. 17261 passed filtering [88.19%]. 0 in different contigs [0.00%].
  closing gaps ...
```

Above, you can see that in iteration 1.1 above 10% of reads aligned in different contigs, while in subsequent iteration, less than 2% reads. This informs that in general it's enough to use 1-2 iterations of scaffolding per library. 

##### Gap closing
For the time-being, only the iteration number is reported.

```bash
 iteration 1.1 ...
 iteration 1.2 ...
 iteration 2.1 ...
 iteration 2.2 ...
```

##### Summary statistics
At the end, the program reports details of each step:
- file name
- fragmentation and size of the assembly at this stage
- %GC content
- fragmentation and size of assembly in contigs >1 Kb
- N50 & N90
- cumulative size of gaps
- the size of the longest contig

fname | contigs | bases | GC [%] | contigs >1kb | bases in contigs >1kb | N50 | N90 | Ns | longest
:----- | -----: | -----: | :-----: | -----: | -----: | -----: | -----: | -----: | -----: 
test/run6/contigs.fa | 245 | 163897 | 40.298 | 24 | 117391 | 3975 | 233 | 0 | 29603
test/run6/contigs.reduced.fa | 23 | 97319 | 39.333 | 17 | 94157 | 7321 | 2195 | 0 | 29603
test/run6/_sspace.1.1.fa | 4 | 97242 | 39.315 | 3 | 96924 | 87727 | 87727 | 473 | 87727
test/run6/_sspace.1.2.fa | 3 | 97321 | 39.315 | 2 | 97003 | 87727 | 87727 | 552 | 87727
test/run6/_sspace.2.1.fa | 1 | 100280 | 39.315 | 1 | 100280 | 100280 | 100280 | 3511 | 100280
test/run6/_sspace.2.1.filled.fa | 1 | 100251 | 39.363 | 1 | 100251 | 100251 | 100251 | 2782 | 100251
test/run6/_sspace.2.2.fa | 1 | 100251 | 39.363 | 1 | 100251 | 100251 | 100251 | 2782 | 100251
test/run6/_sspace.2.2.filled.fa | 1 | 100123 | 39.360 | 1 | 100123 | 100123 | 100123 | 2677 | 100123
test/run6/scaffolds.fa | 1 | 100123 | 39.360 | 1 | 100123 | 100123 | 100123 | 2677 | 100123
test/run6/_gapcloser.1.1.fa | 1 | 100495 | 39.542 | 1 | 100495 | 100495 | 100495 | 1265 | 100495
test/run6/_gapcloser.1.2.fa | 1 | 100498 | 39.579 | 1 | 100498 | 100498 | 100498 | 975 | 100498
test/run6/_gapcloser.2.1.fa | 1 | 100502 | 39.679 | 1 | 100502 | 100502 | 100502 | 4 | 100502
test/run6/_gapcloser.2.2.fa | 1 | 100506 | 39.679 | 1 | 100506 | 100506 | 100506 | 4 | 100506
test/run6/scaffolds.filled.fa | 1 | 100506 | 39.679 | 1 | 100506 | 100506 | 100506 | 4 | 100506


## Accuracy estimation
Accuracy of recovered contigs can be assessed by alignment of final scaffolds (`run1/scaffolds.filled.fa`) back onto reference (`ref.fa`).

### NUCMER
This can be quickly achieved with nucmer (works for small genomes):

```bash
cd run1
# align
nucmer -p nucmer.ref ../ref.fa scaffolds.filled.fa
# generate pairwise plot
mummerplot --png --large nucmer.ref.delta -p nucmer.ref.plot
# open plot
eog nucmer.ref.plot.png
```
[Figure 1: Pairwise alignment of reference sequence and redundans scaffolds](/docs/nucmer.ref.plot.png)

You may want to compare final redundans scaffolds (`run1/scaffolds.filled.fa`) with initial SPAdes contigs (`contigs.fa`). 

```bash
cd run1
# align
nucmer -p nucmer.contigs scaffolds.filled.fa ../contigs.fa
# generate pairwise plot
mummerplot --png --large nucmer.contigs.delta -p nucmer.contigs.plot
# open plot
eog nucmer.contigs.plot.png
```
[Figure 2: Pairwise alignment of redundans scaffolds and SPAdes contigs](/docs/nucmer.contigs.plot.png)

Finally, comparison of consensus contigs from dipSPAdes (`consensus_contigs.fa`) and reference sequnece (`ref.fa`) informs about missing regions in dipSPAdes reconstruction.

```bash
# align
nucmer -p nucmer.ref_dipspades ref.fa consensus_contigs.fa
# generate pairwise plot
mummerplot --png --large nucmer.ref_dipspades.delta -p nucmer.ref_dipspades.plot
# open plot
eog nucmer.ref_dipspades.plot.png
```
[Figure 3: Pairwise alignment of reference sequence and dipSPAdes consensus contigs](/docs/nucmer.ref_dipspades.plot.png)

### LAST
The same can be done using LAST (this will work also for large genomes):

```bash
# index reference genome (need to be done only once per each reference)
lastdb ../ref.fa ../ref.fa

# align & generate dotplot on the fly
lastal -f TAB ../ref.fa scaffolds.filled.fa | last-dotplot - scaffolds.filled.fa.png

# open plot
eog scaffolds.filled.fa.png
```
[Figure 4: Pairwise alignment of reference sequence and Redundans contigs](/docs/scaffolds.filled.fa.png)


All programs/scripts not disclosed in this repository (i.e. fasta2diverged.py), can be found in https://github.com/lpryszcz/bin. 
