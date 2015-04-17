# Test dataset
This directory contains data from simulated hybrid genome with 5% divergence between parentals and 40% of loss of heterozygosity (LOH).

To run the test example, execute: 
```bash
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1 
```
## Generation of test set
Below, you can find detailed description for regeneration of test dataset with any input genome.  

### Simulations
Entire test dataset is based on 100 Kb region from *C. parapsilosis* CDC317 (REF!). 
```bash
samtools faidx CANPA.fa HE605202:100000-200000 > ref.fa
```

#### Heterozygous genome
Heterozygous regions, 5% divergence and 40% LOH, have been simulated with [fasta2diverged.py](https://github.com/lpryszcz/bin). 
```bash
fasta2diverged.py -v -i ref.fa --loh 0.40 --dna -d 0.05 > ref.d05.l40.fa
```

#### Short reads
Two Illumina libraries were simulated using [GemSIM](http://sourceforge.net/projects/gemsim/): 
- paired-end with 600bp insert:
 - 600_1.fq.gz
 - 600_2.fq.gz
- mate pairs with 5kb insert:
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
```

### *De novo* genome assembly
Simulated reads were assembled with SPAdes v.3.5.0. In order 
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
- number of reads that it's going to align (based on `-l/--limit` parameter
- library statistics:
 - insert size: median, mean, stdev
 - mates orientation: FF, FR, RF, RR

```bash
 Aligning 32779 mates per library...
```
###### Table 1: Library statistics

| Insert size statistics |  |  |  | Mates orientation stats |
FastQ files | median | mean | stdev | FF | FR | RF | RR
----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- 
test/5000_1.fq.gz test/5000_2.fq.gz | 5028 | 5031.51 | 600.93 | 5 | 17494 | 159 | 1
test/600_1.fq.gz test/600_2.fq.gz | 598 | 598.28 | 39.26 | 0 | 32778 | 1 | 0

Above, you can see, there are two libraries:
- FR: paired-end with 600 bp insert size +/- 39 bp
- FR: mate-pairs with 5000 bp insert size +/- 600 bp

Based on these statistics, the program will perform scaffolding and gap closing starting from libraries with the smallest insert size.
Here, the user is notified about libraries with poor statistics i.e. large stdev or not consisted orientation. For these cases, the statistics are updated before every scaffolding iteration. 

##### Reduction
In the first step, Redundans will recognise and selectively remove heterozygous contigs from initial assembly.
At the end of this step, the user will be notified about:
- initial assembly size and fragmentation
- what portion of the the initial assembly have been recognised as heterozygous
- identity between heterozygous regions
- size and fragmentation of reduced assembly

###### Table 2: Reduction statistics

file name | genome size | contigs | heterozygous size | [%] | heterozygous contigs | [%] | identity [%] | possible joins | homozygous size | [%] | homozygous contigs | [%]
----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- | ----- 
run1/contigs.fa | 163897 | 245 | 65287 | 39.83 | 217 | 88.57 | 95.243 | 0 | 98610 | 60.17 | 28 | 11.43

Above, you can see that:
- initial SPAdes assembly was fragmented (245 contigs) and larger (163 Kb) than reference sequence (100 Kb)
- Redundans:
 - removed 217 heterozygous contigs, resulting in reduced assembly of 98 Kb (very similar to reference) in 28 contigs
 - estimated 95.2% identity between heterozygous regions, which is very similar to the divergence that was used for simulations (5%)

##### Scaffolding
Redundans repeat scaffolding on each library several times (`--iters`). For each scaffolding iteration, the program reports number of:
- processed pairs
- pairs passing alignement quality
- pairs that aligned in different contigs

```bash
 iteration 1.1 ...
   32780 pairs. 31126 passed filtering [94.95%]. 3484 in different contigs [10.63%].
 iteration 1.2 ...
   32780 pairs. 31249 passed filtering [95.33%]. 626 in different contigs [1.91%].
 iteration 2.1 ...
   20000 pairs. 18581 passed filtering [92.91%]. 1523 in different contigs [7.62%].
 iteration 2.2 ...
   20000 pairs. 18862 passed filtering [94.31%]. 0 in different contigs [0.00%].
```

Above, you can see that in iteration 1.1 above 10% of reads aligned in different contigs, while in subsequent iteration, less than 2% reads. This informs that in general it's enough to use 1-2 iterations of scaffolding per library. 

##### Gap closing
For the time-being, only the iteration number is reported.

```bash
 iteration 1 ...
 iteration 2 ...
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
run1/contigs.fa | 245 | 163897 | 40.298 | 24 | 117391 | 3975 | 233 | 0 | 29603
run1/contigs.reduced.fa | 28 | 98610 | 39.516 | 17 | 94157 | 7321 | 1858 | 0 | 29603
run1/_sspace.1.1.fa | 6 | 98479 | 39.507 | 4 | 97405 | 87549 | 4745 | 584 | 87549
run1/_sspace.1.2.fa | 4 | 98937 | 39.507 | 4 | 98937 | 88627 | 4745 | 1042 | 88627
run1/_sspace.2.1.fa | 1 | 100747 | 39.507 | 1 | 100747 | 100747 | 100747 | 2852 | 100747
run1/_sspace.2.2.fa | 1 | 100747 | 39.507 | 1 | 100747 | 100747 | 100747 | 2852 | 100747
run1/scaffolds.fa | 1 | 100747 | 39.507 | 1 | 100747 | 100747 | 100747 | 2852 | 100747
run1/_gap2seq.1.1.fa | 1 | 100818 | 39.762 | 1 | 100818 | 100818 | 100818 | 38 | 100818
run1/_gap2seq.2.1.fa | 1 | 100818 | 39.762 | 1 | 100818 | 100818 | 100818 | 38 | 100818
run1/scaffolds.filled.fa | 1 | 100818 | 39.762 | 1 | 100818 | 100818 | 100818 | 38 | 100818

## Accuracy estimation
Accuracy of recovered contigs can be assessed by alignment of final scaffolds (`run1/scaffolds.filled.fa`) back onto reference (`ref.fa`).
This can be quickly achieved with nucmer:

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


All programs/scripts not disclosed in this repository (i.e.  fasta2diverged.py), can be found in https://github.com/lpryszcz/bin. 
