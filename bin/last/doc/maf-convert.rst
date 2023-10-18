maf-convert
===========

This script reads alignments in maf_ format, and writes them in
another format.  It can write them in these formats: axt_, blast,
blasttab, chain_, gff, html, psl_, sam, tab.  You can use it like this::

  maf-convert psl my-alignments.maf > my-alignments.psl

It's often convenient to pipe in the input, like this::

  ... | maf-convert psl > my-alignments.psl

This script takes the first (topmost) maf sequence as the "reference"
/ "subject" / "target", and the second sequence as the "query".
(Exception: when converting DNA-to-protein alignments to gff or psl,
the protein becomes the "query" and the DNA becomes the "reference".)

For html: if the input includes probability lines starting with 'p',
then the output will be colored by column probability.  (To get lines
starting with 'p', run lastal with option -j set to 4 or higher.)

.. _maf: http://genome.ucsc.edu/FAQ/FAQformat.html#format5
.. _axt: https://genome.ucsc.edu/goldenPath/help/axt.html
.. _chain: https://genome.ucsc.edu/goldenPath/help/chain.html
.. _psl: https://genome.ucsc.edu/FAQ/FAQformat.html#format2

Options
-------

-h, --help
       Print a help message and exit.

-p, --protein
       Specify that the alignments are of proteins, rather than
       nucleotides.  This affects psl format only (the first 4
       columns), and has no effect for DNA-to-protein alignments.

-j N, --join=N
       Join alignments that are co-linear (align different parts of
       the same sequences and strands, with the parts being in the
       same order in each sequence), are separated by at most N
       letters in each sequence, and are consecutive in the input.
       This affects psl and gff formats only.

-J N, --Join=N
       Join alignments that are co-linear, are separated by at most
       N letters in each sequence, and are nearest in each sequence.
       This affects psl and gff formats only, and reads the whole
       input into memory.

-n, --noheader
       Omit any header lines from the output.  This may be useful if
       you concatenate outputs, e.g. from parallel jobs.

-d, --dictionary
       Include a dictionary of sequence lengths in the sam header
       section (lines starting with @SQ).  This requires reading the
       input twice, so it must be a real file (not a pipe).  This
       affects sam format only.

-f DICTFILE, --dictfile=DICTFILE
       Get a sequence dictionary from DICTFILE.  This affects sam
       format only.  You can create a dict file using
       CreateSequenceDictionary (http://picard.sourceforge.net/).

-r READGROUP, --readgroup=READGROUP
       Specify read group information.  This affects sam format
       only.  Example: -r 'ID:1 PL:ILLUMINA SM:mysample'

-l CHARS, --linesize=CHARS
       Write CHARS characters per line.  This affects blast and html
       formats only.

DNA-versus-protein blast output
-------------------------------

If the input has protein-versus-DNA alignments like this::

  GluTrpThrAlaLeuIleAsnLeuLysAsnArg--AspLeuValIleLysAlaAlaAsp
  GAATAGTCCGGTTGAAAAAATGTACAAAAACATAAGAGAACATTACAAAACTTGCAGTC

then the blast-format output will look like this::

  Glu***SerGly***LysAsnValGlnLysHis  GluAsnIleThrLysLeuAlaVal
  GAATAGTCCGGTTGAAAAAATGTACAAAAACATAAGAGAACATTACAAAACTTGCAGTC
  |||::::::         |||...:::...:::  :::   :::...|||   |||
  GluTrpThrAlaLeuIleAsnLeuLysAsnArg--AspLeuValIleLysAlaAlaAsp

The DNA's translation (assuming the standard genetic code) is shown
above it.  ``|||`` indicates a match, ``:::`` a positive alignment
score, and ``...`` an alignment score of 0.  ``:::`` and ``...`` are
shown only for alignments preceded by a substitution score matrix of
the sort in lastal's output header.

Hints for sam/bam
-----------------

* To run fast on multiple CPUs, and get a correct header at the top,
  this may be the least-awkward way.  First, make a header (perhaps by
  using CreateSequenceDictionary).  Then, concatenate the output of a
  command like this::

    parallel-fastq "... | maf-convert -n sam" < q.fastq

* Here is yet another way to get a sequence dictionary, using samtools
  (http://samtools.sourceforge.net/).  Assume the reference sequences
  are in ref.fa.  These commands convert x.sam to y.bam while adding a
  sequence dictionary::

    samtools faidx ref.fa
    samtools view -bt ref.fa.fai x.sam > y.bam

* If a query name ends in "/1" or "/2", maf-convert interprets it as a
  paired sequence.  (This affects sam format only.)  However, it does
  not calculate all of the sam pairing information (because it's hard
  and better done by specialized sam manipulators).

  Fix the pair information in y.sam, putting the output in z.bam.
  Using picard::

    java -jar FixMateInformation.jar I=y.sam O=z.bam VALIDATION_STRINGENCY=SILENT

  Alternatively, using samtools::

    samtools sort -n y.bam ysorted
    samtools fixmate ysorted.bam z.bam

"Bugs"
------

* For sam: the QUAL field (column 11) is simply copied from the maf q
  line.  The QUAL is supposed to be encoded as ASCII(phred+33),
  whereas maf q lines are encoded differently according to the UCSC
  Genome FAQ.  However, if you run lastal with option -Q1, the maf q
  lines will in fact be ASCII(phred+33).

* The blast format is merely blast-like: it is not identical to NCBI
  BLAST.
