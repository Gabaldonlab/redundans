last-split
==========

last-split finds "split alignments" (typically for DNA) or "spliced
alignments" (typically for RNA).

It reads candidate alignments of query sequences to a genome, and cuts
from them a unique best alignment for each part of each query.  This
is useful for DNA queries that cross rearrangement breakpoints, or RNA
queries that cross splice junctions.

For a detailed explanation of what last-split does, please see `A
survey of localized sequence rearrangements in human DNA
<https://doi.org/10.1093/nar/gkx1266>`_ (update in `A pipeline for
complete characterization of complex germline rearrangements from long
DNA reads <https://doi.org/10.1186/s13073-020-00762-1>`_).  The
algorithm, and application to whole genomes, is in `Split-alignment of
genomes finds orthologies more accurately
<https://doi.org/10.1186/s13059-015-0670-9>`_.

Examples
--------

Split alignment of DNA reads to a genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assume the DNA reads are in a file called "q.fastq" (in fastq-sanger
format), and the genome is in "genome.fasta" (in fasta format).  We
can do the alignment like this::

  lastdb -uNEAR db genome.fasta
  last-train -Q1 db q.fastq > train.out
  lastal -p train.out db q.fastq | last-split > out.maf

Spliced alignment of RNA reads to a genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now we assume that "q.fastq" has reads from RNA forward (sense)
strands.  This time, we provide the genome information to last-split,
which causes it to do spliced instead of split alignment, and also
tells it where the splice signals are (GT, AG, etc)::

  lastal -p train.out -D10 db q.fastq | last-split -g db > out.maf

This will favor splices starting at GT (and to a lesser extent GC and
AT), and ending at AG (and to a lesser extent AC).  The output shows
the donor (``don``) and acceptor (``acc``) dinucleotides.

It also favors splices with introns of typical length, specified by a
log-normal distribution (cis-splices).  However, it allows arbitrary
trans-splices between any two places in the genome.

``-D10`` sets a very loose significance threshold, so that we can find
very short parts of a spliced alignment (e.g. short exons).  Note that
last-split discards the lowest-significance alignments, but it uses
them to estimate the ambiguity of higher-significance alignments.

If your reads are from unknown/mixed RNA strands, add ``-d2`` to the
last-split options.

Alignment of two whole genomes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

See the cookbook_.

Output
------

The output is in MAF(-like) format::

  a score=150 mismap=0.000413
  s chr21  15963638 25 + 48129895 TCAGATGAGGACCTAATTTATTACT
  s query7       50 25 +       75 TCAGATGAGGACCTAATTTATTACT
  q query7                        EBEEC@CE=EEE?FEDAED5?@@D@
  p                               !#$'BBBBBBBBBBBBBBBBBBBBB

The "mismap" is the estimated probability that this part of the query
should be aligned to a different part of the genome.

The line starting with "p" (shown by option ``-fMAF+``) indicates the
probability that each base should be aligned to a different part of
the genome.  It uses a compact code:

======  =================   ======  =================
Symbol  Error probability   Symbol  Error probability
======  =================   ======  =================
``!``   0.79 -- 1           ``0``   0.025 -- 0.032
``"``   0.63 -- 0.79        ``1``   0.02  -- 0.025
``#``   0.5  -- 0.63        ``2``   0.016 -- 0.02
``$``   0.4  -- 0.5         ``3``   0.013 -- 0.016
``%``   0.32 -- 0.4         ``4``   0.01  -- 0.013
``&``   0.25 -- 0.32        ``5``   0.0079 -- 0.01
``'``   0.2  -- 0.25        ``6``   0.0063 -- 0.0079
``(``   0.16 -- 0.2         ``7``   0.005  -- 0.0063
``)``   0.13 -- 0.16        ``8``   0.004  -- 0.005
``*``   0.1  -- 0.13        ``9``   0.0032 -- 0.004
``+``   0.079 -- 0.1        ``:``   0.0025 -- 0.0032
``,``   0.063 -- 0.079      ``;``   0.002  -- 0.0025
``-``   0.05  -- 0.063      ``<``   0.0016 -- 0.002
``.``   0.04  -- 0.05       ``=``   0.0013 -- 0.0016
``/``   0.032 -- 0.04       ``>``   0.001  -- 0.0013
======  =================   ======  =================

Other symbols indicate lower error probabilities, and "~" is the
lowest possible.  In general::

  Error probability <= 10 ^ -((ASCII value - 33) / 10)

The "mismap" is simply the lowest probability from the "p" line.  (If
you run last-split twice, as in the genome alignment recipes, the
mismap is the lowest combined error probability from both "p" lines.)

Tips
----

To omit alignments with mismap probability > ``10^-6`` (say), you can
use the ``-m`` option (see below), or do this::

  awk -F= '/^a/ {i = $3 <= 1e-6} i' out.maf > out2.maf

FAQ
---

:Q: Before aligning RNA, should poly-A tails be trimmed?

:A: It's not essential, but it might make things faster.  Poly-A
    tracts tend to have many matches in the genome.  By trimming, you
    can prevent lastal and last-split from wasting time on such
    matches.

Split versus spliced alignment
------------------------------

Here is a split alignment::

  Query         ttctttgat--gctagtcctgatgttatggtattttttatcgaatgataa
                  |||||||--||||||                |||x||||||||||||
  Genome chrA  ...ctttgatatgctagt...             |||x||||||||||||
  Genome chrB                                 ...tttatatcgaatgata...

And here is a spliced alignment::

  Query        ctagtcgatatt--gctgtacgtctgttagctat-tttttcctctgtttg
                  |||x|||||--|||||||||----|||||||-|||||x|||||
  Genome chrA  ...gtctatattatgctgtacgt... |||||||-|||||x|||||
  Genome chrB                          ...tagctatattttttctctg...

Split alignment allows arbitrarily large unaligned parts in the middle
of the query, whereas spliced alignment applies a standard gap
penalty.  (Both allow arbitrarily large unaligned parts at the edges
of the query.)

Specialized examples
--------------------

Faster spliced alignment
~~~~~~~~~~~~~~~~~~~~~~~~

Spliced alignment can be slow.  It can be sped up, at a small cost in
accuracy, by not favoring cis-splices::

  lastal -p train.out -D10 db q.fastq | last-split -c0 -t0.004 -g db > out.maf

The ``-c0`` turns off cis-splicing, and the ``-t0.004`` specifies a
higher probability of "trans-splicing" (which includes cis-splicing,
just doesn't favor it).

"Spliced" alignment of DNA reads to a genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If we do not wish to allow arbitrarily large unaligned parts in the
middle of the query, we can do "spliced" alignment without considering
splice signals or favoring cis-splices::

  lastal -p train.out db q.fastq | last-split -c0 > out.maf

Options
-------

-h, --help
       Show a help message, with default option values, and exit.

-f, --format=FMT
       Choose the output format: ``MAF`` (without "p" lines), or
       ``MAF+`` (with "p" lines).  The format name is not
       case-sensitive.  The default is ``MAF`` (unless the input
       alignments have "p" lines from ``lastal -j``, in which case
       the default is ``MAF+``).

-r, --reverse
       Reverse the roles of the 2 sequences in each alignment: use the
       1st (top) sequence as the "query" and the 2nd as the
       "reference".

-g, --genome=NAME
       Do spliced alignment, and read splice signals (GT, AG, etc)
       from the named genome.  NAME should be the name of a lastdb
       database.

-d, --direction=D
       Do spliced alignment, and set the strandedness of the
       queries: 0=antisense, 1=sense, 2=unknown/mixed.  This
       determines whether forward and/or reverse-complement splice
       signals are used.

       If you use ``-d2``, the output will have an extra ``sense``
       field, indicating the log-odds that the query is
       sense-stranded::

	   log2[ prob(sense) / prob(antisense) ]

       The donor and acceptor annotations also indicate strandedness.
       This order means that the query sequence is from the forward
       strand of a spliced RNA::

           don
           acc don
           acc don
           acc

       And this order means it is from the reverse strand::

           acc
           don acc
           don acc
           don

       It's possible for ``don`` and ``acc`` to conflict with
       ``sense``.  That happens if one orientation is overall more
       likely (indicated by ``sense``), but any one alignment with
       that orientation is not the most likely.

-c, --cis=PROB
       Do spliced alignment, and set the average probability per
       base of cis-splicing.  The default value roughly fits human
       RNA.

-t, --trans=PROB
       Do spliced alignment, and set the average probability per
       base of trans-splicing.

-M, --mean=MEAN
       Do spliced alignment, and set the mean of ln(intron length).
       The default value fits human RNA.

-S, --sdev=SDEV
       Do spliced alignment, and set the standard deviation of
       ln(intron length).  The default value fits human RNA.

-m, --mismap=PROB
       Don't write alignments with mismap probability > PROB.

-s, --score=INT
       Don't write alignments with score < INT.

       For SPLIT alignment, the default value is e (the lastal score
       threshold).  Alignments with score just above INT will get
       high mismap probabilities.

       For SPLICED alignment, the default value is e + t * ln(100),
       where t is a scale factor that is written in the lastal
       header.  This roughly means that, for every alignment it
       writes, it has considered alternative alignments with
       one-hundredth the probability.  Alignments with score just
       above INT will not necessarily get high mismap probabilities.

-n, --no-split
       Do probability calculations as usual, but write the
       *original* alignments, annotated with "p" lines and mismap
       probabilities.  Note that the mismap and score limits still
       apply.

-b, --bytes=B
       Skip any query sequence that would require more than B bytes
       of memory to process.  (This only limits the size of some
       core data-structures: the total memory use will be greater.)
       A warning is written for each skipped sequence.  You can use
       suffixes such as K (KibiBytes), M (MebiBytes), G (GibiBytes),
       T (TebiBytes), e.g. ``-b20G``.

-v, --verbose
       Show progress information on the screen.

-V, --version
       Show version information and exit.

Details
-------

* The input must be in MAF format, and it must include header lines
  (of the kind produced by lastal) describing the alignment score
  parameters.

* The input must not mix alignments of different query sequences.  In
  other words, all the alignments of one query must be next to each
  other.  If you use ``-r``/``--reverse``, however, there is no such
  requirement, and the whole input gets read into memory.

* lastal can optionally write "p" lines, indicating the probability
  that each base is misaligned due to wrong gap placement.
  last-split, on the other hand, writes "p" lines indicating the
  probability that each base is aligned to the wrong genomic locus.
  You can combine both sources of error (roughly) by taking the
  maximum of the two error probabilities for each base.

The following points matter only if you are doing something unusual
(e.g. bisulfite alignment):

* If the header has more than one score matrix, last-split will use
  the first one.

* It assumes this score matrix applies to all alignments, when the
  alignments are oriented to use the forward strand of the query.

last-split5
-----------

last-split5 is almost identical to last-split.  The only difference is
the ``-g`` option: last-split can only read the output of lastdb,
whereas last-split5 can only read the output of lastdb5_.

Limitations
-----------

last-split does not support:

* DNA-versus-protein alignments.
* Generalized affine gap costs.

To do
-----

* An option to specify splice signals and their strengths.

.. _lastdb5: doc/lastdb.rst
.. _cookbook: doc/last-cookbook.rst
