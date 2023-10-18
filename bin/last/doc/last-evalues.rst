LAST E-values
=============

It is useful to know whether a similarity is significant,
i.e. unlikely to occur by chance between random sequences.  LAST
indicates this by EG2 and E::

  a score=40 EG2=0.031 E=0.025
  s mouse.chrY  2905908 130 +  91744698 CAGTCAC---AAATTTCTATCAAATATAA--CAGCT...
  s human.chrX 57162275 127 + 156040895 CACTCACTGCAAATCCTTATCAAACGTAATTCAGCT...

  a score=40 EG2=0.031 E=0.0093
  s mouse.chrY  2905908 130 + 91744698 CAGTCAC---AAATTTCTATCAAATATAA--CAGCT...
  s human.chrY 24460641 127 + 57227415 CACTCACTGCAAATCCTTATCAAACGTAATTCAGCT...

**EG2** is "expected alignments per square giga".  This is the
expected number of alignments with greater or equal score, between two
randomly-shuffled sequences of length 1 billion each.

**E** is the expected number of alignments with greater or equal
score, between: a random sequence with the same length as the query
sequence, and random sequences with the same length as the database.

In this example, the two alignments are identical, but their "E"
values are not.  That is because the query sequences have different
lengths: human.chrX is longer than human.chrY.

Default threshold
-----------------

By default, LAST reports alignments that are expected by chance at
most once per million query letters (for a given database).

Setting a threshold
-------------------

You can make lastal report alignments that are expected by chance at
most once per (say) thousand query letters, with option ``-D``::

  lastal -D1000 humdb fuguMito.fa > myalns.maf

You can make lastal report alignments with EG2 ≤ (say) 10, with option
``-E``::

  lastal -E10 humdb fuguMito.fa > myalns.maf

Advanced issues
---------------

Letter frequencies
~~~~~~~~~~~~~~~~~~

LAST's E-values refer to randomly-shuffled sequences *with specific
letter frequencies*.  Those frequencies are determined by the
substitution score matrix.  You can see them by running lastal with
option -v (verbose).

If your sequences have very different letter frequencies (e.g. very
AT-rich DNA), the E-values will be misleading.  The best solution is
to use a suitable score matrix, such as AT77 or ATMAP.

Calculation of EG2
~~~~~~~~~~~~~~~~~~

EG2 is calculated from the score like this::

  EG2 = K * exp(-lambda * score) * (1 billion) * (1 billion)

The parameters lambda and K are printed in the header of lastal's
output.

Calculation of E
~~~~~~~~~~~~~~~~

E is calculated by this formula::

  E = K * exp(-lambda * score) * strands * area(m, n, score) * N / n,

where:

* ``strands`` is either 1 or 2 (if forward and reverse DNA strands are
  searched).
* ``m`` is the length of the query sequence.
* ``n`` is the length of the longest sequence in the database.
* ``N`` is the total length of all database sequences.
* ``area(m, n, score)`` is `slightly less than
  <https://doi.org/10.1186/1756-0500-5-286>`_ ``m * n``.

Effect of option -D
~~~~~~~~~~~~~~~~~~~

Option -D sets the minimum alignment score, by this formula::

  K * exp(-lambda * score) * strands * D * length(n, score) * N / n  ≤  1,

where ``length(n, score)`` is slightly less than ``n``.

Bit score
~~~~~~~~~

Some people like to use "bit scores".  If you are one of them, the bit
score can be calculated like this::

  Bit score = (lambda * score - ln[K]) / ln[2]
            = log2( 1e18 / EG2 )

Limitations
~~~~~~~~~~~

* E-values cannot be calculated for scoring schemes with weak mismatch
  or gap costs (e.g. match score = 99, mismatch cost = 1).  In such
  cases, lastal requires you to set a score threshold with option -e.

* There may be a long startup time to calculate the E-value parameters
  (lambda, K, and finite-size correction parameters).  This is
  avoided, for particular scoring schemes, by storing the parameters
  in LAST's source code.  Parameters for other scoring schemes can be
  added on request.  (Or you can do it yourself: run lastal with
  option -v, so it prints the E-value parameters, then copy them into
  LastEvaluer.cc.)

* The E-values do not account for lastal option -c.  If you use this
  option, the E-values will be under-estimates.

* The E-values are for local alignment.  So if you use lastal option
  -T1 (overlap alignment), the E-values will be over-estimates.

Other resources
~~~~~~~~~~~~~~~

For more flexible E-value calculation, try ALP and FALP:
http://www.ncbi.nlm.nih.gov/CBBresearch/Spouge/html_ncbi/html/index/software.html
