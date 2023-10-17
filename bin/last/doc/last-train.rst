last-train
==========

last-train finds the rates (probabilities) of insertion, deletion, and
substitutions between two sets of sequences.  It thereby finds
suitable substitution and gap scores for aligning them.

It (probabilistically) aligns the sequences using some initial score
parameters, then estimates better score parameters based on the
alignments, and repeats this procedure until the parameters stop
changing.

The usage is like this::

  lastdb mydb reference.fasta
  last-train mydb queries.fasta

last-train prints a summary of each alignment step, followed by the
final score parameters, in a format that can be read by `lastal's -p
option <doc/lastal.rst>`_.

last-train can read .gz files, or from pipes::

  bzcat queries.fasta.bz2 | last-train mydb

Options
-------

-h, --help
       Show a help message, with default option values, and exit.
-v, --verbose
       Show more details of intermediate steps.

Training options
~~~~~~~~~~~~~~~~

--revsym
       Force the substitution scores to have reverse-complement
       symmetry, e.g. score(A→G) = score(T→C).  This is often
       appropriate, if neither strand is "special".
--matsym
       Force the substitution scores to have directional symmetry,
       e.g. score(A→G) = score(G→A).
--gapsym
       Force the insertion costs to equal the deletion costs.
--pid=PID
       Ignore alignments with > PID% identity (matches / [matches +
       mismatches]).  This aims to optimize the parameters for
       low-similarity alignments (similarly to the BLOSUM matrices).
--postmask=NUMBER
       By default, last-train ignores alignments of mostly-lowercase
       sequence (by using `last-postmask <doc/last-postmask.rst>`_).
       To turn this off, do ``--postmask=0``.
--sample-number=N
       Use N randomly-chosen chunks of the query sequences.  The
       queries are chopped into fixed-length chunks (as if they were
       first concatenated into one long sequence).  If there are ≤ N
       chunks, all are picked.  Otherwise, if the final chunk is
       shorter, it is never picked.  0 means use everything.
--sample-length=L
       Use randomly-chosen chunks of length L.
--scale=S
       Output scores in units of 1/S bits.  Traditional values
       include 2 for half-bit scores and 3 for 1/3-bit scores.
       (Note that 1/3-bit scores essentially equal Phred scores
       a.k.a. decibans, because log10(2) ≈ 3/10.)  The default is to
       infer a scale from the initial score parameters.
--codon
       Do training for DNA query sequences versus protein reference
       sequences.  These options will be ignored: ``--revsym
       --matsym --gapsym --postmask -q -p -S``.  If ``--pid`` is used,
       matches are defined by the genetic code, which is inferred from
       the substitution rates.  When ``--codon`` is used, the "initial
       parameter options" are initial probabilities, not scores/costs.

All options below this point are passed to lastal to do the
alignments: they are described in more detail at `<doc/lastal.rst>`_.

Initial parameter options
~~~~~~~~~~~~~~~~~~~~~~~~~

-r SCORE   Initial match score.
-q COST    Initial mismatch cost.
-p NAME    Initial match/mismatch score matrix.
-a COST    Initial gap existence cost.
-b COST    Initial gap extension cost.
-A COST    Initial insertion existence cost.
-B COST    Initial insertion extension cost.
-F LIST    Initial frameshift probabilities (only used with ``--codon``).

Alignment options
~~~~~~~~~~~~~~~~~

-D LENGTH  Query letters per random alignment.  (See `here
           <doc/last-evalues.rst>`_.)
-E EG2     Maximum expected alignments per square giga.  (See `here
           <doc/last-evalues.rst>`_.)
-s NUMBER  Which query strand to use: 0=reverse, 1=forward, 2=both.
           If specified, this parameter is written in last-train's
           output, so it will override lastal's default.
-S NUMBER  Specify how to use the substitution score matrix for
           reverse strands.  If you use ``--revsym``, this makes no
           difference.  "0" means that the matrix is used as-is for
           all alignments.  "1" (the default) means that the matrix
           is used as-is for alignments of query sequence forward
           strands, and the complemented matrix is used for query
           sequence reverse strands.

           This parameter is always written in last-train's output,
           so it will override lastal's default.

-C COUNT   Before extending gapped alignments, discard any gapless
           alignment whose query range lies in COUNT other gapless
           alignments with higher score-per-length.  This aims to
           reduce run time.
-T NUMBER  Type of alignment: 0=local, 1=overlap.
-R DIGITS  Lowercase & simple-sequence options.  If specified, this is
           written in last-train's output, so it will override
           lastal's default.
-m COUNT   Maximum number of initial matches per query position.
-k STEP    Look for initial matches starting only at every STEP-th
           position in each query.
-P COUNT   Number of parallel threads.
-X NUMBER  How to score a match/mismatch involving N (for DNA) or X
           (otherwise).  By default, the lowest match/mismatch score
           is used. 0 means the default; 1 means treat reference
           Ns/Xs as fully-ambiguous letters; 2 means treat query
           Ns/Xs as ambiguous; 3 means treat reference and query
           Ns/Xs as ambiguous.

           If specified, this parameter is written in last-train's
           output, so it will override lastal's default.

-Q NAME    How to read the query sequences (the NAME is not
           case-sensitive)::

             Default         fasta
             "0", "fastx"    fasta or fastq: discard per-base quality data
             "1", "sanger"   fastq-sanger

           The ``fastq`` formats are described here:
           `<doc/lastal.rst>`_.  last-train assumes the per-base
           quality codes indicate substitution error probabilities,
           *not* insertion or deletion error probabilities.  If this
           assumption is dubious (e.g. for data with many insertion
           or deletion errors), it may be better to discard the
           quality data.  For ``fastq-sanger``, last-train finds the
           rates of substitutions not explained by the quality data
           (ideally, real substitutions as opposed to errors).

           If specified, this parameter is written in last-train's
           output, so it will override lastal's default.

Details
-------

* last-train (and lastal) uses "Model A", in Figure 5A of btz576_.

* last-train (and lastal) converts between path and alignment
  parameters as in Supplementary Section 3.1 of btz576_.

* last-train uses parameters with "homogeneous letter probabilities"
  and "balanced length probability" (btz576_).

* last-train rounds the scores to integers, which makes them slightly
  inaccurate.  It then finds an adjusted scale factor (without
  changing the scores), which makes the integer-rounded scores
  correspond to homogeneous letter probabilities and balanced length
  probability.  It writes this adjusted scale (in nats, not bits) as a
  "-t" option for lastal, e.g. "-t4.4363".

* In rare cases, it may be impossible to find such an adjusted scale
  factor.  If that happens, last-train doubles the original scale (to
  reduce the inaccuracy of integer rounding), until the problem goes
  away.

.. _btz576: https://doi.org/10.1093/bioinformatics/btz576

Bugs
----

* last-train assumes that gap lengths roughly follow a geometric
  distribution.  If they do not (which is often the case), the results
  may be poor.

* last-train can fail for various reasons, e.g. if the sequences are
  too dissimilar.  If it fails to find any alignments, you could try
  reducing the alignment significance_ threshold with option ``-D``.

.. _significance: doc/last-evalues.rst
