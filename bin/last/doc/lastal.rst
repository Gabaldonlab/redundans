lastal
======

This program finds local alignments between query sequences, and
reference sequences that have been prepared using lastdb_.  You can
use it like this::

  lastdb humanDb humanChromosome*.fasta
  lastal humanDb dna*.fasta > myalns.maf

The lastdb command reads files called ``humanChromosome*.fasta`` and
writes several files whose names begin with ``humanDb``.  The lastal
command reads files called ``dna*.fasta``, compares them to
``humanDb``, and writes alignments to a file called ``myalns.maf``.

You can use gzip (.gz) compressed query files.  You can also pipe
query sequences into lastal, for example::

  bzcat seqs.fasta.bz2 | lastal humanDb > myalns.maf

Steps in lastal
---------------

1) Find initial matches.  For each possible start position in the
   query: find the shortest match with length ≥ l that *either* occurs
   ≤ m times in the reference, *or* has length L.

2) Extend a gapless alignment from each initial match, and keep those
   with score ≥ d.

3) Define cores: find the longest run of identical matches in each
   gapless alignment.

4) Extend a gapped alignment from either side of each core, and keep
   those with score ≥ e.

5) Non-redundantize the alignments: if several alignments share an
   endpoint (same coordinates in both sequences), remove all but one
   highest-scoring one.

6) Estimate the ambiguity of each aligned column (OFF by default).

7) Redo the alignments to minimize column ambiguity, using either
   gamma-centroid or LAMA (OFF by default).

8) Cut the alignments down to a unique best alignment for each part of
   each query sequence (OFF by default).

Options
-------

Cosmetic options
~~~~~~~~~~~~~~~~

-h, --help
    Show all options and their default settings, and exit.

-V, --version
    Show version information, and exit.

-v  Be verbose: write messages about what lastal is doing.

-f NAME
    Choose the output format.  The NAME is not case-sensitive.

    **MAF** format looks like this::

      a score=15 EG2=4.7e+04 E=2.6e-05
      s chr3 9 23 + 939557 TTTGGGAGTTGAAGTTTTCGCCC
      s seqA 2 21 +     25 TTTGGGAGTTGAAGGTT--GCCC

    Lines starting with "s" contain: the sequence name, the start
    coordinate of the alignment, the number of sequence letters
    spanned by the alignment, the strand, the sequence length, and
    the aligned letters.  The start coordinates are zero-based.  If
    the strand is "-", the start coordinate is in the reverse
    strand.

    The same alignment in **TAB** format looks like this::

      15 chr3 9 23 + 939557 seqA 2 21 + 25 17,2:0,4 EG2=4.7e+04 E=2.6e-05

    The ``17,2:0,4`` shows the sizes and offsets of gapless blocks in
    the alignment.  In this case, we have a block of size 17, then an
    offset of size 2 in the upper sequence and 0 in the lower
    sequence, then a block of size 4.

    The same alignment in **BlastTab** format looks like this::

      seqA chr3 86.96 23 1 1 3 23 10 32 2.6e-05 44.3

    The fields are: query name, reference name, percent identity,
    alignment length, mismatches, gap opens, query start, query end,
    reference start, reference end, E-value, bit score.  The start
    coordinates are one-based.  *Warning:* this is a lossy format,
    because it does not show gap positions.  *Warning:* the other
    LAST programs cannot read this format.  *Warning:* `"bit score"
    is not the same as "score" <doc/last-evalues.rst>`_.

    **BlastTab+** format is the same as BlastTab, with 3 extra
    columns at the end: length of query sequence, length of
    reference sequence, and (raw) score.  More columns might be
    added in future.

    For backwards compatibility, a NAME of 0 means TAB and 1 means
    MAF.

E-value options
~~~~~~~~~~~~~~~

-D LENGTH
    Report alignments that are expected by chance at most once per
    LENGTH query letters.  This option only affects the default value
    of ``-E``, so if you specify ``-E`` then ``-D`` has no effect.

-E THRESHOLD
    Maximum EG2 (`expected alignments per square giga
    <doc/last-evalues.rst>`_).  This option only affects the default
    value of ``-e``, so if you specify ``-e`` then ``-E`` has no
    effect.

Score options
~~~~~~~~~~~~~

-r SCORE
    Match score.

-q COST
    Mismatch cost.

-p NAME
    Specify a match/mismatch score matrix.  Options -r and -q will
    be ignored.  The built-in matrices are described in
    `<doc/last-matrices.rst>`_.

    Any other NAME is assumed to be a file name.  For an example of
    the format, see the matrix files in the data directory.
    Asymmetric scores are allowed: query letters correspond to
    columns and reference letters correspond to rows.

    Any letters that aren't in the matrix get default match/mismatch
    scores.  For doubly- and triply-ambiguous bases (such as "W"
    meaning A or T), these default scores are derived from the ACGT
    scores, and are shown in the header of lastal's output.  Any
    other letters get the lowest score in the matrix when aligned to
    anything.

    Other options can be specified on lines starting with ``#last``,
    but command line options override them.

-X NUMBER
    How to score a match/mismatch involving N (for DNA) or X
    (otherwise), if not specified by a score matrix.  By default,
    the lowest match/mismatch score is used.  0 means the default; 1
    means treat reference Ns/Xs as fully-ambiguous letters; 2 means
    treat query Ns/Xs as ambiguous; 3 means treat reference and
    query Ns/Xs as ambiguous.

-a COST
    Gap existence cost.

-b COST
    Gap extension cost.  A gap of size k costs: a + (b × k).

-A COST
    Insertion existence cost.  This refers to insertions in the
    query relative to the reference.  If -A is not used, the
    insertion existence cost will equal the deletion existence cost,
    which is set by -a.

-B COST
    Insertion extension cost.

-c COST
    This option allows use of "generalized affine gap costs" (SF
    Altschul 1998, Proteins 32(1):88-96).  Here, a "gap" may consist
    of unaligned regions of both sequences.  If these unaligned
    regions have sizes j and k, where j ≤ k, the cost is: a +
    b⋅(k-j) + c⋅j.  If c ≥ a + 2b (the default), it reduces to
    standard affine gaps.

-F LIST
    Align DNA queries to protein reference sequences, using the
    specified frameshift cost(s): either one cost (old-style
    frameshifts), or 4 comma-separated costs (new-style
    frameshifts).  As a special case, ``-F0`` means
    DNA-versus-protein alignment without frameshifts, which is
    faster.

    The four new-style frameshift costs are for, in order: deletion
    of length k mod 3 = 1 bases, deletion of k mod 3 = 2 bases,
    insertion of k mod 3 = 1 bases, insertion of k mod 3 = 2 bases.
    (You're expected to get them from last-train_, not set them
    manually.)  New-style frameshifts can only be used with "full
    scores", and old-style frameshifts can only be used with ordinary
    scores.

    The output looks like this::

      a score=108
      s prot 2  40 + 649 FLLQAVKLQDP-STPHQIVPSP-VSDLIATHTLCPRMKYQDD
      s dna  8 117 + 999 FFLQ-IKLWDP\STPH*IVSSP/PSDLISAHTLCPRMKSQDN

    The ``\`` indicates a forward shift by one nucleotide, and the
    ``/`` indicates a reverse shift by one nucleotide.  The ``*``
    indicates a stop codon.  The same alignment in tabular format
    looks like this::

      108 prot 2 40 + 649 dna 8 117 + 999 4,1:0,6,0:1,10,0:-1,19

    The "-1" indicates the reverse frameshift.

-z DROP
    Maximum score drop for gapped alignments.  Gapped alignments are
    forbidden from having any internal region with score < -DROP.
    The default value is e-1, which arguably produces the best
    alignments.  Lower values improve speed, by quitting unpromising
    extensions sooner.  You can specify this parameter in 3 ways:

    * A score (e.g. ``-z20``).

    * A percentage.  For example, ``-z50%`` specifies 50% of the
      default value (rounded down to the nearest integer).

    * A maximum gap length.  For example, ``-z8g`` sets the maximum
      score drop to: min[a+8b, A+8B].  However, this never increases
      the value above the default.

-x DROP
    This option makes lastal extend gapped alignments twice.  First,
    it extends gapped alignments with a maximum score drop of x, and
    discards those with score < e.  The surviving alignments are
    redone with a (presumably higher) maximum score drop of z.  This
    aims to improve speed with minimal effect on the final
    alignments.  You can specify -x in the same ways as -z (with the
    default value of x being z).

-y DROP
    Maximum score drop for gapless alignments.

-d SCORE
    Minimum score for gapless alignments.

-e SCORE
    Minimum alignment score.  (If you do gapless alignment with
    option -j1, then -d and -e mean the same thing.  If you set
    both, -e will prevail.)

Initial-match options
~~~~~~~~~~~~~~~~~~~~~

-m MULTIPLICITY
    Maximum multiplicity for initial matches.  Each initial match is
    lengthened until it occurs at most this many times in the
    reference.

    If the reference was split into volumes by lastdb_, then lastal
    uses one volume at a time.  The maximum multiplicity then applies
    to each volume, not the whole reference.  This is why voluming
    changes the results.

-l LENGTH
    Minimum length for initial matches.  Length means the number of
    letters spanned by the match.

-L LENGTH
    Maximum length for initial matches.

-k STEP
    Look for initial matches starting only at every STEP-th position
    in each query (positions 0, STEP, 2×STEP, etc).  This makes
    lastal faster but less sensitive.

-W SIZE
    Look for initial matches starting only at query positions that
    are "minimum" in any window of SIZE consecutive positions (see
    `<doc/lastdb.rst>`_).  By default, this parameter takes the same
    value as was used for lastdb -W.

Miscellaneous options
~~~~~~~~~~~~~~~~~~~~~

-s STRAND
    Specify which query strand should be used: 0 means reverse only,
    1 means forward only, and 2 means both.

-S NUMBER
    Specify how to use the substitution score matrix for reverse
    strands.  This matters only for unusual matrices that lack
    strand symmetry (e.g. if the a:g score differs from the t:c
    score).  "0" means that the matrix is used as-is for all
    alignments.  "1" means that the matrix is used as-is for
    alignments of query sequence forward strands, and the
    complemented matrix is used for query sequence reverse strands.

-K LIMIT
    Omit any alignment whose query range is contained in LIMIT or more
    other alignments with higher score (and on the same strand).  This
    is a useful way to get just the top few hits to each part of each
    query (P Berman et al. 2000, J Comput Biol 7:293-302).  As a
    special case, a LIMIT of 0 means: omit any alignment whose query
    range overlaps an alignment with higher score (and on the same
    strand).

-C LIMIT
    Before extending gapped alignments, discard any gapless
    alignment whose query range lies in LIMIT or more others (for
    the same strand and volume) with higher score-per-length.  This
    can reduce run time and output size (MC Frith & R Kawaguchi
    2015, Genome Biol 16:106).

-P THREADS
    Divide the work between this number of threads running in
    parallel.  0 means use as many threads as your computer claims it
    can handle simultaneously.  Single query sequences are not divided
    between threads, so you need multiple queries for this option to
    take effect.  With multiple threads, the order of the output is
    not fixed, but there are two guarantees:

    * All the alignments for one query sequence will appear together,
      and in a fixed order.
    * If each query sequence has length <= 2000, then pairs of queries
      stay together, e.g. the output for the 2nd query will be
      immediately after the output for the 1st query.

-i BYTES
    Process the query sequences in batches, of at most this many
    bytes.  If a single sequence exceeds this amount, however, it is
    not split.  You can use suffixes K, M, and G to specify KibiBytes,
    MebiBytes, and GibiBytes.  This option makes ``-P`` less
    efficient, because each batch is separately multi-threaded, but it
    fixes the output order to be the same as the input.

-M  Find minimum-difference alignments, which is faster but cruder.
    This treats all matches the same, and minimizes the number of
    differences (mismatches plus gaps).

    * Any substitution score matrix will be ignored.  The
      substitution scores are set by the match score (r) and the
      mismatch cost (q).
    * The gap cost parameters will be ignored.  The gap existence
      cost will be 0 and the gap extension cost will be q + r/2.
    * The match score (r) must be an even number.
    * Any sequence quality data (e.g. fastq) will be ignored.

-T NUMBER
    Type of alignment: 0 means "local alignment" and 1 means
    "overlap alignment".  Local alignments can end anywhere in the
    middle or at the ends of the sequences.  Overlap alignments must
    extend to the left until they hit the end of a sequence (either
    query or reference), and to the right until they hit the end of
    a sequence.

    **Warning:** it's often a bad idea to use -T1.  This setting
    does not change the maximum score drop allowed inside
    alignments, so if an alignment cannot be extended to the end of
    a sequence without exceeding this drop, it will be discarded.

-n COUNT
    Maximum number of gapless alignments per query position.  When
    lastal extends gapless alignments from initial matches that
    start at one query position, if it gets COUNT successful
    extensions, it skips any remaining initial matches starting at
    that position.

-N COUNT
    Stop after finding COUNT alignments per query strand.  This
    makes lastal faster: it quits gapless and gapped extensions as
    soon as it finds COUNT alignments with score ≥ e.

-R DIGITS
    Specify lowercase-marking of repeats, by two digits (e.g. "-R 01"),
    with the following meanings.

    First digit:

    0. Convert the input sequences to uppercase while reading them.
    1. Keep any lowercase in the input sequences.

    Second digit:

    0. Do not check for simple repeats.
    1. Convert simple repeats (e.g. cacacacacacacacac) to lowercase.
    2. Convert simple repeats, within AT-rich DNA, to lowercase.
    3. Convert simple repeats, including weaker simple repeats, to
       lowercase (with tantan's ``r`` parameter = 0.02).

    The default is to use the same ``-R`` setting as was used by
    lastdb_, except that if lastdb's 2nd ``R`` digit was ``3``, it
    defaults to ``1``.

    Details: Tantan_ is applied separately to forward and reverse
    strands.  For DNA-versus-protein alignment, if you use a codon
    substitution matrix (e.g. from ``last-train --codon``), tantan
    is applied to the DNA before translation, else it is applied
    after translation.

-u NUMBER
    Specify treatment of lowercase letters when extending
    alignments:

    0. Mask them for neither gapless nor gapped extensions.
    1. Mask them for gapless but not gapped extensions.
    2. Mask them for gapless but not gapped extensions, and then
       discard alignments that lack any segment with score ≥ e when
       lowercase is masked.  (For "full scores": mask them for gapless
       and gapped extensions, then recalculate the alignments *but not
       the scores* without masking.)
    3. Mask them for gapless and gapped extensions.

    "Mask" means change their match/mismatch scores to min(unmasked
    score, 0), a.k.a. `gentle masking`_.  (But if you use a codon
    substitution matrix, a lowercase-containing base-triplet will be
    scored as ``nnn``, which defaults to the lowest match/mismatch
    score.)

    This option does not affect treatment of lowercase for initial
    matches.

-w DISTANCE
    This option is a kludge to avoid catastrophic time and memory
    usage when self-comparing a large sequence.  If the sequence
    contains a tandem repeat, we may get a gapless alignment that is
    slightly offset from the main self-alignment.  In that case, the
    gapped extension might "discover" the main self-alignment and
    extend over the entire length of the sequence.

    To avoid this problem, gapped alignments are not triggered from
    any gapless alignment that:

    * is contained, in both sequences, in the "core" of another
      alignment
    * has start coordinates offset by DISTANCE or less relative to
      this core

    Use ``-w0`` to turn this off.

-G GENETIC-CODE
    Specify the genetic code for translating DNA to protein.  Codes
    are specified by numbers (e.g. 1 = standard, 2 = vertebrate
    mitochondrial), listed here:
    https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi.  Any
    other GENETIC-CODE is assumed to be a file name: for an example
    of the format, see vertebrateMito.gc in the examples directory.

-t TEMPERATURE
    Parameter for converting between scores and probability ratios.
    This affects the column ambiguity estimates.  A score is
    converted to a probability ratio by this formula: exp(score /
    TEMPERATURE).  The default value is 1/lambda, where lambda is
    the scale factor of the scoring matrix, which is calculated by
    the method of Yu and Altschul (YK Yu et al. 2003, PNAS
    100(26):15688-93).

-g GAMMA
    This option affects gamma-centroid and LAMA alignment only.

    Gamma-centroid alignments minimize the ambiguity of paired
    letters.  In fact, this method aligns letters whose column error
    probability is less than GAMMA/(GAMMA+1).  When GAMMA is low, it
    aligns confidently-paired letters only, so there tend to be many
    unaligned letters.  When GAMMA is high, it aligns letters more
    liberally.

    LAMA (Local Alignment Metric Accuracy) alignments minimize the
    ambiguity of columns (both paired letters and gap columns).
    When GAMMA is low, this method produces shorter alignments with
    more-confident columns, and when GAMMA is high it produces
    longer alignments including less-confident columns.

    In summary: to get the most accurately paired letters, use
    gamma-centroid.  To get accurately placed gaps, use LAMA.

    Note that the reported alignment score is that of the gapped
    alignment before realigning with gamma-centroid or LAMA.

-j NUMBER
    Output type: 0 means counts of initial matches (of all lengths);
    1 means gapless alignments; 2 means gapped alignments before
    non-redundantization; 3 means gapped alignments after
    non-redundantization; 4 means alignments with ambiguity
    estimates; 5 means gamma-centroid alignments; 6 means LAMA
    alignments; 7 means alignments with expected counts.

    If you use -j0, lastal will count the number of initial matches,
    per length, per query sequence.  Options -l and -L will set the
    minimum and maximum lengths, and -m will be ignored.  If you
    compare a large sequence to itself with -j0, it's wise to set
    option -L.

    If you use -j7, lastal will print an extra MAF line starting
    with "c" for each alignment.  The first 16 numbers on this line
    are the expected counts of matches and mismatches: first the
    count of reference As aligned to query As, then the count of
    reference As aligned to query Cs, and so on.  For proteins there
    will be 400 such numbers.  The next 5 numbers are expected
    counts related to gaps.  They are:

    * The count of matches plus mismatches.  (This may exceed the
      total of the preceding numbers, if the sequences have non-ACGT
      letters.)
    * The count of deleted letters.
    * The count of inserted letters.
    * The count of delete opens (= count of delete closes).
    * The count of insert opens (= count of insert closes).

-J NUMBER
    Score type: 0 means ordinary score, 1 means "full score" (also
    known as "forward score" or "sum-of-paths score").  Both types of
    score are measures of how significant a similarity is.  An
    ordinary score is based on one alignment, whereas a "full score"
    is based on many alternative ways of aligning the similar regions.
    Full scores are expected to be more sensitive, but they are not
    recognized by last-split_.  Full score E-values_ can be calculated
    only for parameters from last-train_.

-Q NAME
    Specify how to read the query sequences (the NAME is not
    case-sensitive)::

      Default           fasta
      "0", "fastx"      fasta or fastq: discard per-base quality data
      "keep"            fasta or fastq: keep but ignore per-base quality data
      "1", "sanger"     fastq-sanger
      "2", "solexa"     fastq-solexa
      "3", "illumina"   fastq-illumina
      "4", "prb"        prb
      "5", "pssm"       PSSM

    *Warning*: Illumina data is not necessarily in fastq-illumina
    format; it is often in fastq-sanger format.

    The fastq formats look like this::

      @mySequenceName
      TTTTTTTTGCCTCGGGCCTGAGTTCTTAGCCGCG
      +
      55555555*&5-/55*5//5(55,5#&$)$)*+$

    The "+" may be followed by text (ignored).  The symbols below
    the "+" are quality codes, one per sequence letter.  The
    sequence and quality codes may wrap onto more than one line.

    lastal assumes the quality codes indicate substitution error
    probabilities, *not* insertion or deletion error probabilities.
    If this assumption is dubious (e.g. for data with many insertion
    or deletion errors), it may be better to discard or ignore them.

    For fastq-sanger, quality scores are obtained by subtracting 33
    from the ASCII values of the quality codes.  For fastq-solexa
    and fastq-illumina, they are obtained by subtracting 64.

    prb format stores four quality scores (A, C, G, T) per position,
    with one sequence per line, like this::

      -40  40 -40 -40      -12   1 -12  -3      -10  10 -40 -40

    Since prb does not store sequence names, lastal uses the line
    number (starting from 1) as the name.

    In fastq-sanger and fastq-illumina format, the quality scores
    are related to error probabilities like this: qScore =
    -10⋅log10[p].  In fastq-solexa and prb, however, qScore =
    -10⋅log10[p/(1-p)].  In lastal's MAF output, the quality scores
    are written on lines starting with "q".  For fastq, they are
    written with the same encoding as the input.  For prb, they are
    written in the fastq-solexa (ASCII-64) encoding.

    Finally, PSSM means "position-specific scoring matrix".  The
    format is::

      myLovelyPSSM
           A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
      1 M -2 -2 -3 -4 -2 -1 -3 -3 -2  1  2 -2  8 -1 -3 -2 -1 -2 -2  0
      2 S  0 -2  0  1  3 -1 -1 -1 -2 -3 -3 -1 -2 -3 -2  5  0 -4 -3 -2
      3 D -1 -2  0  7 -4 -1  1 -2 -2 -4 -4 -2 -4 -4 -2 -1 -2 -5 -4 -4

    The sequence appears in the second column, and columns 3 onwards
    contain the position-specific scores.  Any letters not specified
    by any column will get the lowest score in each row.  This
    format is a simplified version of PSI-BLAST's ASCII format: the
    non-simplified version is allowed too.

    *Warning*: lastal cannot directly calculate E-values for PSSMs.
    The E-values (and the default value of -y) are determined by the
    otherwise-unused match and mismatch scores (options -r -q and
    -p).  There is evidence these E-values will be accurate if the
    PSSM is "constructed to the same scale" as the match/mismatch
    scores (SF Altschul et al. 1997, NAR 25(17):3389-402).

Split options
~~~~~~~~~~~~~

--split
    Cut the alignments down to a unique best alignment for each part
    of each query sequence.  This is useful for DNA queries that cross
    rearrangement breakpoints.  It's the same as running last-split_
    with default options.

--splice
    This is similar to ``--split``, but it favors alignments that
    would be expected due to intron/exon splicing.  So it favors
    alignments where parts of a query sequence are separated by
    intron-like gaps with GT-AG splice signals.  This is the same as
    running last-split_ with option ``-g``.

--split-f, --split-m, --split-s, --split-n, --split-b
    These are equivalent to last-split_ options ``-f``, ``-m``,
    ``-s``, ``-n``, and ``-b``.  If you use one of these, you don't
    need to also specify ``--split``.

--split-d, --split-c, --split-t, --split-M, --split-S
    These are equivalent to last-split_ options ``-d``, ``-c``,
    ``-t``, ``-M``, and ``-S``.  They imply ``--splice``, so you don't
    need to also specify ``--splice``.

Parallel processes and memory sharing
-------------------------------------

If you run several lastal commands (i.e. processes) at the same time
on the same computer, using the same set of reference files prepared
by lastdb, then they will share memory for the reference files.

Multiple volumes
----------------

If lastdb_ creates multiple volumes::

  lastdb hugeDb huge.fasta

You can either run lastal on the whole thing::

  lastal hugeDb queries.fasta > myalns.maf

Or on one volume at a time::

  lastal hugeDb0 queries.fasta > myalns0.maf
  lastal hugeDb1 queries.fasta > myalns1.maf
  lastal hugeDb2 queries.fasta > myalns2.maf

The former method reads the queries in large batches, and aligns each
batch to one volume at a time.  If you run several processes in
parallel, they will not necessarily use the same volume at the same
time.

Therefore, with parallel processes, you should either ensure you have
enough memory to hold several volumes simultaneously, or run lastal on
one volume at a time.  An efficient scheme is to use a different
computer for each volume.

lastal5
-------

lastal5 has identical usage to lastal, and is used with lastdb5_.
lastal cannot read the output of lastdb5, and lastal5 cannot read the
output of lastdb.

.. _lastdb5:
.. _lastdb: doc/lastdb.rst
.. _last-train: doc/last-train.rst
.. _last-split: doc/last-split.rst
.. _E-values: doc/last-evalues.rst
.. _tantan: https://gitlab.com/mcfrith/tantan
.. _gentle masking: https://doi.org/10.1371/journal.pone.0028819
