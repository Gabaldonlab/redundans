lastdb
======

This program prepares sequences for subsequent comparison and
alignment using lastal_.  You can use it like this::

  lastdb humanDb humanChromosome*.fasta

This will read files called ``humanChromosome*.fasta``, and write
several files whose names begin with ``humanDb``.

Input
-----

The input should be one or more files in fasta format, which looks
like this::

  >MyFirstSequence
  ATCGGGATATATGGAGAGCTTAGAG
  TTTGGATATG
  >My2ndSequence
  TTTAGAGGGTTCTTCGGGATT

These files may be compressed in gzip (.gz) format.  You can also pipe
sequences into lastdb, for example::

  bzcat humanChromosome*.fasta.bz2 | lastdb humanDb

Options
-------

Main Options
~~~~~~~~~~~~

-h, --help
    Show all options and their default settings, and exit.

-p  Interpret the sequences as proteins.  The default is to interpret
    them as DNA.

-c  Soft-mask lowercase letters.  This means that, when we compare
    these sequences to some other sequences using lastal, lowercase
    letters will be excluded from initial matches.  This will apply
    to lowercase letters in *both* sets of sequences.

-u NAME
    Specify a seeding scheme.  The -m and -d options will then be
    ignored.  The built-in schemes are described in
    `<doc/last-seeds.rst>`_.

    Any other NAME is assumed to be a file name.  For an example of
    the format, see the seed files in the data directory.  You can
    set other lastdb options on lines starting with ``#lastdb``, but
    command line options override them.  You can also set lastal_
    options on lines starting with ``#lastal``, which are overridden
    by options from a `scoring scheme <doc/last-matrices.rst>`_ or
    the lastal command line.

-P THREADS
    Make ``lastdb`` faster, with no effect on results, by running this
    many threads in parallel.  ``0`` means use as many threads as your
    computer claims it can handle simultaneously.

Advanced Options
~~~~~~~~~~~~~~~~

-R DIGITS
    Specify lowercase usage, by two digits (e.g. ``-R01``), with the
    following meanings.

    First digit:

    0. Convert the input sequences to uppercase while reading them.
    1. Keep any lowercase in the input sequences.

    Second digit:

    0. Do not check for simple repeats.
    1. Convert simple repeats (e.g. ``cacacacacacacacac``) to
       lowercase.  This uses tantan_, which reliably prevents
       non-homologous alignments, unlike other repeat finders.
    2. Convert simple DNA repeats to lowercase, with tantan tuned
       for ~80% AT-rich genomes.
    3. Convert simple repeats, including weaker simple repeats, to
       lowercase (with tantan's ``r`` parameter = 0.02).

    The default is ``-R01`` (unless ``-q`` is specified, in which case
    the default is ``-R03``).

-w STEP
    Allow initial matches to start only at every STEP-th position in
    each of the sequences given to lastdb (positions 0, STEP,
    2×STEP, etc).  This reduces the memory usage of lastdb and
    lastal, and it makes lastdb faster.  Its effect on the speed and
    sensitivity of lastal is not entirely clear.

-W SIZE
    Allow initial matches to start only at positions that are
    "minimum" in any window of SIZE consecutive positions.
    "Minimum" means that the sequence starting here is
    alphabetically earliest.

    The "alphabetical" order depends on the `seed pattern
    <doc/last-seeds.rst>`_.  The letter order is determined by the
    order of the letter groups, and letters in the same group are
    considered equivalent.

    The fraction of positions that are "minimum" is roughly: 2 /
    (SIZE + 1).

-S STRAND
    Specify which strand of the input sequences should be prepared
    for alignment: 0 means reverse only, 1 means forward only, and 2
    means both.

-s BYTES
    Limit memory usage, by splitting the output files into smaller
    "volumes" if necessary.  This will limit the memory usage of
    both lastdb and lastal, but it will make lastal slower.  It is
    also likely to change the exact results found by lastal.

    Each volume will have at most BYTES bytes.  (Roughly.  The more
    masked letters or DNA "N"s, the more it will undershoot.)  You
    can use suffixes K, M, and G to specify KibiBytes, MebiBytes,
    and GibiBytes (e.g. "-s 5G").

    However, the output for one sequence is never split.  Since the
    output files are several-fold bigger than the input (unless you
    use -w/-W/-d/-u), this means that mammalian chromosomes cannot
    be processed using much less than 2G (unless you use -w/-W/-d/-u).

    There is a hard upper limit of about 4 billion sequence letters
    per volume.  Together with the previous point, this means that
    lastdb will refuse to process any single sequence longer than
    about 4 billion.

-Q NAME
    Specify how to read the sequences (the NAME is not case-sensitive)::

      Default           fasta
      "0", "fastx"      fasta or fastq: discard per-base quality data
      "keep"            fasta or fastq: keep but ignore per-base quality data
      "1", "sanger"     fastq-sanger
      "2", "solexa"     fastq-solexa
      "3", "illumina"   fastq-illumina

    The fastq formats are described in `<doc/lastal.rst>`_.

-q  Interpret the sequences as proteins, use a 21-letter alphabet
    with ``*`` meaning STOP, and append ``*`` to each sequence.

-m PATTERN
    Specify a spaced seed pattern, for example "-m 110101".  In this
    example, mismatches will be allowed at every third and fifth
    position out of six in initial matches.

    This option does not constrain the length of initial matches.
    The pattern will get cyclically repeated as often as necessary
    to cover any length.

    Although the 0 positions allow mismatches, they exclude
    non-standard letters (e.g. non-ACGT for DNA).  If option -c is
    used, they also exclude lowercase letters.

    You can also specify transition constraints, e.g "-m 100TT1".
    In this example, transitions (but not transversions) will be
    allowed at every fourth and fifth position out of six.
    Alternatively, you can use Iedera_'s notation, for example
    "-m '#@#--##--#-#'" ('#' for match, '@' for transition, '-' or
    '_' for mismatch).

    You can specify multiple patterns by separating them with commas
    and/or using "-m" multiple times.

-d PATTERN
    Specify DNA seed patterns, for example: ``-d RYrNn@N,YyRn@NN`` .
    The symbols have these meanings::

      N  any match is allowed (a:a, c:c, g:g, t:t)
      n  any match or mismatch is allowed
      R  only purine matches are allowed (a:a, g:g)
      r  purine matches or mismatches are allowed (a:a, g:g, a:g, g:a)
      Y  only pyrimidine matches are allowed (c:c, t:t)
      y  pyrimidine matches or mismatches are allowed (c:c, t:t, c:t, t:c)
      A  only a:a matches are allowed
      C  only c:c matches are allowed
      G  only g:g matches are allowed
      T  only t:t matches are allowed
      @  any match or transition is allowed

    The pattern will get cyclically repeated as often as necessary
    to cover any length.  *However*, in 2nd and subsequent cycles,
    the base-restricted symbols are replaced with unrestricted
    symbols: ``RYACGT`` => ``N``, ``ry`` => ``@``.

-a SYMBOLS
    Specify your own alphabet, e.g. "-a 0123".  The default (DNA)
    alphabet is equivalent to "-a ACGT".  The protein alphabet (-p)
    is equivalent to "-a ACDEFGHIKLMNPQRSTVWY".  Non-alphabet
    letters are allowed in sequences, but by default they are
    excluded from initial matches and get the mismatch score when
    aligned to anything.  As a special case, for the DNA alphabet,
    Us are converted to Ts.  If -a is specified, -p is ignored.

-i MATCHES
    This option makes lastdb faster, at the expense of limiting your
    options with lastal_.  If you use (say) "-i 10", then you cannot
    use lastal with option m < 10.

-b LENGTH
    Specify the maximum length for "buckets", a lookup table for all
    possible sequences of length <= LENGTH.  Larger values increase
    the memory usage of lastdb and lastal, make lastal faster, and
    have no effect on lastal's results.

-B FOLD
    Use the maximum possible bucket length, such that the buckets'
    memory use is at most 1/FOLD times that of the stored sequence
    positions.

-C NUMBER
    Specify the type of "child table" to make: 0 means none, 1 means
    byte-size (uses a little more memory), 2 means short-size (uses
    somewhat more memory), 3 means full (uses a lot more memory).
    Choices > 0 make lastdb slower, and have no effect on lastal's
    results.  This option seems to be mostly obsolete: it used to make
    lastal faster, but now makes it slower (sometimes).

-x  Just count sequences and letters.  This is much faster.  Letter
    counting is never case-sensitive.

-D  Read pre-existing lastdb files, and print all the sequences in them.

-v  Be verbose: write messages about what lastdb is doing.

-V, --version
    Show version information, and exit.

lastdb5
-------

lastdb5 is identical to lastdb, except that it internally uses larger
(5-byte) integers.  This means it can handle more than 4 billion
sequence letters per volume, but it uses more memory.

Memory and disk usage
---------------------

Suppose we give lastdb N letters of sequence data, of which M are
non-masked "real" letters (e.g. excluding N for DNA and X for
proteins).  The output files will include:

* The sequences (N bytes).

* An "index" consisting of:
  positions (4M bytes), and "buckets" (<= M bytes).

* The sequence names (*usually* negligible).

This is modified by several options.

* -C1 adds M bytes to the index, -C2 adds 2M bytes, and -C3 adds 4M
  bytes.

* -w STEP: makes the index STEP times smaller.

* -W SIZE: makes the index about (SIZE+1)/2 times smaller.

* lastdb5: makes the index 25% bigger.

* -u, -m, -d: Multiple patterns multiply the index size.  For example,
  MAM8_ makes it 8 times bigger.

* -u, -d: may reduce the index, e.g. RY32_ makes it 32 times smaller.

* -s: does not change the total size, but splits it into volumes.

* -S2: doubles the size of everything.

Limitations
-----------

lastdb can become catastrophically slow for highly redundant
sequences, e.g. two almost-identical genomes.  It usually processes
several GB per hour, but if it becomes much slower, try option "-i
10", which is likely to solve the problem.

.. _lastal: doc/lastal.rst
.. _RY32:
.. _MAM8: doc/last-seeds.rst
.. _tantan: https://gitlab.com/mcfrith/tantan
.. _Iedera: https://bioinfo.lifl.fr/yass/iedera.php
