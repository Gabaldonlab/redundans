Running LAST in parallel
========================

You can make LAST faster by running it on multiple CPUs / cores.  The
easiest way is with lastal_'s -P option::

  lastal -P8 my-index queries.fasta > out.maf

This will use 8 parallel threads.  If you specify -P0, it will use as
many threads as your computer claims it can handle simultaneously.

This works by aligning different query sequences in different threads
- so if you only have one query you won't get any parallelization!

Long query sequences use more memory
------------------------------------

If the query sequences are very long (e.g. chromosomes), many threads
may use a lot of memory, because each thread holds one query and its
alignments.  So you have to trade thread number against memory use.

Dealing with pipelines
----------------------

If you have a multi-command "pipeline", such as::

  lastal -P4 my-index queries.fasta | last-split > out.maf

then the ``-P`` option may help, because lastal is often the slowest
step, but it would be nice to parallelize the whole thing.

You can use ``parallel-fasta`` and ``parallel-fastq``, which accompany
LAST, but require `GNU Parallel`_ to be installed.

Instead of this::

  lastal -P8 -p my.train mydb seqs.fasta | last-split > out.maf

try this::

  parallel-fasta -j8 "lastal -p my.train mydb | last-split" < seqs.fasta > out.maf

Instead of this::

  lastal -P8 -p my.train mydb seqs.fastq.gz | last-split | gzip > out.maf.gz

try this::

  zcat seqs.fastq.gz | parallel-fastq -j8 "lastal -p my.train mydb | last-split | gzip" > out.maf.gz


Notes:

* ``parallel-fasta`` and ``parallel-fastq`` simply run GNU Parallel
  with a few options for fasta or fastq: you can specify other GNU
  Parallel options (like ``-j8``).

* The ``gzip`` example above runs ``gzip`` several times in parallel,
  and concatenates the results.  Fortunately, that is valid ``gz``
  format.

* ``lastal`` reads the ``mydb`` files into shared memory, so this
  memory use does not get mutiplied 8 times.

* ``parallel-fastq`` assumes that each fastq record is 4 lines, so
  there should be no line wrapping or blank lines.  (Actually, it
  assumes one "record" is 8 lines, to keep paired reads together.)

* In the current version, ``parallel-fasta`` and ``parallel-fastq``
  use GNU Parallel's ``--round`` option.  This avoids the overhead of
  restarting the command for each input chunk, but it prevents keeping
  the order of sequences, and stashes the output in large temporary
  files.  If necessary, use ``$TMPDIR`` or ``--tmpdir`` to choose a
  temp directory with enough space.

.. _lastal: doc/lastal.rst
.. _GNU parallel: http://www.gnu.org/software/parallel/
