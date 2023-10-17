LAST FAQ
========

:Q: Why does LAST sometimes miss strong alignments?

:A: This is because it avoids getting huge numbers of repetitive
    alignments, by just finding a few strongest matches for each part
    of each query sequence.  To find more repetitive matches, increase
    the ``lastal -m`` parameter, e.g. ``-m100``.

..

:Q: Why does LAST sometimes show two separate alignments, and
    sometimes one merged alignment, between the same sequences?

:A: This is because LAST forbids highly dissimilar regions inside
    alignments, but the maximum-allowed dissimilarity depends on the
    alignment score / E-value threshold (which depends on the database
    size).  For more details, please see the discussion of "X-drop" in
    `Parameters for accurate genome alignment`_.

..

:Q: Does it matter which sequence is used as the "reference" (given to
    lastdb) and which is used as the "query" (given to lastal)?

:A: It may do.  In short, LAST tries hard to find alignments for every
    position in the query.  When mapping reads to a genome, you
    probably want the genome to be the reference, and the reads to be
    the query.  That way, for each read, it will search for several
    most-similar locations in the genome.  The other way, for each
    location in the genome, it will search for several most-similar
    reads.  As another example, if you compare a genome to a library
    of repeat sequences, you probably want the genome to be the query
    and the repeat library to be the reference.

..

:Q: How can I get the percent-identity of each alignment?

:A: Use lastal option -fBlastTab, or use maf-convert.

..

:Q: Why is LAST so slow at reading/writing files?

:A: Probably because the files are huge and your disk is slow (or lots
    of people are using it).  Try to use a reasonably fast disk.

    Typically, after reading a large file once, subsequent reads of
    the same file are much faster, because the operating system caches
    it.

    lastal reads files into shared memory using "mmap", and this
    occasionally seems to have trouble, notably on "advanced" file
    systems (Lustre and GlusterFS).

    On one Linux system, the above-mentioned cache occasionally got
    into some kind of bad state, making lastal very slow.  This was
    solved by running the following command (as root)::

      echo 1 > /proc/sys/vm/drop_caches

..

:Q: How does LAST get the sequence names?  How can I get nice, short,
    unique names?

:A: The first whitespace-delimited word in the sequence header line is
    used as the name.  You can arbitrarily customise the names using
    standard Unix tools.  For example, this will replace each FASTA
    name with a unique serial number::

      awk '/>/ {$0 = ">" ++n} 1' queries.fasta | lastal myDb

    This will do the same for FASTQ (assuming 4 lines per record,
    i.e. no line wrapping)::

      awk 'NR % 4 == 1 {$0 = "@" ++n} 1' queries.fastq | lastal myDb

    Sometimes you can make LAST's output significantly smaller by
    shortening the names.

..

:Q: How can I find alignments with > 95% identity?

:A: One way is to use a scoring scheme like this: +5 for a match, and
    -95 for a mismatch or a gap::

      lastal -r5 -q95 -a0 -b95

RAQ (Rarely Asked Questions)
----------------------------

:Q: Why does lastal hang forever on my quad processor system?

:A: This problem was reported once, and was fixed by prefacing the
    command with::

      numactl --interleave=all

.. _Parameters for accurate genome alignment: https://doi.org/10.1186/1471-2105-11-80
