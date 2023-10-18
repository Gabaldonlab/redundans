LAST Performance Tuning
=======================

Here are some ways to trade-off **speed**, **sensitivity**, and
**memory and disk usage**.

Sparsity options
~~~~~~~~~~~~~~~~

It's advisable to use at most one sparsity option: combining them will
likely give poor sensitivity.

lastal -k
---------

By default lastal looks for initial matches starting at every position
in the query sequence(s), but -k2 makes it check every 2nd position,
-k3 every 3rd position, etc.  Compared to the other sparsity options,
this **increases speed** the most while **reducing sensitivity** the
least.

lastdb -w
---------

By default lastdb indexes every position in the reference sequence(s),
but -w2 makes it index every 2nd position, -w3 every 3rd position,
etc.  Compared to the other sparsity options, this **decreases memory
and disk use** the most while reducing sensitivity the least.

This has a complex effect on the speed and sensitivity of lastal.
LAST uses initial matches that are sufficiently rare: this option
makes it lose some matches (because those positions are not indexed),
but gain others (because they are rarer).  In practice, small values
of w (e.g. 2) sometimes make lastal slower and more sensitive, but
very large values will eventually reduce sensitivity.

Among other aligners, MegaBLAST indexes every 5th position, and BLAT
indexes every 11th position.

lastdb -u
---------

``-uRY4`` makes LAST check for initial matches starting at only ~1/4
of positions, in both query and reference.  This **decreases memory
and disk use**, and **increases speed** of lastdb and lastal, but
**reduces sensitivity**.  ``-uRY8`` checks only ~1/8 of positions,
``-uRY16`` checks only ~1/16 of positions, ``-uRY32`` checks only
~1/32 of positions.  These ``-uRY`` options only apply to DNA, not
protein.

lastdb -W
---------

This is another method to check for initial matches starting at only
some positions, in both query and reference.  It's slightly inferior,
but more flexible: you can use any sparsity, and it works for protein.

Specifically, this makes LAST look for initial matches starting only
at positions that are "minimum" in any window of W consecutive
positions.  "Minimum" means that the sequence starting here is
alphabetically earliest.

The fraction of positions that are "minimum" is roughly: 2 / (W + 1).

lastdb5 & lastal5
~~~~~~~~~~~~~~~~~

If your reference has more than about 4 billion letters, 5-byte LAST
may be beneficial.  Ordinary (4-byte) LAST cannot directly handle so
much data, so it splits it into volumes, which is inefficient.  5-byte
LAST can handle such data without voluming, but it uses more memory.

5-byte LAST combines well with the lastdb sparsity options, which
reduce memory usage.  Something like ``lastdb5 -uRY32`` enables rapid,
huge-scale homology search, with moderate memory usage, but low
sensitivity.

Other options
~~~~~~~~~~~~~

lastal -m
---------

This option **trades speed for sensitivity**.  It sets the rareness
limit for initial matches: initial matches are lengthened until they
occur at most this many times in the lastdb volume.  The default is
10.  So -m100 makes it more sensitive but slower, by using more
initial matches.

lastal -l
---------

This option makes lastal **faster** but **less sensitive**.  It sets
the minimum length of initial matches, e.g. -l50 means length 50.
(The default is 1).  This can make it *much* faster, and the
sensitivity is adequate if the alignments contain long, gapless,
high-identity matches.

lastal -C
---------

This option (gapless alignment culling) can make lastal **faster** but
**less sensitive**.  It can also **reduce redundant output**.  For
example, ``-C2`` makes it discard alignments (before gapped extension)
whose query coordinates lie in those of 2 or more stronger alignments.
This works well for aligning long, repeat-rich, indel-poor sequences
(e.g. mammal chromosomes) without repeat-masking.

lastal -z
---------

This option can make lastal **faster** but **less sensitive**.  It
sets the maximum score drop in alignments, in the gapped extension
phase.  Lower values make it faster, by quitting unpromising
extensions sooner.  The default aims at best accuracy.

You can set this option in several ways: perhaps the most intuitive is
via maximum gap length.  For example, ``-z10g`` sets the maximum score
drop such that the longest possible gap length is 10.

lastal -x
---------

This option (preliminary gapped extension) can make lastal **faster**
but **less sensitive**.  For example, ``-x2g`` makes it extend gapped
alignments with a maximum gap length of 2, discard those with score
below the gapped score threshold, then redo the survivors with the
final max score drop (z).

lastal -M
---------

This option requests "minimum-difference" alignment, which is **faster
but cruder** than standard gapped alignment.  This treats all matches
the same, and minimizes the number of differences (mismatches plus
gaps).

lastal -j1
----------

This option requests **gapless** alignment, which is even **faster**.
(You could get the same effect by using very high gap costs, but
``-j1`` is faster because it skips the gapping phase entirely.)

lastal -f
---------

Option ``-fTAB`` **reduces the output size**, which can improve speed.

lastdb -i
---------

This option **makes lastdb faster**, but disables some lastal options.
If lastdb is too slow, try ``-i10``.

lastdb -C2
----------

**Obsolete(?)** This made old versions (< 1253) of lastal faster, but
makes new versions of lastal slower.  If you already used this option,
you can undo it by deleting the ``.chi2`` files (or move/rename them
to test which is faster).

lastdb -B
---------

Lower values (e.g. 1) make lastal **faster**, but use **more memory
and disk**.  This has no effect on the results.

Repeat masking
--------------

This can make LAST **much faster**, produce **less output**, and
reduce memory and disk usage.  Please see `<doc/last-repeats.rst>`_.
