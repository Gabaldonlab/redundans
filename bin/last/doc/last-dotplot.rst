last-dotplot
============

This script makes a dotplot, a.k.a. Oxford Grid, of pair-wise sequence
alignments in MAF_ or LAST tabular format.  It requires the `Python
Imaging Library`_ to be installed.  It can be used like this::

  last-dotplot my-alignments my-plot.png

The output can be in any format supported by the Imaging Library::

  last-dotplot alns alns.gif

It can read alignments from a pipe like this::

  ... | last-dotplot - my-plot.png

Terminology
-----------

last-dotplot shows alignments of one set of sequences against another
set of sequences.  This document calls a "set of sequences" a
"genome", though it need not actually be a genome.

Options
-------

-h, --help
    Show a help message, with default option values, and exit.
-v, --verbose
    Show progress messages & data about the plot.
-x INT, --width=INT
    Maximum width in pixels.
-y INT, --height=INT
    Maximum height in pixels.
-m M, --maxseqs=M
    Maximum number of horizontal or vertical sequences.  If there
    are >M sequences, the smallest ones (after cutting) will be
    discarded.
-1 PATTERN, --seq1=PATTERN
    Which sequences to show from the 1st (horizontal) genome.
-2 PATTERN, --seq2=PATTERN
    Which sequences to show from the 2nd (vertical) genome.
-a FILE
    Read annotations for the 1st (horizontal) genome, and draw them as
    vertical stripes.  For backwards compatibility, these options have
    the same meaning: ``--bed1 --rmsk1 --genePred1 --gap1``.
-b FILE
    Read annotations for the 2nd (vertical) genome, and draw them as
    horizontal stripes.  For backwards compatibility, these options
    have the same meaning: ``--bed2 --rmsk2 --genePred2 --gap2``.
--alignments=FILE
    Read secondary alignments.  For example: we could use primary
    alignment data with one human DNA read aligned to the human
    genome, and secondary alignment data with the whole chimpanzee
    versus human genomes.  last-dotplot will show the parts of the
    secondary alignments that are near the primary alignments.
--sort1=N
    Put the 1st genome's sequences left-to-right in order of: their
    appearance in the input (0), their names (1), their lengths (2),
    the top-to-bottom order of (the midpoints of) their alignments
    (3).  You can use two colon-separated values, e.g. "2:1" to
    specify 2 for primary and 1 for secondary alignments.
--sort2=N
    Put the 2nd genome's sequences top-to-bottom in order of: their
    appearance in the input (0), their names (1), their lengths (2),
    the left-to-right order of (the midpoints of) their alignments
    (3).
--strands1=N
    Put the 1st genome's sequences: in forwards orientation (0), in
    the orientation of most of their aligned bases (1).  In the
    latter case, the labels will be colored (in the same way as the
    alignments) to indicate the sequences' orientations.  You can
    use two colon-separated values for primary and secondary
    alignments.
--strands2=N
    Put the 2nd genome's sequences: in forwards orientation (0), in
    the orientation of most of their aligned bases (1).
--max-gap1=FRAC
    Maximum unaligned gap in the 1st genome.  For example, if two
    parts of one DNA read align to widely-separated parts of one
    chromosome, it's probably best to cut the intervening region
    from the dotplot.  FRAC is a fraction of the length of the
    (primary) alignments.  You can specify "inf" to keep all
    unaligned gaps.  You can use two comma-separated values,
    e.g. "0.5,3" to specify 0.5 for end-gaps (unaligned sequence
    ends) and 3 for mid-gaps (between alignments).  You can use two
    colon-separated values (each of which may be comma-separated)
    for primary and secondary alignments.
--max-gap2=FRAC
    Maximum unaligned gap in the 2nd genome.
--pad=FRAC
    Length of pad to leave when cutting unaligned gaps.
-j N, --join=N
    Draw horizontal or vertical lines joining adjacent alignments.
    0 means don't join, 1 means draw vertical lines joining
    alignments that are adjacent in the 1st (horizontal) genome, 2
    means draw horizontal lines joining alignments that are adjacent
    in the 2nd (vertical) genome.
--border-pixels=INT
    Number of pixels between sequences.

Text options
~~~~~~~~~~~~

-f FILE, --fontfile=FILE
    TrueType or OpenType font file.
-s SIZE, --fontsize=SIZE
    TrueType or OpenType font size.
--labels1=N
    Label the displayed regions of the 1st genome with their:
    sequence name (0), name:length (1), name:start:length (2),
    name:start-end (3).
--labels2=N
    Label the displayed regions of the 2nd genome with their:
    sequence name (0), name:length (1), name:start:length (2),
    name:start-end (3).
--rot1=ROT
    Text rotation for the 1st genome: h(orizontal) or v(ertical).
--rot2=ROT
    Text rotation for the 2nd genome: h(orizontal) or v(ertical).

Color options
~~~~~~~~~~~~~

-c COLOR, --forwardcolor=COLOR
    Color for forward alignments.
-r COLOR, --reversecolor=COLOR
    Color for reverse alignments.
--border-color=COLOR
    Color for pixels between sequences.
--margin-color=COLOR
    Color for the margins.
--exon-color=COLOR
    Color for exons.
--cds-color=COLOR
    Color for protein-coding regions.
--bridged-color=COLOR
    Color for unsequenced gaps with "yes" evidence of linkage.
--unbridged-color=COLOR
    Color for unsequenced gaps with "no" evidence of linkage.

Annotations
-----------

Options ``-a`` and ``-b`` can read annotations in these formats:

* BED_: The color is specified by the RGB field if present, else pale
  red if the strand is "+", pale blue if "-", or pale purple.  BED
  lines with higher score are drawn on top of ones with lower score.

* Repeatmasker_ .out, rmsk.txt: The color is pale purple for "low
  complexity", "simple repeats", and "satellites", else pale red for
  "+" strand and pale blue for "-" strand.

* genePred_, GFF/GTF: Exons are shown in green, with a darker shade
  for protein-coding regions.

* AGP_, gap.txt: Unsequenced gaps are shown, but only if the gap
  covers at least one whole pixel.

You can use these options multiple times, e.g.
``-a stuff.bed -a more.bed -a rmsk.txt``.  Annotations look good only
if reasonably sparse, e.g. you can't sensibly view 20000 gene
annotations in one small dotplot.

Choosing sequences
------------------

For example, you can exclude sequences with names like
"chrUn_random522" like this::

  last-dotplot -1 'chr[!U]*' -2 'chr[!U]*' alns alns.png

Option "-1" selects sequences from the 1st (horizontal) genome, and
"-2" selects sequences from the 2nd (vertical) genome.  'chr[!U]*' is
a *pattern* that specifies names starting with "chr", followed by any
character except U, followed by anything.

==========  =============================
Pattern     Meaning
==========  =============================
``*``       zero or more of any character
``?``       any single character
``[abc]``   any character in abc
``[!abc]``  any character not in abc
==========  =============================

If a sequence name has a dot (e.g. "hg19.chr7"), the pattern is
compared to both the whole name and the part after the dot.

You can specify more than one pattern, e.g. this gets sequences with
names starting in "chr" followed by one or two characters::

  last-dotplot -1 'chr?' -1 'chr??' alns alns.png

You can also specify a sequence range; for example this gets the first
1000 bases of chr9::

  last-dotplot -1 chr9:0-1000 alns alns.png

Text font
---------

You can improve the font quality by increasing its size, e.g. to 20
points::

  last-dotplot -s20 my-alignments my-plot.png

last-dotplot tries to find a nice font on your computer, but may fail
and use an ugly font.  You can specify a font like this::

  last-dotplot -f /usr/share/fonts/liberation/LiberationSans-Regular.ttf alns alns.png

Colors
------

Colors can be specified in `various ways described here
<https://pillow.readthedocs.io/en/stable/reference/ImageColor.html>`_.

.. _Python Imaging Library: https://pillow.readthedocs.io/
.. _MAF: https://genome.ucsc.edu/FAQ/FAQformat.html#format5
.. _BED: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _genePred: https://genome.ucsc.edu/FAQ/FAQformat.html#format9
.. _RepeatMasker: http://www.repeatmasker.org/
.. _AGP: https://www.ncbi.nlm.nih.gov/assembly/agp/
