last-postmask
=============

This script does post-alignment masking.  It reads pair-wise sequence
alignments, and only writes alignments that include a significant
amount of uppercase-to-uppercase alignment.  (Lowercase is often used
to indicate undesirable "simple" sequence.)

The input should be in MAF format, with header lines (of the kind
produced by lastal) describing the alignment score parameters.

The script discards alignments that lack any segment with score >=
threshold, when applying "gentle masking" of lowercase letters: this
means that the score for such a letter is min(unmasked score, 0).

The usage is simple::

  last-postmask in.maf > out.maf

It may be useful to feed the input from a pipe::

  ... | last-postmask > out.maf

Gentle masking is described in:

  Gentle masking of low-complexity sequences improves homology search.
  Frith MC 2011 PLoS One 6:e28819.

Limitations
-----------

last-postmask does not (yet) handle sequence quality data,
frameshifts, or generalized affine gaps.
