LAST and repeats
================

Biomolecule sequences are full of repeats, which have a massive effect
on sequence comparison.  There are two different things called
"repeats":

* Interspersed repeats, such as LINEs, Alus, transposons.

* Simple sequence / low-complexity sequence / tandem repeats, such as
  atatatatatatatatat.

They cause two different problems:

* Simple tracts cause alignments of sequences that are not homologous
  (not descended from a common ancestor), because similar simple
  sequences can evolve independently.

* Both kinds of repeat can make sequence comparison very slow, and
  produce an overwhelming number of alignments.

LAST options
------------

LAST has built-in recognition of simple sequences, but not (yet)
interspersed repeats.  So you may wish to annotate interspersed
repeats with a tool such as RepeatMasker or WindowMasker.  You can
often obtain pre-masked genomes, with repeats indicated by lowercase.

LAST's built-in simple-sequence masking uses `tantan
<http://cbrc3.cbrc.jp/~martin/tantan/>`_, which prevents
non-homologous alignments more reliably than previous methods.

To mask all kinds of repeat, use lastdb -cR11:

  This keeps lowercase in the input sequences, additionally lowercases
  simple tracts found by tantan, and masks lowercase during alignment.

To mask simple tracts only, use lastdb -cR01:

  This does not keep lowercase in the input sequences, lowercases
  simple tracts found by tantan, and masks lowercase during alignment.

In LAST, "mask" means that lowercase is excluded when finding similar
sequences, but included when making the final alignments.

Orthology search and post-masking
---------------------------------

Sometimes we wish to find orthologs but not paralogs: e.g. if we
compare two whole genomes, or align DNA reads to a genome.  Masking
can cause false positives, if an ortholog is masked but a paralog is
not.

We'd like to avoid this problem, but also avoid non-homologous
alignments of simple sequences.  This can be done by post-masking:

1. Align using lastdb -R01: this does not keep lowercase in the input
   sequences, lowercases simple tracts found by tantan, but does not
   mask it (i.e. case has no effect).

2. Predict orthologs, e.g. using `last-split <last-split.html>`_.

3. Discard mostly-lowercase alignments, using `last-postmask
   <last-postmask.html>`_.
