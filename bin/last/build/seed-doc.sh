#! /bin/sh

# This generates text documentation from subset seeds.

cat <<'EOF'
LAST seeding schemes
====================

LAST's critical first step is to find *seeds*, i.e. initial matches
between query and reference sequences.  It can use various seeding
schemes, which allow different kinds of mismatches at different seed
positions.

   A seeding scheme consists of a seed alphabet, such as::

     1  A C G T
     0  ACGT
     T  AG CT

   and one or more patterns, such as this one::

     1T1T10T1101101

   Each symbol in a pattern represents a grouping of sequence letters:
   in this example, ``T`` represents the grouping ``AG CT``.  At each
   position in an initial match, mismatches are allowed between
   letters that are grouped at that position in the pattern.

   Although the patterns have fixed lengths, LAST's initial matches do
   not.  LAST finds shorter matches by using a prefix of the pattern,
   and longer matches by cyclically repeating the pattern.

   A *restricted* symbol omits letters of the main sequence alphabet,
   which are then forbidden at those positions::

     r  AG

   An *exact* symbol groups no letters::

     Y  C T

   In 2nd and subsequent cycles, restricted symbols are made
   unrestricted: if it is exact then the omitted letters are added as
   separate groups, else they are added as one group.

EOF

for i in "$@"
do
    name=$(basename $i .seed)
    abbr=$(grep '^#abbreviation' $i | cut -d' ' -f2)
    test "$abbr" && name="$name (abbreviation: $abbr)"
    echo $name
    echo $name | sed 's/./-/g'  # underline
    echo
    grep '^# ' $i | cut -d' ' -f2-
    echo It uses this seed alphabet::
    echo
    awk '!/^#/ && NF > 1 && length($1) == 1' $i | sed 's/^/  /'
    echo
    if [ $(awk '!/^#/ && length($1) > 1 || NF == 1' $i | wc -w) = 1 ]
	then echo And this pattern::
	else echo And these patterns::
    fi
    echo
    awk '!/^#/ && length($1) > 1 || NF == 1' $i | sed 's/^/  /'
    echo
    grep -q '^#lastdb' $i && {
	echo It sets this lastdb default:
	grep '^#lastdb' $i | cut -d' ' -f2-
	echo
    }
    grep -q '^#lastal' $i && {
	echo It sets this lastal default:
	grep '^#lastal' $i | cut -d' ' -f2-
	echo
    }
done

exit 0
