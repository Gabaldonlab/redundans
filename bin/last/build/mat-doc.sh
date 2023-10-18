#! /bin/sh

# This generates text documentation from substitution score matrices.

cat <<EOF
LAST built-in scoring schemes
=============================

EOF

for i in "$@"
do
    name=$(basename $i .mat)
    nick=$(grep '^#nickname' $i | cut -d' ' -f2)
    test "$nick" && name="$name or $nick"
    echo $name
    echo $name | sed 's/./-/g'  # underline
    echo
    grep '^# ' $i | cut -d' ' -f2-
    echo It uses this matrix::
    echo
    grep -v '^#' $i | awk NF | sed 's/^/  /'
    echo
    grep -q '^#last' $i && {
	echo It sets these default lastal parameter values:
	grep '^#last' $i | cut -d' ' -f2-
	echo
    }
done
