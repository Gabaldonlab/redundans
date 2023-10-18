#! /bin/sh

# This generates source code from substitution score matrices.

cat <<EOF
const struct {
  const char *nickname;
  const char *realname;
} scoreMatrixNicknames [] = {
EOF
for i in "$@"
do
    name=$(basename $i .mat)
    grep '^#nickname' $i | cut -d' ' -f2 | sed 's/.*/{"&", "'$name'"},/'
done
echo "};"

echo

cat <<EOF
const struct {
  const char *name;
  const char *text;
} scoreMatrices[] = {
EOF
for i in "$@"
do
    basename $i .mat | sed 's/.*/{"&", "\\/'
    grep -v '^#[n ]' $i | awk NF | sed 's/$/\\n\\/'
    echo '"},'
    echo
done
echo "};"
