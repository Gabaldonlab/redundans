#! /bin/sh

# This generates source code from subset seeds.

cat <<EOF
const struct {
  const char *nickname;
  const char *realname;
} subsetSeedNicknames [] = {
EOF
for i in "$@"
do
    name=$(basename $i .seed)
    grep '^#abbreviation' $i | cut -d' ' -f2 | sed 's/.*/{"&", "'$name'"},/'
done
echo "};"

echo

cat <<EOF
const struct {
  const char *name;
  const char *text;
} subsetSeeds[] = {
EOF
for i in "$@"
do
    basename $i .seed | sed 's/.*/{"&", "\\/'
    grep -v '^#[a ]' $i | awk NF | sed 's/$/\\n\\/'
    echo '"},'
    echo
done
echo "};"
