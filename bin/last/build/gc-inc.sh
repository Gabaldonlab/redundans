#! /bin/sh

# This generates source code from genetic codes.

cat <<EOF
const struct {
  const char *name;
  const char *text;
} geneticCodes[] = {
EOF

cat "$@" | tr -d '",' |
awk '
$1 == "id" {print "{\"" $2 "\", \"\\"}
$1 == "ncbieaa" {print "  AAs = " $2 "\\n\\"}
/-- Base/ {print $2 " = " $3 "\\n\\"}
/-- Base3/ {print "\"},"}
'

echo "};"
