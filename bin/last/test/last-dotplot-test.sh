#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

png=/tmp/$$.png

trap 'rm -f "$png"' EXIT

plotSum () {
    echo TEST "$@"
    last-dotplot -s0 "$@" "$png"
    cksum < "$png"
    echo
}

{
    plotSum 102.maf
    plotSum -v -x600 -c orange -r teal --sort1=2 --border-color=brown 102.maf
    plotSum -v -y30 --max-gap1=0 --pad=0.08 --border-pixels=2 --sort1=3 102.maf
    plotSum -v --strands1=1 102.maf
    plotSum -v --sort2=3 bs100.maf
    plotSum -v -2'?' --sort2=0 --strands2=1 bs100.maf
} 2>&1 |
grep -v font |
diff -u $(basename $0 .sh).out -
