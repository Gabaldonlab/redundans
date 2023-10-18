#! /bin/sh

# Regression tests

cd $(dirname $0)

PATH=../bin:$PATH

tmp=${TMPDIR-/tmp}/$$
trap 'rm -f $tmp.*' EXIT

lastdb -m1111110 $tmp.x hg19-M.fa

lastal -Q1 -e120 -i1 $tmp.x bs1.fastq > $tmp.maf1
lastal -Q1 -e120 -i1 $tmp.x bs2.fastq > $tmp.maf2

maf-convert tab $tmp.maf1 > $tmp.tab1
maf-convert tab $tmp.maf2 > $tmp.tab2

{
    last-pair-probs -h

    last-pair-probs $tmp.maf1 $tmp.maf2

    last-pair-probs -m0.001 $tmp.tab1 $tmp.tab2

    last-pair-probs -r $tmp.tab1 $tmp.tab2

    last-pair-probs -f248 -s40.7715 $tmp.tab1 $tmp.tab2

    last-pair-probs -d0.001 $tmp.tab1 $tmp.tab2

    last-pair-probs -cX $tmp.tab1 $tmp.tab2

    last-pair-probs -c. $tmp.tab1 $tmp.tab2

    sed 's:/1::' $tmp.maf1 > $tmp.m1
    sed 's:/2::' $tmp.maf2 > $tmp.m2
    last-pair-probs $tmp.m1 $tmp.m2

    sed 's:/1::' $tmp.tab1 > $tmp.t1
    sed 's:/2::' $tmp.tab2 > $tmp.t2
    last-pair-probs $tmp.t1 $tmp.t2
} | diff -u last-pair-test.out -
