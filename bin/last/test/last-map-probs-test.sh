#! /bin/sh

cd $(dirname $0)

PATH=../bin:$PATH

{
    last-map-probs -s400 SRR359290-1k.maf
    maf-convert tab SRR359290-1k.maf | last-map-probs -m0.1
} | diff -u last-map-probs-test.out -
