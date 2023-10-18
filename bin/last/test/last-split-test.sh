#! /bin/sh

# Regression tests

cd $(dirname $0)

PATH=../bin:$PATH

maf=SRR359290-1k.maf

{
    last-split -h

    last-split -m0.01 -fMAF+ $maf

    last-split -m0.01 -c0 $maf

    last-split -m0.01 -t0 $maf

    last-split -m0.01 -M7.5 -S2 $maf

    last-split -m0.001 -s180 $maf

    last-split -m0.01 -n -fMAF+ $maf

    last-split -d0 -m0.001 -s180 -fMAF+ $maf
    last-split -d1 -m0.001 -s180 -fMAF+ $maf
    last-split -d2 -m0.001 -s180 -fMAF+ $maf

    grep -v '^q' $maf | last-split -m0.001 -s180 -fMAF+

    last-split -fMAF+ 102.maf
    last-split -fMAF 102.maf

    last-split -d1 -fMAF+ spliceWithGap.maf

    last-split -r split1.maf

} | diff -u last-split-test.out -
