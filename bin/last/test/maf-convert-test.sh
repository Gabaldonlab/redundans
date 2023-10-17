#! /bin/sh

# Regression tests

cd $(dirname $0)

PATH=../bin:$PATH

r=maf-convert

maf1=SRR359290-1k.maf
maf2=bs100.maf

{
    $r -h
    head -n999 $maf1 | $r axt
    head -n999 $maf1 | $r blast
    head -n999 $maf1 | $r -l100 blast
    $r blast $maf2
    $r -l120 blast frameshift-new.maf
    tail -n8 frameshift-new.maf | maf-swap | $r blast
    $r blasttab $maf2
    $r blasttab frameshift-new.maf
    head -n999 $maf1 | $r chain
    $r gff 102.maf
    $r -J1e9 gff 102.maf
    $r gff frameshift-new.maf
    $r html -l100 $maf2
    head -n999 $maf1 | $r -n html
    head -n999 $maf1 | $r psl
    head -n999 $maf1 | $r -p psl
    $r psl $maf2
    $r -j1e9 psl 102.maf
    $r -J1e9 psl 102.maf
    $r psl 90089.maf
    $r psl frameshift-old.maf
    $r psl frameshift-new.maf
    $r -n sam $maf2
    head -n999 $maf1 | $r -r 'ID:1 PL:ILLUMINA SM:x' sam
    $r -d sam $maf1
    head -n999 $maf1 | $r -n tab
    head -n999 $maf1 | $r tab
    $r -n tab frameshift-new.maf
} | diff -u maf-convert-test.out -
