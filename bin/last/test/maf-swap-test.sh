#! /bin/sh

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

PATH=../bin:$PATH

{
    try maf-swap -h
    try maf-swap bs100.maf
    try maf-swap -n1 90089.maf
    try maf-swap -n3 ../examples/multiMito.maf
    try maf-swap frameshift-new.maf
} 2>&1 | diff -u maf-swap-test.out -
