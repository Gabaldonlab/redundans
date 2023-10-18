#! /bin/sh

try () {
    echo TEST "$@"
    eval "$@"
    echo
}

cd $(dirname $0)

PATH=../bin:$PATH

{
    try last-postmask 102.maf
    try last-postmask 90089.maf

    try "sed 's/humdb/h/' ../examples/myalns.maf | last-postmask"

} 2>&1 |
diff -u $(basename $0 .sh).out -
