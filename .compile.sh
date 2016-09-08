#!/bin/bash
# Compile redundans dependencies

log="/tmp/compiling.log"
if [ ! -z $1 ]; then log=$1; fi

echo `date` "Compiling dependencies..." 
echo " === You can find log in: $log === "
cores=`grep -c ^processor /proc/cpuinfo | awk '{if($1>1){print $1-1} else {print 1}}'`
echo " I'll use $cores threads for compiling"

echo `date` " BWA"
(cd bin/bwa && make clean && make -j $cores) >> $log 2>&1 

echo `date` " LASTal"
(cd bin/last && make clean && make -j $cores) >> $log 2>&1

echo `date` "Done!"
exit 0
