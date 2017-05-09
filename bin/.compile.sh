#!/bin/bash
# Compile redundans dependencies

log="/tmp/compiling.log"
if [ ! -z $1 ]; then log=$1; fi

echo `date` "Updating submodules..."
git submodule update --init --recursive >> $log 2>&1
#git submodule update --recursive  >> $log 2>&1
git submodule foreach git pull origin master >> $log 2>&1 # compatible with git < v1.6.1

echo `date` "Compiling dependencies..." 
echo " === You can find log in: $log === "
cores=`grep -c ^processor /proc/cpuinfo | awk '{if($1>1){print $1-1} else {print 1}}'`
echo " I'll use $cores thread(s) for compiling"

echo `date` " GNU parallel"
#wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2 && tar xpfj parallel-latest.tar.bz2 && rm parallel-latest.tar.bz2 && mv parallel bin/parallel
(cd bin/parallel && ./configure && make clean && make -j $cores)  >> $log 2>&1 
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " BWA"
(cd bin/bwa && make clean && make -j $cores) >> $log 2>&1 
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " LASTal"
(cd bin/last && make clean && make -j $cores) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " SNAP-aligner"
(cd bin/snap && make clean && make -j $cores) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` "Done!"

