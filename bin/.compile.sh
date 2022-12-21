#!/bin/bash
# Compile redundans dependencies

log=`mktemp -t redundans.compiling.XXXXXXXXX.log`
if [ ! -z $1 ]; then log=$1; fi

echo `date` "Updating submodules..."
git submodule update --init --recursive >> $log 2>&1
#git submodule update --recursive  >> $log 2>&1
git submodule foreach git pull origin master >> $log 2>&1 # compatible with git < v1.6.1

echo `date` "Compiling dependencies..." 
echo " === You can find log in: $log === "
cores=`grep -c ^processor /proc/cpuinfo | awk '{if($1>1){print $1-1} else {print 1}}'`
echo " I'll use $cores thread(s) for compiling"

#echo `date` " GNU parallel"
#wget http://ftp.gnu.org/gnu/parallel/parallel-latest.tar.bz2 && tar xpfj parallel-latest.tar.bz2 && rm parallel-latest.tar.bz2 && mv parallel bin/parallel
#(cd bin/parallel && ./configure && make clean && make -j $cores)  >> $log 2>&1 
#retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

#echo `date` " idba"
#(cd bin/idba && ./build.sh && ./configure && make -j $cores) >> $log 2>&1 
#retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

#echo `date` " SPAdes"
#(cd bin/ && wget http://cab.spbu.ru/files/release3.11.1/SPAdes-3.11.1-Linux.tar.gz && tar xpfz SPAdes-3.11.1-Linux.tar.gz && rm SPAdes-3.11.1-Linux.tar.gz && ln -s SPAdes-3.11.1-Linux SPAdes) >> $log 2>&1 
#retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " Platanus"
(wget -O- http://platanus.bio.titech.ac.jp/?ddownload=145 > bin/platanus && chmod +x bin/platanus) >> $log 2>&1 
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " BWA"
(cd bin/bwa && make clean && make -j $cores) >> $log 2>&1 
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " LASTal"
(cd bin/last && make clean && make -j $cores) >> $log 2>&1
#(cd bin && wget http://last.cbrc.jp/last/index.cgi/archive/tip.tar.gz && tar xfpz tip.tar.gz && ln -s last-* last && \
#rm tip.tar.gz && cd last && make -j $cores) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " SNAP-aligner"
(cd bin && wget -nc https://github.com/amplab/snap/releases/download/v2.0.1/snap-aligner && chmod +x snap-aligner)
#(cd bin/snap && make clean && make -j $cores) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` "Done!"

