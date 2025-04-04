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

#echo `date` " Platanus"
#(wget -O- http://platanus.bio.titech.ac.jp/?ddownload=145 > bin/platanus && chmod +x bin/platanus) >> $log 2>&1 
#retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " Miniasm"
(cd bin/miniasm && make clean && make) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " BWA"
(cd bin/bwa && make clean && make -j $cores) >> $log 2>&1 
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " LASTal"
(cd bin/last && make clean && make -j $cores) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " SNAP-aligner"

#Had issues while compiling SNAP-aligner so download directly the binary from git, keep track of the file
(cd bin && wget -nc https://github.com/amplab/snap/releases/download/v2.0.1/snap-aligner && chmod +x snap-aligner)  >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " Minimap2"
(cd bin/minimap2 && make clean && make) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

#echo `date` " k8"
#(cd bin && curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf - && mv -t . ./k8-0.2.4/k8-Linux && rm -r ./k8-0.2.4/) >> $log 2>&1
#retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " Meryl"
(cd bin && wget -nc https://github.com/marbl/meryl/releases/download/v1.3/meryl-1.3.Linux-amd64.tar.xz && tar -xJf meryl-1.3.Linux-amd64.tar.xz && mv -t . ./meryl-1.3/bin/*  && rm -r ./meryl-1.3/ && rm meryl-1.3.Linux-amd64.tar.xz) >> $log 2>&1
retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` " gfastats"
(cd bin/gfastats && make -j) >> $log 2>&1

retcode=$?; if [ $retcode -gt 0 ]; then exit $retcode; fi

echo `date` "Done!"

