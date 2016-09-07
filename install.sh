#!/bin/bash
log="install.log"
echo `date` "Installing dependencies..." 
echo " === You can find log in: $log === "

echo `date` " BWA"
(cd bin/bwa && make clean && make -j 4) >> $log 2>&1 

echo `date` " LASTal"
(cd bin/last && make clean && make -j 4) >> $log 2>&1

echo `date` "Done!"
