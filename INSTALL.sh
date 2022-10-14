#!/bin/bash
###
# Redundans installer for UNIX.
# bash <(curl -Ls http://bit.ly/redundans_installer)
# version 2.00
###

log="install.log"
branch="master" 
if [ ! -z $1 ]; then branch=$1; fi

echo "#####################################################################################"
echo "#                                                                                   #"
echo "#                               Redundans installer                                 #"
echo "#                                                                                   #"
echo "#       version 2.00                   l.p.pryszcz+git AT gmail   & dfuentes@bsc.es #"
echo "#####################################################################################"
echo ""
echo "Redundans and its dependencies will be installed in:" `pwd`/redundans
echo ""
echo "Installation will take 2-3 minutes. "
echo "To track the installation status execute in the new terminal:"
echo "  tail -f `pwd`/redundans/$log"
echo ""

# sleep
echo "I'll proceed with installation in 10 seconds... Press Ctrl-C to cancel."
sleep 10s

echo ""
echo `date` "Checking dependencies..."

exists()
{
  command -v "$1" >/dev/null 2>&1
}

error=""
# check if all programs exists
for cmd in echo awk wget tar git gcc g++ make cd ln date ldconfig unzip perl python; do
    if ! exists $cmd; then
        echo "Install $cmd first (ie. 'sudo apt-get install $cmd')!"
        error=1
    fi
done

# check if all libs present #BWA
for lib in libz; do
    if [ -z "$(ldconfig -p | grep $lib.so)" ] ; then
        echo " Missing library $lib !"
        error=1
    fi
done
# check headers #BWA
for lib in zlib.h; do
    if [ ! -s /usr/include/$lib ] && [ ! -s /usr/lib/$lib ]; then
        echo " Missing headers $lib !"
        error=1
    fi
done

# skip if error
if [ ! -z $error ]; then
    echo -e "\nAborted due to missing dependencies (see above)!"
    return 1;
fi

# check python version 
PyVer=`python --version 2>&1 | cut -f2 -d" " | cut -f-1 -d"."`
if [ $PyVer != 3 ]; then 
    echo ""
    echo "[ERROR] Install Python >= 3.0!"
    echo "If you have Python 3 already installed, you can either "
    echo "make an alias before installation and use of Redundans ('alias python=python3' should do)"
    echo "or use Python virtual environment (https://virtualenv.pypa.io)."
    return 1
fi
echo " Everything looks good :) Let's proceed..."

echo `date` "Downloading Redundans..."
git clone -b $branch --recursive https://github.com/lpryszcz/redundans.git >> /dev/null 2>&1 
cd redundans 

# compile dependencies
sh bin/.compile.sh `pwd`/$log
retcode=$?; 
if [ $retcode -gt 0 ]; then
    echo "  [$retcode] ERROR!"
    tail -n 20 $log
    #return $retcode
fi

echo `date` "Installation finished!"
echo ""
echo "To try redundans, execute:"
echo "cd redundans; ./redundans.py -v -i test/*_?.fq.gz -f test/contigs.fa -o test/run1"
echo ""
echo "To uninstall execute:"
echo "rm -rI `pwd`"
echo ""
echo "#####################################################################################"
echo "# Redundans depends on several programs (http://bit.ly/redundans_dependencies)      #"
echo "# Acknowledge their authors, get familiar with licensing and register if necessary. #"
echo "#####################################################################################"
echo ""
