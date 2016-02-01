#!/bin/bash
###
# Redundans installer for UNIX.
###

# sudo apt-get install curl build-essential libbz2-dev libsqlite3-dev zlib1g-dev libxml2-dev libxslt1-dev libreadline5 libgdbm-dev  libxml2 libssl-dev tk-dev libgdbm-dev libexpat1-dev libncursesw5-dev
# zlib.h


log="/tmp/install.log"
installdir="$HOME/src"
pyversion="2.7.10"

exists()
{
  command -v "$1" >/dev/null 2>&1
}

# check if all programs exists
for cmd in echo wget curl gcc make cd ln date bash; do
    if ! exists $cmd; then
        echo "Install $cmd first!"
        exit 1
    fi
done

echo -e "\n#####################################################################"
echo -e "#                        Redundans installer                        #"
echo -e "#####################################################################\n"
echo "Redundans and its dependencies will be installed in $installdir"
echo "Python $pyversion and all necessary dependencies will be installed in ~/.pythonbrew"
echo -e " Necessary imports will be added to ~/.bashrc automatically\n"
echo "Installation may take several minutes! Installation log can be found in $log."
echo -e '\n! Make sure libssl & zlib.h are installed ie. `sudo apt-get install zlib1g-dev libssl-dev`!\n'

# YES/NO prompt
echo -n "Do you want to proceed with installation (y/n)? "
read answer
if echo "$answer" | grep -viq "^y" ; then
    echo "Aborted!"
    exit 0
fi

echo ""

# make copy of .bashrc
cp ~/.bashrc ~/.bashrc_bak

if [ ! -d $installdir ]; then mkdir -p $installdir; fi
cd $installdir

echo `date` "Installing Python & dependencies..."
# install pythonbrew to ~/.pythonbrew
curl -kLs http://xrl.us/pythonbrewinstall | bash >> $log 2>&1
 
# add to ~/.bashrc to automatically activate pythonbrew
echo -e "\n###\n# redundans imports" >> ~/.bashrc
echo "# python brew activation" >> ~/.bashrc
echo '[[ -s "$HOME/.pythonbrew/etc/bashrc" ]] && source "$HOME/.pythonbrew/etc/bashrc"' >> ~/.bashrc
echo 'export PATH=$PATH:'$installdir/SSPACE:$installdir/bwa:$installdir/last/src:$installdir >> ~/.bashrc
echo "###" >> ~/.bashrc

# export PATH 
# source ~/.bashrc # not working as not interactive shell
source "$HOME/.pythonbrew/etc/bashrc"
export PATH=$PATH:$installdir/SSPACE:$installdir/bwa:$installdir/last/src:$installdir

# install python 
pythonbrew install $pyversion >> $log 2>&1
 
# and enable the new version
pythonbrew switch $pyversion >> $log 2>&1

# biopython, numpy, scpy
pip install -U biopython numpy scipy >> $log 2>&1


echo `date` "Installing redundans dependencies..."
echo `date` " BWA"
# BWA
wget -q http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.12.tar.bz2
tar xpfj bwa-0.7.12.tar.bz2
ln -s bwa-0.7.12 bwa
cd bwa 
make >> $log 2>&1
cd ..

# LAST
echo `date` " LAST"
wget -q http://last.cbrc.jp/last-714.zip
unzip -q last-714.zip
ln -s last-714 last
cd last
make >> $log 2>&1
cd ..

# SSPACE - note tar.gz and dir are different!
echo `date` " SSPACE"
wget -q http://www.baseclear.com/base/download/41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
tar xpfz 41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
ln -s SSPACE-STANDARD-3.0_linux-x86_64 SSPACE

# GapCloser
echo `date` " GapCloser"
wget -q http://downloads.sourceforge.net/project/soapdenovo2/GapCloser/bin/r6/GapCloser-bin-v1.12-r6.tgz
tar xpfz GapCloser-bin-v1.12-r6.tgz

# BLAT
echo `date` " BLAT"
wget -q https://raw.githubusercontent.com/lpryszcz/bin/master/blat
chmod +x blat

# check if installed corretly
for cmd in blat lastal bwa GapCloser SSPACE_Standard_v3.0.pl; do
    if ! exists $cmd; then
        echo "[WARNING] Make sure $cmd installed properly!"
    fi
done


echo `date` " Redundans"
#git clone https://github.com/lpryszcz/redundans
wget -O redundans.tgz https://github.com/lpryszcz/redundans/archive/master.tar.gz
tar xpfz redundans.tgz
mv redundans-master redundans
cd redundans


echo `date` "Trying Redundans..."
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1


echo -e "#####\nMake sure to register as user of:"
echo "- SSPACE http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE"
echo -e "#####\n"

echo "To uninstall execute:"
echo " rm -rI ~/.pythonbrew ~/src/{*SSPACE,bwa,blat,GapCloser,last,redundans}*"
echo " cp ~/.bashrc_bak ~/.bashrc"

exit 0
