#!/bin/bash
###
# Redundans installer for UNIX.
###

# sudo apt-get install wget curl screen make git gcc make bash libc-dev

log="/tmp/install.log"
installdir="$HOME/src"
pyversion="2.7.11"
waiting="30s"

exists()
{
  command -v "$1" >/dev/null 2>&1
}

# check if all programs exists
for cmd in echo wget curl git gcc make cd ln date bash; do
    if ! exists $cmd; then
        echo "Install $cmd first!"
        exit 1
    fi
done

echo -e "\n#######################"
echo -e "# Redundans installer #"
echo -e "#######################\n"
echo "Installation may take several minutes! Installation log can be found in $log."
echo "Redundans and its dependencies will be installed in $installdir"
echo "Python $pyversion and all necessary dependencies will be installed in ~/.pythonbrew"
echo -e " Necessary imports will be added to ~/.bashrc automatically\n"

# YES/NO prompt
echo -n "Do you want to proceed with installation (y/n)?"
read answer
if echo "$answer" | grep -iq "^y" ; then
    echo -n ""
else
    echo "Aborted!"
    exit 0
fi

if [ ! -d $installdir ]; then mkdir -p $installdir; fi
cd $installdir

echo `date` "Installing Python & dependencies..."
# install pythonbrew to ~/.pythonbrew
curl -kLs http://xrl.us/pythonbrewinstall | bash >> $log 2>&1
 
# add to ~/.bashrc to automatically activate pythonbrew
echo -e "###\n# redundans imports" >> ~/.bashrc
echo "# python brew activation" >> ~/.bashrc
echo '[[ -s "$HOME/.pythonbrew/etc/bashrc" ]] && source "$HOME/.pythonbrew/etc/bashrc"' >> ~/.bashrc
echo 'export PATH=$PATH:'$installdir/SSPACE:$installdir/bwa:$installdir/last/src:$installdir >> ~/.bashrc
echo "###" >> ~/.bashrc
# reload 
source ~/.bashrc
 
# install python 
pythonbrew install $pyversion >> $log 2>&1
 
# and enable the new version
pythonbrew switch $pyversion >> $log 2>&1

# install pip
#wget -q wget https://bootstrap.pypa.io/get-pip.py
#python get-pip.py
#rm get-pip.py

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

# SSPACE 
echo `date` " SSPACE"
wget -q http://www.baseclear.com/base/download/41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
tar xpfz 41SSPACE-STANDARD-3.0_linux-x86_64.tar.gz
ln -s 41SSPACE-STANDARD-3.0_linux-x86_64 SSPACE

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
git clone git@github.com:lpryszcz/redundans.git
cd redundans

echo `date` "Trying Redundans..."
./redundans.py -v -i test/*.fq.gz -f test/contigs.fa -o test/run1


echo -e "#####\nMake sure to register as user of:"
echo "- SSPACE http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE"
echo -e "#####\n"

exit 0
