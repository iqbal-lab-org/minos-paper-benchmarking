#!/usr/bin/env bash
set -vexu

install_root=$1

apt-get install -y software-properties-common
apt-add-repository universe
apt-get update

apt-get install -y \
  automake \
  bioperl \
  build-essential \
  cmake \
  curl \
  gcc \
  gdb \
  git \
  openjdk-8-jre \
  libbz2-dev \
  libcurl4-openssl-dev \
  libhts-dev \
  liblzma-dev \
  libncurses5-dev \
  libncursesw5-dev \
  libssl-dev \
  libtabixpp-dev \
  pkg-config \
  python-dev \
  python3 \
  python3-dev \
  python3-pip \
  python3-setuptools \
  r-base-core \
  rsync \
  tabix \
  time \
  wget \
  zlib1g-dev


# bioperl installs mummer v3. We install v4 later, so remove here to avoid
# conflicts. Also, v3 gets put in /usr/bin, which is before where v4 is put
# /usr/local/bin in $PATH.
apt-get remove -y mummer

if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root

# At time of writing, pysam 0.16.0 is latest, but throws error:
# >>> import pysam
# Traceback (most recent call last):
#   File "<stdin>", line 1, in <module>
#   File "/usr/local/lib/python3.6/dist-packages/pysam/__init__.py", line 5, in <module>
#     from pysam.libchtslib import *
#   File "pysam/libchtslib.pyx", line 1, in init pysam.libchtslib
# ImportError: libchtslib.cpython-36m-x86_64-linux-gnu.so: cannot open shared object file: No such file or directory
pip3 install cython
pip3 install pysam==0.15.4



#__________________________ ART _______________________________#
cd $install_root
wget -q https://www.niehs.nih.gov/research/resources/assets/docs/artbinmountrainier20160605linux64tgz.tgz
tar xf artbinmountrainier20160605linux64tgz.tgz
cp -s art_bin_MountRainier/art_illumina .


#________________________ simutator _________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/simutator.git
cd simutator
git checkout 0218c8a5b37fd72eb4e5b2df4cba9f6118f96788
pip3 install .


#_________________________ bcftools _______________________#
cd $install_root
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
tar xf bcftools-1.10.2.tar.bz2
rm bcftools-1.10.2.tar.bz2
cd bcftools-1.10.2/
make
cd ..
cp -s bcftools-1.10.2/bcftools .


#__________________________ BWA____________________________#
cd $install_root
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar xf bwa-0.7.17.tar.bz2
rm bwa-0.7.17.tar.bz2
cd bwa-0.7.17/
make
cd ..
cp -s bwa-0.7.17/bwa .


#________________________ KMC _______________________________#
cd $install_root
wget https://github.com/refresh-bio/KMC/releases/download/v3.1.1/KMC3.1.1.linux.tar.gz
tar xf KMC3.1.1.linux.tar.gz
rm KMC3.1.1.linux.tar.gz


#________________________ mummer ____________________________#
cd $install_root
wget https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
tar xf mummer-4.0.0rc1.tar.gz
rm mummer-4.0.0rc1.tar.gz
cd mummer-4.0.0rc1
./configure
make
make install


#_________________________ samtools ______________________#
cd $install_root
wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
tar xf samtools-1.10.tar.bz2
rm samtools-1.10.tar.bz2
cd samtools-1.10
make
cd ..
cp -s samtools-1.10/samtools .


#________________________ svimmer ___________________________#
cd $install_root
wget https://github.com/DecodeGenetics/svimmer/archive/v0.1.tar.gz
tar xf v0.1.tar.gz
rm v0.1.tar.gz
cp -s svimmer-0.1/svimmer .


#________________________ graphtyper ________________________#
cd $install_root
wget https://github.com/DecodeGenetics/graphtyper/releases/download/v2.7.1/graphtyper
chmod a+x graphtyper


#________________________ bayestyper ________________________#
cd $install_root
wget https://github.com/bioinformatics-centre/BayesTyper/releases/download/v1.5/bayesTyper_v1.5_linux_x86_64.tar.gz
tar xf bayesTyper_v1.5_linux_x86_64.tar.gz
rm bayesTyper_v1.5_linux_x86_64.tar.gz
chmod --recursive a+rx bayesTyper_v1.5_linux_x86_64
cp -s bayesTyper_v1.5_linux_x86_64/bin/* .


#________________________ minimap2 etc ______________________#
cd $install_root
git clone https://github.com/lh3/minimap2.git minimap2_git
cd minimap2_git
git checkout ccb0f7b05df3a17011f0d0f4388ddb301198871b
make
cd ..
cp -s minimap2_git/minimap2 .
curl -L https://github.com/attractivechaos/k8/releases/download/v0.2.4/k8-0.2.4.tar.bz2 | tar -jxf -
cp -s k8-0.2.4/k8-Linux k8
cp -s minimap2_git/misc/paftools.js .


#____________________ cluster_vcf_records _________________#
pip3 install 'cluster_vcf_records==0.13.2'

#________________________ varifier __________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/varifier.git
cd varifier
git checkout ad3d6db3b851eb46357b0308e7496580e191700f
pip3 install .


#________________________ gramtools __________________________#
cd $install_root
# Why six>=1.14.0?
# See https://github.com/pypa/virtualenv/issues/1551
#pip3 install tox "six>=1.14.0"
git clone https://github.com/iqbal-lab-org/gramtools
cd gramtools
git checkout 8af53f6c8c0d72ef95223e89ab82119b717044f2
pip3 install .

#________________________ vt __________________________________#
cd $install_root
git clone https://github.com/atks/vt.git vt-git
cd vt-git
git checkout 2187ff6347086e38f71bd9f8ca622cd7dcfbb40c
make
cd ..
cp -s vt-git/vt .

#______________________vcflib _________________________________#
cd $install_root
git clone --recursive https://github.com/vcflib/vcflib.git
cd vcflib
# Note date is 20210112. There's a bug when running make, when it tries to
# install the man page. We don't care about this, so comment out the line
awk '/^install.*\/man\/man1)/ {$0="#"$0} 1' CMakeLists.txt > z
mv z CMakeLists.txt
make

#________________________ minos ______________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/minos.git
cd minos
git checkout c22f4b00c61be7854336dad5b31b44c0d11e8e81
pip3 install .


#________________________ Trimmomatic ____________________#
cd $install_root
# Website is (permanently?) down. So copy of the zip is now
# included in the repo and when the script runs, that copy
# is assumed to already be in $install_root
#wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip


#_________________________ fqtools ________________________#
cd $install_root
wget https://github.com/alastair-droop/fqtools/archive/986e451.tar.gz
tar xf 986e451.tar.gz
rm 986e451.tar.gz
cd fqtools-986e451/
make
cd ..
cp -s fqtools-986e451/bin/fqtools .


#________________________ stampy _________________________#
cd $install_root
wget http://www.well.ox.ac.uk/~gerton/software/Stampy/stampy-latest.tgz
tar xf stampy-latest.tgz
rm stampy-latest.tgz
cd stampy-*
make
cd ..
cp -s stampy-*/stampy.py .


#________________________ vcftools _______________________#
cd $install_root
wget https://github.com/vcftools/vcftools/releases/download/v0.1.15/vcftools-0.1.15.tar.gz
tar xf vcftools-0.1.15.tar.gz
rm vcftools-0.1.15.tar.gz
cd vcftools-0.1.15
./configure --prefix $PWD/install
make
make install

# cortex needs the perl/ directory. It expects it to be in the vcftools root,
# but somehwere between v0.1.9 and v0.1.15 it moved into src/.
ln -s src/perl/ .


#________________________ cortex _________________________#
cd $install_root
git clone https://github.com/iqbal-lab/cortex.git
cd cortex
# After this commit, cortex changed to use minimap2 instead
# of stampy, but also CLI changed. So pin to this commit,
# otherwise clockwork calls to cortex need changing
git checkout 3a235272e4e0121be64527f01e73f9e066d378d3
bash install.sh
make NUM_COLS=1 cortex_var
make NUM_COLS=2 cortex_var


#________________________ seqtk __________________________#
cd $install_root
wget https://github.com/lh3/seqtk/archive/v1.2.tar.gz
tar xf v1.2.tar.gz
rm v1.2.tar.gz
cd seqtk-1.2/
make
cd ..
cp -s seqtk-1.2/seqtk .


#________________________ clockwork  ______________________________#
cd $install_root
wget https://github.com/iqbal-lab-org/clockwork/archive/v0.9.0.tar.gz
tar xf v0.9.0.tar.gz
rm v0.9.0.tar.gz
cd clockwork-0.9.0/python/
pip3 install .

#______________________ ivcmerge ______________________________#
cd $install_root
git clone https://github.com/iqbal-lab-org/ivcfmerge.git
cd ivcfmerge
git checkout 5819787614a263a9f35fd0c247442f092ab174ff
pip3 install .

