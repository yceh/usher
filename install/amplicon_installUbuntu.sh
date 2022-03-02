# create build directory
startDir=$PWD
cd $(dirname "$0")
mkdir -p build
cd build
buildDir=$PWD

# install minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf -
cp minimap2-2.24_x64-linux/{minimap2,k8,paftools.js} .
export PATH="$PATH:"$PWD                             
chmod +x minimap2

# install fastqToFa
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/fastqToFa .
chmod +x fastqToFa

# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
chmod +x faToVcf

# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bamToFastq .
chmod +x bamToFastq

# install bedtools
#wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.30.0.tar.gz
#tar -zxvf bedtools-2.30.0.tar.gz
#cd bedtools2
#make
#cd $buildDir

# install samtools
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar -vxjf samtools-1.9.tar.bz2
cd samtools-1.9
make
cd $buildDir

# install bcftools
bcftools_VERSION="1.15"
wget https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.15.tar.bz2 -O bcftools.tar.bz2
tar -xjvf bcftools.tar.bz2
cd bcftools-"${bcftools_VERSION}"
make

# install bw2-mem2
#git clone https://github.com/bwa-mem2/bwa-mem2
#cd bwa-mem2
#git submodule init
#git submodule update
#git clone --recursive https://github.com/bwa-mem2/bwa-mem2
#cd bwa-mem2
#make

cd $buildDir
export PATH=$PATH:$PWD
cd $startDir
