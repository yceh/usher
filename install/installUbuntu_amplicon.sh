set -e
sudo -E apt update 
sudo -E apt-get install -yq --no-install-recommends build-essential \
 wget cmake  libboost-filesystem-dev libboost-program-options-dev libboost-iostreams-dev libboost-date-time-dev \
 libprotoc-dev libprotoc-dev protobuf-compiler \
 mafft rsync libtbb-dev openmpi-bin libopenmpi-dev automake libtool autoconf make nasm ninja-build

# create build directory
startDir=$PWD
cd $(dirname "$0")
mkdir -p ../build
cd ../build
wget https://github.com/intel/isa-l/archive/refs/tags/v2.30.0.tar.gz
tar -xvf v2.30.0.tar.gz
cd isa-l-2.30.0
./autogen.sh
./configure
make -j$(nproc)
sudo -E make install
cd ..

#download and install TBB
wget https://github.com/oneapi-src/oneTBB/archive/2019_U9.tar.gz 
tar -xvzf 2019_U9.tar.gz

# install amplicon placement specific tools
# install minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf -
cp minimap2-2.24_x64-linux/{minimap2,k8,paftools.js} .  # copy executables
#export PATH="$PATH:"$PWD                               # put the current directory on PATH
chmod +x minimap2

# install faToVcf
rsync -aP rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faToVcf .
chmod +x faToVcf

# install samtools
#wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
#tar -vxjf samtools-1.9.tar.bz2
#cd samtools-1.9
#make
#mv samtools-1.9 samtools
#cd ..

# Debug mode
cmake -GNinja -DTBB_DIR=${PWD}/oneTBB-2019_U9  -DCMAKE_PREFIX_PATH=${PWD}/oneTBB-2019_U9/cmake -DCMAKE_BUILD_TYPE=Debug ..
#make -j$(nproc) VERBOSE=1
ninja

cd $startDir
export PATH=$PATH:$PWD/build/

