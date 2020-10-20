brew install cmake boost protobuf wget python@3

# TBB
wget https://github.com/oneapi-src/oneTBB/releases/download/2019_U9/tbb2019_20191006oss_mac.tgz
tar -xvzf tbb2019_20191006oss_mac.tgz

# Build UShER
mkdir -p build
cd build
cmake -DTBB_DIR=${PWD}/../tbb2019_20191006oss -DCMAKE_PREFIX_PATH=${PWD}/../tbb2019_20191006oss/cmake ..
make -j
cd ..

# install faToVcf
wget https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/faToVcf
chmod 777 ./faToVcf
mv ./faToVcf ./scripts/faToVcf

# install mafft
if ! command -v mafft &> /dev/null; then 
wget https://mafft.cbrc.jp/alignment/software/mafft-7.471-mac.zip
unzip mafft-7.471-mac.zip
cd mafft-mac/
mv mafft.bat /usr/local/bin/mafft; mv mafftdir /usr/local/bin/
cd ..
rm -rf mafft-mac
fi
