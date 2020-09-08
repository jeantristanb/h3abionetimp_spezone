#scl enable devtoolset-8 bash
#wget -c https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
#tar -xvjf htslib-1.10.2.tar.bz2
#DirInstall=$PWD/htslib/
#mkdir -p $DirInstall
#cd htslib-1.10.2
#./configure --prefix=$DirInstall
#make 
#make install
cd ../
#wget -c https://dl.bintray.com/boostorg/release/1.74.0/source/boost_1_74_0.tar.bz2
#tar -xvjf boost_1_74_0.tar.bz2
#cd boost_1_74_0.tar
#./bootstrap.sh
#./b2
cd ../
wget -c https://curl.haxx.se/download/curl-7.72.0.tar.bz2
tar -xvjf curl-7.72.0.tar.bz2
DirCurl=$PWD/curl/
mkdir -p $DirCurl
cd curl-7.72.0/
./configure --prefix=$DirCurl
make 
make install
cd ..
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"$PWD/curl/lib/":"$PWD/htslib/lib"


