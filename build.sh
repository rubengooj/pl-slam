echo "Building 3rdparty/line_descriptor ... "
cd 3rdparty/line_descriptor
mkdir build
cd build
cmake ..
make -j8
cd ../../../

echo "Building 3rdparty/DBoW2 ... "
cd 3rdparty/DBoW2
mkdir build
cd build
cmake ..
make -j8
cd ../../../

echo "Uncompressing vocabulary ..."
cd vocabulary
tar -xf voc.tar.gz
cd ..

echo "Building PL-SLAM ... "
mkdir build
cd build
cmake ..
make -j8
