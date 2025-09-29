#!/bin/bash

# if [ -d "build" ]; then
#     rm -rf build
# fi
# if [ -d "install" ]; then
#     rm -rf install
# fi
# mkdir build
if [ ! -d "build" ]; then
    mkdir build
fi
cd build
cmake .. \
  -DCMAKE_BUILD_TYPE=Debug \
  -DCMAKE_EXPORT_COMPILE_COMMANDS=1 \
  -DCMAKE_INSTALL_PREFIX=./install \
  -DNAVTOOLS_BUILD_TESTS=ON \

cmake --build . --parallel 8
cmake --install .
cd ..
