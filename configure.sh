#!/bin/bash

cd spoa
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -Dspoa_build_executable=ON ..
make
cd ..
if [ -f ../consensus.cpp ]
then
    mv ../consensus.cpp .
fi
g++ consensus.cpp -std=c++11 -Iinclude/ -Lbuild/lib/ -lspoa -o consensus -static
cp consensus ../TRiCoLOR/REFER/
