#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
CUR_DIR=$PWD

cd $SCRIPT_DIR/..

if [ ! -d "build" ]; then
    mkdir build
    cd build
    cmake ..
else
    cd build
fi

echo "Building the documentation"
cmake ..
cmake --build . --target doxygen

cd $CUR_DIR
