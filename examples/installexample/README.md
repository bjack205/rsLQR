## Install Example
This example shows how to link against the rsLQR library after it is 
installed on your computer. See the installation instructions in the 
documentation for more details on how to install it. If the library 
was installed in a standard location, you shouldn't have to do anything
special for CMake to find it. However, if you installed it in some 
custom location by specifying `CMAKE_INSTALL_PREFIX`, you'll need to 
point CMake to that location. 

### Build Steps
```shell
mkdir build
cd build
cmake -D CMAKE_PREFIX_PATH=<path to install location> ..
cmake --build .
```