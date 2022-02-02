## Import Example
This example demonstrates how to incorporate the rsLQR package into your 
CMake build using the `FetchContent` module. This clones the rsLQR 
repository directly into your `build` folder and adds it to your build. 
This is useful since it doesn't require the end user to download and build 
the rsLQR separately from the library that uses it.

### Build Steps
```shell
mkdir build
cd build
cmake ..  # can pass any of the rsLQR build options here, e.g. -D RSLQR_LINALG_LIBRARY=Eigen
cmake --build .
```