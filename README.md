# Scalable Mixed-Mode MPC

# Installation

### Dependencies

1. `wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py`
2. `python install.py -install -tool -ot`
    1. By default it will build for Release. `-DCMAKE_BUILD_TYPE=[Release|Debug]` option is also available.
    2. No sudo? Change [`CMAKE_INSTALL_PREFIX`](https://cmake.org/cmake/help/v2.8.8/cmake.html#variable%3aCMAKE_INSTALL_PREFIX).
    3. On Mac [homebrew](https://brew.sh/) is needed for installation. 
3. Install openFHE
    ```console
    git clone https://github.com/openfheorg/openfhe-development.git --branch v1.0.4
    cd openfhe-development && mkdir build && cd build
    cmake .. && make -j8 && sudo make install
    ```
### Build this project

1. Clone the repository and build using cmake.

```console
git clone https://github.com/radhika1601/ScalableMixedModeMPC.git
cd ScalableMixedModeMPC
cmake .
make
```

### Question
Please send email to radhikaradhika2028@u.northwestern.edu
