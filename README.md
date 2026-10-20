# Scalable Mixed-Mode MPC

# Installation

### Dependencies

This project builds against emp-tool and emp-ot 1.0. These need C++20, CMake 3.25 and OpenSSL 3.

```console
git clone https://github.com/emp-toolkit/emp-tool.git
git clone https://github.com/emp-toolkit/emp-ot.git
cmake -S emp-tool -B emp-tool/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<prefix>
cmake --build emp-tool/build -j && cmake --install emp-tool/build
cmake -S emp-ot -B emp-ot/build -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<prefix> -DCMAKE_PREFIX_PATH=<prefix>
cmake --build emp-ot/build -j && cmake --install emp-ot/build
```

On macOS, remove any older emp headers from `/usr/local/include`. Apple clang searches that directory before the package include path.

### Build this project

```console
git clone https://github.com/radhika1601/ScalableMixedModeMPC.git
cd ScalableMixedModeMPC
cmake -S . -B build -DCMAKE_PREFIX_PATH=<prefix>
cmake --build build -j
cd build && ctest
```

Tests use a triple buffer of 2^20 (`test/test_config.h`). The library default is about 15M triples, which needs several GB per party.

# Running tests

To run the tests across multiple servers you can pass the ip configuration as a txt file with the following format.

```
<ip1> <port1>
<ip2> <port2>
```

### Question
Please send email to radhikaradhika2028@u.northwestern.edu
