# FiPS

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Build status](https://github.com/ByteHamster/FiPS/actions/workflows/build.yml/badge.svg)

A minimal perfect hash function (MPHF) maps a set S of n keys to the first n integers without collisions.
Perfect hash functions have applications in databases, bioinformatics, and as a building block of various space-efficient data structures.

FiPS (**Fi**nperprint **P**erfect Hashing through **S**orting) is a very simple
and fast implementation of [perfect hashing through fingerprinting](https://doi.org/10.1007/978-3-319-07959-2_12).

The idea of perfect hashing through fingerprinting is to hash each input key to a fingerprint.
An array stores a bit for each possible fingerprint indicating whether keys collided.
Colliding keys are handled in another layer of the same data structure.
There are many implementations of the approach, for example [FMPH (Rust)](https://docs.rs/ph/latest/ph/).

The construction of perfect hashing through fingerprinting is usually implemented as a linear scan,
producing a fault for every key.
Instead, FiPS is based on sorting the keys, which is more cache friendly and faster.
Also, it interleaves its select data structure with the payload data, which enables very fast queries.
Lastly, the FiPS implementation is very simple, with just about 200 lines of code.

### Library usage

Clone this repository (with submodules) and add the following to your `CMakeLists.txt`.

```
add_subdirectory(path/to/FiPS)
target_link_libraries(YourTarget PRIVATE FiPS::fips)
```

You can construct a FiPS perfect hash function as follows.

```cpp
std::vector<std::string> keys = {"abc", "def", "123", "456"};
fips::FiPS<> hashFunc(keys, /* gamma = */ 2.0f);
std::cout << hashFunc("abc") << std::endl;
```

### Licensing
This code is licensed under the [GPLv3](/LICENSE).
