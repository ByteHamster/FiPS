# FiPS

**Fi**nperprint **P**erfect Hashing through **S**orting.
A very simple (125 lines of code) and fast implementation of perfect hashing through fingerprinting.

```
@inproceedings{muller2014retrieval,
  title={Retrieval and perfect hashing using fingerprinting},
  author={M{\"u}ller, Ingo and Sanders, Peter and Schulze, Robert and Zhou, Wei},
  booktitle={Experimental Algorithms: 13th International Symposium, SEA 2014, Copenhagen, Denmark, June 29--July 1, 2014. Proceedings 13},
  pages={138--149},
  year={2014},
  organization={Springer}
}
```

### Library usage

Clone this repository (with submodules) and add the following to your `CMakeLists.txt`.

```
add_subdirectory(path/to/FiPS)
target_link_libraries(YourTarget PRIVATE FiPS)
```

You can construct a PHF as follows.

```cpp
std::vector<std::string> keys = ...;
fips::FiPS hashFunc(keys);
std::cout<<hashFunc(keys.at(0))<<std::endl;
```

### Licensing
This code is licensed under the [GPLv3](/LICENSE).
