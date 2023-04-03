# Instructions

## Install NFLlib with AVX2 flag enabled
## Choose parameters
Parameters can be adjusted by modifying `params.hpp` file. We suggest to use d=2048 and d=4096 for 32-bit and 64-bit modulus, respectively. Other parameters can be generated through sage script provided in `params.py`.
## Compile & Run 
The easiest way to run the program is to open console and run `make Test_vericypt && ./Test_vericrypt` command. `NFL_INCLUDE_PATH` and `NFL_LIBRARY_PATH` should be defined if necessary.