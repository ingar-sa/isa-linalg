#!/bin/bash

BuildFolder="$(pwd)/build_linux"

CommonCCompilerFlags="-mt -g0 -fno-exceptions -O0 -fomit-frame-pointer -Wall -Wextra -Wno-unused-parameter -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -o $BuildFolder -c"

#find . -name "*.cpp" -exec gcc $CommonCPPCompilerFlags {} \;
find . -name "*.c" -exec gcc $CommonCCompilerFlags {} \;
