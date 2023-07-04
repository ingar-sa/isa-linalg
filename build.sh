#!/bin/bash

BuildFolder="$(pwd)/build_linux"
CommonCCompilerFlags="-c -Os -fmessage-length=0 -Wall -Wno-unused-function -fno-exceptions -ffunction-sections -fdata-sections -fno-inline -finline-small-functions -fno-common -DMCU -DFPM_DEFAULT -DNODEBUG -DNO_RELOC='0' -g1 -ggdb"
CommonCLinkerFlags="-lm"

# Create the build folder if it does not exist
mkdir -p "${BuildFolder}"

# Find all .c files recursively and compile them
find . -type f -name "*.c" -print0 | while IFS= read -r -d '' file; do
    # Get the base name of the file without extension
    base_name=$(basename "${file}" .c)

    # Compile the .c file to an object file
    gcc ${CommonCCompilerFlags} "${file}" -o "${BuildFolder}/${base_name}.o"

    # Link the object file to create the executable
    gcc "${BuildFolder}/${base_name}.o" ${CommonCLinkerFlags} -o "${BuildFolder}/${base_name}"
done
