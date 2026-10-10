#!/bin/bash

# Description:
# Run check using diff and diffstat to see if solutions mutated


fields=("U" "T" "H" "p_rgh" "alphas" "divB")

# file paths to solution and reference (i.e., MTHD solution)
fpSource=../../laserbeamFoam/Plate2D
fpReference=../../laserbeamFoam/Plate2D_MTHD_original

## Compressible solver:
# fpSource=./richter2d_rerun2
# # Other paths to check: ./richter2dAutoPtr/
# fpReference=./richter2D_MTHD_original


# Check filepaths exist. Throw error if not found
if [ -d "${fpSource}" ]; then
    echo "${fpSource} exists"
else
    echo "${fpSource} Not found; exitting"
    exit -1
fi

if [ -d "${fpReference}" ]; then
    echo "${fpReference} exists"
else
    echo "${fpReference} Not found; exitting"
    exit -1
fi


# Use `foamListTimes` to get all the time directories
folders=$(foamListTimes -case ${fpSource}) # could pass `-time 0:``


# Quick list
## Richter 2D:
# AutoPtr, Orig = Match
# `==` old, Orig = Match
# `==` Recompile, Orig = 
#
## Plate 2D:
#  AutoPtr (T+U), orig = Match
#  AutoPtr (T only), orig = Match
#  With `==`

# File path to output
fpReport=./solution\_check.log

echo "Checking files" > ${fpReport}
echo "fpSource = ${fpSource}" >> ${fpReport}
echo "fpReference = ${fpReference}" >> ${fpReport}


# Main loop
for folder in ${folders[@]}; do
    echo "(${folder})"
    for field in ${fields[@]}; do
        # checks if directory exists:
        if [ -d "${fpSource}/${folder}" ]; then
            echo "(${folder}, ${field})"
            

            fileRef=${fpReference}/${folder}/${field}
            fileSrc=${fpSource}/${folder}/${field}
            
            # Report the differences, if any
            echo "${folder}/${field}:" >> ${fpReport}
            diff -u ${fileRef} ${fileSrc} | diffstat 2>&1 | tee -a ${fpReport}
        fi
    done
done

echo "[checkSolutions.sh] Done."