#!/bin/bash

make

NR=(0016 0032 0064 0128 0256 0512)

for nr in "${NR[@]}"
do
    if [ "$(uname)" == "Darwin" ]; then
        sed -e "s/^Num_R[[:blank:]].*$/Num_R ${nr}/" -i '' in.par
        sed -e "s/^Num_Checkpoints[[:blank:]].*$/Num_Checkpoints 0/" -i '' in.par
    else
        sed -i "s/^Num_R\s.*$/Num_R ${nr}/" in.par
        sed -i "s/^Num_Checkpoints\s.*$/Num_Checkpoints 0/" in.par
    fi
    if [ "$nr" -lt 100 ]; then
        ./disco
    else
        mpiexec -np 3 ./disco
    fi
    mv output.h5 output.$nr.h5
done

python Python/acousticwaveAnalysis.py output.*.h5
#python Python/alfvenwaveAnalysis.py output.*.h5
