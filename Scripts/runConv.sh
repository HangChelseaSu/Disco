#!/bin/bash

make

# NR=(0016 0032 0064 0128 0256)
# NR=(0032 0048 0064 0096 0128)
NR=(0032 0048 0064 0096 0128 0192 0256)
# NR=(0032 0048 0064 0096 0128 0196 0256 0384 0512)
#NR=(0128 0256 0512 1024 2048)

for nr in "${NR[@]}"
do
    if [ "$(uname)" == "Darwin" ]; then
        sed -e "s/^Num_R[[:blank:]].*$/Num_R ${nr}/" -i '' in.par
        sed -e "s/^Num_Checkpoints[[:blank:]].*$/Num_Checkpoints 0/" -i '' in.par
    else
        sed -i "s/^Num_R\s.*$/Num_R ${nr}/" in.par
        sed -i "s/^Num_Checkpoints\s.*$/Num_Checkpoints 0/" in.par
    fi
    if [ "$nr" -lt 65 ]; then
        ./disco
    else
        mpiexec -np 3 ./disco
    fi
    mv output.h5 output.$nr.h5
done

#python3 Python/shearCartAnalysis.py output.*.h5
python3 Python/acousticwaveAnalysis.py output.*.h5
#python3 Python/alfvenwaveAnalysis.py output.*.h5
#python3 Python/magnetosonicwaveAnalysis.py output.*.h5
# python3 Python/advectionAnalysis.py output.*.h5
# python3 Python/entropywaveAnalysis.py output.*.h5
