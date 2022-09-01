#!/bin/bash

make

# NR=(0016 0032 0064 0128 0256)
# NR=(0032 0048 0064 0096 0128)
# NR=(0128 0256 0512 1024 2048)
# NS=(020 030 040 050 060)
# NS=(050 040 020 030 010)
NS=(050 070 100)
NRa=(3.0 1.0 5.0)
# NRa=(5.0)

for nr in "${NS[@]}"
do
  for nra in "${NRa[@]}"
  do
    if [ "$(uname)" == "Darwin" ]; then
        sed -e "s/^Num_R[[:blank:]].*$/Num_R ${nr}/" -i '' in.par
        sed -e "s/^Num_Checkpoints[[:blank:]].*$/Num_Checkpoints 0/" -i '' in.par
    else
        sed -i "s/^Softening\s.*$/Softening 0.${nr}/" in.par
        sed -i "s/^Sink_Par1\s.*$/Sink_Par1 ${nra}/" in.par
        sed -i "s/^Num_Checkpoints\s.*$/Num_Checkpoints 0/" in.par
    fi
    mpirun -n 32 ./disco
    mv output.h5 output.$nr.$nra.h5
  done
done

#python3 Python/shearCartAnalysis.py output.*.h5
# python3 Python/acousticwaveAnalysis.py output.*.h5
#python3 Python/alfvenwaveAnalysis.py output.*.h5
#python3 Python/magnetosonicwaveAnalysis.py output.*.h5
#python3 Python/advectionAnalysis.py output.*.h5
