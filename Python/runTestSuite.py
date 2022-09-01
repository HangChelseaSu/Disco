import sys
import subprocess
from pathlib import Path
import numpy as np
import discopy as dp


def runAnalysis(size, parFileName, analysisFileName, name):

    analysisDir = Path("test") / "pole_cartInterp_{0:s}_n{1:02d}".format(
        name, size)
    analysisDir.mkdir(exist_ok=True)

    pars = dp.pars.readParfile(parFileName, piMult=False)
    
    hd = Path.cwd()

    nrs = [32, 48, 64, 96, 128, 192, 256, 384, 512]
    Rmax = 1.0

    outputFiles = []

    for nr in nrs:
        pars['Num_Checkpoints'] = 0
        pars['Num_R'] = nr
        pars['Cartesian_Interp_R0'] = size * Rmax / nr
        dp.pars.writeParfile(analysisDir / "in.par", pars)

        if nr <= 64:
            subprocess.run(["../../disco"], cwd=analysisDir)
        else:
            subprocess.run(["mpiexec", "-np", "3", "../../disco"],
                           cwd=analysisDir)

        outfile = "output.{0:04d}.h5".format(nr)
        subprocess.run(["mv", "output.h5", outfile], cwd=analysisDir)
        outputFiles.append(outfile)

    subprocess.run(["python", analysisFileName.resolve(), *outputFiles],
                   cwd=analysisDir)


if __name__ == "__main__":

    size_in_cells = [2, 4, 8, 16]

    parFileName = Path(sys.argv[1])
    analysisFileName = Path(sys.argv[2])
    analysisName = sys.argv[3]

    for n in size_in_cells:
        runAnalysis(n, parFileName, analysisFileName, analysisName)


