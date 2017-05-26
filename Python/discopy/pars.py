import sys

parnames = ['Restart',
            'T_Start',
            'T_End',
            'Num_Reports',
            'Num_Snapshots',
            'Num_Checkpoints',
            'Use_Logtime',
            'Num_R',
            'Num_Z',
            'aspect',
            'Max_Aspect_Short',
            'Max_Aspect_Long',
            'R_Min',
            'R_Max',
            'Z_Min',
            'Z_Max',
            'Z_Periodic',
            'Phi_Max',
            'Log_Zoning',
            'Log_Radius',
            'CFL',
            'PLM',
            'Riemann_Solver',
            'Mesh_Motion',
            'Absorbing_BC',
            'Initial_Regrid',
            'Density_Floor',
            'Pressure_Floor',
            'Constrained_Transport',
            'Adiabatic_Index',
            'Isothermal',
            'Use_Viscosity',
            'Viscosity',
            'Use_As_Alpha',
            'Mass_Ratio',
            'Eccentricity',
            'Drift_Rate',
            'Drift_Exp',
            'Grav2D',
            'Mach_Number',
            'Include_Atmos',
            'Metric_Par0',
            'Metric_Par1',
            'Metric_Par2',
            'Metric_Par3',
            'Metric_Par4',
            'Init_Par0',
            'Init_Par1',
            'Init_Par2',
            'Init_Par3',
            'Init_Par4',
            'Noise_Type',
            'Noise_Abs',
            'Noise_Rel']

def readParfile(filename):
    # Read a parameter file and load the contents into a dict.

    f = open(filename, "r")

    pars = dict()

    for line in f:

        words = line.split()
        if len(words) < 2:
            continue
        if words[0] in parnames:
            key = words[0]
            sval = words[1]
            try:
                val = int(sval)
            except ValueError:
                try:
                    val = float(sval)
                except ValueError:
                    val = None
            pars[key] = val

    f.close()

    return pars

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python pars.py <parfile>")
        print("Loads <parfile> into a dict, prints to screen")
        sys.exit()

    pars = readParfile(sys.argv[1])
    print(pars)
