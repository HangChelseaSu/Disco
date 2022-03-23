import sys
import numpy as np

parnames = ['Restart',
            'T_Start',
            'T_End',
            'T_Times_2pi',
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
            'R_Periodic',
            'Z_Periodic',
            'Phi_Max',
            'P_Times_2pi',
            'Z_Times_pi',
            'Log_Zoning',
            'Log_Radius',
            'NoBC_Rmin',
            'NoBC_Rmax',
            'NoBC_Zmin',
            'NoBC_Zmax',
            'CFL',
            'PLM',
            'Max_DT',
            'Timestep',
            'Riemann_Solver',
            'Mesh_Motion',
            'Exact_Mesh_Omega',
            'Exact_Mesh_Omega_Par',
            'Absorbing_BC',
            'Initial_Regrid',
            'Density_Floor',
            'Pressure_Floor',
            'Constrained_Transport',
            'Cartesian_Interp',
            'Cartesian_Interp_R0',
            'Adiabatic_Index',
            'Isothermal',
            'Cs2_Profile',
            'Cs2_Par',
            'Energy_Omega',
            'Energy_Omega_Par',
            'Use_Viscosity',
            'Viscosity',
            'Use_As_Alpha',
            'Mass_Ratio',
            'Eccentricity',
            'Drift_Rate',
            'Drift_Exp',
            'Grav2D',
            'RotFrame',
            'RotOmega',
            'RotD',
            'Mach_Number',
            'Include_Atmos',
            'Init_Par0',
            'Init_Par1',
            'Init_Par2',
            'Init_Par3',
            'Init_Par4',
            'Init_Par5',
            'Init_Par6',
            'Init_Par7',
            'Init_Par8',
            'Sink_Type',
            'Sink_Par1',
            'Sink_Par2',
            'Sink_Par3',
            'Sink_Par4',
            'Sink_Par5',
            'Cool_Type',
            'Cool_Par1',
            'Cool_Par2',
            'Cool_Par3',
            'Cool_Par4',
            'Softening',
            'Noise_Type',
            'Noise_Abs',
            'Noise_Rel',
            'Nozzle_Type',
            'Nozzle_Par1',
            'Nozzle_Par2',
            'Nozzle_Par3',
            'Nozzle_Par4',
            'Metric_Par0',
            'Metric_Par1',
            'Metric_Par2',
            'Metric_Par3',
            'Metric_Par4']

def readParfile(filename, piMult=True):
    # Read a parameter file and load the contents into a dict.

    f = open(filename, "r")

    pars = dict()

    t2pi = False
    p2pi = False
    z2pi = False

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

        if words[0] == "T_Times_2pi":
            if int(words[1]) != 0:
                t2pi = True
        if words[0] == "P_Times_2pi":
            if int(words[1]) != 0:
                p2pi = True
        if words[0] == "Z_Times_pi":
            if int(words[1]) != 0:
                z2pi = True

    f.close()

    if piMult:
        if t2pi:
            pars['T_Start'] *= 2*np.pi
            pars['T_End'] *= 2*np.pi
        if p2pi:
            pars['Phi_Max'] *= 2*np.pi
        if z2pi:
            pars['Z_Min'] *= np.pi
            pars['Z_Max'] *= np.pi

    return pars

def writeParfile(filename, pars):

    with open(filename, "w") as f:
        for key in pars:
            f.write("{0:s}    {1:s}\n".format(key, str(pars[key])))

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("usage: $ python pars.py <parfile>")
        print("Loads <parfile> into a dict, prints to screen")
        sys.exit()

    pars = readParfile(sys.argv[1])
    print(pars)
