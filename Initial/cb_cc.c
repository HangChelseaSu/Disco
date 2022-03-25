
#include "../paul.h"

static double gam = 0.0;
static double Mach = 0.0;

double get_cs2(double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    Mach = theDomain->theParList.Disk_Mach;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double phi = x[1];

    double cs2 = get_cs2(x);

    double rho0 = 1.0;
    double Rcav = 2.5;
    double d0 = 1.0e-5;
    //double Rout = 10;
    double v0 = 1.0e-3;

    //double f = 1 - 1/(1 + exp(-2*(r-Rout)));
    double f = 1.0;

    double rho = rho0 * ((1-d0)*exp(-pow(Rcav/r, 12)) + d0) * f;

    double om0 = sqrt((1 - 1/(Mach*Mach)) / (r*r*r));
    double n = 4;
    double omega = pow(pow(om0, -n) + 1, -1/n);
    double vr = v0 * sin(phi) * r * exp(-pow(r/3.5, 6));

    double P = rho*cs2/gam;

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = vr;
    prim[UPP] = omega;
    prim[UZZ] = 0.0;

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        prim[q] = 0.0;
}
