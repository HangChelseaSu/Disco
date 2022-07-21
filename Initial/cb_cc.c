
#include "../paul.h"

static double gam = 0.0;
static double Mach = 0.0;
static double rho0 = 0.0;
static double d0 = 0.0;
static double Rcav = 0.0;
static double Rout = 0.0;
static double v0 = 0.0;

double get_cs2(double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    Mach = theDomain->theParList.Disk_Mach;

    rho0 = theDomain->theParList.initPar1;
    d0 = theDomain->theParList.initPar2;
    Rcav = theDomain->theParList.initPar3;
    Rout = theDomain->theParList.initPar4;
    v0 = theDomain->theParList.initPar5;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double phi = x[1];

    double cs2 = get_cs2(x);

    double f = 1.0;
    if(Rout > 0.0)
        f = 1 - 1/(1 + exp(-2*(r-Rout)));

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
