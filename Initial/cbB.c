#include "../paul.h"

static double gam = 0.0;
static double visc = 0.0;
static int alpha_flag = 0;
static double Mdot = 0.0;
static double R0 = 0.0;
static double xi = 0.0;
static double delta = 0.0;
static double rho0 = 0.0;
static struct planet *thePlanets = NULL;
static double Mach = 0.0;
static int Npl = 0;

double get_cs2(double *);

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    alpha_flag = theDomain->theParList.alpha_flag;
    Mdot = theDomain->theParList.initPar1;
    R0 = theDomain->theParList.initPar2;
    delta = theDomain->theParList.initPar5;
    xi = theDomain->theParList.initPar5;
    rho0 = theDomain->theParList.initPar3;
    Mach = theDomain->theParList.Disk_Mach;
    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double R = r/R0;
    double phi = x[1];

    double cs2 = get_cs2(x);

    double rho, fact;

    double sig0 = Mdot*Mach*Mach/(3.0*3.14159265*visc);

    fact = exp(-pow((R),-xi));
    rho = sig0*pow(R,-delta)*fact+rho0;
    double om = 1.0;
    if (r > 1.0){
      om = pow(r,-1.5);
      om = om*(1.0 + 3.0/(r*r*16));
      double dpdr = sig0*fact*pow(R, -delta)*(xi*pow(R, -xi - 1.0) - delta/R);
      dpdr = dpdr*cs2 + rho/(r*r*Mach);
      //om = sqrt(om*om + dpdr/(r*rho));
    }
    
    double nu;
    nu = visc*cs2/om;

    double v = -1.5*nu/(r+0.0001);
    double P = rho*cs2/gam;
 
    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v;
    prim[UPP] = om;
    prim[UZZ] = 0.0;

    int q2;
    for(q2 = 5; q2 < NUM_Q; q2++)
        prim[q2] = 0.0;
}
