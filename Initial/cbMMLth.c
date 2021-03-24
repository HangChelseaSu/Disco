#include "../paul.h"
#include "../omega.h"

static double gam = 0.0;
static double visc = 0.0;
static int alpha_flag = 0;
static struct planet *thePlanets = NULL;
static double Mach = 0.0;
static int Npl = 0;
static double massq = 0.0;
static double xi = 0.0;
static double rin = 0.0;
static double redge = 0.0;
static double rswitch = 0.0;
static double epsfl = 0.0;

void setICparams( struct domain * theDomain )
{
    gam = theDomain->theParList.Adiabatic_Index;
    visc = theDomain->theParList.viscosity;
    alpha_flag = theDomain->theParList.alpha_flag;
    Mach = theDomain->theParList.Disk_Mach;
    massq = theDomain->theParList.Mass_Ratio;
    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;
    xi = theDomain->theParList.initPar1;
    rin = theDomain->theParList.initPar2;
    redge = theDomain->theParList.initPar3;
    rswitch = theDomain->theParList.rmax;
    rswitch = rswitch*0.5;
    epsfl = theDomain->theParList.Density_Floor;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double R = r + 0.05;
    double phi = x[1];

    double cs2 = get_cs2(x);

    double rho, efact, fth, dfth;

    int np;
    double phitot = 0.0;
    double dphitot = 0.0;
    for (np = 0; np<Npl; np++){
      double denom = thePlanets[np].eps*thePlanets[np].eps + thePlanets[np].r*thePlanets[np].r;
      denom += r*r - r*thePlanets[np].r*cos(phi - thePlanets[np].phi);
      double sqdenom = sqrt(denom);
      phitot -= thePlanets[np].M/sqdenom;
      dphitot += thePlanets[np].M*(r - thePlanets[np].r*cos(phi - thePlanets[np].phi))/(denom*sqdenom);
    }

    double nu = visc;
    double dnu = 0.0;
    if (alpha_flag == 1){
      double ihom = 1.0/get_height_om(x);
      dnu = nu*(-dphitot*ihom/(Mach*Mach) + 1.5*cs2*ihom*ihom*pow(r,-2.5)); //approximates dOmega/dr as Keplerian
      nu = nu*cs2*ihom;
    }

    double sig0 = 1.0/(3.0*M_PI*nu);
    double dsig0 = -dnu/(3.0*M_PI*nu*nu);

    efact = exp(-pow((R/redge),-xi));
    rho = sig0*efact + epsfl;
    double drho = sig0*efact*xi*pow((R/redge),-xi)/R + efact*dsig0;

    double v = -1.5*nu/(R);
    double P = rho*cs2/gam;

    double multom = 1.0 + 0.75*massq/(R*R*(1.0 + massq)*(1.0 + massq));
    double addom = rho*dphitot + phitot*drho;
    addom *= 1.0/(Mach*Mach*r*rho);
    double om = sqrt(fabs(om*om*multom + addom));

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v;
    prim[UPP] = om;
    prim[UZZ] = 0.0;

    int q2;
    for(q2 = 5; q2 < NUM_Q; q2++)

        prim[q2] = 0.0;
}
