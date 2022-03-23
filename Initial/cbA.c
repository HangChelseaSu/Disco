#include "../paul.h"

static double gam = 0.0;
static double visc = 0.0;
static int alpha_flag = 0;
static double Mdot = 0.0;
static double R0 = 0.0;
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

    double xi = 4.0;
    double p = 0.0;
    double rho, fact;

    double om = 1.0;
    if (r > 1.0) om = pow(r,-1.5);
    int np;
    double alpha = visc;
    double nu = visc;
    if (alpha_flag == 1){
      if (Npl < 2){
          nu = alpha*cs2/sqrt(om);
      }
      else{
        double omtot = 0;
        double cosp, sinp, px, py, dx, dy, gx, gy, mag;
        gx = r*cos(phi);
        gy = r*sin(phi);
        for(np = 0; np<Npl; np++){
          cosp = cos(thePlanets[np].phi);
          sinp = sin(thePlanets[np].phi);
          px = thePlanets[np].r*cosp;
          py = thePlanets[np].r*sinp;
          dx = gx-px;
          dy = gy-py;
          mag = dx*dx + dy*dy + thePlanets[np].eps*thePlanets[np].eps;
          omtot +=	thePlanets[np].M*pow(mag, -1.5);
        }  	
        nu = alpha*cs2/sqrt(omtot);
        om = sqrt(omtot);
      }
    }

    double sig0 = Mdot;
    fact = exp(-pow((R),-xi));
    rho = sig0*pow(R,-p)*fact + 0.01;
 
    double v = -1.5*nu/r;
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
