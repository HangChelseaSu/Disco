#include "../paul.h"

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

double get_cs2(double *);

//N.B. actuall -phi and -f
double phigrav( double , double , double , int); //int here is type
double fgrav( double , double , double , int);

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
    double R = pow( pow(r, 20.0) + pow(0.01, 20.0) , 0.05);

    double phi = x[1];

    double cs2 = get_cs2(x);

    double rho, efact, fth, dfth;

    //double om = 1.0;
    double om = pow(R,-1.5);
    if (R<0.1) om = pow(0.1, -1.5);
    int np;
    double alpha = visc;
    double nu = visc;

    double cosp, sinp, px, py, dx, dy, gx, gy, mag;
    gx = r*cos(phi);
    gy = r*sin(phi);

    if (alpha_flag == 1){
      if (Npl < 2){
          nu = alpha*cs2/sqrt(om);
      }
      else{
        double omtot = 0;
        for(np = 0; np<Npl; np++){
          cosp = cos(thePlanets[np].phi);
          sinp = sin(thePlanets[np].phi);
          px = thePlanets[np].r*cosp;
          py = thePlanets[np].r*sinp;
          dx = gx-px;
          dy = gy-py;
          //mag = dx*dx + dy*dy + thePlanets[np].eps*thePlanets[np].eps;
          mag = sqrt(dx*dx + dy*dy);
          omtot +=  fgrav( thePlanets[np].M, mag , thePlanets[np].eps, thePlanets[np].type)/mag;
        }  	
        nu = alpha*cs2/sqrt(omtot);
        om = sqrt(omtot);
      }
    }

    double phitot = 0.0;
    double dphitot = 0.0;
    for (np = 0; np<Npl; np++){
      cosp = cos(thePlanets[np].phi);
      sinp = sin(thePlanets[np].phi);
      px = thePlanets[np].r*cosp;
      py = thePlanets[np].r*sinp;
      dx = gx-px;
      dy = gy-py;
      mag = sqrt(dx*dx + dy*dy);
      phitot -= phigrav( thePlanets[np].M, mag , thePlanets[np].eps, thePlanets[np].type);
      dphitot += fgrav( thePlanets[np].M, mag , thePlanets[np].eps, thePlanets[np].type);
    }

    double sig0 = 1.0/(3.0*M_PI*nu);
    efact = exp(-pow((R/redge),-xi));
    
    rho = sig0*efact + epsfl;
    double drho = sig0*efact*xi*pow((R/redge),-xi)/R;
 
    double v = -1.5*nu/(R);
    double P = -rho*phitot/(Mach*Mach);

    double multom = 1.0 + 0.75*massq/(R*R*(1.0 + massq)*(1.0 + massq));
    double addom = rho*dphitot + phitot*drho;
    addom *= 1.0/(Mach*Mach*r*rho);
    om = sqrt(fabs(om*om*multom + addom));

    prim[RHO] = rho;
    prim[PPP] = P;
    prim[URR] = v;
    prim[UPP] = om;
    prim[UZZ] = 0.0;

    int q2;
    for(q2 = 5; q2 < NUM_Q; q2++)

        prim[q2] = 0.0;
}
