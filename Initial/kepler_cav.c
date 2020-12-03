#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;
static int alpha_flag = 0;
static double Npl = 1;
static struct planet *thePlanets = NULL;


double phigrav( double , double , double , int); //int here is type
double fgrav( double , double , double , int);

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
   alpha_flag = theDomain->theParList.alpha_flag;
   thePlanets = theDomain->thePlanets;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double epsfl  = 0.00001;
   double R = 1.0/phigrav( thePlanets[0].M, r , thePlanets[0].eps, thePlanets[0].type);

   double visc = nu;
   int np;

   double phitot = 0.0;
   double dphitot = 0.0;
   for (np = 0; np<Npl; np++){     
     phitot -= phigrav( thePlanets[np].M, r , thePlanets[np].eps, thePlanets[np].type);
     dphitot -= fgrav( thePlanets[np].M, r , thePlanets[np].eps, thePlanets[np].type);
   }

   double redge = 0.15;
   double xi = 10.0;
   double efact = exp(-pow((R/redge),-xi));

   double rho = efact*(1-epsfl) + epsfl;
   double drho = (1-epsfl)*efact*xi*pow((R/redge),-xi)/R;

   double P = -rho*phitot/(Mach*Mach);


   double addom = rho*dphitot + phitot*drho;
   addom *= 1.0/(Mach*Mach*r*rho);
   double omega = 1.0/(R*R*R);
   omega = sqrt(fabs(omega + addom));

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = P;
   prim[URR] = -1.5*visc/R;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
