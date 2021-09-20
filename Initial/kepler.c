#include "../paul.h"
#include "../omega.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;

double get_nu( const double *, const double *);


void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double R = sqrt(r*r + eps*eps);

   double omega02 = 1.0/pow(R,3.);
   double omegaP2 = (1./(Mach*Mach*R*R*R));

   double omega2 = fmax( (omega02 - omegaP2), 0.0 );
   double omega = sqrt(omega2);

   double cs2 = get_cs2(x);
   double visc = get_nu(x, prim);
   double rho = 1.0/visc;
   //if (nu > 0.0) rho = rho/nu;
   double Pp = rho*cs2/gam;

   double X = 0.0;
   if( r*cos(x[1]) > 0.0 ) X = 1.0;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*visc/(R);
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
