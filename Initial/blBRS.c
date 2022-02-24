#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double eps = 0.0;
static int alpha_flag = 0;

static double rbl = 0.0;
static double dbl = 0.0;
static double om0 = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
   alpha_flag = theDomain->theParList.alpha_flag;
   rbl = theDomain->theParList.initPar1;
   dbl = theDomain->theParList.initPar2;
   om0 = theDomain->theParList.initPar3;
}

void initial( double * prim , double * x ){

   double r = x[0];
   double R = sqrt(r*r);

   double cs2 = 1.0/(Mach*Mach);
   double omega = 1.0/sqrt(r*r*r);
   double rho = 1.0;
   double A, B, C;
   double off1, off2, off3;

   A = rbl - dbl*0.5;
   B = (1.0/pow(rbl + dbl*0.5, 1.5) - om0)/dbl;
   C = om0 - B*A;
   off1 = 1./(rbl+dbl*0.5) + C*C*0.5*pow(rbl+dbl*0.5, 2.0) + 2.0*C*B*pow(rbl+dbl*0.5, 3.0)/3.0+ B*B*0.25*pow(rbl+dbl*0.5, 4.0);
   off2 = 1./(rbl-dbl*0.5) + C*C*0.5*pow(rbl-dbl*0.5, 2.0) + 2.0*C*B*pow(rbl-dbl*0.5, 3.0)/3.0+ B*B*0.25*pow(rbl-dbl*0.5, 4.0);
   off3 = 1./(rbl - dbl*0.5) + om0*om0*0.5*pow(rbl-dbl*0.5, 2.0);


   if (r <= rbl + 0.5*dbl) { 
      omega = C + r*B;
      rho = exp(gam*Mach*Mach*(1./r + C*C*0.5*r*r + 2.0*C*B*r*r*r/3.0 + B*B*0.25*r*r*r*r - off1));
   }
   if (r < rbl - 0.5*dbl) {
      omega = om0;
      rho = exp(gam*Mach*Mach*(1./r + om0*om0*0.5*r*r + off2 - off1 - off3)); 
   }

   double visc = nu;
   //if (alpha_flag == 1) visc = nu*cs2/omega;
   //double rho = 1.0;
   //if (nu > 0.0) rho = rho/nu;
   double Pp = rho*cs2/gam;

   double X = 0.0; 
   if( r > rbl - 0.5*dbl ) X = (r-(rbl-dbl*0.5))/dbl; 
   if( r > rbl + 0.5*dbl ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
