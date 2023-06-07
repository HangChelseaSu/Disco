
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double rho = 1.0;
   double P0  = 0.01;
   double GM  = 1.0;

   double r   = x[0];
   double phi = x[1];
   double omega = sqrt(GM / (r * r * r));

   double X = 0.0;
   if( cos(phi) > 0.0 ) X = 1.0;

   prim[RHO] = rho;
   prim[PPP] = P0;
   prim[URR] = 0.0;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
