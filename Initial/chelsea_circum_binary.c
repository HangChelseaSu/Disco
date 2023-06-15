
#include "../paul.h"
#include "../omega.h"

static double nu = 0;
static double gam = 0.0;
static int isothermal_flag = 0;

void setICparams( struct domain * theDomain ){
   nu = theDomain->theParList.viscosity;
   gam = theDomain->theParList.Adiabatic_Index;
   isothermal_flag = theDomain->theParList.isothermal_flag;
}

void initial( double * prim , double * x ){

   double rho = 1.0;
   double P0  = 0.01;
   double GM  = 1.0;
   double M_dot = 0.1;
   double J_dot = 0;
   double vr = 0;

   double r   = x[0];
   double phi = x[1];
   double omega = sqrt(GM / (r * r * r));
   double x_coord = r * cos(phi);
   double y_coord = r * sin(phi);
   double dist = sqrt((x_coord - 1) * (x_coord - 1) + (y_coord) * (y_coord));
   double width = 0.1;
   double amp = 0.1;
   double f = 1 - J_dot / (M_dot * sqrt(GM * r));
   
   rho = M_dot * f / (3 * M_PI * nu);
   vr = 3 * nu / (-2 * r * f);

   if (isothermal_flag == 1){
   P0 = rho * get_cs2(x) / gam;
   }
   // P0 = rho * get_cs2(x) / gamma;

   rho = -1 * (J_dot / (sqrt(GM * r)) - M_dot) / (3 * M_PI * nu);
   rho += exp(- (dist * dist) / (width * width)) * amp;

   if (r < 1){
      vr = 0;
      omega = 1;
   }

   double X = 0.0;
   if( cos(phi) > 0.0 ) X = 1.0;

   prim[RHO] = rho;
   prim[PPP] = P0;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
