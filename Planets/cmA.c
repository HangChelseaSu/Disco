
#include "../paul.h"

static double q_planet = 1.0;
static double Mach = 1.0;
static double eps = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2; 
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;
 
   double q = q_planet;
   double mu = q/(1.+q);

   double om = pow( a , -1.5 );

   thePlanets[0].M     = 1.-mu; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = om; 
   thePlanets[0].r     = a*mu; 
   thePlanets[0].phi   = M_PI; 
   thePlanets[0].eps   = eps;
   thePlanets[0].type  = PLPOINTMASS;
   thePlanets[0].RK_dM = 0.0;
   thePlanets[0].dM = 0.0;

   thePlanets[0].accL = 0.0;
   thePlanets[0].RK_accL = 0.0;
   thePlanets[0].Ls = 0.0;
   thePlanets[0].RK_Ls = 0.0;
   thePlanets[0].gravL = 0.0;
   thePlanets[0].RK_gravL = 0.0;
   thePlanets[0].therm = 0.0;
   thePlanets[0].RK_therm = 0.0;

   thePlanets[0].accE = 0.0;
   thePlanets[0].RK_accE = 0.0;
   thePlanets[0].gravE = 0.0;
   thePlanets[0].RK_gravE = 0.0;


   thePlanets[1].M     = mu;
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om;  
   thePlanets[1].r     = a*(1.-mu); 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = eps;
   thePlanets[1].type  = PLPOINTMASS;
   thePlanets[1].RK_dM = 0.0;
   thePlanets[1].dM = 0.0;

   thePlanets[1].accL = 0.0;
   thePlanets[1].RK_accL = 0.0;
   thePlanets[1].gravL = 0.0;
   thePlanets[1].RK_gravL = 0.0;
   thePlanets[1].Ls = 0.0;
   thePlanets[1].RK_Ls = 0.0;
   thePlanets[1].therm = 0.0;
   thePlanets[1].RK_therm = 0.0;

   thePlanets[1].gravE = 0.0;
   thePlanets[1].RK_gravE = 0.0;
   thePlanets[1].accE = 0.0;
   thePlanets[1].RK_accE = 0.0;
}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[0].phi += thePlanets[0].omega*dt;
   thePlanets[1].phi += thePlanets[1].omega*dt;
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

