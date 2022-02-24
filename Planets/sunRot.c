
#include "../paul.h"

static double q_binary = 1.0;
static double q_bary = 1.0;
static double a_binary = 1.0;
static double a_bary = 1.0;
static double Mach = 1.0;
static double eps = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 3;
   q_bary = theDomain->theParList.Mass_Ratio;
   q_binary = theDomain->theParList.initPar1;
   //a_bary = theDomain->theParList.RotD;
   a_binary = theDomain->theParList.initPar2;

   Mach = theDomain->theParList.Disk_Mach;
   eps = theDomain->theParList.grav_eps;
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double q = q_bary;
   double mu = q/(1.+q);

   double mu2 = q_binary/(1.+q_binary);

   double om = sqrt(mu/pow(a_binary,3.0));

   thePlanets[0].M     = 1.-mu;
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 0.0; 
   thePlanets[0].r     = a_bary;
   thePlanets[0].phi   = 0.0;
   thePlanets[0].eps   = 0.0;
   thePlanets[0].type  = PLPOINTMASS;
   thePlanets[0].RK_dM = 0.0;
   thePlanets[0].dM = 0.0;

   thePlanets[0].accL = 0.0;
   thePlanets[0].RK_accL = 0.0;
   thePlanets[0].Ls = 0.0;
   thePlanets[0].RK_Ls = 0.0;
   thePlanets[0].gravL = 0.0;
   thePlanets[0].RK_gravL = 0.0;
   thePlanets[0].kin = 0.0;
   thePlanets[0].RK_kin = 0.0;
   thePlanets[0].therm = 0.0;
   thePlanets[0].RK_therm = 0.0;

   thePlanets[0].linXmom = 0.0;
   thePlanets[0].RK_linXmom = 0.0;
   thePlanets[0].linYmom = 0.0;
   thePlanets[0].RK_linYmom = 0.0;


   thePlanets[1].M     = mu;
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om; 
   thePlanets[1].r     = 0; 
   thePlanets[1].phi   = M_PI; 
   thePlanets[1].eps   = eps;
   thePlanets[1].type  = PLPOINTMASS;
   thePlanets[1].RK_dM = 0.0;
   thePlanets[1].dM = 0.0;

   thePlanets[1].accL = 0.0;
   thePlanets[1].RK_accL = 0.0;
   thePlanets[1].Ls = 0.0;
   thePlanets[1].RK_Ls = 0.0;
   thePlanets[1].gravL = 0.0;
   thePlanets[1].RK_gravL = 0.0;
   thePlanets[1].kin = 0.0;
   thePlanets[1].RK_kin = 0.0;
   thePlanets[1].therm = 0.0;
   thePlanets[1].RK_therm = 0.0;

   thePlanets[1].linXmom = 0.0;
   thePlanets[1].RK_linXmom = 0.0;
   thePlanets[1].linYmom = 0.0;
   thePlanets[1].RK_linYmom = 0.0;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   //Silence is golden.
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

