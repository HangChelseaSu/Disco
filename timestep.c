#include "paul.h"
#include "analysis.h"
#include "planet.h"

void onestep( struct domain * , double , double , int , int , double );

void timestep( struct domain * theDomain , double dt ){
   
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int stepper = theDomain->theParList.Timestep;

   int jk;

   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      memcpy(theDomain->RKcons[jk], theDomain->cons[jk],
              Np[jk]*NUM_Q*sizeof(double));
      memcpy(theDomain->RK_Phi[jk], theDomain->Phi[jk],
              Np[jk]*NUM_FACES*sizeof(double));
   }

   copy_RK_diag(theDomain);
   copyPlanetsRK(theDomain);

   if(stepper == 1)
   {
      onestep( theDomain , 0.0 ,     dt , 1 , 1 , dt );
      theDomain->t += dt;
   }
   else if(stepper == 3)
   {
      onestep( theDomain ,   0.0,       dt, 1 , 0 , dt );
      theDomain->t += dt;
      onestep( theDomain ,  0.75,  0.25*dt, 0 , 0 , dt );
      theDomain->t -= dt*0.5;
      onestep( theDomain , 1./3., 2.*dt/3., 0 , 1 , dt );
      theDomain->t += dt*0.5;
   }
   else
   {
      onestep( theDomain , 0.0 ,     dt , 1 , 0 , dt );
      // Second RK2 timestep occurs at t^n+1
      theDomain->t += dt;   
      onestep( theDomain , 0.5 , 0.5*dt , 0 , 1 , dt );
   }

   add_diagnostics( theDomain , dt );

   theDomain->count_steps += 1;

}
