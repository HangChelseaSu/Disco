#include "paul.h"

void planet_RK_copy( struct planet * );
void onestep( struct domain * , double , double , int , int , double );
void add_diagnostics( struct domain * , double );

void timestep( struct domain * theDomain , double dt ){
   
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   int stepper = theDomain->theParList.Timestep;

   int i,jk,p;

   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
         memcpy( c->RK_Phi , c->Phi  , NUM_FACES*sizeof(double) );
      }
   }
   for( p=0 ; p<Npl ; ++p ){
      planet_RK_copy( theDomain->thePlanets + p );
   }

   if(stepper == 1)
   {
      onestep( theDomain , 0.0 ,     dt , 1 , 1 , dt );
      theDomain->t += dt;
   }
   if(stepper == 3)
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
