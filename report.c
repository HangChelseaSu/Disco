#include "paul.h"
#include "geometry.h"

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

void report( struct domain * theDomain ){

   double t = theDomain->t;
#if USE_MPI
   MPI_Comm grid_comm = theDomain->theComm;
#endif

   struct planet * thePlanets = theDomain->thePlanets;
   int Npl = theDomain->Npl;
   int rank = theDomain->rank;
   int j=0;

   double q = 0.0;
   if (Npl > 1) q = ( thePlanets[1].M / thePlanets[0].M );


   double * M_acc, * La_pls, * Ls_pls, *therm_pls, *Lg_pls, *Eg_pls, *Eacc_pls, *Esink_pls;
   M_acc = calloc(Npl, sizeof(double) );
   La_pls = calloc(Npl, sizeof(double) );
   Ls_pls = calloc(Npl, sizeof(double) );
   Lg_pls = calloc(Npl, sizeof(double) );
   therm_pls = calloc(Npl, sizeof(double) );
   Eacc_pls = calloc(Npl, sizeof(double) );
   Esink_pls = calloc(Npl, sizeof(double) );
   Eg_pls = calloc(Npl, sizeof(double) );


   for(j=0; j<Npl; ++j){
      M_acc[j] = thePlanets[j].dM;
      La_pls[j] = thePlanets[j].accL;
      Ls_pls[j] = thePlanets[j].Ls;
      therm_pls[j] = thePlanets[j].therm;
      Lg_pls[j] = thePlanets[j].gravL;
      Eg_pls[j] = thePlanets[j].gravE;
      Eacc_pls[j] = thePlanets[j].accE;
      Esink_pls[j] = thePlanets[j].sinkE;

      thePlanets[j].dM = 0.0;
      thePlanets[j].RK_dM = 0.0;
      thePlanets[j].accL = 0.0;
      thePlanets[j].RK_accL = 0.0;
      thePlanets[j].Ls = 0.0;
      thePlanets[j].RK_Ls = 0.0;
      thePlanets[j].therm = 0.0;
      thePlanets[j].RK_therm = 0.0;
      thePlanets[j].gravL = 0.0;
      thePlanets[j].RK_gravL = 0.0;
      thePlanets[j].accE = 0.0;
      thePlanets[j].RK_accE = 0.0;
      thePlanets[j].sinkE = 0.0;
      thePlanets[j].RK_sinkE = 0.0;
      thePlanets[j].gravE = 0.0;
      thePlanets[j].RK_gravE = 0.0;

  }

#if USE_MPI
   //MPI_Allreduce( MPI_IN_PLACE , &Mass    , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   //MPI_Allreduce( MPI_IN_PLACE , &Torque  , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );
   //MPI_Allreduce( MPI_IN_PLACE , &Torque2 , 1 , MPI_DOUBLE , MPI_SUM , grid_comm );

   MPI_Allreduce( MPI_IN_PLACE , M_acc  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , La_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Ls_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Lg_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , therm_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Eg_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Eacc_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Esink_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
#endif

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      fprintf(rFile,"%.15le %.15le ",  t, q);
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", Lg_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", M_acc[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", La_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", Ls_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", therm_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", Eg_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", Eacc_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile,"%.15le ", Esink_pls[j]);
      }
      fprintf(rFile,"\n");

      fclose(rFile);
   }
   free(M_acc);
   free(La_pls);
   free(Ls_pls);
   free(therm_pls);
   free(Lg_pls);
   free(Eacc_pls);
   free(Esink_pls);
   free(Eg_pls);
}
