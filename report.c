#include "paul.h"
#include "hydro.h"
#include "geometry.h"
#include "planet.h"


void report( struct domain * theDomain )
{
   double t = theDomain->t;
#if USE_MPI
   MPI_Comm grid_comm = theDomain->theComm;
#endif

   //double * r_jph = theDomain->r_jph;
   //double * z_kph = theDomain->z_kph;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int *Np = theDomain->Np;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;

   struct cell ** theCells = theDomain->theCells;
   struct planet * thePlanets = theDomain->thePlanets;
   int Npl = theDomain->Npl;

   int jmin = NgRa;
   int jmax = Nr-NgRb;
   int kmin = NgZa;
   int kmax = Nz-NgZb;

   int j,k,i;

   double cons_tot[NUM_Q] = {0};

   for( j=jmin ; j<jmax ; ++j ){
      for( k=kmin ; k<kmax ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            /*
            double phi = c->piph - .5*c->dphi;
            double Pp  = c->prim[PPP];
            double rho = c->prim[RHO];

            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );

            double r = get_centroid( xp[0] , xm[0], 1 );
            */
            
            int l;
            for(l=0; l<NUM_Q; l++)
                cons_tot[l] += c->cons[l];
         }
      }
   }

   int rank = theDomain->rank;

   double q = 0.0;
   if (Npl > 1) q = ( thePlanets[1].M / thePlanets[0].M );

   double * M_acc, * La_pls, * Ls_pls, *therm_pls, *Lg_pls, *Eg_pls, *Eacc_pls;
   M_acc = calloc(Npl, sizeof(double) );
   La_pls = calloc(Npl, sizeof(double) );
   Ls_pls = calloc(Npl, sizeof(double) );
   Lg_pls = calloc(Npl, sizeof(double) );
   therm_pls = calloc(Npl, sizeof(double) );
   Eacc_pls = calloc(Npl, sizeof(double) );
   Eg_pls = calloc(Npl, sizeof(double) );


   for(j=0; j<Npl; ++j){
      M_acc[j] = thePlanets[j].dM;
      La_pls[j] = thePlanets[j].accL;
      Ls_pls[j] = thePlanets[j].Ls;
      therm_pls[j] = thePlanets[j].therm;
      Lg_pls[j] = thePlanets[j].gravL;
      Eg_pls[j] = thePlanets[j].gravE;
      Eacc_pls[j] = thePlanets[j].accE;
  }

   double planet_aux[Npl * NUM_PL_AUX];
   int iq;
   for(j=0; j<Npl; j++)
       for(iq=0; iq<NUM_PL_AUX; iq++)
           planet_aux[j*NUM_PL_AUX+iq] = thePlanets[j].aux[iq];
   
   for(j=0; j<Npl; ++j)
       planet_zero_aux(thePlanets + j);

#if USE_MPI
   MPI_Allreduce( MPI_IN_PLACE , cons_tot    , NUM_Q , MPI_DOUBLE , MPI_SUM , grid_comm );

   MPI_Allreduce( MPI_IN_PLACE , M_acc  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , La_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Ls_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Lg_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , therm_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Eg_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   MPI_Allreduce( MPI_IN_PLACE , Eacc_pls  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   
   MPI_Allreduce( MPI_IN_PLACE , planet_aux  , Npl*NUM_PL_AUX , MPI_DOUBLE ,
                MPI_SUM , grid_comm );
#endif

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      fprintf(rFile,"%.15le",  t);

      for(j=0; j<NUM_Q; j++)
          fprintf(rFile, " %.15le", cons_tot[j]);

      fprintf(rFile," %.15le", q);

      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", Lg_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", M_acc[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", La_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", Ls_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", therm_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", Eg_pls[j]);
      }
      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", Eacc_pls[j]);
      }
      for( j=0; j<Npl; ++j){
        for( iq=0; iq<NUM_PL_AUX; iq++){
         fprintf(rFile," %.15le", planet_aux[j*NUM_PL_AUX+iq]);
        }
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
   free(Eg_pls);
}
