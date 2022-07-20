#include "paul.h"
#include "hydro.h"
#include "geometry.h"
#include "planet.h"

void initializeReport(struct domain *theDomain)
{
    if( theDomain->rank != 0 || theDomain->theParList.restart_flag )
        return;
    FILE * rFile = fopen("report.dat","w");

    //Stamp file with compiler options
    fprintf(rFile, "# GIT_VERSION %s\n", GIT_VERSION);
    fprintf(rFile, "# INITIAL %s\n", INITIAL);
    fprintf(rFile, "# HYDRO %s\n", HYDRO);
    fprintf(rFile, "# GEOMETRY %s\n", GEOMETRY);
    fprintf(rFile, "# BOUNDARY %s\n", BOUNDARY);
    fprintf(rFile, "# OUTPUT %s\n", OUTPUT);
    fprintf(rFile, "# RESTART %s\n", RESTART);
    fprintf(rFile, "# PLANETS %s\n", PLANETS);
    fprintf(rFile, "# HLLD %s\n", HLLD);
    fprintf(rFile, "# ANALYSIS %s\n", ANALYSIS);
    fprintf(rFile, "# METRIC %s\n", METRIC);
    fprintf(rFile, "# FRAME %s\n", FRAME);
    fprintf(rFile, "# CT_MODE %d\n", CT_MODE);

    //Print numerical parameters helpful for parsing
    fprintf(rFile, "# NUM_C %d\n", NUM_C);
    fprintf(rFile, "# NUM_N %d\n", NUM_N);
    fprintf(rFile, "# Npl %d\n", theDomain->Npl);
    fprintf(rFile, "# NUM_PL_KIN %d\n", NUM_PL_KIN);
    fprintf(rFile, "# NUM_PL_AUX %d\n", NUM_PL_AUX);

    fclose(rFile);
}


void report( struct domain * theDomain )
{
   double t = theDomain->t;
#if USE_MPI
   MPI_Comm grid_comm = theDomain->theComm;
#endif

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
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

   int p;
   for(p=0; p<theDomain->Npl; p++)
        theDomain->thePlanets[p].Uf = 0.0;

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
            double phip = c->piph;
            double phim = phip-c->dphi;
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};

            double x[3];
            double rpz[3];
            get_centroid_arr(xp, xm, x);
            get_rpz(x, rpz);
            
            for(p=0; p<theDomain->Npl; p++)
            {
                double Phi = planetaryPotential(theDomain->thePlanets+p,
                        rpz[0], rpz[1], rpz[2]);
                theDomain->thePlanets[p].Uf += c->cons[DDD] * Phi;
            }
         }
      }
   }

   int rank = theDomain->rank;

   double q = 0.0;
   if (Npl > 1) q = ( thePlanets[1].M / thePlanets[0].M );

   double * M_acc, * La_pls, * Ls_pls, *therm_pls, *Lg_pls, *Eg_pls,
          *Eacc_pls, *Uf;
   M_acc = calloc(Npl, sizeof(double) );
   La_pls = calloc(Npl, sizeof(double) );
   Ls_pls = calloc(Npl, sizeof(double) );
   Lg_pls = calloc(Npl, sizeof(double) );
   therm_pls = calloc(Npl, sizeof(double) );
   Eacc_pls = calloc(Npl, sizeof(double) );
   Eg_pls = calloc(Npl, sizeof(double) );
   Uf = calloc(Npl, sizeof(double) );
   
   for(p=0; p<Npl; ++p){
      M_acc[p] = thePlanets[p].dM;
      La_pls[p] = thePlanets[p].accL;
      Ls_pls[p] = thePlanets[p].Ls;
      therm_pls[p] = thePlanets[p].therm;
      Lg_pls[p] = thePlanets[p].gravL;
      Eg_pls[p] = thePlanets[p].gravE;
      Eacc_pls[p] = thePlanets[p].accE;
      Uf[p] = thePlanets[p].Uf;
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
   MPI_Allreduce( MPI_IN_PLACE , Uf  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   
   if(!theDomain->planet_gas_track_synced)
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
         fprintf(rFile," %.15le", Uf[j]);
      }
      for( j=0; j<Npl; ++j){
        for( iq=0; iq<NUM_PL_KIN; iq++){
         fprintf(rFile," %.15le", theDomain->thePlanets[j].kin[iq]);
        }

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
