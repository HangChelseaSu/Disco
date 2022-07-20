#include "paul.h"
#include "hydro.h"
#include "geometry.h"
#include "planet.h"

void initializeReport(struct domain *theDomain)
{
    if( theDomain->rank != 0)
        return;

    // Open read-only to check if it exists
    FILE *rFile = fopen("report.dat", "r");    

    // if report.dat exists and we're restarting, then leave!
    if(rFile && theDomain->theParList.restart_flag )
    {
        fclose(rFile);
        return;
    }

    // Clean up.
    if(rFile)
        fclose(rFile);


    // If we're still here, then either we're restarting but report.dat
    // doesn't exist (and so we need to write the header)
    // or
    // we are not restarting, and need to write report.dat anyways.

    rFile = fopen("report.dat","w");

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

   double *Uf;
   Uf = malloc(Npl * sizeof(double) );
   
   for(p=0; p<Npl; ++p){
      Uf[p] = thePlanets[p].Uf;
   }

#if USE_MPI
   MPI_Allreduce( MPI_IN_PLACE , cons_tot    , NUM_Q , MPI_DOUBLE , MPI_SUM , grid_comm );

   MPI_Allreduce( MPI_IN_PLACE , Uf  , Npl , MPI_DOUBLE , MPI_SUM , grid_comm );
   
   if(!theDomain->planet_gas_track_synced)
      MPI_Allreduce( MPI_IN_PLACE , theDomain->pl_aux, Npl*NUM_PL_AUX,
                    MPI_DOUBLE , MPI_SUM , grid_comm );
#endif

   int iq;

   if( rank==0 ){
      FILE * rFile = fopen("report.dat","a");
      fprintf(rFile,"%.15le",  t);

      for(j=0; j<NUM_Q; j++)
          fprintf(rFile, " %.15le", cons_tot[j]);

      for( j=0; j<Npl; ++j){
         fprintf(rFile," %.15le", Uf[j]);
      }
      for( j=0; j<Npl; ++j){
        for( iq=0; iq<NUM_PL_KIN; iq++){
         fprintf(rFile," %.15le", theDomain->pl_kin[j*NUM_PL_KIN+iq]);
        }

        for( iq=0; iq<NUM_PL_AUX; iq++){
         fprintf(rFile," %.15le", theDomain->pl_aux[j*NUM_PL_AUX+iq]);
        }
      }
      fprintf(rFile,"\n");

      fclose(rFile);
   }

   zeroAuxPlanets(theDomain);

   free(Uf);
}
