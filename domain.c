
#include "paul.h"
#include "boundary.h"
#include "geometry.h"
#include "hydro.h"
#include "omega.h"
#include "analysis.h"
#include "planet.h"
#include "report.h"



void setICparams( struct domain * );
void setRiemannParams( struct domain * );
void setGravParams( struct domain * );
void setHlldParams( struct domain * );
void setRotFrameParams( struct domain * );
void setMetricParams( struct domain * );
void setFrameParams(struct domain * );
void setNoiseParams( struct domain * );
void setSinkParams( struct domain * );

int get_num_rzFaces( int , int , int );

void setupDomain( struct domain * theDomain ){

   srand(314159);
   rand();

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   // Allocate all the grid data!

   // Total number of fields
   int N_data = 16;

   // Master arrays
   theDomain->N_data = N_data;
   theDomain->data = (double ***)malloc(N_data * sizeof(double **));
   theDomain->data_len = (int *)malloc(N_data * sizeof(int));

   // Set the internal lengths of each grid field.
   theDomain->data_len[0] = NUM_Q;          //prim
   theDomain->data_len[1] = NUM_Q;          //cons
   theDomain->data_len[2] = NUM_Q;          //RKcons
   theDomain->data_len[3] = NUM_Q;          //gradr
   theDomain->data_len[4] = NUM_Q;          //gradp
   theDomain->data_len[5] = NUM_Q;          //gradz
   theDomain->data_len[6] = 1;              //piph
   theDomain->data_len[7] = 1;              //dphi
   theDomain->data_len[8] = 1;              //wiph
   theDomain->data_len[9] = 3;              //xyz
   theDomain->data_len[10] = NUM_EDGES;      //E
   theDomain->data_len[11] = NUM_EDGES;     //B
   theDomain->data_len[12] = NUM_AZ_EDGES;  //E_phi
   theDomain->data_len[13] = NUM_FACES;     //Phi
   theDomain->data_len[14] = NUM_FACES;     //RK_Phi
   theDomain->data_len[15] = 1;             //tempDoub
   
   // Allocate each of the data fields!
   int field;
   for(field = 0; field < N_data; field++)
   {
       int size = theDomain->data_len[field];

       if(size == 0)
       {
           theDomain->data[field] = NULL;
           continue;
       }

       theDomain->data[field] = (double **)malloc(Nr*Nz*sizeof(double *));

       int jk;
       for(jk=0; jk < Nr*Nz; jk++)
           theDomain->data[field][jk] = (double *)malloc(Np[jk] * size
                                                            * sizeof(double));
   }

   // Set up the named aliases
   theDomain->prim = theDomain->data[0];
   theDomain->cons = theDomain->data[1];
   theDomain->RKcons = theDomain->data[2];
   theDomain->gradr = theDomain->data[3];
   theDomain->gradp = theDomain->data[4];
   theDomain->gradz = theDomain->data[5];
   theDomain->piph = theDomain->data[6];
   theDomain->dphi = theDomain->data[7];
   theDomain->wiph = theDomain->data[8];
   theDomain->xyz = theDomain->data[9];
   theDomain->E = theDomain->data[10];
   theDomain->B = theDomain->data[11];
   theDomain->E_phi = theDomain->data[12];
   theDomain->Phi = theDomain->data[13];
   theDomain->RK_Phi = theDomain->data[14];
   theDomain->tempDoub = theDomain->data[15];

   setGravParams( theDomain );
   setPlanetParams( theDomain );
   int Npl = theDomain->Npl;

   theDomain->thePlanets = NULL;
   theDomain->pl_gas_track = NULL;
   theDomain->pl_kin = NULL;
   theDomain->pl_RK_kin = NULL;
   theDomain->pl_aux = NULL;
   theDomain->pl_RK_aux = NULL;

   if(Npl > 0)
   {
       theDomain->thePlanets = (struct planet *) malloc(
                                    Npl * sizeof(struct planet) );

       theDomain->pl_gas_track = (double *) malloc(
                                    Npl * NUM_PL_INTEGRALS * sizeof(double));
       theDomain->pl_kin = (double *) malloc(
                                    Npl * NUM_PL_KIN * sizeof(double));
       theDomain->pl_RK_kin = (double *) malloc(
                                    Npl * NUM_PL_KIN * sizeof(double));
       theDomain->pl_aux = (double *) malloc(
                                    Npl * NUM_PL_AUX * sizeof(double));
       theDomain->pl_RK_aux = (double *) malloc(
                                    Npl * NUM_PL_AUX * sizeof(double));
   }
   setupPlanets(theDomain);

   int num_tools = num_diagnostics();
   theDomain->num_tools = num_tools;
   theDomain->theTools.t_avg = 0.0;
   theDomain->theTools.Qrz = (double *) malloc( Nr*Nz*num_tools*sizeof(double) );
   theDomain->theTools.F_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                               * sizeof(double));
   theDomain->theTools.Fvisc_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                                   * sizeof(double));
   theDomain->theTools.RK_F_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                                  * sizeof(double));
   theDomain->theTools.RK_Fvisc_r = (double *) malloc((Nr-1) * Nz * NUM_Q
                                                      * sizeof(double));
   theDomain->theTools.F_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                               * sizeof(double));
   theDomain->theTools.Fvisc_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                                   * sizeof(double));
   theDomain->theTools.RK_F_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                                  * sizeof(double));
   theDomain->theTools.RK_Fvisc_z = (double *) malloc(Nr * (Nz-1) * NUM_Q
                                                      * sizeof(double));
   
   theDomain->theTools.S = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
   theDomain->theTools.Sgrav = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
   theDomain->theTools.Svisc = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
   theDomain->theTools.Ssink = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
   theDomain->theTools.Scool = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
   theDomain->theTools.Sdamp = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
   theDomain->theTools.RK_S = (double *) malloc( Nr*Nz*NUM_Q*sizeof(double) );
   theDomain->theTools.RK_Sgrav = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
   theDomain->theTools.RK_Svisc = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
   theDomain->theTools.RK_Ssink = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
   theDomain->theTools.RK_Scool = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));
   theDomain->theTools.RK_Sdamp = (double *) malloc(Nr*Nz*NUM_Q*sizeof(double));

   zero_diagnostics(theDomain);

    //Setup independent of node layout: pick the right rand()'s
    double Pmax = theDomain->phi_max;
    int j, k;

    //Discard everything from lower (global) z-layers.
    for(k=theDomain->N0z_glob; k<theDomain->N0z; k++)
        for(j=0; j<theDomain->Nr_glob; j++)
            rand();

    int i;

    for( k=0 ; k<Nz ; ++k )
    {
        //Discard randoms from inner (global) annuli
        for(j=theDomain->N0r_glob; j<theDomain->N0r; j++)
            rand();

        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];

        //DO the work
        for( j=0 ; j<Nr ; ++j )
        {
            int jk = k*Nr + j;

            double rm = theDomain->r_jph[j-1];
            double rp = theDomain->r_jph[j];

            double p0 = Pmax*(double)rand()/(double)RAND_MAX;
            double dp = Pmax/(double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i )
            {
                double phi = p0+dp*(double)i;
                if( phi > Pmax )
                    phi -= Pmax;
                
                theDomain->piph[jk][i] = phi;
                theDomain->dphi[jk][i] = dp;

                double xp[3] = {rp, phi, zp};
                double xm[3] = {rm, phi-dp, zm};
                double x[3];
                get_centroid_arr(xp, xm, x);
                get_xyz(x, theDomain->xyz[jk] + 3*i);
            }
        }
        //Discard randoms from outer (global) annuli
        for(j=theDomain->N0r+Nr; j<theDomain->N0r_glob + theDomain->Nr_glob;
                j++)
        {
            rand();
        }

    }

   theDomain->t       = theDomain->theParList.t_min;
   theDomain->t_init  = theDomain->theParList.t_min;
   theDomain->t_fin   = theDomain->theParList.t_max;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->count_steps = 0;
   theDomain->final_step  = 0;
   theDomain->check_plz   = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   theDomain->theFaces_1 = NULL;
   theDomain->theFaces_2 = NULL;
   theDomain->N_ftracks_r = get_num_rzFaces( Nr , Nz , 1 );
   theDomain->N_ftracks_z = get_num_rzFaces( Nr , Nz , 2 );

   theDomain->fIndex_r = (int *) malloc( (theDomain->N_ftracks_r+1)*sizeof(int) );
   theDomain->fIndex_z = (int *) malloc( (theDomain->N_ftracks_z+1)*sizeof(int) );

   setICparams( theDomain );
   setHydroParams( theDomain );
   setRiemannParams( theDomain );
   setHlldParams( theDomain );
   setOmegaParams( theDomain );
   setRotFrameParams( theDomain );
   setMetricParams( theDomain );
   setFrameParams( theDomain );
   setDiagParams( theDomain );
   setReportParams(theDomain);
   setNoiseParams( theDomain );
   setBCParams( theDomain );
   setSinkParams( theDomain );
}

void initial( double * , double * ); 
void restart( struct domain * ); 
void calc_dp( struct domain * );
void set_wcell( struct domain * );
void adjust_gas( struct planet * , double * , double * , double );
void set_B_fields( struct domain * );
void subtract_omega( double * );
void addNoise(double *prim, double *x);
void exchangeData(struct domain *, int);

void setupCells( struct domain * theDomain ){

   int restart_flag = theDomain->theParList.restart_flag;
   if( restart_flag ) restart( theDomain );

   int noiseType = theDomain->theParList.noiseType;

   calc_dp( theDomain );

   int i,j,k;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int Npl = theDomain->Npl;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   
   double **prim = theDomain->prim;
   double **cons = theDomain->cons;
   double **gradr = theDomain->gradr;
   double **gradp = theDomain->gradp;
   double **gradz = theDomain->gradz;
   double **wiph = theDomain->wiph;
   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;

   int atmos = theDomain->theParList.include_atmos;

   //Null setup for all cells
   for(k=0; k<Nz; k++){
      for(j=0; j<Nr; j++){
         int jk = Nr*k + j;
         memset(wiph[jk], 0, Np[jk] * sizeof(double));
         memset(gradr[jk], 0, Np[jk] * NUM_Q * sizeof(double));
         memset(gradp[jk], 0, Np[jk] * NUM_Q * sizeof(double));
         memset(gradz[jk], 0, Np[jk] * NUM_Q * sizeof(double));
      }
   }

   //Setup real cells.
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double z = get_centroid( z_kph[k], z_kph[k-1], 2);
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         double r = get_centroid( r_jph[j], r_jph[j-1], 1);
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            double phip = piph[jk][i];
            double phim = phip - dphi[jk][i];
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );
            double phi = piph[jk][i] - 0.5*dphi[jk][i];
            double x[3] = {r, phi, z};

            int iq = NUM_Q*i;

            if( !restart_flag )
            {

               initial( &(prim[jk][iq]) , x );
               subtract_omega( &(prim[jk][iq]) ); 
               if( atmos ){
                  int p;
                  for( p=0 ; p<Npl ; ++p ){
                     double gam = theDomain->theParList.Adiabatic_Index;
                     adjust_gas( theDomain->thePlanets+p , x ,
                                &(prim[jk][iq]) , gam );
                  }
               }
            }
            if(noiseType != 0)
                addNoise(&(prim[jk][iq]), x);
            prim2cons( &(prim[jk][iq]) , &(cons[jk][iq]) , x , dV, xp, xm);
         }    
      }    
   }

#if NUM_FACES > 0
   if(!restart_flag && set_B_flag() && theDomain->theParList.CT)
   {
      // Communicate piph values to ghost zones.
      // TODO: WHY is this obly for CT??
      exchangeData(theDomain, 0);
      if( Nz > 1 )
         exchangeData(theDomain, 1);

      set_B_fields(theDomain);
   }
#endif

   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double z = get_centroid( z_kph[k], z_kph[k-1], 2);
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         double r = get_centroid( r_jph[j], r_jph[j-1], 1);
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            double phip = piph[jk][i];
            double phim = phip - dphi[jk][i];
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV( xp , xm );
            double phi = piph[jk][i] - 0.5*dphi[jk][i];
            double x[3] = {r, phi, z};
            cons2prim( &(cons[jk][i*NUM_Q]) , &(prim[jk][i*NUM_Q]) ,
                        x , dV, xp, xm);
         }
      }
   }

   set_wcell( theDomain );
}


/*
void clear_cell( struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      c->prim[q]   = 0.0;
      c->cons[q]   = 0.0;
      c->RKcons[q] = 0.0;
      c->gradr[q]   = 0.0;
      c->gradp[q]  = 0.0;
      c->gradz[q]  = 0.0;
   }
   c->riph = 0.0;
   c->RKriph = 0.0;
   c->dr = 0.0;
   c->wiph = 0.0;
}
*/

void freeDomain( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   int field;
    for(field = 0; field < theDomain->N_data; field++)
    {
        if(theDomain->data[field] != NULL)
        {
            for( jk=0 ; jk<Nr*Nz ; ++jk )
                if(theDomain->data[field][jk] != NULL)
                    free(theDomain->data[field][jk]);

            free(theDomain->data[field]);
        }
    }
    theDomain->prim = NULL;
    theDomain->cons = NULL;
    theDomain->RKcons = NULL;
    theDomain->gradr = NULL;
    theDomain->gradp = NULL;
    theDomain->gradz = NULL;
    theDomain->piph = NULL;
    theDomain->dphi = NULL;
    theDomain->wiph = NULL;
    theDomain->xyz = NULL;
    theDomain->E = NULL;
    theDomain->B = NULL;
    theDomain->E_phi = NULL;
    theDomain->Phi = NULL;
    theDomain->RK_Phi = NULL;
    theDomain->tempDoub = NULL;

    free(theDomain->data);
    free(theDomain->data_len);

   free( theDomain->Np );
   theDomain->r_jph--;
   free( theDomain->r_jph );
   theDomain->z_kph--;
   free( theDomain->z_kph );

   if(theDomain->thePlanets != NULL)
      free( theDomain->thePlanets );
   if(theDomain->pl_gas_track != NULL)
      free( theDomain->pl_gas_track);
   if(theDomain->pl_kin != NULL)
      free( theDomain->pl_kin);
   if(theDomain->pl_RK_kin != NULL)
      free( theDomain->pl_RK_kin);
   if(theDomain->pl_aux != NULL)
      free( theDomain->pl_aux);
   if(theDomain->pl_RK_aux != NULL)
      free( theDomain->pl_RK_aux);

   free( theDomain->theTools.Qrz );
   free( theDomain->theTools.F_r );
   free( theDomain->theTools.Fvisc_r );
   free( theDomain->theTools.F_z );
   free( theDomain->theTools.Fvisc_z );
   free( theDomain->theTools.RK_F_r );
   free( theDomain->theTools.RK_Fvisc_r );
   free( theDomain->theTools.RK_F_z );
   free( theDomain->theTools.RK_Fvisc_z );
   
   free( theDomain->theTools.S );
   free( theDomain->theTools.Sgrav );
   free( theDomain->theTools.Svisc );
   free( theDomain->theTools.Ssink );
   free( theDomain->theTools.Scool );
   free( theDomain->theTools.Sdamp );
   free( theDomain->theTools.RK_S );
   free( theDomain->theTools.RK_Sgrav );
   free( theDomain->theTools.RK_Svisc );
   free( theDomain->theTools.RK_Ssink );
   free( theDomain->theTools.RK_Scool );
   free( theDomain->theTools.RK_Sdamp );

   free( theDomain->fIndex_r );
   free( theDomain->fIndex_z );

}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   int check=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   if( theDomain->rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
      FILE * latest = NULL;
      latest = fopen("latest","r");
      if( latest ){ check = 1; fclose(latest); remove("latest");}
   }

#if USE_MPI
   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   MPI_Allreduce( MPI_IN_PLACE , &check , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
#endif
   if( final ) theDomain->final_step = 1;
   if( check ) theDomain->check_plz = 1;

}

void report( struct domain * );
void snapshot( struct domain * , char * );
void output( struct domain * , char * );

void possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nsnp = theDomain->N_snp;
   double Nchk = theDomain->N_chk;
   int LogOut = theDomain->theParList.Out_LogTime;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      //longandshort( &theDomain , &L , &S , &iL , &iS , theDomain.theCells[0] , 0 , 0 );
      report( theDomain );
      //if( theDomain->rank==0 ) printf("t = %.5e\n",t);
   }
   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override || theDomain->check_plz ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if( !theDomain->check_plz ){
            if(theDomain->rank==0) printf("Creating Checkpoint #%04d...\n",n0);
            sprintf(filename,"checkpoint_%04d",n0);
         }else{
            if(theDomain->rank==0) printf("Creating Requested Checkpoint...\n");
            sprintf(filename,"checkpoint_latest");
            theDomain->check_plz = 0;
         }
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }
   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override) sprintf( filename , "snapshot_%04d" , n0 );
      else sprintf( filename , "snapshot" );
      //snapshot( theDomain , filename );
   }

}


