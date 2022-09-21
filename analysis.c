
#include "paul.h"
#include "geometry.h"
#include "analysis.h"

void zero_diagnostics( struct domain * theDomain )
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int Nq = theDomain->num_tools;
    struct diagnostic_avg * theTools = &(theDomain->theTools);

    if(Nq > 0)
    {
        memset(theTools->Qrz, 0, Nz*Nr*Nq * sizeof(double));
    }

    if(Nr > 1)
    {
        memset(theTools->F_r,        0, Nz*(Nr-1)*NUM_Q * sizeof(double));
        memset(theTools->Fvisc_r,    0, Nz*(Nr-1)*NUM_Q * sizeof(double));
        memset(theTools->RK_F_r,     0, Nz*(Nr-1)*NUM_Q * sizeof(double));
        memset(theTools->RK_Fvisc_r, 0, Nz*(Nr-1)*NUM_Q * sizeof(double));
    }

    if(Nz > 1)
    {
        memset(theTools->F_z,        0, (Nz-1)*Nr*NUM_Q * sizeof(double));
        memset(theTools->Fvisc_z,    0, (Nz-1)*Nr*NUM_Q * sizeof(double));
        memset(theTools->RK_F_z,     0, (Nz-1)*Nr*NUM_Q * sizeof(double));
        memset(theTools->RK_Fvisc_z, 0, (Nz-1)*Nr*NUM_Q * sizeof(double));
    }


    memset(theTools->S,     0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Sgrav, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Svisc, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Ssink, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Scool, 0, Nz*Nr*NUM_Q * sizeof(double));
    memset(theTools->Sdamp, 0, Nz*Nr*NUM_Q * sizeof(double));

    theTools->t_avg = 0.0;
}

void avg_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   double dt = theTools->t_avg;

   int j,k,q; 
   for(k=0; k<Nz; k++){
      for( j=0 ; j<Nr ; ++j ){
         for( q=0 ; q<Nq ; ++q ){
            int iq = k*Nq*Nr + j*Nq + q; 
            theTools->Qrz[iq] /= dt; 
         }
      }
   }

   for(k=0; k<Nz; k++){
      for( j=0 ; j<Nr-1 ; ++j ){
         int jk = k*(Nr-1) + j;
         for( q=0 ; q<NUM_Q ; ++q ){
            int iq = NUM_Q*jk + q; 
            theTools->F_r[iq] /= dt; 
            theTools->Fvisc_r[iq] /= dt; 
         }
      }
   }

   for(k=0; k<Nz-1; k++){
      for( j=0 ; j<Nr ; ++j ){
         int jk = k*Nr + j;
         for( q=0 ; q<NUM_Q ; ++q ){
            int iq = NUM_Q*jk + q; 
            theTools->F_z[iq] /= dt; 
            theTools->Fvisc_z[iq] /= dt; 
         }
      }
   }

   for(k=0; k<Nz; k++){
       for(j=0; j<Nr; j++){
           int jk = k*Nr + j;
           for(q=0; q<NUM_Q; q++) {
               int iq = NUM_Q*jk + q;
               theTools->S[iq] /= dt;
               theTools->Sgrav[iq] /= dt;
               theTools->Svisc[iq] /= dt;
               theTools->Ssink[iq] /= dt;
               theTools->Scool[iq] /= dt;
               theTools->Sdamp[iq] /= dt;
           }
       }
   }
       
   theTools->t_avg = 0.0; 

}

void adjust_RK_diag(struct domain *theDomain, double RK)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    struct diagnostic_avg * theTools = &(theDomain->theTools);
    int j,k,q; 

    for(k=0; k<Nz; k++){
        for( j=0 ; j<Nr-1 ; ++j ){
            int jk = k*(Nr-1) + j;
            for( q=0 ; q<NUM_Q ; ++q ){
                int iq = NUM_Q*jk + q; 
                theTools->F_r[iq] = (1-RK)*theTools->F_r[iq]
                                    + RK*theTools->RK_F_r[iq]; 
                theTools->Fvisc_r[iq] = (1-RK)*theTools->Fvisc_r[iq]
                                        + RK*theTools->RK_Fvisc_r[iq]; 
            }
        }
    }

    for(k=0; k<Nz-1; k++){
        for( j=0 ; j<Nr; ++j ){
            int jk = k*Nr + j;
            for( q=0 ; q<NUM_Q ; ++q ){
                int iq = NUM_Q*jk + q; 
                theTools->F_z[iq] = (1-RK)*theTools->F_z[iq]
                                    + RK*theTools->RK_F_z[iq]; 
                theTools->Fvisc_z[iq] = (1-RK)*theTools->Fvisc_z[iq]
                                        + RK*theTools->RK_Fvisc_z[iq]; 
            }
        }
    }
    for(k=0; k<Nz; k++) {
        for(j=0; j<Nr; j++) {
            int jk = k*Nr + j;
            for( q=0 ; q<NUM_Q ; ++q ){
                int iq = NUM_Q*jk + q; 
                theTools->S[iq] = (1-RK) * theTools->S[iq]
                                    + RK * theTools->RK_S[iq];
                theTools->Sgrav[iq] = (1-RK) * theTools->Sgrav[iq]
                                        + RK * theTools->RK_Sgrav[iq];
                theTools->Svisc[iq] = (1-RK) * theTools->Svisc[iq]
                                        + RK * theTools->RK_Svisc[iq];
                theTools->Ssink[iq] = (1-RK) * theTools->Ssink[iq]
                                        + RK * theTools->RK_Ssink[iq];
                theTools->Scool[iq] = (1-RK) * theTools->Scool[iq]
                                        + RK * theTools->RK_Scool[iq];
                theTools->Sdamp[iq] = (1-RK) * theTools->Sdamp[iq]
                                        + RK * theTools->RK_Sdamp[iq];
            }
        }
    }
}

void copy_RK_diag(struct domain *theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    struct diagnostic_avg * theTools = &(theDomain->theTools);

    if(Nr > 1)
    {
        memcpy(theTools->RK_F_r, theTools->F_r,
               Nz*(Nr-1)*NUM_Q*sizeof(double));
        memcpy(theTools->RK_Fvisc_r, theTools->Fvisc_r,
               Nz*(Nr-1)*NUM_Q*sizeof(double));
    }

    if(Nz > 1)
    {
        memcpy(theTools->RK_F_z, theTools->F_z,
               (Nz-1)*Nr*NUM_Q*sizeof(double));
        memcpy(theTools->RK_Fvisc_z, theTools->Fvisc_z,
               (Nz-1)*Nr*NUM_Q*sizeof(double));
    }

    memcpy(theTools->RK_S,     theTools->S,     Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Sgrav, theTools->Sgrav, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Svisc, theTools->Svisc, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Ssink, theTools->Ssink, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Scool, theTools->Scool, Nr*Nz*NUM_Q*sizeof(double));
    memcpy(theTools->RK_Sdamp, theTools->Sdamp, Nr*Nz*NUM_Q*sizeof(double));
}


void add_diagnostics( struct domain * theDomain , double dt ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   struct cell ** theCells = theDomain->theCells;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double temp_sum[Nr*Nz*Nq];
   double temp_vol[Nr*Nz];
   memset( temp_sum, 0 , Nr*Nz*Nq*sizeof(double) );
   memset( temp_vol, 0 , Nr*Nz*sizeof(double) );
   int i,j,k,q;
   int kmin = 0;
   int kmax = Nz;
   int jmin = 0;
   int jmax = Nr;

    for(k=kmin; k<kmax; k++)
    {
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
            for(i=0; i<Np[jk]; i++)
            {
                struct cell * c = theCells[jk]+i;
                double phip = c->piph;
                double phim = phip - c->dphi;
                double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};  
                double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
                double xc[3];
                get_centroid_arr(xp, xm, xc);  
                double dV = get_dV(xp,xm);
                double Qrz[Nq];
                get_diagnostics( xc , c->prim , Qrz , theDomain );
                for( q=0 ; q<Nq ; ++q ) 
                    temp_sum[ jk*Nq + q ] += Qrz[q]*dV;
                temp_vol[jk] += dV;
            }
        }
    }

    for(k=kmin; k<kmax; k++)
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
            for(q=0; q<Nq; q++)
                theTools->Qrz[jk*Nq+q] += temp_sum[jk*Nq+q]*dt / temp_vol[jk];
        }

   theTools->t_avg += dt;
}


void run_inst_diagnostics( struct domain * theDomain ,
                            struct diagnostic_inst *theSlowTools )
{
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Nq = num_inst_diagnostics();
   
   theSlowTools->Ntools = Nq;
   if(Nq == 0)
   {
       theSlowTools->Qrz = NULL;
       return;
   }

   theSlowTools->Qrz = (double *)malloc(Nr*Nz*Nq*sizeof(double));
   memset( theSlowTools->Qrz, 0 , Nr*Nz*Nq*sizeof(double) );
   
   struct cell ** theCells = theDomain->theCells;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int i,j,k,q;
   int kmin = 0;
   int kmax = Nz;
   int jmin = 0;
   int jmax = Nr;

    for(k=kmin; k<kmax; k++)
    {
        for(j=jmin; j<jmax; j++)
        {
            int jk = k*Nr + j;
            for(i=0; i<Np[jk]; i++)
            {
                struct cell * c = theCells[jk]+i;
                double phip = c->piph;
                double phim = phip - c->dphi;
                double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};  
                double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
                double xc[3];
                get_centroid_arr(xp, xm, xc);  
                double dV = get_dV(xp,xm);
                double Qrz[Nq];
                get_inst_diagnostics( xc , c->prim , Qrz , theDomain );
                for( q=0 ; q<Nq ; ++q ) 
                    theSlowTools->Qrz[ jk*Nq + q ] += Qrz[q]*dV;
            }
        }
    }
}

void free_inst_diagnostics(struct diagnostic_inst *theSlowTools)
{
    if(theSlowTools->Qrz != NULL)
        free(theSlowTools->Qrz);
}


