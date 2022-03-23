
#include "paul.h"
#include "geometry.h"

int num_diagnostics(void);
void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain );

void zero_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);

   int i,k,q;
   for( k=0 ; k<Nz ; ++k ){
      for( i=0 ; i<Nr ; ++i ){
         for( q=0 ; q<Nq ; ++q ){
            int iq = k*Nr*Nq + i*Nq + q;
               theTools->Qrz[iq] = 0.0;
         }
      }
   }
   theTools->t_avg = 0.0;

}

void avg_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   double dt = theTools->t_avg;

   int i,k,q; 
   for(k=0; k<Nz; k++){
      for( i=0 ; i<Nr ; ++i ){
         for( q=0 ; q<Nq ; ++q ){
            int iq = k*Nq*Nr + i*Nq + q; 
            theTools->Qrz[iq] /= dt; 
         }
      }
   }
       
   theTools->t_avg = 0.0; 

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

