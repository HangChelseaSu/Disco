
#include "../paul.h"
#include <string.h>

void initial( double * , double * );
double get_dV( double * , double * );
double get_centroid( double , double , int );
void cons2prim( double * , double * , double * , double );
void prim2cons( double * , double * , double * , double );
void subtract_omega( double * );

void boundary_trans( struct domain * theDomain , int dim ){

   struct cell ** theCells = theDomain->theCells;
   int * fIndex = theDomain->fIndex_r;

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   if( dim==1 && dim_rank[0] == dim_size[0]-1 ){
      int i,j,k;
      for( j=Nr-1 ; j>Nr-1-NgRb ; --j ){
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { .5*(r_jph[j]+r_jph[j-1]) , phi , .5*(z_kph[k]+z_kph[k-1]) };
               initial( c->prim , x );
               subtract_omega( c->prim );
            }
         }
      }
   }
   if( dim==1 && dim_rank[0] == 0 ){
      double rho_in = 0.0;
      double P_in   = 0.0;
//      double v_in   = 0.0;
      int n_tot = 0;
      int i,k;
      int j=Ng;
      double rp = r_jph[j];
      double rm = r_jph[j-1];
      double rg = .5*(rp+rm);

      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            rho_in += c->prim[RHO];
            P_in   += c->prim[PPP];
//            v_in   += c->prim[RHO]*c->prim[URR];
            ++n_tot;
         }
      }
      rho_in /= (double)n_tot;
      P_in   /= (double)n_tot;

//      v_in   /= (double)n_tot;
//      v_in   /= rho_in;
//      if( v_in > 0.0 ) v_in = 0.0;
      double k0 = 0.5;
      double nu0 = 0.5;

      for( j=0 ; j<NgRa ; ++j ){
         double r = get_centroid(r_jph[j], r_jph[j-1], 1);
         for( k=0 ; k<Nz ; ++k ){
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = {r , phi , z};
               initial( c->prim , x ); 
               subtract_omega( c->prim );
               c->prim[RHO] = rho_in*pow(r/rg,-k0);
               c->prim[PPP] = P_in*pow(r/rg,-k0-2.*nu0);
//               c->tempDoub  = 0.0;
//               c->prim[URR] = v_in;
            }    
         }  
      }
   }

}


