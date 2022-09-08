#include "paul.h"
#include "geometry.h"
#include "hydro.h"
#include "omega.h"
#include "planet.h"
#include "riemann.h"
#include "sink.h"

#include <string.h>

void clean_pi( struct domain * theDomain ){
   
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   double phi_max = theDomain->phi_max;

   double **piph = theDomain->piph;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         double phi = piph[jk][i];
         while( phi > phi_max ){ phi -= phi_max; }
         while( phi < 0.0     ){ phi += phi_max; }
         piph[jk][i] = phi; 
      }    
   }
   int p;
   for( p=0 ; p<Npl ; ++p ){
      struct planet * pl = theDomain->thePlanets + p;
      double phi = pl->phi;
      while( phi > phi_max ){ phi -= phi_max;}
      while( phi < 0.0     ){ phi += phi_max;}
      pl->phi = phi; 

      phi = theDomain->pl_kin[p*NUM_PL_KIN + PL_PHI];
      double RKphi = theDomain->pl_RK_kin[p*NUM_PL_KIN + PL_PHI];
      while( phi > phi_max ){ phi -= phi_max; RKphi -= phi_max; }
      while( phi < 0.0     ){ phi += phi_max; RKphi += phi_max; }
      theDomain->pl_kin[p*NUM_PL_KIN + PL_PHI] = phi;
      theDomain->pl_RK_kin[p*NUM_PL_KIN + PL_PHI] = RKphi;
   }

}


double getmindt( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;
   double **wiph = theDomain->wiph;
   double **prim = theDomain->prim;

   double dt = theDomain->theParList.maxDT / theDomain->theParList.CFL;
   if(dt <= 0.0)
       dt = 1.0e100; //HUGE_VAL

   int i,j,k;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            double phip = piph[jk][i];
            double phim = phip - dphi[jk][i];
            double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};
            double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
            int im = i-1;
            if( i==0 ) im = Np[jk]-1; 
            double wm = wiph[jk][im];
            double wp = wiph[jk][i];
            double w = .5*(wm+wp);
            double dt_temp = mindt( &(prim[jk][i*NUM_Q]) , w , xp , xm );
            if( dt > dt_temp ) dt = dt_temp;
         }
      }
   }
   dt *= theDomain->theParList.CFL; 
#if USE_MPI
   MPI_Allreduce( MPI_IN_PLACE , &dt , 1 , MPI_DOUBLE , MPI_MIN , theDomain->theComm );
#endif

   return( dt );
}

void initial( double * , double * );
void restart( struct domain * );
/*
void clear_w( struct domain * theDomain ){
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         theDomain->theCells[jk][i].wiph = 0.0;
      }
   }
}*/

void set_wcell( struct domain * theDomain ){
   int mesh_motion = theDomain->theParList.Mesh_Motion;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;
   double **wiph = theDomain->wiph;
   double **prim = theDomain->prim;

   int i,j,k;
   for( k=0 ; k<Nz ; ++k ){
      double z = get_centroid(z_kph[k], z_kph[k-1], 2);
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
         double r = get_centroid(r_jph[j], r_jph[j-1], 1);
         for( i=0 ; i<Np[jk] ; ++i ){
            double w = 0.0;
            if( mesh_motion ){
               double phip = piph[jk][i];
               double phim = phip - dphi[jk][i];
               double x[3] = {r, 0.5*(phim+phip), z};
               double wL = get_omega( &(prim[jk][i*NUM_Q]) , x );

               int ip = (i+1)%Np[jk];
               phip = piph[jk][ip];
               phim = phip - dphi[jk][ip];
               x[1] = 0.5*(phim+phip);
               double wR = get_omega( &(prim[jk][ip*NUM_Q]) , x );

               w = .5*(wL + wR); 
            }
            wiph[jk][i] = w;
         }
      }    
   }

   if( mesh_motion == 3 ){
      for( k=0 ; k<Nz ; ++k ){
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double w = 0.0;
            for( i=0 ; i<Np[jk] ; ++i ){
               w += wiph[jk][i]; 
            }    
            w /= (double)Np[jk];
            for( i=0 ; i<Np[jk] ; ++i ){
               wiph[jk][i] = w; 
            }    
         }    
      } 
   }
   if( mesh_motion == 4 ){
      for( k=0 ; k<Nz ; ++k ){
         double z = get_centroid(z_kph[k], z_kph[k-1], 2);
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double r = get_centroid(r_jph[j], r_jph[j-1], 1);
            for( i=0 ; i<Np[jk] ; ++i ){
               double phip = piph[jk][i];
               double phim = phip-dphi[jk][i];
               double x[3] = {r, 0.5*(phim+phip), z};
               wiph[jk][i] = mesh_om(x); 
            }    
         }    
      } 
   }
}

void clear_cell( struct cell * );

void regrid( struct domain * theDomain ){
/*
   struct cell ** theCells = theDomain->theCells;
   int Nt = theDomain->Nt;
   int Np = theDomain->Np;
   int * Nr = theDomain->Nr;
   double * t_jph = theDomain->t_jph;
   double * p_kph = theDomain->p_kph;
   double jthresh = 0.1;
   int i,j,k;
   for( j=0 ; j<Nt ; ++j ){
   for( k=0 ; k<Np ; ++k ){
      int jk = j+Nt*k;
      for( i=2 ; i<Nr[jk]-1 ; ++i ){
         double jump = log( theCells[jk][i-1].prim[RHO]/theCells[jk][i+1].prim[RHO] );
         if( fabs(jump) > jthresh ){
            struct cell * c = theCells[jk]+i;
            double rp = c->riph;
            double rm = (c-1)->riph;
            int Nsplit = (int)(fabs(jump/jthresh));
            if(Nsplit>5) Nsplit=5;
            if(Nsplit>3) printf("r=%.2e jump=%.2e split = %d (No Cause For Alarm)\n",theCells[jk][i].riph,jump,Nsplit);
            int blocksize = (Nr[jk]-1) - i;
            Nr[jk] += Nsplit;
            theCells[jk] = (struct cell *) realloc( theCells[jk] , Nr[jk]*sizeof(struct cell) );
            memmove( theCells[jk]+i+1+Nsplit , theCells[jk]+i+1 , blocksize*sizeof(struct cell) );
            int l;
            for( l=0 ; l<Nsplit+1 ; ++l ){
               int m = l+i;
               struct cell * cnew = theCells[jk]+m;
               clear_cell( cnew );
               double rlp = rm + (rp-rm)*( (double)(l+1)/(double)(Nsplit+1) );
               double rlm = rm + (rp-rm)*( (double)(l  )/(double)(Nsplit+1) );
               double xp[3] = {rlp,t_jph[j]  ,p_kph[k]  };
               double xm[3] = {rlm,t_jph[j-1],p_kph[k-1]};
               double x[3];
               int d;
               for( d=0 ; d<3 ; ++d ) x[d] = .5*(xp[d]+xm[d]);
               initial( cnew->prim , x );
               cnew->riph = rlp;
               cnew->dr = rlp-rlm;
               double dV = get_dV(xp,xm);
               double rr = (2./3.)*(rlp*rlp*rlp-rlm*rlm*rlm)/(rlp*rlp-rlm*rlm);
               prim2cons( cnew->prim , cnew->cons , rr , dV );
            }
            i += Nsplit;
         }
      }
   }
   }
*/
}


void adjust_RK_cons( struct domain * theDomain , double RK ){
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   double **cons = theDomain->cons;
   double **RKcons = theDomain->RKcons;
   double **Phi = theDomain->Phi;
   double **RK_Phi = theDomain->RK_Phi;

   int i,jk,q;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         int iq = NUM_Q*i;
         for( q=0 ; q<NUM_Q ; ++q ){
            cons[jk][iq+q] = (1.-RK)*cons[jk][iq+q] + RK*RKcons[jk][iq+q];
         }
         int ip = NUM_FACES*i;
         for( q=0 ; q<NUM_FACES ; ++q ){
            Phi[jk][ip+q] = (1.-RK)*Phi[jk][ip+q] + RK*RK_Phi[jk][ip+q];
         }
      }
   }
}

void move_cells( struct domain * theDomain , double dt){
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;

   double **piph = theDomain->piph;
   double **wiph = theDomain->wiph;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         piph[jk][i] += wiph[jk][i]*dt;
      }
   }
}

void set_cell_xyz(struct domain * theDomain)
{
    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int * Np = theDomain->Np;
    int i,j,k;
   
    for(k=0; k<Nz; k++)
    {
        double zm = theDomain->z_kph[k-1];
        double zp = theDomain->z_kph[k];
        
        for(j=0; j<Nr; j++)
        {
            int jk = k*Nr + j;
            double rm = theDomain->r_jph[j-1];
            double rp = theDomain->r_jph[j];

            double *piph = theDomain->piph[jk];
            double *dphi = theDomain->dphi[jk];
            double *xyz = theDomain->xyz[jk];

            for( i=0 ; i<Np[jk] ; ++i )
            {
                double x[3];
                double xp[3] = {rp, piph[i], zp};
                double xm[3] = {rm, piph[i] - dphi[i], zm};
                get_centroid_arr(xp, xm, x);
                get_xyz(x, xyz + 3*i);
            }
        }
   }
}


void calc_dp( struct domain * theDomain ){
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   
   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;

   int i,jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         int im = i-1;
         if( i == 0 ) im = Np[jk]-1;
         double phim = piph[jk][im];
         double phip = piph[jk][i ];
         dphi[jk][i] = get_dp(phip,phim);
      }
   }
}

void calc_prim( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;
   double **prim = theDomain->prim;
   double **cons = theDomain->cons;

   int i,j,k;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double zm = z_kph[k-1];
      double zp = z_kph[k];
      double z = get_centroid(zp, zm, 2);
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double r = get_centroid(rp, rm, 1);
         for( i=0 ; i<Np[jk] ; ++i ){
            double phip = piph[jk][i];
            double phim = phip - dphi[jk][i];
            double xp[3] = {rp, phip, zp};
            double xm[3] = {rm, phim, zm};
            double dV = get_dV( xp , xm );
            double x[3] = {r, 0.5*(phim+phip), z};
            cons2prim( &(cons[jk][i*NUM_Q]) , &(prim[jk][i*NUM_Q]) ,
                        x , dV, xp, xm );
         }
      }
   }
}

void calc_cons( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;
   double **prim = theDomain->prim;
   double **cons = theDomain->cons;

   int i,j,k;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      double zm = z_kph[k-1];
      double zp = z_kph[k];
      double z = get_centroid(zp, zm, 2);

      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         double rm = r_jph[j-1];
         double rp = r_jph[j];
         double r = get_centroid(rp, rm, 1);

         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            double phip = piph[jk][i];
            double phim = phip - dphi[jk][i];
            double xp[3] = {rp, phip, zp};
            double xm[3] = {rm, phim, zm};
            double dV = get_dV( xp , xm );
            double x[3] = {r, 0.5*(phim+phip), z};
            prim2cons( &(prim[jk][i*NUM_Q] ), &(cons[jk][i*NUM_Q]) ,
                        x , dV , xp, xm);
         }
      }
   }
}

void phi_flux( struct domain * theDomain , double dt ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   int bflag = set_B_flag();
   int CT = theDomain->theParList.CT;
   int jmin, jmax, kmin, kmax;
   if(bflag && CT)
   {
       jmin = NgRa==0 ? 0  : NgRa-1;
       jmax = NgRb==0 ? Nr : Nr-NgRb+1;
       kmin = NgZa==0 ? 0  : NgZa-1;
       kmax = NgZb==0 ? Nz : Nz-NgZb+1;
   }
   else
   {
       jmin = NgRa;
       jmax = Nr-NgRb;
       kmin = NgZa;
       kmax = Nz-NgZb;
   }

   double **E = theDomain->E;
   double **B = theDomain->B;


   int i,j,k;
   for( k=kmin ; k<kmax ; ++k ){
      double zp = z_kph[k];
      double zm = z_kph[k-1];
      double z = get_centroid(zp, zm, 2);

      for( j=jmin ; j<jmax ; ++j ){
         double rp = r_jph[j];
         double rm = r_jph[j-1];
         double r = get_centroid(rp, rm, 1);
         
         int jk = j+Nr*k;

         double *prim = theDomain->prim[jk];
         double *cons = theDomain->cons[jk];
         double *piph = theDomain->piph[jk];
         double *dphi = theDomain->dphi[jk];
         double *wiph = theDomain->wiph[jk];
         double *gradr = theDomain->gradr[jk];
         double *gradp = theDomain->gradp[jk];
         double *gradz = theDomain->gradz[jk];

         for( i=0 ; i<Np[jk] ; ++i ){
            int ip = i<Np[jk]-1 ? i+1 : 0;
            double phi = piph[i];
            double xp[3] = {rp, phi, zp};
            double xm[3] = {rm, phi, zm};
            double dA = get_dA(xp,xm,0); 
            double x[3] = {r, phi, z};

            int iL = i;
            int iR = ip;
            int iqL = NUM_Q*iL;
            int iqR = NUM_Q*iR;
            int ieL = NUM_EDGES*iL;
            int ieR = NUM_EDGES*iR;

            double *EL = (E == NULL) ? NULL : &(E[jk][ieL]);
            double *ER = (E == NULL) ? NULL : &(E[jk][ieR]);
            double *BL = (B == NULL) ? NULL : &(B[jk][ieL]);
            double *BR = (B == NULL) ? NULL : &(B[jk][ieR]);

            riemann_phi( &(prim[iqL]), &(prim[iqR]),
                        &(cons[iqL]), &(cons[iqR]),
                        &(gradr[iqL]), &(gradp[iqL]), &(gradz[iqL]),
                        &(gradr[iqR]), &(gradp[iqR]), &(gradz[iqR]),
                        piph[iL], dphi[iL], piph[iR], dphi[iR], wiph[iL],
                        x, xp, xm, dA*dt, EL, BL, ER, BR );
         }
      }
   }

}

void buildfaces( struct domain * , int , int );

void trans_flux( struct domain * theDomain , double dt , int dim ){

    /*
     * Need to calculate fluxes for all interior faces.
     * For CT also need one level of fluxes outside the interior
     * ie. need the z-fluxes directly outside the r-boundary and vice-versa
     */ 

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int NgRa = theDomain->NgRa;
    int NgRb = theDomain->NgRb;
    int NgZa = theDomain->NgZa;
    int NgZb = theDomain->NgZb;
    
    int bflag = set_B_flag();
    int CT = theDomain->theParList.CT;

    int jmin, jmax, kmin, kmax, Nfr;

    int *fI;
    struct face * theFaces;

    double *fhydro_diag;
    double *fvisc_diag;

    if( dim==1 )
    {
        fI = theDomain->fIndex_r;
        theFaces = theDomain->theFaces_1;
        Nfr = Nr-1;
        jmin = NgRa==0 ? 0 : NgRa-1;
        jmax = NgRb==0 ? Nr-1 : Nr-NgRb;

        if(bflag && CT)
        {
            kmin = NgZa==0 ? 0 : NgZa-1;
            kmax = NgZb==0 ? Nz : Nz-NgZb+1;
        }
        else
        {
            kmin = NgZa;
            kmax = Nz-NgZb;
        }

        fhydro_diag = theDomain->theTools.F_r;
        fvisc_diag = theDomain->theTools.Fvisc_r;
    }
    else
    {
        fI = theDomain->fIndex_z;
        theFaces = theDomain->theFaces_2;
        Nfr = Nr;
        if(bflag && CT)
        {
            jmin = NgRa==0 ? 0  : NgRa-1;
            jmax = NgRb==0 ? Nr : Nr-NgRb+1;
        }
        else
        {
            jmin = NgRa;
            jmax = Nr-NgRb;
        }
        kmin = NgZa==0 ? 0  : NgZa-1;
        kmax = NgZb==0 ? Nz-1 : Nz-NgZb;
        
        fhydro_diag = theDomain->theTools.F_z;
        fvisc_diag = theDomain->theTools.Fvisc_z;
    }

    double **E = theDomain->E;
    double **B = theDomain->B;
    double **E_phi = theDomain->E_phi;

    int j, k, q;
    for(k=kmin; k<kmax; k++)
    {
        double zm, zp;
        if(dim == 1)
            zm = theDomain->z_kph[k-1];
        else
            zm = theDomain->z_kph[k];
        zp = theDomain->z_kph[k];

        for(j=jmin; j<jmax; j++)
        {
            double rm, rp;
            int jkL, jkR;
            if(dim == 1)
            {
                rm = theDomain->r_jph[j];
                jkL = Nr*k + j;
                jkR = Nr*k + j+1;
            }
            else
            {
                rm = theDomain->r_jph[j-1];
                jkL = Nr*k + j;
                jkR = Nr*(k+1) + j;
            }
            rp = theDomain->r_jph[j];

            int JK = j + Nfr*k;

            double *primL = theDomain->prim[jkL];
            double *consL = theDomain->cons[jkL];
            double *gradrL = theDomain->gradr[jkL];
            double *gradpL = theDomain->gradp[jkL];
            double *gradzL = theDomain->gradz[jkL];
            double *piphL = theDomain->piph[jkL];
            double *dphiL = theDomain->dphi[jkL];
            double *primR = theDomain->prim[jkR];
            double *consR = theDomain->cons[jkR];
            double *gradrR = theDomain->gradr[jkR];
            double *gradpR = theDomain->gradp[jkR];
            double *gradzR = theDomain->gradz[jkR];
            double *piphR = theDomain->piph[jkR];
            double *dphiR = theDomain->dphi[jkR];

            double *ELa = (E == NULL) ? NULL : E[jkL];
            double *BLa = (B == NULL) ? NULL : B[jkL];
            double *EphiLa = (E_phi == NULL) ? NULL : E_phi[jkL];
            double *ERa = (E == NULL) ? NULL : E[jkR];
            double *BRa = (B == NULL) ? NULL : B[jkR];
            double *EphiRa = (E_phi == NULL) ? NULL : E_phi[jkR];

            int f;
            for(f=fI[JK]; f<fI[JK+1]; f++)
            {
                double fdAdt_hydro[NUM_Q];
                double fdAdt_visc[NUM_Q];
                int iL = (theFaces + f)->iL;
                int iR = (theFaces + f)->iR;
                int iqL = iL * NUM_Q;
                int iqR = iR * NUM_Q;

                double *EL = (ELa == NULL) ? NULL : ELa + iL*NUM_EDGES;
                double *BL = (BLa == NULL) ? NULL : BLa + iL*NUM_EDGES;
                double *EphiL = (EphiLa == NULL) ? NULL : EphiLa + iL*NUM_EDGES;
                double *ER = (ERa == NULL) ? NULL : ERa + iR*NUM_EDGES;
                double *BR = (BRa == NULL) ? NULL : BRa + iR*NUM_EDGES;
                double *EphiR = (EphiRa == NULL) ? NULL : EphiRa + iR*NUM_EDGES;

                riemann_trans(theFaces + f, primL+iqL, primR+iqR,
                                consL+iqL, consR+iqR,
                                gradrL+iqL, gradpL+iqL, gradzL+iqL,
                                gradrR+iqR, gradpR+iqR, gradzR+iqR,
                                piphL[iL], dphiL[iL], piphR[iR], dphiR[iR],
                                dt, dim, rp, rm, zp, zm,
                                EL, BL, EphiL, ER, BR, EphiR,
                              fdAdt_hydro, fdAdt_visc);

                for(q=0; q<NUM_Q; q++)
                {
                    fhydro_diag[NUM_Q*JK+q] += fdAdt_hydro[q];
                    fvisc_diag[NUM_Q*JK+q] += fdAdt_visc[q];
                }
            }
        }
    }
}


void setup_faces( struct domain * theDomain , int dim ){

   struct face ** theFaces;
   int NN;
   int * nn;
   if( dim==1 ){ 
      theFaces = &(theDomain->theFaces_1);
      nn = theDomain->fIndex_r;
      NN = theDomain->N_ftracks_r;
   }else{
      theFaces = &(theDomain->theFaces_2);
      nn = theDomain->fIndex_z;
      NN = theDomain->N_ftracks_z;
   }

   buildfaces( theDomain , dim , 0 );
   int Nf = nn[NN];
   *theFaces = (struct face *) malloc( Nf*sizeof(struct face) );
   buildfaces( theDomain , dim , 1 );

}

void omega_src( double * , double * , double * , double * , double );

void add_source( struct domain * theDomain , double dt ){

   struct planet * thePlanets = theDomain->thePlanets;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;

   int visc_flag = theDomain->theParList.visc_flag;

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double **prim = theDomain->prim;
   double **cons = theDomain->cons;
   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;
   double **gradr = theDomain->gradr;
   double **gradp = theDomain->gradp;
   double **gradz = theDomain->gradz;

   int i,j,k,p;
   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;

         double *xyz_a = theDomain->xyz[jk];

         for( i=0 ; i<Np[jk] ; ++i ){
            double phip = piph[jk][i];
            double phim = piph[jk][i] - dphi[jk][i];
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double dV = get_dV(xp,xm);
            double *xyz = xyz_a + 3*i;

            double sdVdt_hydro[NUM_Q] = {0};
            double sdVdt_grav[NUM_Q] = {0};
            double sdVdt_visc[NUM_Q] = {0};
            double sdVdt_sink[NUM_Q] = {0};
            double sdVdt_cool[NUM_Q] = {0};
            double sdVdt_damp[NUM_Q] = {0};

            int iq = i*NUM_Q;
            
            source( &(prim[jk][iq]) , sdVdt_hydro, xp , xm , dV*dt  );
            
            for( p=0 ; p<Npl ; ++p ){
               planet_src( thePlanets+p, &(prim[jk][iq]), sdVdt_grav, xp, xm,
                            xyz, dV, dt,
                          theDomain->pl_gas_track + p*NUM_PL_INTEGRALS);
            }
            if(visc_flag)
                visc_source( &(prim[jk][iq]), &(gradr[jk][iq]),
                            &(gradp[jk][iq]), &(gradz[jk][iq]), sdVdt_visc,
                            xp, xm, dV*dt);
            omega_src( &(prim[jk][iq]) , sdVdt_hydro , xp , xm , dV*dt );
            sink_src( &(prim[jk][iq]) , sdVdt_sink , xp , xm , xyz, dV, dt,
                      theDomain->pl_gas_track);
            cooling( &(prim[jk][iq]) , sdVdt_cool , xp , xm , xyz, dV, dt );
            damping( &(prim[jk][iq]) , sdVdt_damp , xp , xm , xyz, dV, dt );

            int q;
            for(q=0; q<NUM_Q; q++)
            {
                cons[jk][iq+q] += sdVdt_hydro[q] + sdVdt_grav[q] + sdVdt_visc[q]
                              + sdVdt_sink[q] + sdVdt_cool[q] + sdVdt_damp[q];
                int iq = NUM_Q*jk + q;
                theDomain->theTools.S[iq] += sdVdt_hydro[q];
                theDomain->theTools.Sgrav[iq] += sdVdt_grav[q];
                theDomain->theTools.Svisc[iq] += sdVdt_visc[q];
                theDomain->theTools.Ssink[iq] += sdVdt_sink[q];
                theDomain->theTools.Scool[iq] += sdVdt_cool[q];
                theDomain->theTools.Sdamp[iq] += sdVdt_damp[q];
            }
         }    
      }    
   }   
}
/*
void longandshort( struct domain * theDomain , double * L , double * S , int * iL , int * iS , struct cell * sweep , int j , int k ){

   int Nr = theDomain->Nr;
   int jk = j + Nr*k;
   int Np = theDomain->Np[jk];
   double * r_jph = theDomain->r_jph;
   double dr = r_jph[j]-r_jph[j-1];
   double r  = r_jph[j];
   double Long = 0.0;
   double Short = 0.0;
   int iLong = 0;
   int iShort = 0;
   int i;
   for( i=0 ; i<Np ; ++i ){
      double dx = dr;
      double dy = r*sweep[i].dphi;
      double l = dy/dx;
      double s = dx/dy;
      if( Long  < l ){ Long  = l; iLong  = i; } 
      if( Short < s ){ Short = s; iShort = i; } 
   }
   *L  = Long;
   *iL = iLong;
   *S  = Short;
   *iS = iShort;

}

void AMRsweep( struct domain * theDomain , struct cell ** swptr , int jk ){

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int j = jk%Nr;
   int k = jk/Nr;

   double MaxShort = theDomain->theParList.MaxShort;
   double MaxLong  = theDomain->theParList.MaxLong;

   struct cell * sweep = *swptr;
   double L,S;
   int iL=0;
   int iS=0;
   longandshort( theDomain , &L , &S , &iL , &iS , sweep , j , k );
   //printf("Long = %e, Short = %e\n",L,S);

   if( S>MaxShort ){
      int iSp = (iS+1)%Np[jk];

      //Possibly shift iS backwards by 1
      int iSm = iS-1;
      if( iSm == -1 ) iSm = Np[jk]-1;
      double dpL = sweep[iSm].dphi;
      double dpR = sweep[iSp].dphi;
      if( dpL < dpR ){
         --iS;
         --iSm;
         --iSp;
         if( iS  == -1 ) iS  = Np[jk]-1;
         if( iSm == -1 ) iSm = Np[jk]-1;
         if( iSp == -1 ) iSp = Np[jk]-1;
      }

      //Remove Zone at iS+1
      sweep[iS].dphi  += sweep[iSp].dphi;
      sweep[iS].piph   = sweep[iSp].piph;
      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iS].cons[q]   += sweep[iSp].cons[q];
         sweep[iS].RKcons[q] += sweep[iSp].RKcons[q];
      }
      double phip = sweep[iS].piph;
      double phim = phip - sweep[iS].dphi;
      double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double r  = get_centroid( xp[0] , xm[0], 1 );
      double dV = get_dV( xp , xm );
      double x[3] = {r, 0.5*(phim+phip), 0.5*(z_kph[k-1]+z_kph[k])};
      cons2prim( sweep[iS].cons , sweep[iS].prim , x , dV, xp, xm );
      //Shift Memory
      int blocksize = Np[jk]-iSp-1;
      if( iSp != Np[jk]-1 ) memmove( sweep+iSp , sweep+iSp+1 , blocksize*sizeof(struct cell) );
      Np[jk] -= 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      if( iS < iL ) iL--;

   }

   if( L>MaxLong ){
      Np[jk] += 1;
      *swptr = (struct cell *) realloc( sweep , Np[jk]*sizeof(struct cell) );
      sweep = *swptr;
      int blocksize = Np[jk]-iL-1;
      memmove( sweep+iL+1 , sweep+iL , blocksize*sizeof(struct cell) );

      double dphi = sweep[iL].dphi;
      double phip = sweep[iL].piph;
      double phim = phip - dphi;
      double phi0 = .5*(phip+phim);

      sweep[iL].piph   = phi0;
      sweep[iL].dphi   = .5*dphi;
      sweep[iL+1].dphi = .5*dphi;

      int q;
      for( q=0 ; q<NUM_Q ; ++q ){
         sweep[iL].cons[q]     *= .5;
         sweep[iL].RKcons[q]   *= .5;
         sweep[iL+1].cons[q]   *= .5;
         sweep[iL+1].RKcons[q] *= .5;
      }

      double xp[3] = {r_jph[j]  ,phi0,z_kph[k]  };
      double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
      double dV = get_dV( xp , xm );
      double r  = get_centroid( xp[0] , xm[0], 1 );
      double x[3] = {r, 0.5*(xp[1]+xm[1]), 0.5*(xp[2]+xm[2])};
      cons2prim( sweep[iL].cons , sweep[iL].prim , x , dV, xp, xm);

      xp[1] = phip;
      xm[1] = phi0;
      dV = get_dV( xp , xm );
      r  = get_centroid( xp[0] , xm[0], 1 );
      x[0] = r;
      x[1] = 0.5*(xp[1]+xm[1]);
      cons2prim( sweep[iL+1].cons , sweep[iL+1].prim , x , dV , xp, xm);

   }
}

void AMR( struct domain * theDomain ){
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      AMRsweep( theDomain , theCells+jk , jk );
   }
}
*/

void print_welcome()
{
    printf("\nDisco!\n\n");
    printf("Git Version: %s\n\n", GIT_VERSION);
    printf("*Compile-time Options*\n");
    printf("HYDRO: %s\n", HYDRO);
    printf("GEOMETRY: %s\n", GEOMETRY);
    printf("INITIAL: %s\n", INITIAL);
    printf("BOUNDARY: %s\n", BOUNDARY);
    printf("OUTPUT: %s\n", OUTPUT);
    printf("RESTART: %s\n", RESTART);
    printf("PLANETS: %s\n", PLANETS);
    printf("HLLD: %s\n", HLLD);
    printf("ANALYSIS: %s\n", ANALYSIS);
    printf("METRIC: %s\n", METRIC);
    printf("FRAME: %s\n", FRAME);
    printf("NUM_C: %d\n", NUM_C);
    printf("NUM_N: %d\n", NUM_N);
    printf("CT_MODE: %d\n", CT_MODE);
    printf("\n");
}


void dump_grid(struct domain *theDomain, char filename[])
{
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int NgRa = theDomain->NgRa;
   int NgRb = theDomain->NgRb;
   int NgZa = theDomain->NgZa;
   int NgZb = theDomain->NgZb;
   int * Np = theDomain->Np;

   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double **prim = theDomain->prim;
   double **cons = theDomain->cons;
   double **RKcons = theDomain->RKcons;
   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;
   double **gradr = theDomain->gradr;
   double **gradp = theDomain->gradp;
   double **gradz = theDomain->gradz;

   int i,j,k, q;

   char fullName[256];
   sprintf(fullName, "%s.%05d.txt", filename, theDomain->count_steps);

   FILE *f = fopen(fullName, "w");

   for( k=NgZa ; k<Nz-NgZb ; ++k ){
      for( j=NgRa ; j<Nr-NgRb ; ++j ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            double phip = piph[jk][i];
            double phim = piph[jk][i] - dphi[jk][i];
            double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double x[3];
            get_centroid_arr(xp, xm, x);

            int iq = NUM_Q*i;
            fprintf(f, "%d %d %d %.12le %.12le %.12le",
                    k, j, i, x[0], x[1], x[2]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", prim[jk][iq+q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", cons[jk][iq+q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", RKcons[jk][iq+q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", gradr[jk][iq+q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", gradp[jk][iq+q]);
            for(q=0; q<NUM_Q; q++)
                fprintf(f, " %.12le", gradz[jk][iq+q]);
            fprintf(f, "\n");
         }    
      }    
   }
   fclose(f);

}
