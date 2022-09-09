
#include "paul.h"
#include "geometry.h"
#include "hydro.h"


double minmod( double a , double b , double c ){
   double m = a;
   if( a*b < 0.0 ) m = 0.0;
   if( fabs(b) < fabs(m) ) m = b;
   if( b*c < 0.0 ) m = 0.0;
   if( fabs(c) < fabs(m) ) m = c;
   return(m);
}

void plm_phi( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double PLM = theDomain->theParList.PLM;
#if ENABLE_CART_INTERP
   int Cartesian_Interp = theDomain->theParList.Cartesian_Interp;
#endif

   double **dphi = theDomain->dphi;
   double **prim = theDomain->prim;
   double **gradp = theDomain->gradp;

   int i,j,k,q;


   for( k=0 ; k<Nz ; ++k ){
      for( j=0 ; j<Nr ; ++j ){
         int jk = j+Nr*k;
#if ENABLE_CART_INTERP
         double z = get_centroid(theDomain->z_kph[k],theDomain->z_kph[k-1],2);
         double r = get_centroid(theDomain->r_jph[j],theDomain->r_jph[j-1],1);
         double x[3] = {r, 0.0, z};
         double w = 0.0;
         if(Cartesian_Interp)
            w = getCartInterpWeight(x);
#endif
         for( i=0 ; i<Np[jk] ; ++i ){
            int im = (i == 0) ? Np[jk]-1 : i-1;
            int ip = (i == Np[jk]-1) ? 0 : i+1;
            double dpL = dphi[jk][im];
            double dpC = dphi[jk][i];
            double dpR = dphi[jk][ip];

            double *primL = &(prim[jk][NUM_Q*im]);
            double *primC = &(prim[jk][NUM_Q*i]);
            double *primR = &(prim[jk][NUM_Q*ip]);
            double *gradpC = &(gradp[jk][NUM_Q*i]);

            for( q=0 ; q<NUM_Q ; ++q ){
               double pL = primL[q];
               double pC = primC[q];
               double pR = primR[q];
               double sL = pC - pL;
               sL /= .5*( dpC + dpL );
               double sR = pR - pC;
               sR /= .5*( dpR + dpC );
               double sC = pR - pL;
               sC /= .5*( dpL + dpR ) + dpC;
               gradpC[q] = minmod( PLM*sL , sC , PLM*sR );
            }
#if ENABLE_CART_INTERP
            if(w > 0.0)
            {
               double xL[3] = {r, -0.5*(dpL+dpC), z}; 
               double xC[3] = {r, 0.0, z}; 
               double xR[3] = {r, 0.5*(dpC+dpR), z}; 
               double cartPrimL[NUM_Q];
               double cartPrimC[NUM_Q];
               double cartPrimR[NUM_Q];
               double cartGrad[NUM_Q];
               double grad[NUM_Q];
               geom_rebase_to_cart(primL, xL, cartPrimL);
               geom_rebase_to_cart(primC, xC, cartPrimC);
               geom_rebase_to_cart(primR, xR, cartPrimR);
               for( q=0 ; q<NUM_Q ; ++q ){
                  double pL = cartPrimL[q];
                  double pC = cartPrimC[q];
                  double pR = cartPrimR[q];
                  double sL = pC - pL;
                  double sR = pR - pC;
                  double sC = pR - pL;
                  if(q == URR || q == UPP)
                  {
                      sL /= sin(.5*( dpC + dpL ));
                      sR /= sin(.5*( dpR + dpC ));
                      sC /= 2*sin(0.5*(.5*( dpL + dpR ) + dpC));
                  }
                  else
                  {
                      sL /= .5*( dpC + dpL );
                      sR /= .5*( dpR + dpC );
                      sC /= .5*( dpL + dpR ) + dpC;
                  }
                  cartGrad[q] = minmod( PLM*sL , sC , PLM*sR );
               }
               geom_gradCart_to_grad(cartGrad, primC, xC, grad, 0);
               for( q=0 ; q<NUM_Q ; ++q )
                  gradpC[q] = (1-w) * gradpC[q] + w * grad[q];
            }
#endif
         }
      }
   }
}

void plm_geom_boundary(struct domain *theDomain, int jmin, int jmax, int kmin,
                        int kmax, int dim, int LR);

void plm_trans( struct domain * theDomain , struct face* theFaces , int Nf , int dim ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   double PLM = theDomain->theParList.PLM;
#if ENABLE_CART_INTERP
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int Cartesian_Interp = theDomain->theParList.Cartesian_Interp;
#endif

   double **piph = theDomain->piph;
   double **dphi = theDomain->dphi;
   double **prim = theDomain->prim;
   double **gradp = theDomain->gradp;
   double **grad = (dim == 1) ? theDomain->gradr : theDomain->gradz;
   double **tempDoub = theDomain->tempDoub;

   const struct face_strip *theStrips;
   int Nfs;
   double **dAf;
   double **xf;

   if(dim == 1)
   {
       theStrips = theDomain->theFaceStripsR;
       Nfs = theDomain->N_fsR;
       dAf = theDomain->dA_fr;
       xf = theDomain->x_fr;
   }
   else
   {
       theStrips = theDomain->theFaceStripsZ;
       Nfs = theDomain->N_fsZ;
       dAf = theDomain->dA_fz;
       xf = theDomain->x_fz;
   }
   
   int j,k,q;

   //Clear gradients
   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         memset(grad[jk], 0, NUM_Q*Np[jk]*sizeof(double));
      }
   }
#if ENABLE_CART_INTERP
   if(Cartesian_Interp && dim == 1)
   {
      for( k=0 ; k<Nz ; ++k ){
         double z = get_centroid(z_kph[k], z_kph[k-1], 2);
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            double r = get_centroid(r_jph[j], r_jph[j-1], 1);
            double x[3] = {r, 0.0, z};
            double w = getCartInterpWeight(x);
            for( i=0 ; i<Np[jk] ; ++i ){
                int iq = NUM_Q*i;
                grad[jk][iq+UPP] = -w * prim[jk][iq+UPP] / r;
             }
          }
       }
   }
#endif

   if(PLM == 0.0)
       return;

   //printf("Doing grads!\n");
   
    //Calc total weight
    for(j = 0; j < Nr; j++)
        for(k = 0; k < Nz; k++)
        {
            int jk = j+Nr*k;
            memset(tempDoub[jk], 0, Np[jk]*sizeof(double));
        }
    
    int n;
    for(n=0; n < Nfs; n++)
    {
        const struct face_strip *fs = theStrips + n;

        int jkL = (fs->kL)*Nr + fs->jL;
        int jkR = (fs->kR)*Nr + fs->jR;

        int iL = fs->iL0;
        int iR = fs->iR0;

        double *piphL = piph[jkL];
        double *piphR = piph[jkR];
        double *dA = dAf[fs->jkf];
        
        int f;
        for(f=fs->fa; f < fs->fb; f++)
        {
            tempDoub[jkL][iL] += dA[f];
            tempDoub[jkR][iR] += dA[f];

            if(get_signed_dp(piphL[iL], piphR[iR]) > 0)
                iR = (iR < Np[jkR]-1) ? iR+1 : 0;
            else
                iL = (iL < Np[jkL]-1) ? iL+1 : 0;
        }
    }
   
    //Add weighted slopes
    
    for(n=0; n < Nfs; n++)
    {
        const struct face_strip *fs = theStrips + n;

        int jkL = (fs->kL)*Nr + fs->jL;
        int jkR = (fs->kR)*Nr + fs->jR;

        int iL = fs->iL0;
        int iR = fs->iR0;

        double *piphL = piph[jkL];
        double *piphR = piph[jkR];
        double *dphiL = dphi[jkL];
        double *dphiR = dphi[jkR];

        double *dA = dAf[fs->jkf];
        double *x = xf[fs->jkf];

        double dxL, dxR;
        if(dim == 1)
        {
            double rL = get_centroid(theDomain->r_jph[fs->jL],
                                     theDomain->r_jph[fs->jL-1], 1);
            double rR = get_centroid(theDomain->r_jph[fs->jR],
                                     theDomain->r_jph[fs->jR-1], 1);
            dxL = fs->rm - rL;
            dxR = rR - fs->rm;
        }
        else
        {
            double zL = get_centroid(theDomain->z_kph[fs->kL],
                                     theDomain->z_kph[fs->kL-1], 2);
            double zR = get_centroid(theDomain->z_kph[fs->kR],
                                     theDomain->z_kph[fs->kR-1], 2);
            dxL = fs->zm - zL;
            dxR = zR - fs->zm;
        }

        int f;
        for(f=fs->fa; f < fs->fb; f++)
        {
            double phi = x[3*f+1];
            double phiL = piphL[iL] - 0.5*dphiL[iL];
            double phiR = piphR[iR] - 0.5*dphiR[iR];
            double dpL = get_signed_dp(phi,phiL);
            double dpR = get_signed_dp(phi,phiR);

            double dA_fL = dA[f] / tempDoub[jkL][iL];
            double dA_fR = dA[f] / tempDoub[jkR][iR];

            int iqL = iL*NUM_Q;
            int iqR = iR*NUM_Q;

            double *primL = &(prim[jkL][iqL]);
            double *primR = &(prim[jkR][iqR]);
            
            double *gradpL = &(gradp[jkL][iqL]);
            double *gradpR = &(gradp[jkR][iqR]);

            double *gradL = &(grad[jkL][iqL]);
            double *gradR = &(grad[jkR][iqR]);
            
            for( q=0 ; q<NUM_Q ; ++q )
            {
                double WL = primL[q] + dpL*gradpL[q];
                double WR = primR[q] + dpR*gradpR[q];
                double S = (WR-WL)/(dxR+dxL);
            
                gradL[q] += S*dA_fL;
                gradR[q] += S*dA_fR;
            }

            if(get_signed_dp(piphL[iL], piphR[iR]) > 0)
                iR = (iR < Np[jkR]-1) ? iR+1 : 0;
            else
                iL = (iL < Np[jkL]-1) ? iL+1 : 0;
        }
    }
    /*
    for(n=0; n<Nf; ++n)
    {
        struct face * f  = &( theFaces[n] );
        double phi = f->cm[1];
        int jkL = f->jkL;
        int iL = f->iL;
        int jkR = f->jkR;
        int iR = f->iR;
        double dxL = f->dxL;
        double dxR = f->dxR;
        double phiL = piph[jkL][iL] - 0.5*dphi[jkL][iL];
        double phiR = piph[jkR][iR] - 0.5*dphi[jkR][iR];
        double dpL = get_signed_dp(phi,phiL);
        double dpR = get_signed_dp(phi,phiR);
        double dA  = f->dA;

        double dA_fL = dA / tempDoub[jkL][iL];
        double dA_fR = dA / tempDoub[jkR][iR];

        int iqL = iL*NUM_Q;
        int iqR = iR*NUM_Q;

        double *primL = &(prim[jkL][iqL]);
        double *primR = &(prim[jkR][iqR]);
        
        double *gradpL = &(gradp[jkL][iqL]);
        double *gradpR = &(gradp[jkR][iqR]);

        double *gradL = &(grad[jkL][iqL]);
        double *gradR = &(grad[jkR][iqR]);

#if ENABLE_CART_INTERP
        double w = 0.0;
        if(theDomain->theParList.Cartesian_Interp)
            w = getCartInterpWeight(f->cm);

        if(w > 0.0)
        {
            double gradCIL[NUM_Q], gradCIR[NUM_Q];
            geom_cart_interp_grad_trans(primL, primR,
                                        gradpL, gradpR,
                                        dpL, dpR, f->cm, -dxL, dxR,
                                        gradCIL, gradCIR, dim);
            for(q=0; q<NUM_Q; q++)
            {
                double WL = primL[q] + dpL*gradpL[q];
                double WR = primR[q] + dpR*gradpR[q];
                double S = (WR-WL)/(dxR+dxL);

                gradL[q] += ((1-w) * S + w * gradCIL[q]) * dA_fL;
                gradR[q] += ((1-w) * S + w * gradCIR[q]) * dA_fR;
            }
        }
        else
        {
#endif
            for( q=0 ; q<NUM_Q ; ++q )
            {
                double WL = primL[q] + dpL*gradpL[q];
                double WR = primR[q] + dpR*gradpR[q];
                double S = (WR-WL)/(dxR+dxL);
            
                gradL[q] += S*dA_fL;
                gradR[q] += S*dA_fR;
            }
#if ENABLE_CART_INTERP
        }
#endif
    }
    */

    //Slope Limiting
    for(n=0; n < Nfs; n++)
    {
        const struct face_strip *fs = theStrips + n;

        int jkL = (fs->kL)*Nr + fs->jL;
        int jkR = (fs->kR)*Nr + fs->jR;

        int iL = fs->iL0;
        int iR = fs->iR0;

        double *piphL = piph[jkL];
        double *piphR = piph[jkR];
        double *dphiL = dphi[jkL];
        double *dphiR = dphi[jkR];

        double *x = xf[fs->jkf];

        double dxL, dxR;
        if(dim == 1)
        {
            double rL = get_centroid(theDomain->r_jph[fs->jL],
                                     theDomain->r_jph[fs->jL-1], 1);
            double rR = get_centroid(theDomain->r_jph[fs->jR],
                                     theDomain->r_jph[fs->jR-1], 1);
            dxL = fs->rm - rL;
            dxR = rR - fs->rm;
        }
        else
        {
            double zL = get_centroid(theDomain->z_kph[fs->kL],
                                     theDomain->z_kph[fs->kL-1], 2);
            double zR = get_centroid(theDomain->z_kph[fs->kR],
                                     theDomain->z_kph[fs->kR-1], 2);
            dxL = fs->zm - zL;
            dxR = zR - fs->zm;
        }

        int f;
        for(f=fs->fa; f < fs->fb; f++)
        {
            double phi = x[3*f+1];
            double phiL = piphL[iL] - 0.5*dphiL[iL];
            double phiR = piphR[iR] - 0.5*dphiR[iR];
            double dpL = get_signed_dp(phi,phiL);
            double dpR = get_signed_dp(phi,phiR);

            int iqL = iL*NUM_Q;
            int iqR = iR*NUM_Q;

            double *primL = &(prim[jkL][iqL]);
            double *primR = &(prim[jkR][iqR]);
            
            double *gradpL = &(gradp[jkL][iqL]);
            double *gradpR = &(gradp[jkR][iqR]);

            double *gradL = &(grad[jkL][iqL]);
            double *gradR = &(grad[jkR][iqR]);
            
            for( q=0 ; q<NUM_Q ; ++q )
            {
                double WL = primL[q] + dpL*gradpL[q];
                double WR = primR[q] + dpR*gradpR[q];
            
                double S = (WR-WL)/(dxR+dxL);
                double SL = gradL[q];
                double SR = gradR[q];

                if( S*SL < 0.0 )
                    gradL[q] = 0.0; 
                else if( fabs(PLM*S) < fabs(SL) )
                    gradL[q] = PLM*S;

                if( S*SR < 0.0 )
                    gradR[q] = 0.0; 
                else if( fabs(PLM*S) < fabs(SR) )
                    gradR[q] = PLM*S;
            }

            if(get_signed_dp(piphL[iL], piphR[iR]) > 0)
                iR = (iR < Np[jkR]-1) ? iR+1 : 0;
            else
                iL = (iL < Np[jkL]-1) ? iL+1 : 0;
        }
    }
   
    /*
    for( n=0 ; n<Nf ; ++n )
    {
        struct face * f  = &( theFaces[n] );
        double phi = f->cm[1];
        int jkL = f->jkL;
        int iL = f->iL;
        int jkR = f->jkR;
        int iR = f->iR;
        double dxL = f->dxL;
        double dxR = f->dxR;
        double phiL = piph[jkL][iL] - 0.5*dphi[jkL][iL];
        double phiR = piph[jkR][iR] - 0.5*dphi[jkR][iR];
        double dpL = get_signed_dp(phi,phiL);
        double dpR = get_signed_dp(phi,phiR);
     
        int iqL = NUM_Q*iL;
        int iqR = NUM_Q*iR;

        double *primL = &(prim[jkL][iqL]);
        double *primR = &(prim[jkR][iqR]);
        
        double *gradpL = &(gradp[jkL][iqL]);
        double *gradpR = &(gradp[jkR][iqR]);
        
        double *gradL = &(grad[jkL][iqL]);
        double *gradR = &(grad[jkR][iqR]);
      
#if ENABLE_CART_INTERP
        double w = 0.0;
        if(Cartesian_Interp)
            w = getCartInterpWeight(f->cm);
        if(w > 0.0)
        {
            double gradCIL[NUM_Q], gradCIR[NUM_Q];
            geom_cart_interp_grad_trans(primL, primR,
                                        gradpL, gradpR,
                                        dpL, dpR, f->cm, -dxL, dxR,
                                        gradCIL, gradCIR, dim);

            for(q=0; q<NUM_Q; q++)
            {
                double WL = primL[q] + dpL*gradpL[q];
                double WR = primR[q] + dpR*gradpR[q];
                double S = (WR-WL)/(dxR+dxL);
                double SL = gradL[q];
                double SR = gradR[q];
                double SfL = (1-w) * S + w * gradCIL[q];
                double SfR = (1-w) * S + w * gradCIR[q];

                if(dim == 1 && q == UPP)
                {
                    double SL0 = -w*primL[UPP] / (f->cm[0] - dxL);
                    double SR0 = -w*primR[UPP] / (f->cm[0] + dxR);

                    SL -= SL0;
                    SR -= SR0;
                    
                    if( SfL*SL < 0.0 )
                        gradL[q] = SL0; 
                    else if( fabs(PLM*SfL) < fabs(SL) )
                        gradL[q] = PLM*SfL + SL0;

                    if( SfR*SR < 0.0 )
                        gradR[q] = SR0; 
                    else if( fabs(PLM*SfR) < fabs(SR) )
                        gradR[q] = PLM*SfR + SR0;
                }
                else
                {
                    if( SfL*SL < 0.0 )
                        gradL[q] = 0.0; 
                    else if( fabs(PLM*SfL) < fabs(SL) )
                        gradL[q] = PLM*SfL;
                    

                    if( SfR*SR < 0.0 )
                        gradR[q] = 0.0; 
                    else if( fabs(PLM*SfR) < fabs(SR) )
                        gradR[q] = PLM*SfR;
                }
            }
        }
        else
        {
#endif
            for( q=0 ; q<NUM_Q ; ++q )
            {
                double WL = primL[q] + dpL*gradpL[q];
                double WR = primR[q] + dpR*gradpR[q];

                double S = (WR-WL)/(dxR+dxL);
                double SL = gradL[q];
                double SR = gradR[q];
                if( S*SL < 0.0 )
                    gradL[q] = 0.0; 
                else if( fabs(PLM*S) < fabs(SL) )
                    gradL[q] = PLM*S;

                if( S*SR < 0.0 )
                    gradR[q] = 0.0; 
                else if( fabs(PLM*S) < fabs(SR) )
                    gradR[q] = PLM*S;
            }
#if ENABLE_CART_INTERP
        }
#endif
    }
*/
   
   // Geometric boundaries don't have ghost zones, so the gradients there need
   // to be fixed.  
 

   if(dim == 1 && theDomain->NgRa == 0)
      plm_geom_boundary(theDomain, 0, 0, 0, Nz-1, dim, 0);

   if(dim == 1 && theDomain->NgRb == 0)
      plm_geom_boundary(theDomain, Nr-1, Nr-1, 0, Nz-1, dim, 1);

   if(dim == 2 && theDomain->NgZa == 0)
      plm_geom_boundary(theDomain, 0, Nr-1, 0, 0, dim, 0);

   if(dim == 2 && theDomain->NgZb == 0)
      plm_geom_boundary(theDomain, 0, Nr-1, Nz-1, Nz-1, dim, 1);
}

void plm_geom_boundary(struct domain *theDomain, int jmin, int jmax, 
                        int kmin, int kmax, int dim, int LR)
{
    int Nr = theDomain->Nr;
    int * Np = theDomain->Np;
    double * r_jph = theDomain->r_jph;
    double * z_kph = theDomain->z_kph;
    double PLM = theDomain->theParList.PLM;

    int i,j,k;

    double **grad = (dim == 1) ? theDomain->gradr : theDomain->gradz;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;
    double **prim = theDomain->prim;

    double xp[3], xm[3];
    for(k=kmin; k<=kmax; k++)
    {
       xp[2] = z_kph[k];
       xm[2] = z_kph[k-1];
       for(j=jmin; j<=jmax; j++)
       {
           int jk = j+Nr*k;
           xp[0] = r_jph[j];
           xm[0] = r_jph[j-1];
           for(i=0; i < Np[jk]; i++)
           {
               xp[1] = piph[jk][i];
               xm[1] = piph[jk][i] - dphi[jk][i];

               int iq = NUM_Q*i;

               geom_grad(&(prim[jk][iq]), &(grad[jk][iq]),
                         xp, xm, PLM, dim, LR); 
           }
       }
    }
}
