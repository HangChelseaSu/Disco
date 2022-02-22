#include "paul.h"
#include "geometry.h"

double PHI_ORDER = 2.0;
static int grav2D = 0;
static int polar_sources_r = 0;
static int polar_sources_p = 0;
static int polar_sources_z = 0;

void setGravParams( struct domain * theDomain ){

   grav2D = theDomain->theParList.grav2D; 
   if(strcmp(GEOMETRY, "cylindrical") == 0
             && theDomain->theParList.NoBC_Rmin == 1)
   {
       polar_sources_r = 1;
       polar_sources_p = 1;
   }
   if(strcmp(GEOMETRY, "spherical") == 0)
   {
       //if(theDomain->theParList.NoBC_Rmin == 1)
       if(theDomain->theParList.NoBC_Rmin == 0)
       {
           polar_sources_r = 1;
           polar_sources_p = 1;
       }
       if(theDomain->theParList.NoBC_Zmin == 1
          || theDomain->theParList.NoBC_Zmax == 1)
       {
           polar_sources_p = 1;
           polar_sources_z = 1;
       }
   }

}

double phigrav( double M , double r , double eps , int type)
{
    if(type == PLPOINTMASS)
    {
        double n = PHI_ORDER;
        return( M/pow( pow(r,n) + pow(eps,n) , 1./n ) ) ;
    }
    else if(type == PLPW)
    {
        return M / (r - 2*M);
    }
    else if(type == PLWEGGC)
    {
        //potential "C" from Wegg 2012, ApJ 749
        double sq6 = sqrt(6);
        double Alpha = -4*(2.0+sq6)/3;
        double Rx = M*(4.0*sq6 - 9);
        double Ry = -4*M*(2*sq6 - 3.0)/3.0;
        return -Alpha*M/r - (1-Alpha)*M/(r-Rx) - M*Ry/(r*r);
    }
    else if(type == PLSURFACEGRAV)
    {
        return M*r; // M is gravitational acceleration
                    // only makes sense if grav2D is on
    }
    else if(type == PLSPLINE)
    {
        eps = eps*2.8;
        double u = r/eps;
        double val;
        if (u<0.5) val = 16.*u*u/3. - 48.*u*u*u*u/5. + 32.*u*u*u*u*u/5. - 14./5.;
        else if (u < 1.0) val = 1./(15.*u) + 32.*u*u/3. - 16.*u*u*u + 48.*u*u*u*u/5. - 32.*pow(u, 5.0)/15. - 3.2;
        else val = -1./u ;
        return -1*M*val/eps;
    }
    return 0.0;
}

double fgrav( double M , double r , double eps , int type)
{
    if(type == PLPOINTMASS)
    {
        double n = PHI_ORDER;
        return( M*pow(r,n-1.)/pow( pow(r,n) + pow(eps,n) ,1.+1./n) );
    }
    else if(type == PLPW)
    {
        return M / ((r-2*M)*(r-2*M));
    }
    else if(type == PLSURFACEGRAV)
    {
        return M; // M is gravitational acceleration
                  // only makes sense if grav2D is on
    }
    else if(type == PLWEGGC)
    {
        //potential "C" from Wegg 2012, ApJ 749
        double sq6 = sqrt(6);
        double Alpha = -4*(2.0+sq6)/3;
        double Rx = M*(4.0*sq6 - 9);
        double Ry = -4*M*(2*sq6 - 3.0)/3.0;
        return -Alpha*M/(r*r) - (1-Alpha)*M/((r-Rx)*(r-Rx)) - 2*M*Ry/(r*r*r);
    }
    else if(type == PLSPLINE)
    {
        eps = eps*2.8;
        double u = r/eps;
        double val;
        if (u<0.5) val = 32.*u/3. - 192.*u*u*u/5. + 32.*u*u*u*u;
        else if (u < 1.0) val = -1./(15.*u*u) + 64.*u/3. - 48.*u*u + 192.*u*u*u/5. - 32*u*u*u*u/3.;
        else val = 1./u/u;
        return M*val/eps/eps;
    }
    return 0.0;
    
}

void adjust_gas( struct planet * pl , double * x , double * prim , double gam ){

   double r   = x[0];
   double phi = x[1];

   double rp = pl->r;
   double pp = pl->phi;
   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy);

   //double z = M_PI*0.5; //x[2];
   double pot = phigrav( pl->M , script_r , pl->eps , pl->type);

   double c2 = gam*prim[PPP]/prim[RHO];
   double factor = 1. + (gam-1.)*pot/c2;

   prim[RHO] *= factor;
   prim[PPP] *= pow( factor , gam );

}

void planetaryForce( struct planet * pl , double r , double phi , double z , double * fr , double * fp , double * fz , int mode){

   double rp = pl->r;
   double pp = pl->phi;
   if(grav2D == 1)
      z = 0.0;
   else if(grav2D == 2)
   {
       rp = r;
       pp = phi;
       z = 1.0;
   }

   double cosp = cos(phi);
   double sinp = sin(phi);
   double dx = r*cosp-rp*cos(pp);
   double dy = r*sinp-rp*sin(pp);
   double script_r = sqrt(dx*dx+dy*dy+z*z);
   double script_r_perp = sqrt(dx*dx+dy*dy);

   double f1 = -fgrav( pl->M , script_r , pl->eps , pl->type);

   double cosa = dx/script_r_perp;
   double sina = dy/script_r_perp;

   double cosap = cosa*cosp+sina*sinp;
   double sinap = sina*cosp-cosa*sinp;

   if( mode==1 ){
      cosap = cosa*cos(pp)+sina*sin(pp);
      sinap = sina*cos(pp)-cosa*sin(pp);
   }

   if(script_r_perp <= 0.0)
   {
       cosap = 0.0;
       sinap = 0.0;
   }
   /*
   double rH = rp*pow( pl->M/3.,1./3.);
   double pd = 0.8; 
   double fd = 1./(1.+exp(-( script_r/rH-pd)/(pd/10.)));
   */

   double sint = script_r_perp/script_r;
   double cost = z/script_r;

   *fr = cosap*f1*sint; //*fd;
   *fp = sinap*f1*sint; //*fd;
   *fz = f1*cost;
}

void planet_src( struct planet * pl , double * prim , double * cons , double * xp , double * xm , double dVdt ){

   double rho = prim[RHO];
   double vr  = prim[URR];
   double vz  = prim[UZZ];
   double omega = prim[UPP];
   
   //double rp = xp[0];
   //double rm = xm[0];
   //double r = 0.5*(rp+rm);
   
   double dphi = get_dp(xp[1],xm[1]);
   double phi = xm[1] + 0.5*dphi;
   double x[3] = {get_centroid(xp[0],xm[0],1), phi, 
                    get_centroid(xp[2],xm[2],2)};
   double xcyl[3];
   get_rpz(x, xcyl);
   
   double Fcyl[3], F[3];
   planetaryForce( pl, xcyl[0], xcyl[1], xcyl[2],
                    &(Fcyl[0]), &(Fcyl[1]), &(Fcyl[2]), 0);
   get_vec_from_rpz(x, Fcyl, F);

   double hr = get_scale_factor(x, 1);
   double hp = get_scale_factor(x, 0);
   double hz = get_scale_factor(x, 2);

   if(polar_sources_r || polar_sources_p || polar_sources_z)
   {
       double adjust[3];
       geom_polar_vec_adjust(xp, xm, adjust);
       if(polar_sources_r)
           F[0] *= adjust[0];
       if(polar_sources_p)
           F[1] *= adjust[1];
       if(polar_sources_z)
           F[2] *= adjust[2];
   }

   cons[SRR] += rho*hr*F[0]*dVdt;
   cons[LLL] += rho*hp*F[1]*dVdt;
   cons[SZZ] += rho*hz*F[2]*dVdt;
   cons[TAU] += rho*( hr*F[0]*vr + hz*F[2]*vz + hp*F[1]*omega )*dVdt;

}

void planet_RK_copy( struct planet * pl ){
   pl->RK_r     = pl->r;
   pl->RK_phi   = pl->phi;
   pl->RK_M     = pl->M;
   pl->RK_omega = pl->omega;
   pl->RK_vr    = pl->vr;
   pl->RK_dM    = pl->dM;
   pl->RK_accL  = pl->accL;
   pl->RK_Ls    = pl->Ls;
   pl->RK_therm = pl->therm;
   pl->RK_gravL = pl->gravL;
   pl->RK_accE = pl->accE;
   pl->RK_gravE = pl->gravE;

}

void planet_RK_adjust_kin( struct planet * pl , double RK ){
   pl->r     = (1.-RK)*pl->r     + RK*pl->RK_r;
   pl->phi   = (1.-RK)*pl->phi   + RK*pl->RK_phi;
   pl->M     = (1.-RK)*pl->M     + RK*pl->RK_M;
   pl->omega = (1.-RK)*pl->omega + RK*pl->RK_omega;
   pl->vr    = (1.-RK)*pl->vr    + RK*pl->RK_vr;
}

void planet_RK_adjust_aux( struct planet * pl , double RK ){
   pl->dM    = (1.-RK)*pl->dM    + RK*pl->RK_dM;
   pl->accL  = (1.-RK)*pl->accL  + RK*pl->RK_accL;
   pl->Ls    = (1.-RK)*pl->Ls    + RK*pl->RK_Ls;
   pl->therm = (1.-RK)*pl->therm + RK*pl->RK_therm;
   pl->gravL = (1.-RK)*pl->gravL + RK*pl->RK_gravL;
   pl-> accE = (1.-RK)*pl->accE  + RK*pl->RK_accE;
   pl->sinkE = (1.-RK)*pl->sinkE + RK*pl->RK_sinkE;
   pl->gravE = (1.-RK)*pl->gravE + RK*pl->RK_gravE;
}
