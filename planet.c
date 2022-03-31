#include "paul.h"
#include "geometry.h"
#include "planet.h"

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
       if(theDomain->theParList.NoBC_Rmin == 1)
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
    else if(type == PLQUAD)
    {
        return 0.5*M*r*r;
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
    else if(type == PLQUAD)
    {
        return M * r;
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

void planetaryForce( struct planet * pl , double r , double phi , double z , double * fr , double * fp , double * fz , int mode )
{
    /*
     * Calculates the specific gravitational force (ie. acceleration) between
     * the gas and planet.
     * The force is returned in an orthonormal basis in cylindrical coords:
     *    |F| = sqrt(fr*fr + fp*fp + fz*fz)
     * mode == 0 computes the force components At The Gas Position
     * mode == 1 computes the force components At The Planet Position
     */

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

   //TODO: WAY TOO MANY sincos calls!
   double Fp_xyz[3];
   get_vec_xyz(x, F, Fp_xyz);
   double cosp = cos(pl->phi);
   double sinp = sin(pl->phi);
   double Fp[3] = {cosp*Fp_xyz[0] + sinp*Fp_xyz[1],
                   -sinp*Fp_xyz[0] + cosp*Fp_xyz[1], Fp_xyz[2]};

   pl->gravL -= rho*(pl->r)*Fp[1]*dVdt;
   pl->gravE -= rho*(pl->vr*Fp[0] + pl->omega*pl->r*Fp[1])*dVdt;

   pl->gas_track[PL_GRV_PX] -= rho*Fp_xyz[0]*dVdt;
   pl->gas_track[PL_GRV_PY] -= rho*Fp_xyz[1]*dVdt;
   pl->gas_track[PL_GRV_PZ] -= rho*Fp_xyz[2]*dVdt;
   pl->gas_track[PL_GRV_JZ] -= rho*(pl->r)*Fp[1]*dVdt;
   pl->gas_track[PL_GRV_EGAS] -= rho*(hr*F[0]*vr + hz*F[2]*vz
                                      + hp*F[1]*omega )*dVdt;
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
    
    int q;
    for(q=0; q<NUM_PL_AUX; q++)
        pl->RK_aux[q] = pl->aux[q];
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
   pl->gravE = (1.-RK)*pl->gravE + RK*pl->RK_gravE;

    int q;
    for(q=0; q<NUM_PL_AUX; q++)
        pl->aux[q] = (1-RK)*pl->aux[q] + RK*pl->RK_aux[q];
}

void planet_zero_aux(struct planet *pl)
{
    pl->dM = 0.0;
    pl->Ls = 0.0;
    pl->accL = 0.0;
    pl->gravL = 0.0;
    pl->therm = 0.0;
    pl->accE = 0.0;
    pl->gravE = 0.0;
    
    pl->RK_dM = 0.0;
    pl->RK_Ls = 0.0;
    pl->RK_accL = 0.0;
    pl->RK_gravL = 0.0;
    pl->RK_therm = 0.0;
    pl->RK_accE = 0.0;
    pl->RK_gravE = 0.0;

    int q;
    for(q=0; q<NUM_PL_AUX; q++)
    {
        pl->aux[q] = 0.0;
        pl->RK_aux[q] = 0.0;
    }
}

void setupPlanets(struct domain *theDomain)
{
   initializePlanets( theDomain->thePlanets );
  
   int Npl = theDomain->Npl;
   int p;
   for(p=0; p<Npl; p++)
       planet_zero_aux(theDomain->thePlanets + p);

   initializePlanetTracking(theDomain);
}

void initializePlanetTracking(struct domain *theDomain)
{
    int Npl = theDomain->Npl;
    int p, q;
    for(p=0; p<Npl; p++)
    {
        struct planet *pl = theDomain->thePlanets + p;

        for(q=0; q<NUM_PL_INTEGRALS; q++)
            pl->gas_track[q] = 0.0;
    }
}

void updatePlanetTracking(struct domain *theDomain)
{
    int p;
    int Npl = theDomain->Npl;
    for(p=0; p<Npl; p++)
    {
        struct planet *pl = theDomain->thePlanets + p;

        //Direct values, easy!
        pl->aux[PL_SNK_M] += pl->gas_track[PL_SNK_M];
        pl->aux[PL_GRV_PX] += pl->gas_track[PL_GRV_PX];
        pl->aux[PL_GRV_PY] += pl->gas_track[PL_GRV_PY];
        pl->aux[PL_GRV_PZ] += pl->gas_track[PL_GRV_PZ];
        pl->aux[PL_GRV_JZ] += pl->gas_track[PL_GRV_JZ];
        pl->aux[PL_SNK_PX] += pl->gas_track[PL_SNK_PX];
        pl->aux[PL_SNK_PY] += pl->gas_track[PL_SNK_PY];
        pl->aux[PL_SNK_PZ] += pl->gas_track[PL_SNK_PZ];
        pl->aux[PL_SNK_JZ] += pl->gas_track[PL_SNK_JZ];
        pl->aux[PL_SNK_SZ] += pl->gas_track[PL_SNK_SZ];
        pl->aux[PL_SNK_X] += pl->gas_track[PL_SNK_X] / pl->M;
        pl->aux[PL_SNK_Y] += pl->gas_track[PL_SNK_Y] / pl->M;
        pl->aux[PL_SNK_Z] += pl->gas_track[PL_SNK_Z] / pl->M;
        pl->aux[PL_GRV_EGAS] += pl->gas_track[PL_GRV_EGAS];
        pl->aux[PL_SNK_EGAS] += pl->gas_track[PL_SNK_EGAS];

        // Some more complicated things.

        double r = pl->r;
        double cosp = cos(pl->phi);
        double sinp = sin(pl->phi);
        double x = r * cosp;
        double y = r * sinp;
        double vx = pl->vr * cosp - r * pl->omega * sinp;
        double vy = pl->vr * sinp + r * pl->omega * cosp;
        pl->aux[PL_SNK_LZ] += x * pl->gas_track[PL_SNK_PY]
                            - y * pl->gas_track[PL_SNK_PX]
                            + pl->gas_track[PL_SNK_X] * vy
                            - pl->gas_track[PL_SNK_Y] * vx;

        pl->aux[PL_GRV_K] += vx * pl->gas_track[PL_GRV_PX]
                            + vy * pl->gas_track[PL_GRV_PY];
        pl->aux[PL_SNK_K] += vx * pl->gas_track[PL_SNK_PX]
                            + vy * pl->gas_track[PL_SNK_PY]
                            - 0.5*(vx*vx + vy*vy)*pl->gas_track[PL_SNK_M];
        
        pl->aux[PL_GRV_U] += pl->gas_track[PL_GRV_EGAS]
                            - vx * pl->gas_track[PL_GRV_PX]
                            - vy * pl->gas_track[PL_GRV_PY];

        double xyz[3] = {x, y, 0};
        double X[3];
        get_coords_from_xyz(X, xyz);

        double Phi_ext = 0.0;
        double gx_ext = 0;
        double gy_ext = 0;
        double gz_ext = 0;
        int p2;
        for(p2=0; p2<Npl; p2++)
        {
            if(p == p2)
                continue;
        
            struct planet *pl2 = theDomain->thePlanets + p2;
            double r2 = pl2->r;
            double cosp2 = cos(pl2->phi);
            double sinp2 = sin(pl2->phi);
            double x2 = r2*cosp2;
            double y2 = r2*sinp2;

            double rsep = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2));

            Phi_ext += -phigrav(pl2->M, rsep, pl2->eps, pl2->type);
            double g = fgrav(pl2->M, rsep, pl2->eps, pl2->type);

            gx_ext += g * (x2-x) / rsep;
            gy_ext += g * (y2-y) / rsep;
            gz_ext += 0;
        }
        
        pl->aux[PL_SNK_U] += Phi_ext * pl->gas_track[PL_SNK_M]
                            - gx_ext * pl->gas_track[PL_SNK_X]
                            - gy_ext * pl->gas_track[PL_SNK_Y]
                            - gz_ext * pl->gas_track[PL_SNK_Z];
    }
}
