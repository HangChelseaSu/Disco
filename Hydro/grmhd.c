
#include "../paul.h"
#include "metric.h"
#include "frame.h"

#define DEBUG 1
#define DEBUG2 0
#define DEBUG3 0
#define DEBUG_RMAX 5.5
#define DEBUG_ZMAX 3.5
#define ND 3

//Global Functions
double get_cs2( double );
double get_dp( double , double );
double get_dL( double * , double * , int );

//Local Functions
void cons2prim_prep(double *cons, double *x);
void cons2prim_solve_isothermal(double *cons, double *prim, double *x);
void cons2prim_solve_adiabatic(double *cons, double *prim, double *x);
void cons2prim_finalize(double *prim, double *x);

static double gamma_law = 0.0; 
static double RHO_FLOOR = 0.0; 
static double PRE_FLOOR = 0.0; 
static int isothermal = 0;

void setHydroParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
   isothermal = theDomain->theParList.isothermal_flag;
   RHO_FLOOR = theDomain->theParList.Density_Floor;
   PRE_FLOOR = theDomain->theParList.Pressure_Floor;
}

int set_B_flag(void){
   return(1);
}

double get_omega(double *prim, double *x)
{
    int i,j;
    double r = x[0];
    double l[3] = {prim[URR], r*r*prim[UPP], prim[UZZ]};
    
    double lapse;
    double shift[3];
    double igam[9];

    lapse = metric_lapse(x);
    metric_shift(x, shift);
    metric_igam(x, igam);

    double u2 = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            u2 += igam[3*i+j]*l[i]*l[j];
    double w = sqrt(1.0 + u2);
    double vp = lapse*(igam[3]*l[0]+igam[4]*l[1]+igam[5]*l[2])/w - shift[1];
    
    return vp;
}

void prim2cons( double *prim, double *cons, double *x, double dV)
{
    double r = x[0];
    double al, be[3], gam[9], igam[9], jac;
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    double sqrtgam = jac / al;

    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[3] = {prim[URR], r*r*prim[UPP], prim[UZZ]};
    double B[3] = {prim[BRR]/sqrtgam,prim[BPP]/(r*sqrtgam),prim[BZZ]/sqrtgam};

    double U[4];
    frame_U(x, U);

    double w, u0, u2;
    double u[3];
    double igaml[3];
    int i,j;
    for(i=0; i<3; i++)
    {
        igaml[i] = 0.0;
        for(j=0; j<3; j++)
            igaml[i] += igam[3*i+j]*l[j];
    }
    u2 = l[0]*igaml[0] + l[1]*igaml[1] + l[2]*igaml[2];
    w = sqrt(1.0 + u2);
    u0 = w/al;
    for(i=0; i<3; i++)
        u[i] = igaml[i] - be[i]*u0;

    double l0 = -al*w + be[0]*l[0] + be[1]*l[1] + be[2]*l[2];
    double uU = U[0]*l0 + U[1]*l[0] + U[2]*l[1] + U[3]*l[2];

    double uB, b0, bd[3], B2, b2, UB, Ub;
    uB = l[0]*B[0] + l[1]*B[1] + l[2]*B[2];
    b0 = uB / al;
    bd[0] = (gam[0]*B[0] + gam[1]*B[1] + gam[2]*B[2] + l[0]*uB) / w;
    bd[1] = (gam[3]*B[0] + gam[4]*B[1] + gam[5]*B[2] + l[1]*uB) / w;
    bd[2] = (gam[6]*B[0] + gam[7]*B[1] + gam[8]*B[2] + l[2]*uB) / w;
    B2 = 0.0;
    UB = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {    
            B2 += gam[3*i+j]*B[i]*B[j];
            UB += gam[3*i+j] * B[i] * (U[0]*be[j] + U[j+1]);
        }
    b2 = (B2 + uB*uB) / (w*w);
    Ub = (UB + uU*uB) / w;
    
    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp + b2;
    double rhoe = Pp / (gamma_law-1.0);

    if(DEBUG2)
    {
        printf("B2 = %.12lg, uB = %.12lg, UB = %.12lg\n", B2, uB, UB);
        printf("bd[0] = %.12lg, bd[1] = %.12lg, bd[2] = %.12lg\n", 
                    bd[0], bd[1], bd[2]);
        printf("b2 = %.12lg, Ub = %.12lg, b0 = %.12lg\n", b2, Ub, b0);
    }

    cons[DDD] = jac * rho*u0 * dV;
    cons[SRR] = jac * (rhoh*u0*l[0] - b0*bd[0]) * dV;
    cons[LLL] = jac * (rhoh*u0*l[1] - b0*bd[1]) * dV;
    cons[SZZ] = jac * (rhoh*u0*l[2] - b0*bd[2]) * dV;
    cons[TAU] = jac * (-(rhoe+0.5*b2)*uU*u0 - (Pp+0.5*b2)*(uU*u0+U[0])
                            - rho*(uU+1.0)*u0 + Ub*b0) * dV;

    cons[BRR] = sqrtgam * B[0]/r * dV;
    cons[BPP] = sqrtgam * B[1] * dV;
    cons[BZZ] = sqrtgam * B[2] * dV;

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        cons[q] = prim[q]*cons[DDD];
  
    if(DEBUG3)
    {
        FILE *f = fopen("p2c.out", "a");
        fprintf(f, "%.10lg %.10lg %.10lg %.10lg %.10lg %.10lg\n",
                    x[0], x[1], prim[URR], prim[UPP], cons[SRR], cons[LLL]);
        fprintf(f, "    %.10lg %.10lg %.10lg\n",
                        B2, prim[BRR], prim[BPP]);
        fclose(f);
    }
}

void getUstar(double *prim, double *Ustar, double *x, double Sk, double Ss, 
                double *n, double *Bpack)
{
    double r = x[0];
    double rho = prim[RHO];
    double Pp = prim[PPP];
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};

    double al, be[3], igam[9], jac, U[4];
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    frame_U(x, U);

    double bn = n[0]*be[0] + n[1]*be[1] + n[2]*be[2];
    double ign = n[0]*igam[0] + n[1]*igam[4] + n[2]*igam[8];
    double Un = n[0]*U[1] + n[1]*U[2] + n[2]*U[3];
    double hn = n[0] + r*n[1] + n[2];

    double ss = Ss/hn;
    double sk = Sk/hn;

    double uS[3];
    double u2 = 0.0;
    int i,j;
    for(i=0; i<3; i++)
    {
        uS[i] = 0.0;;
        for(j=0; j<3; j++)
            uS[i] += igam[3*i+j]*l[j];
        u2 += l[i]*uS[i];
    }

    double w = sqrt(1+u2);
    double u0 = w/al;
    double l0 = -al*w + be[0]*l[0] + be[1]*l[1] + be[2]*l[2];
    double vn = al * (n[0]*uS[0] + n[1]*uS[1] + n[2]*uS[2]) / w - bn;
    double uU = l0*U[0] + l[0]*U[1] + l[1]*U[2] + l[2]*U[3];

    // q == F - s * U
    double rhoh = rho + gamma_law/(gamma_law-1.0) * Pp;
    double qE = rhoh*w*w*(vn-sk) + Pp*(sk+bn);

    double mn = rhoh*w * (n[0]*uS[0] + n[1]*uS[1] + n[2]*uS[2]);
    double qMn = mn*(vn-sk) + al*Pp*ign;

    double ssS = (ss+bn)/al;
    double skS = (sk+bn)/al;

    // P star!
    double Pstar = (qMn - ssS*qE) / (al*(ign - ssS*skS));

    double kappa = (vn - sk) / (ss - sk);
    double alpha1 = (Pp - Pstar) / (ss - sk);
    double alpha2 = (vn*Pp - ss*Pstar) / (ss - sk);

    double rhoe = Pp / (gamma_law - 1.0);
    double tau = -rhoe*uU*u0 - Pp*(u0*uU+U[0]) - rho*(uU+1.0)*u0;

    Ustar[DDD] = jac * rho*u0 * kappa;
    Ustar[SRR] = jac * (rhoh*u0*l[0] * kappa + alpha1 * n[0]);
    Ustar[LLL] = jac * (rhoh*u0*l[1] * kappa + alpha1 * n[1]);
    Ustar[SZZ] = jac * (rhoh*u0*l[2] * kappa + alpha1 * n[2]);
    Ustar[TAU] = jac * (tau * kappa - Un*alpha1 + U[0]*alpha2);

    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        Ustar[q] = prim[q]*Ustar[DDD];
}

void cons2prim(double *cons, double *prim, double *x, double dV)
{
    int q;
    double cons1[NUM_Q];
    for(q=0; q<NUM_Q; q++)
        cons1[q] = cons[q]/dV;

    cons2prim_prep(cons1, x);
    if(isothermal)
        cons2prim_solve_isothermal(cons1, prim, x);
    else
        cons2prim_solve_adiabatic(cons1, prim, x);
    cons2prim_finalize(prim, x);
}

void flux(double *prim, double *flux, double *x, double *n)
{
    double r = x[0];
    double al, be[3], gam[9], igam[9], jac, sqrtgam;
    double U[4];
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    sqrtgam = jac / al;
    frame_U(x, U);

    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[3] = {prim[URR], r*r*prim[UPP], prim[UZZ]};
    double B[3] = {prim[BRR]/sqrtgam,prim[BPP]/(r*sqrtgam),prim[BZZ]/sqrtgam};

    double w, u0, u2, u[3], v[3];
    double igaml[3];
    int i,j;
    for(i=0; i<3; i++)
    {
        igaml[i] = 0.0;
        for(j=0; j<3; j++)
            igaml[i] += igam[3*i+j]*l[j];
    }
    u2 = l[0]*igaml[0] + l[1]*igaml[1] + l[2]*igaml[2];
    w = sqrt(1.0 + u2);
    u0 = w/al;
    for(i=0; i<3; i++)
        u[i] = igaml[i] - be[i]*u0;
    for(i=0; i<3; i++)
        v[i] = u[i]/u0;
    
    double l0 = -al*w + be[0]*l[0] + be[1]*l[1] + be[2]*l[2];
    double uU = U[0]*l0 + U[1]*l[0] + U[2]*l[1] + U[3]*l[2];

    double uB, b0, bd[3], B2, b2, UB, Ub;
    uB = l[0]*B[0] + l[1]*B[1] + l[2]*B[2];
    b0 = uB / al;
    bd[0] = (gam[0]*B[0] + gam[1]*B[1] + gam[2]*B[2] + l[0]*uB) / w;
    bd[1] = (gam[3]*B[0] + gam[4]*B[1] + gam[5]*B[2] + l[1]*uB) / w;
    bd[2] = (gam[6]*B[0] + gam[7]*B[1] + gam[8]*B[2] + l[2]*uB) / w;
    B2 = 0.0;
    UB = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {    
            B2 += gam[3*i+j]*B[i]*B[j];
            UB += gam[3*i+j] * B[i] * (U[0]*be[j] + U[j+1]);
        }
    b2 = (B2 + uB*uB) / (w*w);
    Ub = (UB + uU*uB) / w;

    double un = u[0]*n[0] + u[1]*n[1] + u[2]*n[2];
    double Bn = B[0]*n[0] + B[1]*n[1] + B[2]*n[2];
    double Un = U[1]*n[0] + U[2]*n[1] + U[3]*n[2];
    double vn = v[0]*n[0] + v[1]*n[1] + v[2]*n[2];
    double hn = n[0] + r*n[1] + n[2];
    double bn = (Bn + uB*un)/w;
    
    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp + b2;
    double rhoe = Pp / (gamma_law-1.0);
    double Pt = Pp + 0.5*b2;

    flux[DDD] = jac * hn * rho*un;
    flux[SRR] = jac * hn * (rhoh*un*l[0] + Pt*n[0] - bd[0]*bn);
    flux[LLL] = jac * hn * (rhoh*un*l[1] + Pt*n[1] - bd[1]*bn);
    flux[SZZ] = jac * hn * (rhoh*un*l[2] + Pt*n[2] - bd[2]*bn);
    flux[TAU] = jac * hn * (-(rhoe+0.5*b2)*uU*un - Pt*(uU*un+Un)
                            - rho*(uU+1.0)*un + Ub*bn);
    
    flux[BRR] = sqrtgam * hn * (B[0]*vn - Bn*v[0])/r;
    flux[BPP] = sqrtgam * hn * (B[1]*vn - Bn*v[1]);
    flux[BZZ] = sqrtgam * hn * (B[2]*vn - Bn*v[2]);
    
    int q;
    for(q = NUM_C; q < NUM_Q; q++)
        flux[q] = prim[q]*flux[DDD];
   
    if(DEBUG3)
    {
        FILE *f = fopen("flux.out", "a");
        fprintf(f, "%.1lf %.1lf %.10lg %.10lg %.10lg %.10lg\n",
                n[0], n[1], x[0], x[1], flux[SRR], flux[LLL]);
        fprintf(f, "    %.10lg %.10lg %.10lg %.10lg %.10lg\n",
                        B2, bn, prim[BRR], prim[BPP], w);
        fclose(f);
    }
}

void source(double *prim, double *cons, double *xp, double *xm, double dVdt)
{
    double x[3] = {0.5*(xm[0]+xp[0]), 0.5*(xm[1]+xp[1]), 0.5*(xm[2]+xp[2])};
    int i,j,mu,nu;
    double al, be[3], gam[9], igam[9], ig[16], jac, sqrtgam;
    double U[4], dU[16];
    double r = x[0];
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    frame_U(x, U);
    frame_der_U(x, dU);
    sqrtgam = jac / al;

    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[4] = {0.0, prim[URR], r*r*prim[UPP], prim[UZZ]};
    double B[3] = {prim[BRR]/sqrtgam,prim[BPP]/(r*sqrtgam),prim[BZZ]/sqrtgam};

    double ia2 = 1.0/(al*al);
    ig[0] = -ia2;
    for(mu=0; mu<3; mu++)
    {
        ig[mu+1] = be[mu]*ia2;
        ig[4*(mu+1)] = ig[mu+1];
        for(nu=0; nu<3; nu++)
            ig[4*(mu+1)+nu+1] = igam[3*mu+nu]-be[mu]*be[nu]*ia2;
    }

    double w, u[4], u2;
    double igaml[3];
    igaml[0] = igam[0]*l[1] + igam[1]*l[2] + igam[2]*l[3];
    igaml[1] = igam[3]*l[1] + igam[4]*l[2] + igam[5]*l[3];
    igaml[2] = igam[6]*l[1] + igam[7]*l[2] + igam[8]*l[3];
    u2 = l[1]*igaml[0] + l[2]*igaml[1] + l[3]*igaml[2];
    w = sqrt(1.0 + u2);
    
    u[0] = w/al;
    for(i=0; i<3; i++)
        u[i+1] = igaml[i] - be[i]*u[0];
    l[0] = -al*w + be[0]*l[1] + be[1]*l[2] + be[2]*l[3];
    
    double uB, b[4], bd[4], B2, b2, beB;
    uB = l[1]*B[0] + l[2]*B[1] + l[3]*B[2];
    b[0] = uB / al;
    b[1] = (B[0]+u[1]*uB)/w;
    b[2] = (B[1]+u[2]*uB)/w;
    b[3] = (B[2]+u[3]*uB)/w;
    bd[1] = (gam[0]*B[0] + gam[1]*B[1] + gam[2]*B[2] + l[1]*uB) / w;
    bd[2] = (gam[3]*B[0] + gam[4]*B[1] + gam[5]*B[2] + l[2]*uB) / w;
    bd[3] = (gam[6]*B[0] + gam[7]*B[1] + gam[8]*B[2] + l[3]*uB) / w;
    B2 = 0.0;
    beB = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {
            beB += gam[3*i+j]*be[i]*B[j];
            B2 += gam[3*i+j]*B[i]*B[j];
        }
    bd[0] = (beB + l[0]*uB) / w;
    b2 = (B2 + uB*uB) / (w*w);

    double rhoh = rho + gamma_law/(gamma_law-1.0)*Pp + b2;
    double Pt = Pp + 0.5*b2;

    double S0, Sk[3];

    double rp = xp[0];
    double rm = xm[0];
    double r2_3 = (rp*rp + rp*rm + rm*rm)/3.0;
    double dphi = get_dp(xp[1],xm[1]);
    double rcorr = r2_3/(r*r);
    double pcorr = sin(0.5*dphi)/(0.5*dphi);

    for(i=0; i<3; i++)
    {
        Sk[i] = 0.0;
        if(metric_killing(i+1))
            continue;

        double dg[16];
        metric_der_g(x, i+1, dg);
        for(mu=0; mu<4; mu++)
            for(nu=0; nu<4; nu++)
                Sk[i] += (rhoh*u[mu]*u[nu]*rcorr*pcorr + ig[4*mu+nu]*Pt - b[mu]*b[nu]*pcorr)
                            * dg[4*mu+nu];
        Sk[i] *= 0.5;
    }
    S0 = -U[1]*Sk[0] - U[2]*Sk[1] - U[3]*Sk[2];

    for(mu=0; mu<=ND; mu++)
        for(nu=0; nu<=ND; nu++)
        {
            if(mu == nu)
                S0 += -(rhoh*u[mu]*l[nu] + Pt - b[mu]*bd[nu]) * dU[4*mu+nu];
            else
                S0 += -(rhoh*u[mu]*l[nu] - b[mu]*bd[nu]) * dU[4*mu+nu];
        }

    cons[SRR] += jac * Sk[0] * dVdt;
    cons[LLL] += jac * Sk[1] * dVdt;
    cons[SZZ] += jac * Sk[2] * dVdt;
    cons[TAU] += jac * S0 * dVdt;

    if(DEBUG3)
    {
        FILE *f = fopen("source.out", "a");
        fprintf(f, "%.12lg %.12lg %.12lg %.12lg\n",
                x[0], x[1], jac*Sk[0]*dVdt, jac*S0*dVdt);
        fprintf(f, "    %.10lg %.10lg %.10lg %.10lg\n",
                        B2, prim[BRR], prim[BPP], w);
        fclose(f);
    }
}

void visc_flux(double *prim, double *gprim, double *flux, double *x, 
                double *n){}

void flux_to_E(double *Flux, double *Ustr, double *x, double *E1_riemann, 
                double *B1_riemann, double *E2_riemann, double *B2_riemann, 
                int dim)
{
   double r = x[0];

   if( dim==0 )  //PHI
   {
      *E1_riemann = Flux[BRR]*r;   //Ez 
      *B1_riemann = Ustr[BRR]*r*r; // r*Br
      *E2_riemann = Flux[BZZ];    //Er 
      *B2_riemann = Ustr[BZZ]*r;  //-r*Bz
   }
   else if( dim==1 ) //RRR
   {
      *E1_riemann = -Flux[BPP]*r;  //Ez 
      *B1_riemann = Ustr[BRR]*r*r; // r*Br
      *E2_riemann = 1.0*Flux[BZZ];     //Ephi
   }
   else //ZZZ
   {
      *E1_riemann = -Flux[BPP]*r;   //Er 
      *B1_riemann = Ustr[BZZ]*r;  //-r*Bz
      *E2_riemann = 1.0*-Flux[BRR]*r;  //Ephi
   }
}
void vel(double *prim1, double *prim2, double *Sl, double *Sr, double *Ss, 
            double *n, double *x, double *Bpack)
{
    double r = x[0];
    double al, be[3], gam[9], igam[9], sqrtgam;
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    sqrtgam = metric_jacobian(x) / (al * r);

    double rho1 = prim1[RHO];
    double P1   = prim1[PPP];
    double l1[3] = {prim1[URR], r*r*prim1[UPP], prim1[UZZ]};
    double B1[3] = {prim1[BRR]/sqrtgam,prim1[BPP]/(r*sqrtgam),
                    prim1[BZZ]/sqrtgam};

    double cs21 = gamma_law*P1/(rho1+gamma_law/(gamma_law-1.0)*P1);

    double rho2 = prim2[RHO];
    double P2   = prim2[PPP];
    double l2[3] = {prim2[URR], r*r*prim2[UPP], prim2[UZZ]};
    double B2[3] = {prim2[BRR]/sqrtgam,prim2[BPP]/(r*sqrtgam),
                    prim2[BZZ]/sqrtgam};

    double cs22 = gamma_law*P2/(rho2+gamma_law/(gamma_law-1.0)*P2);

    int i,j;
    double u21 = 0.0;
    double u22 = 0.0;
    double uS1[3], uS2[3];
    for(i=0; i<3; i++)
    {
        uS1[i] = 0.0;
        uS2[i] = 0.0;
        for(j=0; j<3; j++)
        {
            uS1[i] += igam[3*i+j]*l1[j];
            uS2[i] += igam[3*i+j]*l2[j];
        }
        u21 += uS1[i]*l1[i];
        u22 += uS2[i]*l2[i];
    }

    double w1 = sqrt(1.0+u21);
    double w2 = sqrt(1.0+u22);
    double v21 = u21/(w1*w1);
    double v22 = u22/(w2*w2);
    
    double uB1, uB2, B21, B22, b21, b22;
    uB1 = l1[0]*B1[0] + l1[1]*B1[1] + l1[2]*B1[2];
    uB2 = l2[0]*B2[0] + l2[1]*B2[1] + l2[2]*B2[2];
    B21 = 0.0;
    B22 = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
        {
            B21 += gam[3*i+j]*B1[i]*B1[j];
            B22 += gam[3*i+j]*B2[i]*B2[j];
        }
    b21 = (B21 + uB1*uB1) / (w1*w1);
    b22 = (B22 + uB2*uB2) / (w2*w2);

    double va21 = b21 / (rho1 + gamma_law/(gamma_law-1.0)*P1 + b22);
    double va22 = b22 / (rho2 + gamma_law/(gamma_law-1.0)*P2 + b22);

    double cf21 = cs21 + va21 - cs21*va21;
    double cf22 = cs22 + va22 - cs22*va22;

    //TODO: Use n[] PROPERLY.  This only works if n = (1,0,0) or some 
    //      permutation.
    double vSn1 = (uS1[0]*n[0]+uS1[1]*n[1]+uS1[2]*n[2]) / w1;
    double vSn2 = (uS2[0]*n[0]+uS2[1]*n[1]+uS2[2]*n[2]) / w2;
    double bn = (be[0]*n[0]+be[1]*n[1]+be[2]*n[2]);
    double ign = igam[3*0+0]*n[0] + igam[3*1+1]*n[1] + igam[3*2+2]*n[2];

    double dv1 = sqrt(cf21*(ign - vSn1*vSn1 - cf21*(ign*v21-vSn1*vSn1))) / w1;
    double dv2 = sqrt(cf22*(ign - vSn2*vSn2 - cf22*(ign*v22-vSn2*vSn2))) / w2;
    double hn = n[0] + r*n[1] + n[2];

    double sl1 = hn * (al * (vSn1*(1.0-cf21) - dv1) / (1.0-v21*cf21) - bn);
    double sr1 = hn * (al * (vSn1*(1.0-cf21) + dv1) / (1.0-v21*cf21) - bn);
    double sl2 = hn * (al * (vSn2*(1.0-cf22) - dv2) / (1.0-v22*cf22) - bn);
    double sr2 = hn * (al * (vSn2*(1.0-cf22) + dv2) / (1.0-v22*cf22) - bn);

/*
    printf("cs2(L/R): %.12lg %.12lg\n", cs21, cs22);
    printf("vSn(L/R): %.12lg %.12lg\n", vSn1, vSn2);
    printf("w  (L/R): %.12lg %.12lg\n", w1, w2);
    printf("u2 (L/R): %.12lg %.12lg\n", u21, u22);
    printf("v2 (L/R): %.12lg %.12lg\n", v21, v22);
    printf("sl (L/R): %.12lg %.12lg\n", sl1, sl2);
    printf("sr (L/R): %.12lg %.12lg\n", sr1, sr2);
*/

    *Sr = sr1 > sr2 ? sr1 : sr2;
    *Sl = sl1 < sl2 ? sl1 : sl2;

    //double maxv = fabs(*Sr) > fabs(*Sl) ? fabs(*Sr) : fabs(*Sl);
    //*Sr = maxv;
    //*Sl = -maxv;

    //Now for the contact wave speed.
    double sL = *Sl / hn;
    double sR = *Sr / hn;

    double rhohL = rho1 + gamma_law/(gamma_law-1)*P1;
    double rhohR = rho2 + gamma_law/(gamma_law-1)*P2;
    double vnL = al*vSn1 - bn;
    double vnR = al*vSn2 - bn;

    double ML[3] = {rhohL*w1*l1[0], rhohL*w1*l1[1], rhohL*w1*l1[2]};
    double MR[3] = {rhohR*w2*l2[0], rhohR*w2*l2[1], rhohR*w2*l2[2]};
    double EL = rhohL*w1*w1-P1;
    double ER = rhohR*w2*w2-P2;

    double FML[3] = {ML[0]*vnL+n[0]*al*P1, ML[1]*vnL+n[1]*al*P1, 
                        ML[2]*vnL+n[2]*al*P1};
    double FMR[3] = {MR[0]*vnR+n[0]*al*P2, MR[1]*vnR+n[1]*al*P2, 
                        MR[2]*vnR+n[2]*al*P2};
    double FEL = EL*vnL + P1*(vnL+bn);
    double FER = ER*vnR + P2*(vnR+bn);

    double UE = (sR*ER - sL*EL - FER + FEL) / (sR - sL);
    double FE = ((sR*FEL - sL*FER + sL*sR*(ER-EL)) / (sR - sL) - bn*UE) / al;

    double UM_hll[3], FM_hll[3];
    for(i=0; i<3; i++)
    {
        UM_hll[i] = (sR*MR[i]-sL*ML[i]-FMR[i]+FML[i]) / (sR-sL);
        FM_hll[i] = ((sR*FML[i]-sL*FMR[i]+sL*sR*(MR[i]-ML[i])) / (sR-sL)
                        - bn*UM_hll[i]) / al;
    }
    double UM = 0.0;
    double FM = 0.0;
    for(i=0; i<3; i++)
    {
        double igi = n[0]*igam[3*i] + n[1]*igam[3*i+1] + n[2]*igam[3*i+2];
        UM += igi * UM_hll[i];
        FM += igi * FM_hll[i];
    }

    double A = FE;
    double B = -FM-ign*UE;
    double C = ign*UM;

    double sS;
    if(fabs(4*A*C/(B*B)) < 1.0e-7)
        sS = -C/B * (1.0 + A*C/(B*B) + 2*A*A*C*C/(B*B*B*B));
    else
        sS = (-B - sqrt(B*B-4*A*C)) / (2*A);

    *Ss = hn * (al*sS - bn);
}

double mindt(double *prim, double wc, double *xp, double *xm)
{
    double x[3] = {0.5*(xm[0]+xp[0]), 0.5*(xm[1]+xp[1]), 0.5*(xm[2]+xp[2])};
    double r = x[0];
    double al, be[3], gam[9], igam[9], sqrtgam;
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    sqrtgam = metric_jacobian(x) / (r*al);

    double rho = prim[RHO];
    double Pp  = prim[PPP];
    double l[3] = {prim[URR], r*r*prim[UPP], prim[UZZ]};
    double B[3] = {prim[BRR]/sqrtgam,prim[BPP]/(r*sqrtgam),prim[BZZ]/sqrtgam};
    double cs  = sqrt(gamma_law*Pp/(rho+gamma_law/(gamma_law-1)*Pp));

    int i,j;
    double uS[3];
    double u2 = 0.0;

    for(i=0; i<3; i++)
    {
        uS[i] = 0.0;
        for(j=0; j<3; j++)
            uS[i] += igam[3*i+j]*l[j];
        u2 += uS[i]*l[i];
    }
    double w = sqrt(1.0+u2);

    double v2, vS[3];
    for(i=0; i<3; i++)
        vS[i] = uS[i]/w;
    v2 = u2/(w*w);

    double uB, B2, b2;
    uB = l[0]*B[0] + l[1]*B[1] + l[2]*B[2];
    B2 = 0.0;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            B2 += gam[3*i+j]*B[i]*B[j];
    b2 = (B2 + uB*uB) / (w*w);

    double va2 = b2 / (rho + gamma_law/(gamma_law-1.0)*Pp + b2);

    double cf = sqrt(cs*cs + va2 - cs*cs*va2);
    double sig = 1-cf*cf;

    double dvr = cf * sqrt(igam[0]*(1-cf*cf*v2) - sig*vS[0]*vS[0]) / w;
    double vrl = fabs(al * (vS[0]*sig - dvr) / (1-v2*cf*cf) - be[0]);
    double vrr = fabs(al * (vS[0]*sig + dvr) / (1-v2*cf*cf) - be[0]);

    double dvp = cf * sqrt(igam[4]*(1-cf*cf*v2) - sig*vS[1]*vS[1]) / w;
    double vpl = fabs(r * (al * (vS[1]*sig - dvp) / (1-v2*cf*cf) - be[1] -wc));
    double vpr = fabs(r * (al * (vS[1]*sig + dvp) / (1-v2*cf*cf) - be[1] -wc));
    
    double dvz = cf * sqrt(igam[8]*(1-cf*cf*v2) - sig*vS[2]*vS[2]) / w;
    double vzl = fabs(al * (vS[2]*sig - dvz) / (1-v2*cf*cf) - be[2]);
    double vzr = fabs(al * (vS[2]*sig + dvz) / (1-v2*cf*cf) - be[2]);

    double maxvr = vrr > vrl ? vrr : vrl;
    double maxvp = vpr > vpl ? vpr : vpl;
    double maxvz = vzr > vzl ? vzr : vzl;

    double dtr = get_dL(xp,xm,1)/maxvr;
    double dtp = get_dL(xp,xm,0)/maxvp;
    double dtz = get_dL(xp,xm,2)/maxvz;

    double dt = dtr;
    dt = dt < dtp ? dt : dtp;
    dt = dt < dtz ? dt : dtz;

    return dt;
}

double getReynolds(double *prim, double w, double *x, double dx)
{
    return 0.0;
}

void cons2prim_prep(double *cons, double *x)
{
    //TODO: complete this.
}

void cons2prim_solve_isothermal(double *cons, double *prim, double *x)
{
    //TODO: complete this.
    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q )
        prim[q] = cons[q]/cons[DDD];
}

void cons2prim_solve_adiabatic(double *cons, double *prim, double *x)
{
    double prec = 1.0e-12;
    double max_iter = 30;
    double Nextra = 10;

    double r = x[0];
    double z = x[2];

    double D = cons[DDD];
    double S[3] = {cons[SRR], cons[LLL], cons[SZZ]};
    double tau = cons[TAU];

    double al, be[3], gam[9], igam[9], jac, sqrtgam;
    double U[4];
    al = metric_lapse(x);
    metric_shift(x, be);
    metric_gam(x, gam);
    metric_igam(x, igam);
    jac = metric_jacobian(x) / r;
    sqrtgam = jac / al;
    frame_U(x, U);
    
    double B[3] = {r*cons[BRR]/sqrtgam, cons[BPP]/sqrtgam, cons[BZZ]/sqrtgam};

    double s2 = 0.0;
    double Us = 0.0;
    double BS = 0.0;
    double B2 = 0.0;

    int i,j;
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            s2 += igam[3*i+j]*S[i]*S[j];
    s2 /= D*D;

    for(i=0; i<3; i++)
        Us += S[i]*(be[i]*U[0] + U[i+1]);
    Us /= D;

    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            B2 += gam[3*i+j]*B[i]*B[j];
    double Q = sqrtgam * B2 / D;

    for(i=0; i<3; i++)
        BS += B[i]*S[i];
    double psi = sqrtgam * BS*BS / (D*D*D);

    
    double e = (tau/D + Us + 1.0) / (al*U[0]);
    double n = (gamma_law-1.0)/gamma_law;

    if(e*e < s2 && DEBUG && r < DEBUG_RMAX && fabs(z)<DEBUG_ZMAX)
    {
        printf("Not enough thermal energy (r=%.12lg, e2=%.12lg, s2=%.12lg)\n",
                r, e*e, s2);

        double cons0[NUM_Q];
        prim2cons(prim, cons0, x, 1.0);

        printf("prim: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                prim[RHO], prim[PPP], prim[URR], prim[UPP], prim[UZZ]);
        printf("cons0: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                cons0[DDD], cons0[TAU], cons0[SRR], cons0[LLL], cons0[SZZ]);
        printf("cons: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                cons[DDD], cons[TAU], cons[SRR], cons[LLL], cons[SZZ]);
    }

    //Newton-Raphson using the 2D scheme of Noble et al.
    //TODO: Rearrange to optimize for cold flows.

    //Initial guess: previous v2, eta
    double u2 = 0.0;
    double l[3] = {prim[URR], prim[UPP], prim[UZZ]};
    for(i=0; i<3; i++)
        for(j=0; j<3; j++)
            u2 += igam[3*i+j]*l[i]*l[j];
    double w = sqrt(1.0 + u2);
    double v20 = u2 / (1.0+u2); // sqrt(1+u2)-1
    double eta0 = w * (1.0 + gamma_law/(gamma_law-1.0)*prim[PPP]/prim[RHO]);

    //Run Newton-Raphson
    double v2, v21, eta, eta1;
    v21 = v20;
    eta1 = eta0;
    
    i = 0;
    int clean = -1;
    
    if(DEBUG2)
    {
        printf("s2 = %.12lg, e = %.12lg, Q = %.12lg, psi = %.12lg\n",
                    s2,e,Q,psi);
        printf("0: (%.12lg, %.12lg)\n", v21, eta1);
    }

    while(1)
    {
        v2 = v21;
        eta = eta1;
        
        w = 1.0/sqrt(1.0-v2);
        double fa = (eta+Q)*(eta+Q)*v2 - psi*(2*eta+Q)/(eta*eta) - s2;
        double fb = eta - n*(1-v2)*(eta-w) + 0.5*Q*(1+v2) - 0.5*psi/(eta*eta)
                    - e;
        double dfadv2 = (eta+Q)*(eta+Q);
        double dfadet = 2*(eta+Q)*(v2 + psi/(eta*eta*eta));
        double dfbdv2 = n * (eta-0.5*w) + 0.5*Q;
        double dfbdet = 1 - n*(1-v2) + psi/(eta*eta*eta);
        double detj = dfadv2*dfbdet - dfbdv2*dfadet;


        v21  = v2  - ( dfbdet*fa - dfadet*fb)/detj;
        eta1 = eta - (-dfbdv2*fa + dfadv2*fb)/detj;

        i++;

        if(v21 > 1.0)
            v21 = 0.5*(v2 + 1.0);
        if(v21 < 0.0)
            v21 = 0.5*v2;
        if(eta1 < 1.0)
            eta1 = 0.5*(eta+1.0);


        double err = (eta1-eta)/eta;
        //if(err != err)
        //    printf("WHAT: v2=%.12lg, eta=%.12lg\n");

        if(DEBUG2)
        {
            printf("%d: (%.12lg, %.12lg) (%.12lg, %.12lg) %.12lg\n", 
                    i, v21, eta1, fa, fb, err);
            printf("    %.12lg %.12lg %.12lg %.12lg\n", dfadv2, dfadet, 
                        dfbdv2, dfbdet);
        }

        if(fabs(err) < prec && clean < 0)
            clean = Nextra+1;
        if(clean >= 0)
            clean--;
        if(clean == 0 || i == max_iter)
            break;
    }

    if(i == max_iter && (DEBUG || DEBUG2) && r < DEBUG_RMAX
                && fabs(z)<DEBUG_ZMAX)
    {
        printf("ERROR: NR failed to converge. x=(%g,%g,%g)  err = %.12lg\n", 
                x[0], x[1], x[2], fabs(eta1-eta)/eta);
        printf("    s2 = %.12lg, e = %.12lg, Q = %.12lg, psi = %.12lg\n",
                s2, e, Q, psi);
        printf("    v20 = %.12lg, et0 = %.12lg, v21 = %.12lg, et1 = %.12lg\n",
                v20, eta0, v21, eta1);
    }

    v2 = v21;
    eta = eta1;

    //Prim recovery
    w = 1.0/sqrt(1-v2);
    double u0 = w/al;
    double h = eta/w;

    double rho = D / (jac*u0);
    if(rho < RHO_FLOOR)
        rho = RHO_FLOOR;
    double Pp = n * rho * (h-1.0);
    if(Pp < PRE_FLOOR*rho)
        Pp = PRE_FLOOR*rho;
    
    h = 1.0 + gamma_law/(gamma_law-1.0) * Pp/rho;
    double Bd[3];
    Bd[0] = gam[0]*B[0] + gam[1]*B[1] + gam[2]*B[2];
    Bd[1] = gam[3]*B[0] + gam[4]*B[1] + gam[5]*B[2];
    Bd[2] = gam[6]*B[0] + gam[7]*B[1] + gam[8]*B[2];

    l[0] = (S[0] + sqrtgam*BS*Bd[0]/(D*h*w)) / (D*h + sqrtgam*B2/w);
    l[1] = (S[1] + sqrtgam*BS*Bd[1]/(D*h*w)) / (D*h + sqrtgam*B2/w);
    l[2] = (S[2] + sqrtgam*BS*Bd[2]/(D*h*w)) / (D*h + sqrtgam*B2/w);

    prim[RHO] = rho;
    prim[URR] = l[0];
    prim[UPP] = l[1]/(r*r);
    prim[UZZ] = l[2];
    prim[PPP] = Pp;

    prim[BRR] = sqrtgam   * B[0];
    prim[BPP] = sqrtgam*r * B[1];
    prim[BZZ] = sqrtgam   * B[2];

    int q;
    for( q=NUM_C ; q<NUM_Q ; ++q )
        prim[q] = cons[q]/cons[DDD];
    
    if(e*e < s2 && DEBUG && r < DEBUG_RMAX && fabs(z) < DEBUG_ZMAX)
    {
        double cons1[NUM_Q];
        prim2cons(prim, cons1, x, 1.0);

        printf("prim1: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                prim[RHO], prim[PPP], prim[URR], prim[UPP], prim[UZZ]);
        printf("cons1: %.16lg %.16lg %.16lg %.16lg %.16lg\n",
                cons1[DDD], cons1[TAU], cons1[SRR], cons1[LLL], cons1[SZZ]);
    }
    
    if(DEBUG3)
    {
        FILE *f = fopen("c2p.out", "a");
        fprintf(f, "%.10lg %.10lg %.10lg %.10lg %.10lg %.10lg\n",
                    x[0], x[1], prim[URR], prim[UPP], cons[SRR], cons[LLL]);
        fprintf(f, "    %.10lg %.10lg %.10lg\n",
                        B2, prim[BRR], prim[BPP]);
        fclose(f);
    }
}

void cons2prim_finalize(double *prim, double *x)
{

}
