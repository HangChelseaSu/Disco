#include "../paul.h"
#include "../geometry.h"
#include "../planet.h"

#define N_AUX_PER_PLANET 12

static double gamma_law = 0.0;
static int Npl = 0;


/*
 * 
 */

void setReportParams(struct domain *theDomain)
{
    gamma_law = theDomain->theParList.Adiabatic_Index;
    Npl = theDomain->Npl;

    if(Npl < 2)
    {
        printf("This report needs AT LEAST 2 planets to work!\n");
        exit(1);
    }
}

int num_shared_reports()
{
    return 6;
}

int num_distributed_aux_reports()
{
    return 2 * N_AUX_PER_PLANET;
}

int num_distributed_integral_reports()
{
    return 29;
}

void get_shared_reports(double *Q, struct domain *theDomain)
{
    Q[0] = theDomain->thePlanets[0].M;
    Q[1] = theDomain->thePlanets[1].M;
    Q[2] = theDomain->thePlanets[0].r;
    Q[3] = theDomain->thePlanets[1].r;
    Q[4] = theDomain->thePlanets[0].phi;
    Q[5] = theDomain->thePlanets[1].phi;
}

void get_distributed_aux_reports(double *Q, struct domain *theDomain)
{
    double *pl_aux1 = theDomain->pl_aux;
    double *pl_aux2 = theDomain->pl_aux + NUM_PL_AUX;

    int idx_aux[N_AUX_PER_PLANET] = {
                PL_SNK_M, PL_GRV_JZ, PL_SNK_JZ, PL_GRV_PX, PL_GRV_PY,
                PL_SNK_PX, PL_SNK_PY, PL_GRV_K, PL_SNK_K,
                PL_SNK_MX, PL_SNK_MY, PL_SNK_SZ};

    int q;
    for(q=0; q<N_AUX_PER_PLANET; q++)
    {
        Q[2*q + 0] = pl_aux1[idx_aux[q]];
        Q[2*q + 1] = pl_aux2[idx_aux[q]];
    }
}

void get_distributed_integral_reports(double *x, double *prim, double *Q,
                                      struct domain *theDomain)
{
    double rpz[3];
    get_rpz(x, rpz);

    double r = rpz[0];
    double phi = rpz[1];
    double z = rpz[1];

    double rho = prim[RHO];
    double vr = prim[URR];
    double vp = r*prim[UPP];

    double cosp = cos(phi);
    double sinp = sin(phi);
    double cos2p = (cosp-sinp)*(cosp+sinp);
    double sin2p = 2*cosp*sinp;

    double xyz[3] = {r*cosp, r*sinp, z};

    double Fxyz[3];
    planetaryForce(theDomain->thePlanets+0, xyz, Fxyz);
    double Fp1 = cosp * Fxyz[1] - sinp * Fxyz[0];

    planetaryForce(theDomain->thePlanets+1, xyz, Fxyz);
    double Fp2 = cosp * Fxyz[1] - sinp * Fxyz[0];

    if(r > 1.0)
    {
        Q[0] = prim[RHO] * r * Fp1;
        Q[1] = prim[RHO] * r * Fp2;
        Q[2] = 0.0;
        Q[3] = 0.0;
    }
    else
    {
        Q[0] = 0.0;
        Q[1] = 0.0;
        Q[2] = prim[RHO] * r * Fp1;
        Q[3] = prim[RHO] * r * Fp2;
    }

    Q[4] = rho;
    Q[5] = rho * r*cosp;
    Q[6] = rho * r*sinp;
    Q[7] = rho * r*r*cos2p;
    Q[8] = rho * r*r*sin2p;


    int a_start = 9;

    int nvals = 4;
    int nsplit = 5;

    double ra = 1.0;
    double rb = 6.0;
    
    double dr = (rb-ra) / nsplit;

    int s, q;
    for(s=0; s<nsplit; s++)
    {
        if((r > ra + s*dr) && (r <= ra + (s+1)*dr))
        {
            double v2 = vr*vr + vp*vp;
            double rdv = r*vr;
            double vx = vr*cosp - vp*sinp;
            double vy = vr*sinp + vp*cosp;

            double ex = (r*v2-1)*cosp - rdv*vx;
            double ey = (r*v2-1)*sinp - rdv*vy;

            Q[a_start + s*nvals + 0] = 1.0;
            Q[a_start + s*nvals + 1] = rho;
            Q[a_start + s*nvals + 2] = rho * ex;
            Q[a_start + s*nvals + 3] = rho * ey;
        }
        else
        {
            for(q=0; q<nvals; q++)
                Q[a_start + s*nvals + q] = 0.0;
        }
    }
}

