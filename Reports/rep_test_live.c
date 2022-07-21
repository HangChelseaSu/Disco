#include "../paul.h"
#include "../geometry.h"
#include "../planet.h"

static double gamma_law = 0.0;
static int Npl = 0;

void setReportParams(struct domain *theDomain)
{
    gamma_law = theDomain->theParList.Adiabatic_Index;
    Npl = theDomain->Npl;
}

int num_shared_reports()
{
    return Npl * (NUM_PL_KIN + NUM_PL_AUX);
}

int num_distributed_aux_reports()
{
    return 0;
}

int num_distributed_integral_reports()
{
    return Npl;
}

void get_shared_reports(double *Q, struct domain *theDomain)
{
    int p, q;

    int stride = NUM_PL_KIN + NUM_PL_AUX;

    for(p=0; p < Npl; p++)
    {
        for(q=0; q < NUM_PL_KIN; q++)
            Q[p*stride + q] = theDomain->pl_kin[p*NUM_PL_KIN + q];
        for(q=0; q < NUM_PL_AUX; q++)
            Q[p*stride + NUM_PL_KIN + q] = theDomain->pl_aux[p*NUM_PL_AUX + q];
    }
}

void get_distributed_aux_reports(double *Q, struct domain *theDomain)
{
    //Silence is golden
}

void get_distributed_integral_reports(double *x, double *prim, double *Q,
                                      struct domain *theDomain)
{
    int p;

    for(p=0; p<Npl; p++)
    {
        double rpz[3];
        get_rpz(x, rpz);
        Q[p] = prim[RHO] * planetaryPotential(theDomain->thePlanets+p,
                                                rpz[0], rpz[1], rpz[2]);
    }
}

