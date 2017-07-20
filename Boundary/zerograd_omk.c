
#include "../paul.h"
#include <string.h>

void boundary_fixed_rinn( struct domain *theDomain);
void boundary_fixed_rout( struct domain *theDomain);
void boundary_fixed_zbot( struct domain *theDomain);
void boundary_fixed_ztop( struct domain *theDomain);
void boundary_zerograd_rinn( struct domain *theDomain, int diode);
void boundary_zerograd_rout( struct domain *theDomain, int diode);
void boundary_zerograd_zbot( struct domain *theDomain, int diode);
void boundary_zerograd_ztop( struct domain *theDomain, int diode);
void boundary_reflect_rinn( struct domain *theDomain);
void boundary_reflect_rout( struct domain *theDomain);
void boundary_reflect_zbot( struct domain *theDomain);
void boundary_reflect_ztop( struct domain *theDomain);
void boundary_fixed_horizon( struct domain *theDomain);
void boundary_fixed_q_rinn( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_rout( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_zbot( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_ztop( struct domain *theDomain, int *q, int nq);

void boundary_trans( struct domain * theDomain , int dim )
{
    int nq = 1;
    int q[1];
    q[0] = UPP;

    if(dim == 1)
    {
        boundary_zerograd_rinn(theDomain, 1);
        boundary_zerograd_rout(theDomain, 1);
        boundary_fixed_q_rinn(theDomain, q, nq);
        boundary_fixed_q_rout(theDomain, q, nq);
    }
    else if(dim == 2)
    {
        boundary_zerograd_zbot(theDomain, 1);
        boundary_zerograd_ztop(theDomain, 1);
    }
}

