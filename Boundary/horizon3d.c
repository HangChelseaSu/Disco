
#include "../paul.h"
#include "../boundary.h"


void boundary_trans( struct domain * theDomain , int dim )
{
    if(dim == 1)
    {
        boundary_fixed_horizon(theDomain);
        boundary_fixed_rout(theDomain);
    }
    else if(dim == 2)
    {
        boundary_fixed_zbot(theDomain);
        boundary_fixed_ztop(theDomain);
    }
}

