
#include "paul.h"
#include <string.h>
#include "geometry.h"
#include "hydro.h"
#include "boundary.h"

static double phi_max = 0.0;

void initial( double * , double * );
void subtract_omega( double * );

void setBCParams(struct domain *theDomain)
{
    phi_max = theDomain->phi_max;
}


void set_cell_init(double *prim, double piph, double dphi,
                    const double *r_jph, const double *z_kph, int j, int k)
{
    double xm[3] = {r_jph[j-1], piph - dphi, z_kph[k-1]};
    double xp[3] = {r_jph[j  ], piph          , z_kph[k  ]};
    double x[3];
    get_centroid_arr(xp, xm, x);

    initial(prim, x);
    subtract_omega(prim);
}

void set_cell_init_q(double *prim, double piph, double dphi,
                     const double *r_jph, const double *z_kph, 
                        int j, int k, int *qarr, int nq)
{
    double r = get_centroid(r_jph[j], r_jph[j-1], 1);
    double z = get_centroid(z_kph[k], z_kph[k-1], 2);
    double phi = piph - 0.5*dphi;
    double x[3] = {r, phi, z};
    double temp_prim[NUM_Q];
    initial(temp_prim, x);
    subtract_omega(temp_prim);

    int iq;
    for(iq=0; iq<nq; iq++)
        prim[qarr[iq]] = temp_prim[qarr[iq]];
}

void set_cells_copy_distant(double *prim, double *piph, double *dphi, int Np,
                        double *prim1, double *piph1, double *dphi1, int Np1)
{
    // Copies annulus c1 to c.

    int i, q;

    // Clear annulus.
    memset(prim, 0, Np*NUM_Q*sizeof(double));

    // Find first cell in c1 that intersects with c[0].
    int i1 = 0;
    double phi0 = piph[0] - dphi[0];
    for(i1=0; i1 < Np1; i1++)
    {
        double dphip = get_signed_dp(piph1[i1], phi0);
        double dphim = dphip - dphi1[i1];

        if(dphip > 0 && dphim <= 0.0)
            break;
    }

    if(i1 >= Np1)
        i1 = 0;

    //Loop through all cells c, add contributions from all neighbours in c1.
    for(i=0; i<Np; i++)
    {
        double phip = piph[i];
        double phim = piph[i] - dphi[i];
        double phip1, phim1;
        do
        {
            phip1 = piph1[i1];
            phim1 = piph1[i1] - dphi1[i1];

            while(phip1 < phim)
            {
                phip1 += phi_max;
                phim1 += phi_max;
            }
            while(phim1 > phip)
            {
                phip1 -= phi_max;
                phim1 -= phi_max;
            }
            double phi1 = phip < phip1 ? phip : phip1;
            double phi2 = phim > phim1 ? phim : phim1;
            double dphi = phi1-phi2;

            if(dphi< 0.0)
                printf("WHOA %d %d %.12lg\n", i, i1, dphi);
            for(q=0; q<NUM_Q; q++)
                prim[NUM_Q*i+q] += prim1[NUM_Q*i1+q] * dphi;

            i1 = i1 < Np1-1 ? i1+1 : 0;
            if(i1 == Np1)
                i1 = 0;
        }
        while(phip1 < phip);
        i1 = i1 > 0 ? i1-1 : Np1-1;

    }

    // Divide by total area.
    for(i=0; i<Np; i++)
        for(q=0; q<NUM_Q; q++)
            prim[NUM_Q*i+q] /= dphi[i];
}

void boundary_fixed_rinn( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<NgRa; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_rout( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k);
            }  
    }
}

void boundary_fixed_zbot( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<NgZa; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k);
            }  
    }
}
void boundary_fixed_ztop( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k);
            }  
    }
}

void boundary_zerograd_rinn( struct domain *theDomain, int diode)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;

    int *dim_rank = theDomain->dim_rank;

    int j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=NgRa-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                int jk1 = j+1+Nr*k;
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);

                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        int iurr = NUM_Q*i + URR;
                        if(prim[jk][iurr] > 0)
                            prim[jk][iurr] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_rout( struct domain *theDomain, int diode)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j + Nr*k;
                int jk1 = j-1+Nr*k;
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);

                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        int iurr = NUM_Q*i + URR;
                        if(prim[jk][iurr] < 0)
                            prim[jk][iurr] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_zbot( struct domain *theDomain, int diode)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;

    int *dim_rank = theDomain->dim_rank;

    int j,k;

    if(dim_rank[1] == 0)
    {
        for(j=0; j<Nr; j++)
            for(k=NgZa-1; k>=0; k--)
            {
                int jk = j+Nr*k;
                int jk1 = j+Nr*(k+1);
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);
                
                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        int iurr = NUM_Q*i + UZZ;
                        if(prim[jk][iuzz] > 0)
                            prim[jk][iuzz] = 0.0;
                    }
                }
            }
    }
}

void boundary_zerograd_ztop( struct domain *theDomain, int diode)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(j=0; j<Nr; j++)
            for(k=Nz-NgZb; k<Nz; k++)
            {
                int jk = j+Nr*k;
                int jk1 = j+Nr*(k-1);
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);
                
                if(diode)
                {
                    int i;
                    for(i=0; i<Np[jk]; i++)
                    {
                        int iurr = NUM_Q*i + UZZ;
                        if(prim[jk][iuzz] < 0)
                            prim[jk][iuzz] = 0.0;
                    }
                }
            }
    }
}

void boundary_reflect_rinn( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);
            for(j=NgRa-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*NgRa-j-1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);
                    
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 0);
                }
            }
        }
    }
}

void boundary_reflect_rout( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*(Nr-NgRb) - j - 1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);
                
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 0);
                }
            }
        }
    }
}

void boundary_reflect_zbot( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=NgZa-1; k>=0; k--)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*NgZa-k-1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 2);
                }
            }
        }
    }
}

void boundary_reflect_ztop( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*(Nz-NgZb) - k - 1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 2);
                }
            }
        }
    }
}

void boundary_fixed_q_rinn( struct domain *theDomain, int *q, int nq)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<NgRa; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k, q, nq);
            }  
    }
}

void boundary_fixed_q_rout( struct domain *theDomain, int *q, int nq)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k, q, nq);
            }  
    }
}

void boundary_fixed_q_zbot( struct domain *theDomain, int *q, int nq)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<NgZa; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k, q, nq);
            }  
    }
}
void boundary_fixed_q_ztop( struct domain *theDomain, int *q, int nq)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                    set_cell_init_q(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                    dphi[jk][i], r_jph, z_kph, j, k, q, nq);
            }  
    }
}

void boundary_fixed_phi_rinn( struct domain *theDomain, double phia,
                              double phib)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
            for(j=0; j<NgRa; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*piph[jk][i];
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                        dphi[jk][i], r_jph, z_kph, j, k);
                }
            }  
    }
}

void boundary_fixed_phi_rout( struct domain *theDomain, double phia,
                            double phib)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*piph[jk][i];
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                        dphi[jk][i], r_jph, z_kph, j, k);
                }
            }  
    }
}

void boundary_fixed_phi_zbot( struct domain *theDomain, double phia,
                              double phib)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=0; k<NgZa; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*piph[jk][i];
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                        dphi[jk][i], r_jph, z_kph, j, k);
                }
            }  
    }
}

void boundary_fixed_phi_ztop( struct domain *theDomain, double phia,
                              double phib)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*piph[jk][i];
                    while(phi < phia)
                        phi += phi_max;
                    if(phi < phib)
                        set_cell_init(&(prim[jk][NUM_Q*i]), piph[jk][i],
                                        dphi[jk][i], r_jph, z_kph, j, k);
                }
            }  
    }
}

void boundary_noslip_rinn( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRa = theDomain->NgRa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[0] == 0 )
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);
            for(j=NgRa-1; j>=0; j--)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*NgRa-j-1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);
                    
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 0);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 1);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 2);
                }
            }
        }
    }
}

void boundary_noslip_rout( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgRb = theDomain->NgRb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[0] == dim_size[0]-1)
    {
        for(k=0; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=Nr-NgRb; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int j1 = 2*(Nr-NgRb) - j - 1;
                int jk1 = j1 + Nr*k;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);
                
                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 0);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 1);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 2);
                }
            }
        }
    }
}

void boundary_noslip_zbot( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int *Np = theDomain->Np;
    int NgZa = theDomain->NgZa;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;

    int i,j,k;

    if(dim_rank[1] == 0)
    {
        for(k=NgZa-1; k>=0; k--)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*NgZa-k-1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 0);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 1);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 2);
                }
            }
        }
    }
}

void boundary_noslip_ztop( struct domain *theDomain)
{
    double **prim = theDomain->prim;
    double **piph = theDomain->piph;
    double **dphi = theDomain->dphi;

    int Nr = theDomain->Nr;
    int Nz = theDomain->Nz;
    int *Np = theDomain->Np;
    int NgZb = theDomain->NgZb;
    double *r_jph = theDomain->r_jph;
    double *z_kph = theDomain->z_kph;

    int *dim_rank = theDomain->dim_rank;
    int *dim_size = theDomain->dim_size;

    int i,j,k;

    if(dim_rank[1] == dim_size[1]-1)
    {
        for(k=Nz-NgZb; k<Nz; k++)
        {
            double z = get_centroid(z_kph[k], z_kph[k-1], 2);

            for(j=0; j<Nr; j++)
            {
                int jk = j+Nr*k;
                
                int k1 = 2*(Nz-NgZb) - k - 1;
                int jk1 = j + Nr*k1;
                
                set_cells_copy_distant(prim[jk], piph[jk], dphi[jk], Np[jk],
                                    prim[jk1], piph[jk1], dphi[jk1], Np[jk1]);

                double r = get_centroid(r_jph[j], r_jph[j-1], 1);

                for(i=0; i<Np[jk]; i++)
                {
                    double phi = piph[jk][i] - 0.5*dphi[jk][i];
                    double x[3] = {r, phi, z};
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 0);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 1);
                    reflect_prims(&(prim[jk][NUM_Q*i]), x, 2);
                }
            }
        }
    }
}

