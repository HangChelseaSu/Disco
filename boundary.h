#ifndef DISCO_BOUNDARY_H
#define DISCO_BOUNDARY_H

#include "paul.h"

// For initialization
void setBCParams(struct domain *theDomain);

//The boundary function, defined in Boundary/
void boundary_trans( struct domain * theDomain , int dim );

//Boundary options
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

void boundary_fixed_q_rinn( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_rout( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_zbot( struct domain *theDomain, int *q, int nq);
void boundary_fixed_q_ztop( struct domain *theDomain, int *q, int nq);

void boundary_fixed_horizon( struct domain *theDomain);

void boundary_fixed_phi_rinn( struct domain *theDomain, double phia,
                                double phib);
void boundary_fixed_phi_rout( struct domain *theDomain, double phia,
                                double phib);
void boundary_fixed_phi_zbot( struct domain *theDomain, double phia,
                                double phib);
void boundary_fixed_phi_ztop( struct domain *theDomain, double phia,
                                double phib);

void boundary_noslip_rinn( struct domain *theDomain);
void boundary_noslip_rout( struct domain *theDomain);
void boundary_noslip_zbot( struct domain *theDomain);
void boundary_noslip_ztop( struct domain *theDomain);

//Low-level routines for boundary
void set_cell_init(double *prim, double piph, double dphi,
                    const double *r_jph, const double *z_kph, int j, int k);
void set_cell_init_q(double *prim, double piph, double dphi,
                     const double *r_jph, const double *z_kph, 
                        int j, int k, int *qarr, int nq);
void set_cells_copy_distant(double *prim, double *piph, double *dphi, int Np,
                        double *prim1, double *piph1, double *dphi1, int Np1);

#endif
