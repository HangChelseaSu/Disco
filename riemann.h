#ifndef DISCO_RIEMANN_H
#define DISCO_RIEMANN_H

#include "paul.h"

void setRiemannParams( struct domain * theDomain );

void riemann_phi( struct cell * cL , struct cell * cR, double * x ,
                const double *xp, const double *xm, double dAdt );
void riemann_trans( struct face * F , double dt , int dim , double rp,
                    double rm, double zp, double zm,
                    double *fdAdt_hydro, double *fdAdt_visc);

void solve_riemann(const double *primL, const double *primR,
                   const double *x, const double *n, 
                   const double *xp, const double *xm,
                   double w, double dAdt, 
                   int dim,
                   double *E1_riemann, double *B1_riemann,
                   double *E2_riemann, double *B2_riemann,
                   double *fdAdt_hydro);
void solve_visc(const double *primL, const double *primR,
                   const double *gradLr, const double *gradLp,
                   const double *gradLz,
                   const double *gradRr, const double *gradRp, 
                   const double *gradRz,
                   const double *x, const double *n, 
                   double dAdt,  double *fdAdt_visc);

#endif
