#ifndef DISCO_RIEMANN_H
#define DISCO_RIEMANN_H

void riemann_phi(   const double *primCL, const double *primCR,
                    double *consL, double *consR,
            const double *gradrL, const double *gradpL, const double *gradzL,
            const double *gradrR, const double *gradpR, const double *gradzR,
            double piphL, double dphiL, double piphR, double dphiR, double wf,
            const double *x, const double *xp, const double *xm, double dAdt,
            double *EL, double *BL, double *ER, double *BR);

void solve_riemann(const double *primL, const double *primR,
                   double *consL, double *consR,
            const double *gradLr, const double *gradLp, const double *gradLz,
            const double *gradRr, const double *gradRp, const double *gradRz,
                   const double *x, const double *n, 
                   const double *xp, const double *xm,
                   double w, double dAdt, 
                   int dim,
                   double *E1_riemann, double *B1_riemann,
                   double *E2_riemann, double *B2_riemann,
                   double *fdAdt_hydro, double *fdAdt_visc);
#endif
