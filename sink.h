#ifndef DISCO_SINK_H
#define DISCO_SINK_H

void setSinkParams(struct domain *theDomain);
void sink_src(const double *prim, double *cons, const double *xp,
                const double *xm, const double *xyz, double dV,
                double dt, double *pl_gas_track);
void cooling(const double *prim, double *cons,
             const double *xp, const double *xm, const double *xyz,
             double dV, double dt );
void damping(const double *prim, double *cons, const double *xp,
            const double *xm, const double *xyz, double dV, double dt );

#endif
