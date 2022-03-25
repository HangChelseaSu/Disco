#ifndef DISCO_ANALYSIS_H
#define DISCO_ANALYSIS_H

#include "paul.h"

int num_diagnostics(void);
void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain );

void zero_diagnostics( struct domain * theDomain );
void avg_diagnostics( struct domain * theDomain );
void adjust_RK_diag(struct domain *theDomain, double RK);
void copy_RK_diag(struct domain *theDomain);
void add_diagnostics( struct domain * theDomain , double dt );

#endif
