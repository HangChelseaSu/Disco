#ifndef DISCO_PLANET_H
#define DISCO_PLANET_H

#include "paul.h"

void setGravParams( struct domain * theDomain );

double phigrav( double M , double r , double eps , int type);
double fgrav( double M , double r , double eps , int type);
void adjust_gas( struct planet * pl , double * x , double * prim , double gam);
void planetaryForce( struct planet * pl , double r , double phi , double z ,
                    double * fr , double * fp , double * fz , int mode );
void planet_src( struct planet * pl , double * prim , double * cons ,
                double * xp , double * xm , double dVdt );

void planet_init_kin(struct planet *pl);
void planet_zero_aux(struct planet *pl);

void setupPlanets(struct domain *theDomain);
void initializePlanetTracking(struct domain *theDomain);
void updatePlanetsKin(struct domain *theDomain, double dt);
void updatePlanetsAux(struct domain *theDomain);
void movePlanetsLive(struct domain *theDomain);
void copyPlanetsRK( struct domain *theDomain);
void adjustPlanetsRKkin( struct domain *theDomain, double RK);
void adjustPlanetsRKaux( struct domain *theDomain, double RK);

// Functions in Planet/ setup
void setPlanetParams( struct domain * theDomain );
int planet_motion_analytic( void );
void initializePlanets( struct planet * thePlanets );
void movePlanets( struct planet * thePlanets , double t , double dt );
void forcePlanets( struct planet * thePlanets , double dt );

#endif
