#ifndef DISCO_METRIC_H
#define DISCO_METRIC_H

void setMetricParams(struct domain *theDomain);

double metric_lapse(const double x[3]);
void metric_shift(const double x[3], double b[3]);
void metric_gam(const double x[3], double gam[9]);
void metric_igam(const double x[3], double igam[9]);
double metric_jacobian(const double x[3]);
void metric_der_g(const double x[3], int i, double dg[16]);
void metric_der_lapse(const double x[3], double da[4]);
void metric_der_shift(const double x[3], double db[12]);
int metric_killing(int mu);

#endif
