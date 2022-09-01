#ifndef DISCO_FRAME_H
#define DISCO_FRAME_H

void setFrameParams(struct domain *theDomain);

void frame_U(const double x[3], double U[4]);
void frame_der_U(const double x[3], double dU[16]);

#endif
