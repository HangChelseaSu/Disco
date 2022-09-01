#include "paul.h"
#include "omega.h"
#include "geometry.h"
#include "planet.h"

static int sinkType = 0;
static int sinkNumber = 0;
static double sinkPar1 = 0.0;
static double sinkPar2 = 0.0;
static double sinkPar3 = 0.0;
static double sinkPar4 = 0.0;
static double sinkPar5 = 0.0;
static int nozzleType;
static double nozzlePar1 = 0.0;
static double nozzlePar2 = 0.0;
static double nozzlePar3 = 0.0;
static double nozzlePar4 = 0.0;

static int coolType = 0;
static double coolPar1 = 0.0;
static double coolPar2 = 0.0;
static double coolPar3 = 0.0;
static double coolPar4 = 0.0;

static double gamma_law = 0.0;
static int twoD = 0;
static int Npl = 0;
static struct planet *thePlanets = NULL;

static double visc = 0.0;
static double Mach = 1.0;

static double rmax;
static double rmin;
static double zmax;
static double zmin;

static int DAMP_OUTER = 0;
static int DAMP_INNER = 0;
static int DAMP_LOWER = 0;
static int DAMP_UPPER = 0;
static double dampTimeInner = 0.0;
static double dampTimeOuter = 0.0;
static double dampTimeLower = 0.0;
static double dampTimeUpper = 0.0;
static double dampLenInner = 0.0;
static double dampLenOuter = 0.0;
static double dampLenUpper = 0.0;
static double dampLenLower = 0.0;

void initial( double * , double * );
void prim2cons( double * , double * , double * , double );

void setSinkParams(struct domain *theDomain)
{
    sinkType = theDomain->theParList.sinkType;
    sinkNumber = theDomain->theParList.sinkNumber;
    sinkPar1 = theDomain->theParList.sinkPar1;
    sinkPar2 = theDomain->theParList.sinkPar2;
    sinkPar3 = theDomain->theParList.sinkPar3;
    sinkPar4 = theDomain->theParList.sinkPar4;
    sinkPar5 = theDomain->theParList.sinkPar5;
    nozzleType = theDomain->theParList.nozzleType;
    nozzlePar1 = theDomain->theParList.nozzlePar1;
    nozzlePar2 = theDomain->theParList.nozzlePar2;
    nozzlePar3 = theDomain->theParList.nozzlePar3;
    nozzlePar4 = theDomain->theParList.nozzlePar4;
    coolType = theDomain->theParList.coolType;
    coolPar1 = theDomain->theParList.coolPar1;
    coolPar2 = theDomain->theParList.coolPar2;
    coolPar3 = theDomain->theParList.coolPar3;
    coolPar4 = theDomain->theParList.coolPar4;

    rmax = theDomain->theParList.rmax;
    rmin = theDomain->theParList.rmin;
    zmax = theDomain->theParList.zmax;
    zmin = theDomain->theParList.zmin;

    gamma_law = theDomain->theParList.Adiabatic_Index;
    if(theDomain->Nz == 1)
        twoD = 1;

    visc = theDomain->theParList.viscosity;
    Mach = theDomain->theParList.Disk_Mach;

    thePlanets = theDomain->thePlanets;
    Npl = theDomain->Npl;

    DAMP_INNER = theDomain->theParList.dampInnerType;
    DAMP_OUTER = theDomain->theParList.dampOuterType;
    DAMP_UPPER = theDomain->theParList.dampUpperType;
    DAMP_LOWER = theDomain->theParList.dampLowerType;
    dampTimeInner = theDomain->theParList.dampTimeInner;
    dampLenInner = theDomain->theParList.dampLenInner;
    dampTimeOuter = theDomain->theParList.dampTimeOuter;
    dampLenOuter = theDomain->theParList.dampLenOuter;
    dampTimeUpper = theDomain->theParList.dampTimeUpper;
    dampLenUpper = theDomain->theParList.dampLenUpper;
    dampTimeLower = theDomain->theParList.dampTimeLower;
    dampLenLower = theDomain->theParList.dampLenLower;
}

void sink_src(double *prim, double *cons, double *xp, double *xm, double dV,
              double dt, double *pl_gas_track)
{
    if(nozzleType == 1)
    {
        double R0 = nozzlePar1;
        double dR = nozzlePar2;
        double v = nozzlePar3;
        double mach = nozzlePar4;

        double r = 0.5*(xp[0]+xm[0]);
        double phi = 0.5*(xp[1]+xm[1]);
        double z = 0.5*(xp[2]+xm[2]);
        if(phi > M_PI) phi -= 2*M_PI;

        if((r-R0)*(r-R0) + r*r*phi*phi + z*z > dR*dR)
            return;

        double Mdot = 1.0;
        double rhodot = Mdot * 3.0/(4.0*M_PI*dR*dR*dR);
        if(twoD)
            rhodot *= 2*sqrt(dR*dR - (r-R0)*(r-R0) - r*r*phi*phi);
        if(rhodot * dV*dt > cons[RHO])
            rhodot = cons[RHO] / dV*dt;

        double x[3] = {r, phi, z};
        double om = get_om(x);
        double cs2 = v*v / (mach*mach);
        double eps = cs2 / (gamma_law*(gamma_law-1));

        cons[RHO] += rhodot * dV*dt;
        cons[LLL] += rhodot * r*v * dV*dt;
        cons[TAU] += (0.5* rhodot * (v-r*om)*(v-r*om) + rhodot*eps) * dV*dt;
        if(NUM_Q > NUM_C)
            cons[NUM_C] += rhodot*dV*dt;
    }



    if(sinkType != 0){

      double x[3];
      get_centroid_arr(xp, xm, x);
      double r = x[0];
      double phi = x[1];
      double z = x[2];

      double rho = prim[RHO];
      double vr  = prim[URR];
      double vp  = prim[UPP]*r;
      double vz  = prim[UZZ];
      double press  = prim[PPP];
      double specenth = press*(1.0 + 1.0/(gamma_law - 1.0))/rho;

      double cosg = cos(phi);
      double sing = sin(phi);
      double gx = r*cosg;
      double gy = r*sing;

      double px, py, dx, dy, mag, eps, epsfactor;
      double rate, surfdiff;
      int pi;
      int numSinks = Npl;

      if ((sinkNumber>0) && (sinkNumber<numSinks)) numSinks = sinkNumber;
      for (pi=0; pi<numSinks; pi++){


          double cosp = cos(thePlanets[pi].phi);
          double sinp = sin(thePlanets[pi].phi);
          px = thePlanets[pi].r*cosp;
          py = thePlanets[pi].r*sinp;

          dx = gx-px;
          dy = gy-py;
          mag = dx*dx + dy*dy;
          mag = sqrt(mag);

          //the part that depends on sinkType
          double f_acc = 0.0;
          if(sinkType == 3){	//constant
            if (mag <= sinkPar3) f_acc = 1.0;
          }
          else if(sinkType == 2){	//exponential
            eps = sinkPar3;
            eps = pow(eps, sinkPar4);
            epsfactor = sinkPar5;
            if(epsfactor <= 0.0) epsfactor = 1.0;
            double magPow = pow(mag, sinkPar4);
            f_acc = exp(-magPow/(eps*epsfactor));
          }
          else if(sinkType == 1){	//polynomial, compact support
            eps = sinkPar3;
            double pwrM = sinkPar4;
            double pwrN = sinkPar5;
            f_acc = 1.0 - pow((mag/eps),pwrM);
            f_acc = pow(f_acc, pwrN);
            if (mag >= eps){
              f_acc = 0.0;
            }
          }

          rate = sinkPar1*thePlanets[pi].omega;
          surfdiff = rho*rate*f_acc;

          if(rate == 0.0 || f_acc == 0.0)
              continue;

          //delta clamped to [0, 1]
          double delta = fmax(fmin(sinkPar2, 1.0), 0.0);
          double rp, omp, vxp, vyp, vxg, vyg, vxr, vyr, vp_p, vp_r, vxn, vyn, cphi, sphi, vg_r, vg_p;
          rp = thePlanets[pi].r;
          omp = thePlanets[pi].omega;
          vp_p = rp*omp;
          vp_r = thePlanets[pi].vr;
          vxp = vp_r*cosp - vp_p*sinp;
          vyp = vp_r*sinp + vp_p*cosp;

          double vxg1 = vr*cosg - vp*sing;
          double vyg1 = vr*sing + vp*cosg;

          // cartesian components of relative velocity
          vxr = vxg1 - vxp;
          vyr = vyg1 - vyp;

          // cos & sin of angular position of gas relative to planet
          cphi = dx/mag;
          sphi = dy/mag;

          //Amount of mass to accrete
          double dM = dV*dt*surfdiff;

          // Velocity of gas to accrete in planet frame
          vxn = (cphi*cphi + (1.0-delta)*sphi*sphi)*vxr + delta*sphi*cphi*vyr;
          vyn = delta*cphi*sphi*vxr + (sphi*sphi + (1.0-delta)*cphi*cphi)*vyr;

          vxg = vxn + vxp;
          vyg = vyn + vyp;
          vg_r =  vxg*cosg + vyg*sing;
          vg_p = -vxg*sing + vyg*cosg;

          cons[DDD] -= dM;
          cons[SRR] -= vg_r*dM;
          cons[LLL] -= r*vg_p*dM;
          cons[SZZ] -= vz*dM;
          double v2 = vxg*vxg + vyg*vyg;
          cons[TAU] -= dM*(specenth + 0.5*v2
                    - 0.5*((vxg-vxg1)*(vxg-vxg1) + (vyg-vyg1)*(vyg-vyg1)));

          double *my_gas_track = pl_gas_track + pi*NUM_PL_INTEGRALS;
          my_gas_track[PL_SNK_M] += dM;
          my_gas_track[PL_SNK_PX] += dM * vxg;
          my_gas_track[PL_SNK_PY] += dM * vyg;
          my_gas_track[PL_SNK_PZ] += dM * vz;
          my_gas_track[PL_SNK_JZ] += dM * r*vg_p;
          my_gas_track[PL_SNK_SZ] += dM * (dx * vyn - dy * vxn);
          my_gas_track[PL_SNK_MX] += dM * dx;
          my_gas_track[PL_SNK_MY] += dM * dy;
          my_gas_track[PL_SNK_MZ] += dM * z;
          my_gas_track[PL_SNK_EGAS] += dM * (specenth + 0.5*v2
                  - 0.5*((vxg-vxg1)*(vxg-vxg1) + (vyg-vyg1)*(vyg-vyg1)));

          int pi2;
          for(pi2=0; pi2<Npl; pi2++)
          {
              double Phi = planetaryPotential(thePlanets+pi2, r, phi, z);
              pl_gas_track[NUM_PL_INTEGRALS*pi2 + PL_SNK_UGAS] += dM * Phi;
          }
      }
    }
}




void cooling(double *prim, double *cons, double *xp, double *xm, double dV, double dt )
{
  if(coolType == COOL_BETA || coolType == COOL_BETA_RELAX)
  {
    //Beta-cooling
    double press = prim[PPP];
    double rho = prim[RHO];
    double gm1 = gamma_law-1.0;
    double beta = coolPar1;
    double x[3];
    get_centroid_arr(xp, xm, x);
    double enTarget = get_cs2(x)/(gamma_law*gm1);
    double enCurrent = press/(rho*gm1);
    double omtot = 0.0;

    if(coolType == COOL_BETA_RELAX || enCurrent > enTarget)
    {
      double r = x[0];
      double phi = x[1];
      double z = x[2];

      double cosg = cos(phi);
      double sing = sin(phi);
      double gx = r*cosg;
      double gy = r*sing;

      int pi;
      double fr,fp,fz, cosp, sinp, px, py, dx, dy, mag;
      for (pi=0; pi<Npl; pi++)
      {
        cosp = cos(thePlanets[pi].phi);
        sinp = sin(thePlanets[pi].phi);
        px = thePlanets[pi].r*cosp;
        py = thePlanets[pi].r*sinp;

        dx = gx-px;
        dy = gy-py;
        mag = dx*dx + dy*dy;
        mag = sqrt(mag);
        planetaryForce( thePlanets + pi, r, phi, z, &fr, &fp, &fz, 1);
        omtot += sqrt(fr*fr + fp*fp + fz*fz)/mag;
      }
      omtot = sqrt(omtot);

      ////direct implementation of source term
      //cons[TAU] += rho*(enCurrent - enTarget)*dt*dV*beta*omtot;

      //integrate source term over timestep
      double Tm1 = expm1(-dt*omtot/beta);
      cons[TAU] += rho*dV*( enCurrent - enTarget)*Tm1;	//N.B. Tm1 is in [-1 and 0]
    }
  }
}


void damping(double *prim, double *cons, double *xp, double *xm, double dV, double dt )
{
  if (DAMP_INNER + DAMP_OUTER + DAMP_LOWER + DAMP_UPPER > 0){
    double x[3];
    get_centroid_arr(xp, xm, x);

    double omtot = 1.0;
    if (DAMP_INNER==2 || DAMP_OUTER==2 || DAMP_LOWER==2 || DAMP_UPPER==2){
      double r = x[0];
      double phi = x[1];
      double z = x[2];

      double cosg = cos(phi);
      double sing = sin(phi);
      double gx = r*cosg;
      double gy = r*sing;

      int pi;
      omtot = 0.0;
      double fr,fp,fz, cosp, sinp, px, py, dx, dy, mag;
      for (pi=0; pi<Npl; pi++)
      {
        cosp = cos(thePlanets[pi].phi);
        sinp = sin(thePlanets[pi].phi);
        px = thePlanets[pi].r*cosp;
        py = thePlanets[pi].r*sinp;

        dx = gx-px;
        dy = gy-py;
        mag = dx*dx + dy*dy;
        mag = sqrt(mag);
        planetaryForce( thePlanets + pi, r, phi, z, &fr, &fp, &fz, 1);
        omtot += sqrt(fr*fr + fp*fp + fz*fz)/mag;
      }
      omtot = sqrt(omtot);
    }
    double ratetot = 0.0;
    double count = 0.0;
    double dampTime;
    double dampFactor = 0.0;
    double dampLen;
    double theta;
    if (DAMP_INNER > 0){
      dampTime = dampTimeInner;
      if (DAMP_INNER == 2) dampTime = dampTime/omtot;
      dampLen = dampLenInner;
      theta = (x[0]-rmin)/dampLen;
      if (theta > 1.0) dampFactor = 0.0;
      else dampFactor = 1.0 - pow(1.0 - pow(theta,2.0), 2.0);
      ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
      count = count + 1.0;
    }
    if (DAMP_OUTER > 0){
      dampTime = dampTimeOuter;
      if (DAMP_OUTER == 2) dampTime = dampTime/omtot;
      dampLen = dampLenOuter;
      theta = (rmax-x[0])/dampLen;
      if (theta > 1.0) dampFactor = 0.0;
      else dampFactor = pow(1.0 - pow(theta,2.0), 2.0);
      ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
      count = count + 1.0;
    }
    if (twoD == 0){
      if (DAMP_UPPER > 0){
        dampTime = dampTimeUpper;
        if (DAMP_UPPER == 2) dampTime = dampTime/omtot;
        dampLen = dampLenUpper;
        theta = (zmax-x[2])/dampLen;
        if (theta > 1.0) dampFactor = 0.0;
        else dampFactor = pow(1.0 - pow(theta,2.0), 2.0);
        ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
        count = count + 1.0;
      }
      if (DAMP_LOWER > 0){
        dampTime = dampTimeLower;
        if (DAMP_LOWER == 2) dampTime = dampTime/omtot;
        dampLen = dampLenLower;
        theta = (zmax-x[2])/dampLen;
        if (theta > 1.0) dampFactor = 0.0;
        else dampFactor = 1.0-pow(1.0 - pow(theta,2.0), 2.0);
        ratetot = (count*ratetot + dampFactor/dampTime)/(1.0+count);
        count = count + 1.0;
      }
    }
    dampFactor = expm1(-dt*dampFactor);
    //dampFactor = -1.0*dt*dampFactor;
    double prims0[NUM_Q];
    double cons0[NUM_Q];
    double cons1[NUM_Q];
    initial(prims0, x);
    prim2cons(prims0, cons0, x, dV);
    prim2cons(prim, cons1, x, dV);
    cons[DDD] += (cons1[DDD] - cons0[DDD])*dampFactor;
    cons[SRR] += (cons1[SRR] - cons0[SRR])*dampFactor;
    cons[LLL] += (cons1[LLL] - cons0[LLL])*dampFactor;
    cons[SZZ] += (cons1[SZZ] - cons0[SZZ])*dampFactor;
    cons[TAU] += (cons1[TAU] - cons0[TAU])*dampFactor;
  }
}
