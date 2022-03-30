
#include "../paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(11);
}

void planetaryForce( struct planet * , double , double , double , double * , double * , double * , int );

/* Generic Diagnostics for Euler*/

void get_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
   double r = x[0];
   double phi = x[1];
   double z = x[2];

   double rho = prim[RHO];
   double vr = prim[URR];
   double omega = prim[UPP];
   //double vz = prim[UZZ];
   //double Pp = prim[PPP];
   double vp = r*omega;

   Qrz[0] = rho;	//Sigma
   Qrz[1] = rho*vr; //Mdot
   Qrz[6] = rho*fabs(vr); //Mdot

   double cosp = cos(phi);
   double sinp = sin(phi);

   double vx = vr*cosp - vp*sinp;
   double vy = vr*sinp + vp*cosp;
   double v2 =  vr*vr + vp*vp;
   double rdv = r * vr;
   Qrz[2] = rho * ((r*v2 - 1)*cosp - rdv*vx);	//e_x
   Qrz[3] = rho * ((r*v2 - 1)*sinp - rdv*vy);	//e_y

   Qrz[4] = 0.0;
   Qrz[5] = rho*vp;

   double Fr,Fp,Fz,rp;
   Fp = 0.0;

   int pi;
   for (pi=0; pi<2; pi++) {
     struct planet * pl = theDomain->thePlanets + pi;
     planetaryForce( pl , r , phi , z , &Fr , &Fp , &Fz , 1 );
     rp = theDomain->thePlanets[pi].r;
     Qrz[4] += rho*rp*Fp;				//Torque density
   }

   double cos2p = (cosp - sinp) * (cosp + sinp);
   double sin2p = 2*sinp*cosp;
   Qrz[7]  = rho * cosp;
   Qrz[8]  = rho * sinp;
   Qrz[9]  = rho * cos2p;
   Qrz[10] = rho * sin2p;
}

