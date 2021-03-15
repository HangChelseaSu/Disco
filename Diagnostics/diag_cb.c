
#include "../paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(5);
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
   Qrz[1] = 2.*M_PI*r*rho*vr; //Mdot

   double gx = r*cos(phi);
   double gy = r*sin(phi);
   double vx = vr*cos(phi) - vp*sin(phi);
   double vy = vr*sin(phi) + vp*cos(phi);
   double vm2 =  vx*vx + vy*vy;
   double rdv = vx*gx + vy*gy;
   Qrz[2] = (vm2 - 1./r)*gx - rdv*vx;	//e_x
   QrZ[2] *= rho;
   Qrz[3] = (vm2 - 1./r)*gy - rdv*vy;	//e_y
   QrZ[3] *= rho;


   double Fr,Fp,Fz,rp;
   Fp = 0.0;

   int pi;
   for (pi=0; pi<2; pi++) {
     struct planet * pl = theDomain->thePlanets + pi;
     planetaryForce( pl , r , phi , z , &Fr , &Fp , &Fz , 0 );
     rp = theDomain->thePlanets[pi].r;
     Qrz[4] += rho*rp*Fp;				//Torque density
   }


}

