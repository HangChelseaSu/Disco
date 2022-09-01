
#include "../paul.h"

static double gamma_law = 0.0;

void setDiagParams( struct domain * theDomain ){
   gamma_law = theDomain->theParList.Adiabatic_Index;
}

int num_diagnostics(void){
   return(15);
}

int num_inst_diagnostics(void){
   return(0);
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
    double Pp = prim[PPP];
    double vp = r*omega;

    Qrz[0] = rho;           // Sigma
    Qrz[1] = rho*vr;        // Mdot
    Qrz[2] = rho*fabs(vr);  // Mdot |vr|
    Qrz[3] = rho * r*vp;    // angular momentum
    Qrz[4] = rho*r*vp*vr;	// (advective) angular momentum flux
    Qrz[5] = rho*r*vp*fabs(vr); 
    Qrz[6] = Pp;

    double cosp = cos(phi);
    double sinp = sin(phi);

    double vx = vr*cosp - vp*sinp;
    double vy = vr*sinp + vp*cosp;
    double v2 =  vr*vr + vp*vp;
    double rdv = r * vr;
    Qrz[7] = rho * ((r*v2 - 1)*cosp - rdv*vx);	//e_x
    Qrz[8] = rho * ((r*v2 - 1)*sinp - rdv*vy);	//e_y

    double Fr, Fp, Fz;

    planetaryForce( theDomain->thePlanets + 0, r, phi, z, &Fr, &Fp, &Fz, 0);
    Qrz[9] = rho * r * Fp;				//Torque density from pl 0
    planetaryForce( theDomain->thePlanets + 1, r, phi, z, &Fr, &Fp, &Fz, 0);
    Qrz[10] = rho * r * Fp;				//Torque density from pl 1

    double cos2p = (cosp - sinp) * (cosp + sinp);
    double sin2p = 2*sinp*cosp;
    Qrz[11]  = rho * r * cosp;
    Qrz[12]  = rho * r * sinp;
    Qrz[13]  = rho * r*r * cos2p;
    Qrz[14] = rho * r*r * sin2p;
}


void get_inst_diagnostics( double * x , double * prim , double * Qrz, 
                        struct domain * theDomain )
{
    //Silence is golden.
}
