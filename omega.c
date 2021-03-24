
#include "paul.h"
#include "omega.h"

static int meshOmChoice = 0;
static double meshOmPar = 0.0;
static int enOmChoice = 0;
static double enOmPar = 0.0;
static int cs2Choice = 0;
static double cs2Par = 0.0;
static double Omega0 = 0.0;
static int Npl = 0;

static double Mach = 0.0;
static double r0 = 0.0;
static double r1 = 0.0;
static double r2 = 0.0;
static double H0 = 0.0;
static double M = 0.0;
static double eps = 0.0;

static struct planet *thePlanets = NULL;


double phigrav( double , double , double , int); //int here is type
double fgrav( double , double , double , int); //int here is type


void setOmegaParams( struct domain * theDomain ){
   meshOmChoice = theDomain->theParList.Exact_Mesh_Omega;
   meshOmPar    = theDomain->theParList.Exact_Mesh_Omega_Par;
   enOmChoice = theDomain->theParList.Energy_Omega;
   enOmPar = theDomain->theParList.Energy_Omega_Par;
   cs2Choice = theDomain->theParList.Cs2_Profile;
   cs2Par = theDomain->theParList.Cs2_Par;
   Omega0 = theDomain->theParList.RotOmega;

   Mach = theDomain->theParList.Disk_Mach;
   r0 = theDomain->theParList.initPar1; // Inner edge
   r1 = theDomain->theParList.initPar2; // Fiducial radius
   r2 = theDomain->theParList.initPar3; // Outer Edge
   H0 = theDomain->theParList.initPar4; // Scale Height
   M = theDomain->theParList.metricPar2;
   Npl = theDomain->Npl;
   eps = theDomain->theParList.grav_eps;

   thePlanets = theDomain->thePlanets;

   if(strcmp(PLANETS, "bin_rot") == 0)
   {
       double q = theDomain->theParList.Mass_Ratio;
       M = q*1.0/(1.0+q);
   }
}

double mesh_om( const double *x)
{
    double r = x[0];
    r = sqrt(r*r + eps*eps);
    double omega;
    if(meshOmChoice == 1)
        omega = 1.0;
    else if(meshOmChoice == 2)
        omega = pow(r,-1.5);
    else if(meshOmChoice == 3)
    {
        double n = 4.0;
        omega = 1./pow( pow( r , 1.5*n ) + 1. , 1./n );
    }
    else
        omega = 0.0;
        
   return( omega );
}

double get_om( const double *x ){
    double r = x[0];
    r = sqrt(r*r + eps*eps);
    double om;

    if(enOmChoice == 1)
        om = 1.0-Omega0;

    else if(enOmChoice == 2)
        om = 1.0/pow(r,1.5)-Omega0;

    else if(enOmChoice == 3)
    {
        double n = 8.0;
        om = 1./pow( pow( r , 1.5*n ) + 1. , 1./n )-Omega0;
    }

    else if(enOmChoice == 4)
        om =  10.*exp(-.5*r*r)-Omega0;

    else
        om = 0.0-Omega0;

    return om;
}

double get_om1( const double *x){
    double r = x[0];
    r = sqrt(r*r + eps*eps);
    double om1;

    if(enOmChoice == 1)
        om1 = 0.0;

    else if(enOmChoice == 2)
        om1 = -1.5/pow(r,2.5);

    else if(enOmChoice == 3)
    {
        double n = 8.0;
        om1 = -1.5 * pow(r,-1+1.5*n) / pow( pow(r,1.5*n) + 1. , 1.0+1./n );
    }

    else if(enOmChoice == 4)
        om1 =  10.*r*exp(-.5*r*r);

    else
        om1 = 0.0;

    return om1;
}
  
double get_om2( const double *x){
    return 0.0;
}

double get_height_om( const double *x){
    double omtot2 = 0.0;
    int pi;

    double r = x[0];
    double phi = x[1];

    double cosp = cos(phi);
    double sinp = sin(phi);
    double gx = r*cosp;
    double gy = r*sinp;

    double px, py, script_r;
    for (pi=0; pi<Npl; pi++){
      cosp = cos(thePlanets[pi].phi);
      sinp = sin(thePlanets[pi].phi);
      px = thePlanets[pi].r*cosp;
      py = thePlanets[pi].r*sinp;
      script_r = sqrt((px-gx)*(px-gx) + (py-gy)*(py-gy));
      double f1 = fgrav( thePlanets[pi].M , script_r , thePlanets[pi].eps , thePlanets[pi].type);
      omtot2 += f1/(script_r + 0.0625*thePlanets[pi].eps);	//Technically, should not include any extra softening, but avoid dividing by zero
    }
    return sqrt(omtot2);
}

double get_cs2( const double *x ){
    double r = x[0];
    double cs2;

    if(cs2Choice == 1)
        cs2 = 1./(Mach*Mach);
    else if(cs2Choice == 2)
    {
        double nu = .5;
        cs2 = .5/Mach/Mach/pow(r,2.*nu);
    }
    else if(cs2Choice == 3)
    {
        cs2 = M*H0*H0 / (2*r1*r1*r1);
    }
    else if(cs2Choice == 4)
    {
        r = sqrt(r*r + eps*eps);
        double v2 = M/r;
        cs2 = v2/(Mach*Mach);
    }
    else if(cs2Choice == 5) 
    {
      double r = x[0];
      double phi = x[1];

      double cosp = cos(phi);
      double sinp = sin(phi);
      double gx = r*cosp;
      double gy = r*sinp;
      int pi;
      double px, py, pr;
      double n = 2.0;
      double phip = 0.0;
 
      for (pi = 0; pi<Npl; pi++)
      {
        cosp = cos(thePlanets[pi].phi);
        sinp = sin(thePlanets[pi].phi);
        px = thePlanets[pi].r*cosp;
        py = thePlanets[pi].r*sinp;
        pr = (px-gx)*(px-gx) + (py-gy)*(py-gy);
        //phip += thePlanets[pi].M/pow( pr + pow(thePlanets[pi].eps,n) , 1./n );
        phip += phigrav( thePlanets[pi].M , sqrt(pr) , thePlanets[pi].eps , thePlanets[pi].type );
      }
      cs2 = phip/(Mach*Mach);        
    }
    else
        cs2 = 1.0;

    return cs2;
}

