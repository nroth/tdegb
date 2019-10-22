#include <math.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "galaxy.h"


using std::vector;

//***************************************************************
// Constructors
//***************************************************************

void Individual_Disruption::Initialize(double m, double T, double b)
{

  mbh = m;
  temperature = T;
  beta = b;
}

double Individual_Disruption::Get_Temperature()
{

  return temperature;
}

// This could be improved
double Individual_Disruption::Main_Sequence_Radius(double mstar)
{

  if (mstar <= 1.)
    return R_SUN * pow(mstar,0.8);
  else 
    return R_SUN * pow(mstar,0.57);
  
}

// ignoring BH spin
// See first paragraph section 2.2 of Stone & Metzger 2016.
// The Newtonian formula in Stone & Metzger requires relativistic correction, as in Beloborodov 1992 (involves IBSCO)
// Final result is consistent with what is quoted in Leloudas et al 2016
// mstar in solar mass, returns a mass in solar masses
double Individual_Disruption::Hills_Mass(double mstar)
{
  double rstar = Main_Sequence_Radius(mstar);
    
  double relativistic_correction = sqrt(5.) / 8. / pow(2.,-1.5); // beloborodov 1992
  return (relativistic_correction * pow(rstar * C_LIGHT * C_LIGHT /(2. * NEWTON_G * pow(mstar * M_SUN,1./3.)) , 3./2.)/M_SUN);
}

// assumes mbh in solar masses
double Individual_Disruption::Eddington_Luminosity()
{
  
  return 4. * PI * NEWTON_G * mbh * M_SUN * M_PROTON * C_LIGHT/(THOMSON_CS);
  
}

// mstar and mbh in solar masses
// mdot in cgs
// Need to think harder about how to handle polytropic index dependence on stellar mass
double Individual_Disruption::Peak_Mdot(double mstar)
{

  if (mstar <= 1.)
    {
      // gamma = 5./3.
      double guillochon_A = exp( (10.253 - 17.380 * beta + 5.9988 * pow(beta,2. ) )/ (1. - 0.46573 * beta - 4.5066 * pow(beta,2. )));

      double rstar = Main_Sequence_Radius(mstar);

      //printf("beta is %f, rstar in solar units is %f\n",beta,rstar/R_SUN);

      // mstar already in solar units, but rstar in cgs
      return guillochon_A * pow(mbh/1.e6,-0.5) * pow(mstar,2.) * pow(rstar/R_SUN,-1.5) * M_SUN / YEAR_TO_SEC;
    }

  else
    {
      printf("ERROR: HAVE NOT DEFINED MDOT YET FOR MASSIVE STARS\n");
      return 0.;
    }
  
}

// at some point might want to think about other ways beta might enter here, aside from setting maximum mass fallback rate
double Individual_Disruption::Max_Luminosity(double mstar)
{

  double Ledd = Eddington_Luminosity();
    
  double mdot_fallback_peak = Peak_Mdot(mstar);
    
  double L_fallback_peak = MAX_RADEFF_STREAMS * mdot_fallback_peak * C_LIGHT * C_LIGHT;

    
  if (L_fallback_peak >  MAX_EDD_RATIO  * Ledd)
    return MAX_EDD_RATIO  * Ledd;
    
  else
    return L_fallback_peak;
    
  
}






