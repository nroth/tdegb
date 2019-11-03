#include <math.h>
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"


Disruption::Disruption(Galaxy gal)
{

  host_gal = gal;
  mbh = gal.Get_Mbh();
  L_Edd = Eddington_Luminosity();

  T_opt = 3.e4;
  beta = 1.;
  mstar = 1.;
  max_L = 0.;
  A_V = 0.;
  R_V = 3.;
  
}

// Rejection sampling to avoid messiness of integrating piecewise function
// could this be made faster by using a better "proposal distribution?".
void Disruption::Rejection_Sample_Mstar(gsl_rng *rangen)
{

    double mstar_min = MSTAR_MIN; //really will want to allow this to be different in different galaxies
    double mstar_max = MSTAR_MAX;
  
    double imf_norm = host_gal.Get_imf_norm();

    // careful about this if you end up using different IMF that is not monotonically decreasing with stellar mass
    double envelope = host_gal.Kroupa_IMF_for_value(mstar_min, imf_norm);

  while (true)
  {
    // pick a random mstar between minimum and maximum allowed values, in units of solar mass
    // this is the "proposal distribution"
    double mstar_e = mstar_min + (mstar_max - mstar_min)*gsl_rng_uniform(rangen);

    // calculate the probability from the present day mass function
    double P = 1./envelope * host_gal.Kroupa_IMF_for_value(mstar_e, imf_norm);
    if (gsl_rng_uniform(rangen) < P)
      {
	mstar = mstar_e;
	return;
      }
  }
    
}

double Disruption::Get_Topt()
{
  return T_opt;
}


double Disruption::Get_beta()
{
  return beta;
}


double Disruption::Get_Mstar()
{
  return mstar;
}

double Disruption::Get_Max_L()
{
  return max_L;
}


// This could be improved
double Disruption::Main_Sequence_Radius()
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
double Disruption::Get_Hills_Mass()
{
  double rstar = Main_Sequence_Radius();
    
  double relativistic_correction = sqrt(5.) / 8. / pow(2.,-1.5); // beloborodov 1992
  return (relativistic_correction * pow(rstar * C_LIGHT * C_LIGHT /(2. * NEWTON_G * pow(mstar * M_SUN,1./3.)) , 3./2.)/M_SUN);
}

// assumes mbh in solar masses
double Disruption::Eddington_Luminosity()
{
  
  return 4. * PI * NEWTON_G * mbh * M_SUN * M_PROTON * C_LIGHT/(THOMSON_CS);
  
}

// mstar and mbh in solar masses
// mdot in cgs
// Need to think harder about how to handle polytropic index dependence on stellar mass
double Disruption::Peak_Mdot()
{

  if (mstar <= 1.)
    {
      // gamma = 5./3.
      double guillochon_A = exp( (10.253 - 17.380 * beta + 5.9988 * pow(beta,2. ) )/ (1. - 0.46573 * beta - 4.5066 * pow(beta,2. )));

      double rstar = Main_Sequence_Radius();

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
void Disruption::Determine_Max_L()
{

  double mdot_fallback_peak = Peak_Mdot();
    
  double L_fallback_peak = MAX_RADEFF_STREAMS * mdot_fallback_peak * C_LIGHT * C_LIGHT;
    
  if (L_fallback_peak >  MAX_EDD_RATIO  * L_Edd)
    max_L =  MAX_EDD_RATIO  * L_Edd;
    
  else
    max_L =  L_fallback_peak;

  return;
  
}



//Not rejection sampling, inverse transform sampling
double Disruption::Sample_Peak_L(gsl_rng *rangen)
{

  double ymin = pow( pow(10., MIN_LOG_LBOL)/max_L ,LF_LOG_POWERLAW);

  double u = gsl_rng_uniform(rangen);

  double this_y = ymin/(1. - u * (1. - ymin));

  return max_L * pow(this_y,1./LF_LOG_POWERLAW);
}



