#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"
#include "magnitudes.h"


// might end up changing the scope of some of these variables, like mbh, if all this gets moved to a galaxy object
double Fraction_Observed(double mstar, double L_c, Galaxy gal, Disruption disrupt)
{

  double L_max = disrupt.Max_Luminosity(mstar);

  if (L_max <= L_c) return 0;

  double logLmax = log10(L_max);
  if (logLmax <= MIN_LOG_LBOL) return 0;
  
  double x_c = log10(L_c) - logLmax;
  double x_min = MIN_LOG_LBOL - logLmax;

  if (x_c < x_min) return 1;

  return ((pow(10., -1. * LF_LOG_POWERLAW * x_c) - 1.)/(pow(10., -1. * LF_LOG_POWERLAW * x_min) - 1.));


}


// could this be made faster by using a better "proposal distribution?". 
double Rejection_Sample_Mstar(gsl_rng *rangen, Galaxy gal)
{
    double mstar_min = MSTAR_MIN; //really will want to allow this to be different in different galaxies
    double mstar_max = MSTAR_MAX;
  
    double imf_norm = gal.Get_imf_norm();

    // careful about this if you end up using different IMF that is not monotonically decreasing with stellar mass
    double envelope = gal.Kroupa_IMF_for_value(mstar_min, imf_norm);

  while (true)
  {
    // pick a random mstar between minimum and maximum allowed values, in units of solar mass
    // this is the "proposal distribution"
    double mstar_e = mstar_min + (mstar_max - mstar_min)*gsl_rng_uniform(rangen);

    // calculate the probability from the present day mass function
    double P = 1./envelope * gal.Kroupa_IMF_for_value(mstar_e, imf_norm);
    if (gsl_rng_uniform(rangen) < P) return mstar_e;
  }


}


//Not rejection sampling ... inverse transform sampling
double Sample_Peak_L(gsl_rng *rangen, double mstar, double L_c, Galaxy gal, Disruption disrupt)
{

  // work in terms of x = log10(L) - log10(Lmax)

  double L_max = disrupt.Max_Luminosity(mstar);  // will also be used to convert the sampled x to a phsyical L

  double ymin = pow( pow(10., MIN_LOG_LBOL)/L_max ,LF_LOG_POWERLAW);

  double u = gsl_rng_uniform(rangen);

  double this_y = ymin/(1. - u * (1. - ymin));

  return L_max * pow(this_y,1./LF_LOG_POWERLAW);
}


void Sample_Disruption_Parameters(gsl_rng *rangen, Galaxy gal,  double& vol_rate_accumulator, double& detected_rate_accumulator)
{

  double m_limit_contrast = find_host_contrast_magnitude(gal); // in the future, maybe specify which band this is for

  double z = gal.Get_z();
  double mbh = gal.Get_Mbh();

  // later will need to move this into the event generation loop. 
  //  double L_c = LCriticalRband(gal,T, z, m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int num_trials = 100;
  vol_rate_accumulator = 0;
  detected_rate_accumulator = 0;

  for (int i = 0; i < num_trials; ++i)
    {

      Disruption disrupt(gal);
      
      // sample mstar
      double mstar = Rejection_Sample_Mstar(rangen, gal);

      double mhills = disrupt.Hills_Mass(mstar);

      if (mbh > mhills) continue;
      else
	{
	  // will then have to do the flare observability criteria for each flare when doing detected disrupt. For volumteric disrupt, accept all of these

	  vol_rate_accumulator += 1.;

	  double T = disrupt.Get_T();
	  double L_c = LCriticalRband(gal,T,z,m_limit_contrast); // careful with the band here
	  double this_peak_L = Sample_Peak_L(rangen,mstar, L_c, gal, disrupt);

	  if (this_peak_L > L_c)
	    {
	      detected_rate_accumulator += 1.;
	    }

	}

    }

  // this can be done at the end to save time
  // normalize by number of trials
  vol_rate_accumulator /= (double) num_trials;
  vol_rate_accumulator *= RATE_NORMALIZATION_COMBINED * pow(gal.Get_nuker_gammaprime()/0.4,RATE_POWERLAW_NUKER); // rate is still in galaxy proper time
  vol_rate_accumulator *=  1./(1. + z);


    // this can be done at the end to save time
  // normalize by number of trials
  detected_rate_accumulator /= (double) num_trials;
  detected_rate_accumulator *= RATE_NORMALIZATION_COMBINED * pow(gal.Get_nuker_gammaprime()/0.4,RATE_POWERLAW_NUKER); // rate is still in galaxy proper time
  detected_rate_accumulator *=  1./(1. + z);


}

