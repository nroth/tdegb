#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_rng.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"
#include "magnitudes.h"
#include "histogramNd.h"


// Rejection sampling to avoid messiness of integrating piecewise function
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


//Not rejection sampling, inverse transform sampling
double Sample_Peak_L(gsl_rng *rangen, double L_max)
{

  double ymin = pow( pow(10., MIN_LOG_LBOL)/L_max ,LF_LOG_POWERLAW);

  double u = gsl_rng_uniform(rangen);

  double this_y = ymin/(1. - u * (1. - ymin));

  return L_max * pow(this_y,1./LF_LOG_POWERLAW);
}


void Sample_Disruption_Parameters(gsl_rng *rangen, Galaxy gal, double& vol_rate_accumulator, double& detected_rate_accumulator, HistogramNd& hist_detected_flares)
{

  double m_limit_contrast = find_host_contrast_magnitude(gal); // in the future, maybe specify which band this is for

  double operating_m_limit = std::min(m_limit_contrast,M_APPARENT_THRESHHOLD);


  double mbh = gal.Get_Mbh();
  double z = gal.Get_z();

  double cosmo_factor = (1. + z)/(4. * PI * pow(LuminosityDistance(z),2.));
  double nu_g_emit = (1. + z) * NU_GBAND;
  double nu_r_emit = (1. + z) * NU_RBAND;  

  // later will need to move this into the event generation loop. 
  //  double L_c = LCriticalRband(gal,T, z, m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int num_trials = 500;
  vol_rate_accumulator = 0.;
  detected_rate_accumulator = 0.;

  double rate_normalization = 1./( (double) num_trials) * RATE_NORMALIZATION_COMBINED * pow(gal.Get_nuker_gammaprime()/0.4,RATE_POWERLAW_NUKER) * 1./(1. + z);

  vector<double> flare_properties(hist_detected_flares.Get_Dimension());

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

      	  double L_max = disrupt.Max_Luminosity(mstar);  // will also be used to convert the sampled x to a phsyical L
	  if (log10(L_max) < MIN_LOG_LBOL) continue;

	  double T = disrupt.Get_T(); // randomly generate
	  double R_V = 3.; // randomly generate? or assing as function of galaxy properties?
	  double A_V = 0.; // will want to randomly generate

	  double this_peak_L = Sample_Peak_L(rangen,L_max);
	    
	  double r_mag_observed = mABFromFnu(ExtinctedFluxObserved(nu_r_emit,cosmo_factor,this_peak_L,T,A_V,R_V));

	  if (r_mag_observed < operating_m_limit)
	    {
	      detected_rate_accumulator += 1.;

	      double g_mag_observed = mABFromFnu(ExtinctedFluxObserved(nu_g_emit,cosmo_factor,this_peak_L,T,A_V,R_V));

	      //	      double flare_properties[5] = {m_g,m_r,z,Lbol_peak,mbh};
	      flare_properties[0] = g_mag_observed;
	      flare_properties[1] = g_mag_observed - r_mag_observed;
	      flare_properties[2] = z;
	      flare_properties[3] = log10(this_peak_L);
	      flare_properties[4] = log10(mbh);

	      hist_detected_flares.Count(flare_properties,rate_normalization);
	    }



	}

    }

  // this can be done at the end to save time
  // normalize by number of trials
  vol_rate_accumulator *=  rate_normalization;


    // this can be done at the end to save time
  // normalize by number of trials
  detected_rate_accumulator *= rate_normalization;


}

