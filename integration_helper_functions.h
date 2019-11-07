#include <math.h>
#include <gsl/gsl_rng.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"
#include "histogramNd.h"
#include "survey.h"


void Sample_Disruption_Parameters(gsl_rng *rangen, Survey surv, Galaxy gal, double& vol_rate_accumulator, double& detected_rate_accumulator, HistogramNd& hist_detected_flares)
{

  double m_r_limit_contrast = surv.Find_Host_Contrast_Magnitude(gal,'r');
  double m_g_limit_contrast = surv.Find_Host_Contrast_Magnitude(gal,'g'); 

  // for all flares in this galaxy
  double operating_m_r_limit = std::min(m_r_limit_contrast,surv.Get_m_r_Threshhold());
  double operating_m_g_limit = std::min(m_g_limit_contrast,surv.Get_m_g_Threshhold());

  double mbh = gal.Get_Mbh();
  double z = gal.Get_z();

  double cosmo_factor = (1. + z)/(4. * PI * pow(gal.Get_Luminosity_Distance(),2.));
  double nu_g_emit = (1. + z) * surv.Get_Nu_Gband();
  double nu_r_emit = (1. + z) * surv.Get_Nu_Rband();  

  int num_trials = 1000;
  vol_rate_accumulator = 0.;
  detected_rate_accumulator = 0.;

  double rate_normalization = 1./( (double) num_trials) * gal.Get_Disruption_Rate_Normalization_Combined() * pow(gal.Get_nuker_gammaprime()/0.4,gal.Get_Disruption_Rate_Powerlaw_Nuker()) * 1./(1. + z);

  vector<double> flare_properties(hist_detected_flares.Get_Dimension());

  Disruption disrupt(gal); // default values filled in now

  for (int i = 0; i < num_trials; ++i)
    {

      // sample mstar
      disrupt.Rejection_Sample_Mstar(rangen);
      double mhills = disrupt.Get_Hills_Mass();

      if (mbh > mhills) continue;
      else
	{
	  // will then have to do the flare observability criteria for each flare when doing detected disrupt. For volumteric disrupt, accept all of these

	  vol_rate_accumulator += 1.;

	  disrupt.Sample_Beta(rangen);
	  
	  // need to haave sampled mstar, beta (and mbh) already 
	  disrupt.Determine_Max_L();
      	  double max_L = disrupt.Get_Max_L();  // will also be used to convert the sampled x to a phsyical L
	  if (log10(max_L) < disrupt.Get_Min_Log_Lbol()) continue;
	  //	  double T_opt = disrupt.Get_Topt(); // randomly generate
	  //	  double R_V = 3.; // randomly generate? or assing as function of galaxy properties?
	  //	  double A_V = 0.; // will want to randomly generate


	  disrupt.Sample_Peak_L(rangen);
	  disrupt.Sample_Topt(rangen);
	  disrupt.Sample_R_V(rangen);
	  disrupt.Sample_A_V(rangen);
	    
	  double r_mag_observed = surv.mAB_From_Fnu(disrupt.Extincted_Flux_Observed(nu_r_emit,cosmo_factor));
	  double g_mag_observed = surv.mAB_From_Fnu(disrupt.Extincted_Flux_Observed(nu_g_emit,cosmo_factor));

	  if (r_mag_observed < operating_m_r_limit && g_mag_observed < operating_m_g_limit)
	    {
	      detected_rate_accumulator += 1.;

	      flare_properties[0] = r_mag_observed;
	      flare_properties[1] = g_mag_observed - r_mag_observed;
	      flare_properties[2] = z;
	      flare_properties[3] = log10(disrupt.Get_Peak_L());
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

