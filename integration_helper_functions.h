#include <math.h>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"
#include "magnitudes.h"

// both mbh and mstar in solar masses
double Stellar_Disruption_Rate(double mstar, Galaxy gal, Disruption disrupt)
{

  double mbh = gal.Get_Mbh();
  //consider using the exp(-M^2) suppression here as suggested by Kesden 2012, instead of a sharp and total cutoff
  double mhills = disrupt.Hills_Mass(mstar);
  if (mbh > mhills) return 0;

  //  return RATE_NORMALIZATION_COMBINED * pow(mbh/1.e8,RATE_POWERLAW); // rate is still in galaxy proper time
  return RATE_NORMALIZATION_COMBINED * pow(gal.Get_nuker_gammaprime()/0.4,RATE_POWERLAW_NUKER); // rate is still in galaxy proper time
  
}


double Stellar_Rate_Integrand(double mstar, void * p)
{

  //  Galaxy gal = * (Galaxy * gal) p;
    struct stellar_rate_params {Galaxy g; Disruption d;};
    struct stellar_rate_params * params = (struct stellar_rate_params *)p;

    Galaxy gal = params->g ;
    Disruption disrupt = params -> d;
    
  double imf_norm = gal.Get_imf_norm();
  //  if (imf_norm < 0.27 || imf_norm > 0.29) printf("problem with norm\n");
  //  if (imf_norm < 0.29) printf("problem with norm\n");
    
  
  return Stellar_Disruption_Rate(mstar, gal, disrupt) * gal.Kroupa_IMF_for_value(mstar, imf_norm); // Stellar disruption rate in terms of galaxy proper time
  
}


double Total_Disruption_Rate(Galaxy gal, Disruption disrupt)
{
 
  int grid_size = 128;
  double relative_error = 1.e-6;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Stellar_Rate_Integrand;

  struct stellar_rate_params {Galaxy g; Disruption d;};
  stellar_rate_params params = {gal,disrupt};
  F.params = &params;
  //  F.params = &gal;

  double result, error;

  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  double z = gal.Get_z();

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units
}


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


// might end up changing the scope of some of these variables, like mbh and z, if all this gets moved to a galaxy object
double Stellar_Rate_Integrand_RbandCut(double mstar, void * p)
{

  struct stellar_rate_params {Galaxy gal; double L_c; Disruption d;};
  struct stellar_rate_params * params = (struct stellar_rate_params *)p;

  Galaxy gal = params->gal;
  double L_c = params->L_c;
  Disruption disrupt = params->d;

  double imf_norm = gal.Get_imf_norm();

  return Stellar_Disruption_Rate(mstar, gal, disrupt) * gal.Kroupa_IMF_for_value(mstar, imf_norm) * Fraction_Observed(mstar, L_c, gal,disrupt); // rate of observable disruptions in terms of galaxy proper time
  
}


double Total_Disruption_Rate_Observed_Rband(Galaxy gal, Disruption disrupt)
{
 
  // Keeping m_g here for testing purposes, really should be m_r
  double m_limit_contrast = find_host_contrast_magnitude(gal);

  double L_c = LCriticalRband(gal,disrupt,m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int grid_size = 128;
  double relative_error = 1.e-5;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Stellar_Rate_Integrand_RbandCut;

  struct stellar_rate_params {Galaxy gal; double this_L_c; Disruption d;};
  stellar_rate_params params = {gal, L_c,disrupt};
  F.params = &params;

  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  double z = gal.Get_z();

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units



}



