#include <math.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "galaxy.h"
//#include "magnitudes.h"


using std::vector;

//***************************************************************
// Constructors
//***************************************************************

Galaxy::Galaxy()
{
 
  imf_normalization = 1.;
 
  disruptions.Initialize(1.e6,3.e4, 1.);

}



Galaxy::Galaxy(double* galaxy_info, double imf_norm, double temperature, double beta)
{
  total_stellar_mass = galaxy_info[0];
  mbh_sigma = log10(galaxy_info[1]);
  mbh_bulge = log10(galaxy_info[2]);
  z = galaxy_info[3];
  sersic_n = galaxy_info[4];
  r50_kpc = galaxy_info[5];
  nuker_gammaprime = Find_Nuker_Gammaprime_From_Sersic(sersic_n, r50_kpc, z);

  imf_normalization = imf_norm;

  mbh = mbh_sigma; // making this choice for now

  disruptions.Initialize(mbh,temperature, beta);

}

double Galaxy::Get_Mbh()
{
  return mbh;
}



// See equation 7 of https://arxiv.org/pdf/1212.0939.pdf
// the params has the overall normalization, so as not to recompute it each time
double Galaxy::Kroupa_IMF(double mstar, void * params)
{
  Galaxy * this_gal = (Galaxy *) params;
  double norm = this_gal -> imf_normalization;
  
  double m1 = 0.08;
  double m2 = 0.5;
  double m3 = MSTAR_MAX;

  double k1 = pow(m1,-0.3 + 1.3);
  double k2 = k1 * pow(m2,-1.3 + 2.3);
  //  double k3 = k2 * pow(m3, -2.3 + 2.3);

  if (mstar >= m1 && mstar < m2)
    return 1./norm * k1 * pow(mstar,-1.3);
  else if (mstar >= m2 && mstar <= m3)
      return 1./norm * k2 * pow(mstar,-2.3);
  else
    return 0.;
  
}

// both mbh and mstar in solar masses
double Galaxy::Stellar_Disruption_Rate(double mstar, Galaxy *gal)
{

  double mbh = gal->Get_Mbh();
  //consider using the exp(-M^2) suppression here as suggested by Kesden 2012, instead of a sharp and total cutoff
  double mhills = gal->disruptions.Hills_Mass(mstar);
  if (mbh > mhills) return 0;

  //  return RATE_NORMALIZATION_COMBINED * pow(mbh/1.e8,RATE_POWERLAW); // rate is still in galaxy proper time
  return RATE_NORMALIZATION_COMBINED * pow(gal->nuker_gammaprime/0.4,RATE_POWERLAW_NUKER); // rate is still in galaxy proper time
  
}


double Galaxy::Stellar_Rate_Integrand(double mstar, void * p)
{

  Galaxy * this_gal = (Galaxy *) p;
  return Stellar_Disruption_Rate(mstar, this_gal) * Kroupa_IMF(mstar, this_gal); // Stellar disruption rate in terms of galaxy proper time
  
}



double Galaxy::Determine_IMF_Normalization()
{


  int grid_size = 128;
  double relative_error = 1.e-6;

  gsl_function F;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  F.function = Galaxy::Kroupa_IMF;
  F.params = this;

  double result, error;
  gsl_integration_qags(&F,MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);
  gsl_integration_workspace_free(workspace);
  
  imf_normalization = result;
  return imf_normalization;

}

/*
void Initialize_IMF()
{ 
   Determine_IMF_Normalization(); // have to be careful about where this is included
}
*/

double Galaxy::Total_Disruption_Rate()
{
 
  int grid_size = 128;
  double relative_error = 1.e-6;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Galaxy::Stellar_Rate_Integrand;
  F.params = this;

  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units


}


// might end up changing the scope of some of these variables, like mbh and z, if all this gets moved to a galaxy object
double Galaxy::Stellar_Rate_Integrand_GbandCut(double mstar, void * p)
{

  struct stellar_rate_params {Galaxy * gal; double L_c;};
  struct stellar_rate_params * params = (struct stellar_rate_params *)p;

  Galaxy * this_gal = params->gal;
  double L_c = params->L_c;

  return Stellar_Disruption_Rate(mstar, this_gal) * Kroupa_IMF(mstar, this_gal) * Fraction_Observed(mstar, L_c, this_gal); // rate of observable disruptions in terms of galaxy proper time
  
}

// might end up changing the scope of some of these variables, like mbh and z, if all this gets moved to a galaxy object
double Galaxy::Stellar_Rate_Integrand_RbandCut(double mstar, void * p)
{

  struct stellar_rate_params {Galaxy * gal; double L_c;};
  struct stellar_rate_params * params = (struct stellar_rate_params *)p;

  Galaxy * this_gal = params->gal;
  double L_c = params->L_c;

  // I don't think Kroupa_IMF currently cares what param pointer is passed to it
  return Stellar_Disruption_Rate(mstar, this_gal) * Kroupa_IMF(mstar, this_gal) * Fraction_Observed(mstar, L_c, this_gal); // rate of observable disruptions in terms of galaxy proper time
  
}



double Galaxy::Total_Disruption_Rate_Observed_Gband(double m_limit_contrast)
{
 

  //  printf("IMF norm is %e\n",IMF_normalization);

  double L_c = LCriticalGband(z,disruptions.Get_Temperature(),m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int grid_size = 128;
  double relative_error = 1.e-5;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Galaxy::Stellar_Rate_Integrand_GbandCut;

  struct stellar_rate_params {Galaxy * gal; double this_L_c;};

  stellar_rate_params params = {this, L_c};
  F.params = &params;

  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units



}

double Galaxy::Total_Disruption_Rate_Observed_Rband(double m_limit_contrast)
{
 

  //  printf("IMF norm is %e\n",IMF_normalization);

  double L_c = LCriticalRband(z,disruptions.Get_Temperature(),m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int grid_size = 128;
  double relative_error = 1.e-5;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Galaxy::Stellar_Rate_Integrand_RbandCut;

  struct stellar_rate_params {Galaxy * gal; double this_L_c;};

  stellar_rate_params params = {this, L_c};
  F.params = &params;

  //  printf("params are %e %f %f %f\n",mbh,beta,z,T);

  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units



}


// might end up changing the scope of some of these variables, like mbh, if all this gets moved to a galaxy object
double Galaxy::Fraction_Observed(double mstar, double L_c, Galaxy * gal)
{

  double L_max = gal->disruptions.Max_Luminosity(mstar);

  if (L_max <= L_c) return 0;

  double logLmax = log10(L_max);
  if (logLmax <= MIN_LOG_LBOL) return 0;
  
  double x_c = log10(L_c) - logLmax;
  double x_min = MIN_LOG_LBOL - logLmax;

  if (x_c < x_min) return 1;

  return ((pow(10., -1. * LF_LOG_POWERLAW * x_c) - 1.)/(pow(10., -1. * LF_LOG_POWERLAW * x_min) - 1.));


}





