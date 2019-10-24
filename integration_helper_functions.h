#include <math.h>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "magnitudes.h"


#define LF_LOG_POWERLAW 1.5

#define MIN_LOG_LBOL 43.0
#define MAX_RADEFF_STREAMS 0.1 // maximum radiative efficiency of stream collisions.

#define MAX_EDD_RATIO 2.



// This could be improved
double Main_Sequence_Radius(double mstar)
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
double Hills_Mass(double mstar)
{
  double rstar = Main_Sequence_Radius(mstar);
    
  double relativistic_correction = sqrt(5.) / 8. / pow(2.,-1.5); // beloborodov 1992
  return (relativistic_correction * pow(rstar * C_LIGHT * C_LIGHT /(2. * NEWTON_G * pow(mstar * M_SUN,1./3.)) , 3./2.)/M_SUN);
}

// assumes mbh in solar masses
double Eddington_Luminosity(double mbh)
{
  
  return 4. * PI * NEWTON_G * mbh * M_SUN * M_PROTON * C_LIGHT/(THOMSON_CS);
  
}

// mstar and mbh in solar masses
// mdot in cgs
// Need to think harder about how to handle polytropic index dependence on stellar mass
double Peak_Mdot(double mstar, double mbh, double beta)
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
double Max_Luminosity(double mstar, double beta, Galaxy gal)
{

  double mbh = gal.Get_Mbh();

  double Ledd = Eddington_Luminosity(mbh);
    
  double mdot_fallback_peak = Peak_Mdot(mstar,mbh, beta);
    
  double L_fallback_peak = MAX_RADEFF_STREAMS * mdot_fallback_peak * C_LIGHT * C_LIGHT;

    
  if (L_fallback_peak >  MAX_EDD_RATIO  * Ledd)
    return MAX_EDD_RATIO  * Ledd;
    
  else
    return L_fallback_peak;
    
  
}



// both mbh and mstar in solar masses
double Stellar_Disruption_Rate(double mstar, Galaxy gal)
{

  double mbh = gal.Get_Mbh();
  //consider using the exp(-M^2) suppression here as suggested by Kesden 2012, instead of a sharp and total cutoff
  double mhills = Hills_Mass(mstar);
  if (mbh > mhills) return 0;

  //  return RATE_NORMALIZATION_COMBINED * pow(mbh/1.e8,RATE_POWERLAW); // rate is still in galaxy proper time
  return RATE_NORMALIZATION_COMBINED * pow(gal.Get_nuker_gammaprime()/0.4,RATE_POWERLAW_NUKER); // rate is still in galaxy proper time
  
}


double Stellar_Rate_Integrand(double mstar, void * p)
{

  Galaxy gal = *(Galaxy *) p;
  double imf_norm = gal.Get_imf_norm();
  //  if (imf_norm < 0.27 || imf_norm > 0.29) printf("problem with norm\n");
  //  if (imf_norm < 0.29) printf("problem with norm\n");
    
  
  return Stellar_Disruption_Rate(mstar, gal) * gal.Kroupa_IMF_for_value(mstar, imf_norm); // Stellar disruption rate in terms of galaxy proper time
  
}


double Total_Disruption_Rate(Galaxy gal)
{
 
  int grid_size = 128;
  double relative_error = 1.e-6;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Stellar_Rate_Integrand;

  //  struct stellar_rate_params {Galaxy gal;};
  //  stellar_rate_params params = {gal};
  //  F.params = &params;
  F.params = &gal;

  double result, error;

  double mbh = gal.Get_Mbh();
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  double z = gal.Get_z();

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units
}


// might end up changing the scope of some of these variables, like mbh, if all this gets moved to a galaxy object
double Fraction_Observed(double mstar, double L_c, double beta, Galaxy gal)
{

  double L_max = Max_Luminosity(mstar,beta, gal);

  if (L_max <= L_c) return 0;

  double logLmax = log10(L_max);
  if (logLmax <= MIN_LOG_LBOL) return 0;
  
  double x_c = log10(L_c) - logLmax;
  double x_min = MIN_LOG_LBOL - logLmax;

  if (x_c < x_min) return 1;

  return ((pow(10., -1. * LF_LOG_POWERLAW * x_c) - 1.)/(pow(10., -1. * LF_LOG_POWERLAW * x_min) - 1.));


}


// might end up changing the scope of some of these variables, like mbh and z, if all this gets moved to a galaxy object
double Stellar_Rate_Integrand_GbandCut(double mstar, void * p)
{

  struct stellar_rate_params {Galaxy gal; double L_c; double T; double b;};
  struct stellar_rate_params * params = (struct stellar_rate_params *)p;

  Galaxy gal = params->gal;
  double L_c = params->L_c;
  double beta = params->b;
  double imf_norm = gal.Get_imf_norm();

  return Stellar_Disruption_Rate(mstar, gal) * gal.Kroupa_IMF_for_value(mstar, imf_norm) * Fraction_Observed(mstar, L_c, beta, gal); // rate of observable disruptions in terms of galaxy proper time
  
}

// might end up changing the scope of some of these variables, like mbh and z, if all this gets moved to a galaxy object
double Stellar_Rate_Integrand_RbandCut(double mstar, void * p)
{

  struct stellar_rate_params {Galaxy gal; double L_c; double T; double b;};
  struct stellar_rate_params * params = (struct stellar_rate_params *)p;

  Galaxy gal = params->gal;
  double L_c = params->L_c;
  double beta = params->b;
  double imf_norm = gal.Get_imf_norm();

  // I don't think Kroupa_IMF currently cares what param pointer is passed to it
  return Stellar_Disruption_Rate(mstar, gal) * gal.Kroupa_IMF_for_value(mstar, imf_norm) * Fraction_Observed(mstar, L_c, beta, gal); // rate of observable disruptions in terms of galaxy proper time
  
}



double Total_Disruption_Rate_Observed_Gband(Galaxy gal, double T, double b)
{
 
  
  //  double m_limit_contrast = find_host_contrast_magnitude(gal.m_g, gal.sersic_n,gal.r50_kpc,gal.z);
  double m_limit_contrast = find_host_contrast_magnitude(gal);
  
  double L_c = LCriticalGband(gal.Get_z(),T,m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int grid_size = 128;
  double relative_error = 1.e-5;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Stellar_Rate_Integrand_GbandCut;

  struct stellar_rate_params {Galaxy gal; double this_L_c; double temperature; double beta;};
  stellar_rate_params params = {gal, L_c, T, b};
  F.params = &params;

  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  double z = gal.Get_z();

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units



}



double Total_Disruption_Rate_Observed_Rband(Galaxy gal, double T, double b)
{
 
  // Keeping m_g here for testing purposes, really should be m_r
  double m_limit_contrast = find_host_contrast_magnitude(gal);
  
  double L_c = LCriticalRband(gal.Get_z(),T,m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int grid_size = 128;
  double relative_error = 1.e-5;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = Stellar_Rate_Integrand_RbandCut;

  struct stellar_rate_params {Galaxy gal; double this_L_c; double temperature; double beta;};
  stellar_rate_params params = {gal, L_c, T, b};
  F.params = &params;

  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  double z = gal.Get_z();

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units



}



