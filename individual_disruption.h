#ifndef INDIVIDUAL_DISRUPTION_H
#define INDIVIDUAL_DISRUPTION_H

#include <math.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "magnitudes.h"

#define MSTAR_MAX 1.0 // maximum stellar mass, in solar units. Really, should consider making this a function of galaxy properties
#define MSTAR_MIN 0.08 // Only considering main sequence stars

// In the future these will depend on galaxy properties
#define RATE_NORMALIZATION 2.9e-5 
#define RATE_POWERLAW -0.404

#define LF_LOG_POWERLAW 1.5

#define MIN_LOG_LBOL 43.0
#define MAX_RADEFF_STREAMS 0.1 // maximum radiative efficiency of stream collisions.

#define MAX_EDD_RATIO 2.

// mstar in solar mass throughout

// should maybe make this a class private variable
// maybe make all of this part of a galaxy class ... that way things like mstar_max, rate noramlization, etc, are actually tied to the galaxy's properties
double IMF_normalization = 1.; // just the initial value, will be modified

// See equation 7 of https://arxiv.org/pdf/1212.0939.pdf
double Kroupa_IMF(double mstar, void * params)
{
  double m1 = 0.08;
  double m2 = 0.5;
  double m3 = MSTAR_MAX;

  double k1 = pow(m1,-0.3 + 1.3);
  double k2 = k1 * pow(m2,-1.3 + 2.3);
  //  double k3 = k2 * pow(m3, -2.3 + 2.3);

  if (mstar >= m1 && mstar < m2)
    return 1./IMF_normalization * k1 * pow(mstar,-1.3);
  else if (mstar >= m2 && mstar <= m3)
      return 1./IMF_normalization * k2 * pow(mstar,-2.3);
  else
    return 0.;
  
}


// This could be improved
double main_sequence_radius(double mstar)
{

  if (mstar <= 1.)
    return R_SUN * pow(mstar,0.8);
  else 
    return R_SUN * pow(mstar,0.57);
  
}

// See first paragraph section 2.2 of Stone & Metzger 2016.
// The Newtonian formula in Stone & Metzger requires relativistic correction, as in Beloborodov 1992 (involves IBSCO)
// Final result is consistent with what is quoted in Leloudas et al 2016
// mstar in solar mass, returns a mass in solar masses
double Hills_Mass(double mstar)
{
  double rstar = main_sequence_radius(mstar);
    
  double relativistic_correction = sqrt(5.) / 8. / pow(2.,-1.5); // beloborodov 1992
  return (relativistic_correction * pow(rstar * C_LIGHT * C_LIGHT /(2. * NEWTON_G * pow(mstar * M_SUN,1./3.)) , 3./2.)/M_SUN);

}

// both mbh and mstar in solar masses
double Stellar_Disruption_Rate(double mbh, double mstar)
{
  //consier using the exp(-M^2) suppression here as suggested by Kesden 2012, or a sharp and total cutoff
  double mhills = Hills_Mass(mstar);
  if (mbh > mhills) return 0;

  return RATE_NORMALIZATION * pow(mbh/1.e8,RATE_POWERLAW); // rate is still in galaxy proper time
  
}


double Stellar_Rate_Integrand(double mstar, void * params)
{

  double mbh = *(double *) params;

  return Stellar_Disruption_Rate(mbh, mstar) * Kroupa_IMF(mstar, params); // Stellar disruption rate in terms of galaxy proper time
  
}



void Determine_IMF_Normalization()
{


  int grid_size = 128;
  double relative_error = 1.e-6;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = &Kroupa_IMF;

  double result, error;
  gsl_integration_qags(&F,MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);
  gsl_integration_workspace_free(workspace);
  
  IMF_normalization = result;

}


void Initialize_IMF()
{ 
 //   Initialize_Stellar_Masses(); // have to be careful about where this is included
   Determine_IMF_Normalization(); // have to be careful about where this is included

}

double Total_Disruption_Rate(double mbh, double z)
{
 

  //  printf("norm is %e\n",IMF_normalization);

  int grid_size = 128;
  double relative_error = 1.e-6;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = &Stellar_Rate_Integrand;

  F.params = &mbh;

  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units


}

// assumes mbh in solar masses
double Eddington_Luminosity(double mbh)
{
  
  return 4. * PI * NEWTON_G * mbh * M_SUN * M_PROTON * C_LIGHT/(THOMSON_CS);
  
}

// mstar and mbh in solar masses
// mdot in cgs
// Need to think harder about how to handle polytropic index dependence on stellar mass
double Peak_Mdot(double mbh, double mstar, double beta)
{

  if (mstar <= 1.)
    {
      // gamma = 5./3.
      double guillochon_A = exp( (10.253 - 17.380 * beta + 5.9988 * pow(beta,2. ) )/ (1. - 0.46573 * beta - 4.5066 * pow(beta,2. )));

      double rstar = main_sequence_radius(mstar);

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
double Max_Luminosity(double mbh, double mstar, double beta)
{

  double Ledd = Eddington_Luminosity(mbh);
    
  double mdot_fallback_peak = Peak_Mdot(mbh,mstar,beta);
    
  double L_fallback_peak = MAX_RADEFF_STREAMS * mdot_fallback_peak * C_LIGHT * C_LIGHT;

    
  if (L_fallback_peak >  MAX_EDD_RATIO  * Ledd)
    return MAX_EDD_RATIO  * Ledd;
    
  else
    return L_fallback_peak;
    
  
}

// might end up changing the scope of some of these variables, like mbh, if all this gets moved to a galaxy object
double Fraction_Observed(double mbh, double mstar, double beta, double L_c)
{

  double L_max = Max_Luminosity(mbh, mstar, beta);

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

  struct stellar_rate_params {double mbh; double beta; double L_c;};
  struct stellar_rate_params * params = (struct stellar_rate_params *)p;

  double mbh = params->mbh;
  double beta = params->beta;
  double L_c = params->L_c;

  // I don't think Kroupa_IMF currently cares what param pointer is passed to it
  return Stellar_Disruption_Rate(mbh, mstar) * Kroupa_IMF(mstar, p) * Fraction_Observed(mbh,mstar, beta,L_c); // rate of observable disruptions in terms of galaxy proper time
  
}

// might end up changing the scope of some of these variables, like mbh and z, if all this gets moved to a galaxy object
double Stellar_Rate_Integrand_RbandCut(double mstar, void * p)
{

  struct stellar_rate_params {double mbh; double beta; double L_c;};
  struct stellar_rate_params * params = (struct stellar_rate_params *)p;

  double mbh = params->mbh;
  double beta = params->beta;
  double L_c = params->L_c;

  // I don't think Kroupa_IMF currently cares what param pointer is passed to it
  return Stellar_Disruption_Rate(mbh, mstar) * Kroupa_IMF(mstar, p) * Fraction_Observed(mbh,mstar, beta,L_c); // rate of observable disruptions in terms of galaxy proper time
  
}



double Total_Disruption_Rate_Observed_Gband(double mbh, double beta, double z, double T, double m_limit_contrast )
{
 

  //  printf("IMF norm is %e\n",IMF_normalization);

  double L_c = LCriticalGband(z,T,m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int grid_size = 128;
  double relative_error = 1.e-5;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = &Stellar_Rate_Integrand_GbandCut;

  struct stellar_rate_params {double this_mbh; double this_beta; double this_L_c;};

  stellar_rate_params params = {mbh, beta, L_c};
  F.params = &params;

  //  printf("params are %e %f %f %f\n",mbh,beta,z,T);


  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units



}

double Total_Disruption_Rate_Observed_Rband(double mbh, double beta, double z, double T, double m_limit_contrast )
{
 

  //  printf("IMF norm is %e\n",IMF_normalization);

  double L_c = LCriticalRband(z,T,m_limit_contrast);
  //  printf("L_c is %e\n", L_c);

  int grid_size = 128;
  double relative_error = 1.e-5;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = &Stellar_Rate_Integrand_RbandCut;

  struct stellar_rate_params {double this_mbh; double this_beta; double this_L_c;};

  stellar_rate_params params = {mbh, beta, L_c};
  F.params = &params;

  //  printf("params are %e %f %f %f\n",mbh,beta,z,T);


  double result, error; 
  
  gsl_integration_qags(&F, MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return 1./(1. + z) * result; // the 1/(1+z) factor converts from galaxy proper time units to observer time units



}





#endif
