#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"


//default. Useful if you only want to use certian member functions
Disruption::Disruption()
{

  printf("Warning: using default disruption constructor. Members are not initialized\n");
  
}

Disruption::Disruption(Galaxy* gal)
{

  lf_log_powerlaw = 1.5;
  min_log_lbol = 43.0;
  max_radeff_streams = 0.1;
  max_edd_ratio = 1.;

  host_gal = gal;
  mbh = host_gal->Get_Mbh();
  L_Edd = Eddington_Luminosity();

  
  // move these to galaxy?
  T_opt_min = 1.e4;
  T_opt_max = 5.e4;
  T_opt_log_powerlaw = -1.;
  beta_mean = 1.;
  beta_sigma = 1.e-6;

  T_opt = 0.;
  beta = 0.;
  A_V = 0.;
  R_V = 0.;
  mstar = 0.;
  peak_L = 0.;

  max_L = 0.;
  
}

void Disruption::Set_Topt(double T)
{
  T_opt = T;
}

void Disruption::Set_Beta(double b)
{
  beta = b;
}

void Disruption::Set_Peak_L(double L)
{
  peak_L = L;
}

void Disruption::Set_A_V(double av)
{
  A_V = av;
}

void Disruption::Set_R_V(double rv)
{
  R_V = rv;
}

// Rejection sampling to avoid messiness of integrating piecewise function
// could this be made faster by using a better "proposal distribution?".
void Disruption::Rejection_Sample_Mstar(gsl_rng *rangen)
{

  double mstar_min = host_gal->Get_Mstar_Min(); //really will want to allow this to be different in different galaxies
  double mstar_max = host_gal->Get_Mstar_Max();
  
  double imf_norm = host_gal->Get_imf_norm();

  // careful about this if you end up using different IMF that is not monotonically decreasing with stellar mass
  double envelope = host_gal->Kroupa_IMF_for_value(mstar_min, imf_norm);

  while (true)
  {
    // pick a random mstar between minimum and maximum allowed values, in units of solar mass
    // this is the "proposal distribution"
    double mstar_e = mstar_min + (mstar_max - mstar_min)*gsl_rng_uniform(rangen);

    // calculate the probability from the present day mass function
    double P = 1./envelope * host_gal->Kroupa_IMF_for_value(mstar_e, imf_norm);
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

double Disruption::Get_Peak_L()
{
  return peak_L;
}

double Disruption::Get_Min_Log_Lbol()
{
  return min_log_lbol;
}

double Disruption::Get_A_V()
{
  return A_V;
}

double Disruption::Get_R_V()
{
  return R_V;
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

double Disruption::Get_Eddington_Luminosity()
{
  
  return Eddington_Luminosity();
  
}

// mstar and mbh in solar masses
// mdot in cgs
// Need to think harder about how to handle polytropic index dependence on stellar mass
double Disruption::Peak_Mdot()
{

  if (mstar <= 1.)
    {
      // gamma = 5./3.
      double guillochon_A = exp( (10.253 - 17.380 * beta + 5.9988 * beta * beta)/ (1. - 0.46573 * beta - 4.5066 * beta * beta));

      double rstar = Main_Sequence_Radius();

      //printf("beta is %f, rstar in solar units is %f\n",beta,rstar/R_SUN);

      // mstar already in solar units, but rstar in cgs
      return guillochon_A * pow(mbh/1.e6,-0.5) * mstar * mstar* pow(rstar/R_SUN,-1.5) * M_SUN / YEAR_TO_SEC;
    }

  else
    {
      
            // gamma = 4./3.
      double guillochon_A = exp( (27.261 - 27.516* beta + 3.8716 * beta * beta )/ (1. - 3.2605 * beta - 1.3865 * beta * beta));

      double rstar = Main_Sequence_Radius();

      //printf("beta is %f, rstar in solar units is %f\n",beta,rstar/R_SUN);

      // mstar already in solar units, but rstar in cgs
      return guillochon_A * pow(mbh/1.e6,-0.5) * mstar * mstar * pow(rstar/R_SUN,-1.5) * M_SUN / YEAR_TO_SEC;

    }
  
}

// at some point might want to think about other ways beta might enter here, aside from setting maximum mass fallback rate
void Disruption::Determine_Max_L()
{

  double mdot_fallback_peak = Peak_Mdot();
    
  double L_fallback_peak = max_radeff_streams * mdot_fallback_peak * C_LIGHT * C_LIGHT;
    
  if (L_fallback_peak >  max_edd_ratio  * L_Edd)
    max_L =  max_edd_ratio  * L_Edd;
    
  else
    max_L =  L_fallback_peak;

  return;
  
}



//Inverse transform sampling (not rejection sampling here)
void Disruption::Sample_Peak_L(gsl_rng *rangen)
{

  double u, y0, this_y;

  if (lf_log_powerlaw == 0)
    {

      y0 = max_L / pow(10., min_log_lbol);

      u = gsl_rng_uniform(rangen);

      this_y = exp(u * log(y0)); // log is natural log                                                                                                                      

      peak_L =  max_L / this_y;

    }

  else
    {
  
      y0 = pow( pow(10., min_log_lbol)/max_L, lf_log_powerlaw);

      u = gsl_rng_uniform(rangen);

      this_y = y0/(1. - u * (1. - y0));

      peak_L = max_L * pow(this_y,1./lf_log_powerlaw);

    }

  //  return;
}

void Disruption::Sample_Beta(gsl_rng *rangen)
{
  beta = beta_mean + gsl_ran_gaussian(rangen, beta_sigma);

}

void Disruption::Sample_Topt(gsl_rng *rangen)
{

  double y0, u, this_y;
  
  if (T_opt_log_powerlaw == 0)
    {

      y0 = T_opt_max / T_opt_min;

      u = gsl_rng_uniform(rangen);

      this_y = exp(u * log(y0)); // log is natural log

      T_opt = T_opt_max / this_y;

    }

  
  else
    {
      y0 = pow(T_opt_min/T_opt_max,T_opt_log_powerlaw);

      u = gsl_rng_uniform(rangen);

      this_y = y0/(1. - u * (1. - y0));

      T_opt =  T_opt_max * pow(this_y,1./T_opt_log_powerlaw);
    }
}

void Disruption::Sample_A_V(gsl_rng *rangen)
{
  
  Sample_R_V(rangen);

  if (log10(host_gal->Get_ssfr()) < -11.3)
    {

      A_V = std::max(host_gal->Get_median_A_V() + gsl_ran_gaussian(rangen, host_gal->Get_sigma_A_V()), 0.);
    }
    else
      {
	double this_Aha = std::max(host_gal->Get_median_A_Ha() + gsl_ran_gaussian(rangen, host_gal->Get_sigma_A_Ha()),0.);
	A_V = this_Aha / host_gal->Calzetti_Extinction(1./(0.65645),R_V);
      }

}

void Disruption::Sample_R_V(gsl_rng *rangen)
{
  
  //  R_V = R_V_mean + gsl_ran_gaussian(rangen, R_V_sigma);
  R_V = host_gal->Get_median_R_V();
  
}

// trying to keep this general, not only for optical emission
// "nu_emit" is to remind you it's in galaxy rest frame, not observer frame
double Disruption::Planck_Function_Frequency(double nu_emit, double T)
{

  return 2. * H_PLANCK * pow(nu_emit,3.)/( pow(C_LIGHT,2.) * (exp(H_PLANCK * nu_emit /(K_BOLTZ * T)) - 1.) );

}

// trying to keep this general, not only for optical emission
// "nu_emit" is to remind you it's in galaxy rest frame, not observer frame
double Disruption::Unobscured_Lnu(double nu_emit, double T, double Lbol)
{

  return PI * Planck_Function_Frequency(nu_emit, T) * Lbol/(STEF_BOLTZ * pow(T,4.));
  
}

double Disruption::Dust_Flux_Factor_Reduction(double nu_emit)
{

  double lambda_cm = C_LIGHT/nu_emit; // assumes nu_emit in Hz (rest frame)
  double lambda_micron = lambda_cm * 1.e4;
  double x = 1./lambda_micron;
  
  //  double cardelli_factor = Cardelli_Extinction(x);
  //  return pow(10., -1. * A_V * cardelli_factor/(2.5));

  double calzetti_factor = host_gal->Calzetti_Extinction(x,R_V);
  return pow(10., -1. * A_V * calzetti_factor/(2.5));

}

// can save time by precomputing the (1. + z)/(4. * PI * d_L * d_L)  for each galaxy and passing it heree
// can't precompute the dust flux factor reduction because A_V might be randomly generated
double Disruption::Extincted_Flux_Observed(double nu_emit, double cosmo_factor)
{

  return Unobscured_Lnu(nu_emit, T_opt,peak_L) * cosmo_factor * Dust_Flux_Factor_Reduction(nu_emit);

}





