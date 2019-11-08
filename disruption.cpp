#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"


Disruption::Disruption(Galaxy gal)
{

  lf_log_powerlaw = 1.5;
  min_log_lbol = 43.0;
  max_radeff_streams = 0.1;
  max_edd_ratio = 2.;

  host_gal = gal;
  mbh = host_gal.Get_Mbh();
  L_Edd = Eddington_Luminosity();



  // move these to galaxy?
  T_opt_mean = 3.e4;
  T_opt_sigma = 1.5e4;
  beta_mean = 1.;
  beta_sigma = 1.e-6;
  A_V_mean = 1.e-6; // should depend on galaxy properties
  A_V_sigma = 0.; // should depend on galaxy properties
  R_V_mean = 3.; // should depend on galaxy properties
  R_V_sigma = 1.e-6; // should depend on galaxy properties

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

  double mstar_min = host_gal.Get_Mstar_Min(); //really will want to allow this to be different in different galaxies
  double mstar_max = host_gal.Get_Mstar_Max();
  
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

double Disruption::Get_Peak_L()
{
  return peak_L;
}

double Disruption::Get_Min_Log_Lbol()
{
  return min_log_lbol;
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
    
  double L_fallback_peak = max_radeff_streams * mdot_fallback_peak * C_LIGHT * C_LIGHT;
    
  if (L_fallback_peak >  max_edd_ratio  * L_Edd)
    max_L =  max_edd_ratio  * L_Edd;
    
  else
    max_L =  L_fallback_peak;

  return;
  
}



//Not rejection sampling, inverse transform sampling
void Disruption::Sample_Peak_L(gsl_rng *rangen)
{

  double ymin = pow( pow(10., min_log_lbol)/max_L, lf_log_powerlaw);

  double u = gsl_rng_uniform(rangen);

  double this_y = ymin/(1. - u * (1. - ymin));

  peak_L = max_L * pow(this_y,1./lf_log_powerlaw);

  return;
}

void Disruption::Sample_Beta(gsl_rng *rangen)
{
  beta = beta_mean + gsl_ran_gaussian(rangen, beta_sigma);

}

void Disruption::Sample_Topt(gsl_rng *rangen)
{
  double u = gsl_rng_uniform(rangen);  
  T_opt = T_opt_mean + T_opt_sigma* (2. * u - 1.);
}

void Disruption::Sample_A_V(gsl_rng *rangen)
{

  A_V = A_V_mean + gsl_ran_gaussian(rangen, A_V_sigma);

}

void Disruption::Sample_R_V(gsl_rng *rangen)
{
  
  R_V = R_V_mean + gsl_ran_gaussian(rangen, R_V_sigma);
  
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
  
  double cardelli_factor = Cardelli_Extinction(x);

  return pow(10., -1. * A_V * cardelli_factor/(2.5));

}

// can save time by precomputing the (1. + z)/(4. * PI * d_L * d_L)  for each galaxy and passing it heree
// can't precompute the dust flux factor reduction because A_V might be randomly generated
double Disruption::Extincted_Flux_Observed(double nu_emit, double cosmo_factor)
{

  return Unobscured_Lnu(nu_emit, T_opt,peak_L) * cosmo_factor * Dust_Flux_Factor_Reduction(nu_emit);

}

double Disruption::Cardelli_Extinction(double x)
{

  double a = 0.;
  double b = 0.;

  if (x < 0.3 || x > 10.)
    {
      printf("ERROR: Cardelli extinction not defined for x < 0.3 or x > 10.\n");
      return 0.;
    }
  if (x >= 0.3 && x < 1.1)
    {
      a = 0.574 * pow(x,1.61);
      b = -0.527 * pow(x,1.61);
    }
  if (x >= 1.1 and x < 3.3)
    {
      double y = x - 1.82;
      a = 1. + 0.17699 * y - 0.50447 * pow(y,2.) -0.02427 * pow(y,3.) + 0.72085 * pow(y,4.) + 0.01979 * pow(y,5.) - 0.77530 * pow(y,6.) + 0.32999 * pow(y,7.);
      b = 1.41338 * y + 2.28305 * pow(y,2.) + 1.07233 * pow(y,3.) - 5.38434 * pow(y,4.) - 0.62251 * pow(y,5.) + 5.3026 * pow(y,6) - 2.09002 * pow(y,7.);
    }
    if (x >= 3.3 and x < 8.)
      {
	double Fa = 0.;
	double Fb = 0.;
	if (x >= 5.9)
	  {
            Fa = -0.04473 * pow(x - 5.9,2.) - 0.009779 *pow(x - 5.9,3.);
            Fb = 0.213 * pow(x - 5.9,2.) + 0.1207 *pow(x - 5.9,3.);
	  }
        a = 1.752 - 0.316 * x - 0.104/(pow(x - 4.67,2.) + 0.341) + Fa;
        b = -3.09 + 1.825 * x + 1.206/(pow(x - 4.62,2.) + 0.263) + Fb;
      }
    if (x >= 8. && x <= 10.)
      {
        a = -1.073 - 0.628 * (x - 8.) + 0.137 * pow(x - 8.,2.) - 0.070 * pow(x - 8.,3.);
        b = 13.67 + 4.257 * (x - 8.) - 0.42 * pow(x - 8.,2.) + 0.374 * pow(x - 8.,3.);
      }
	  /*
    #if (x > 10. and x < 10.96):
    #    #print "ERROR, Cardelli extinction not defined"
    #    a = -1.073 - 0.628 * (10. - 8.) + 0.137 * pow(10. - 8.,2.) - 0.070 * pow(10. - 8.,3.)
    #    b = 13.67 + 4.257 * (10. - 8.) - 0.42 * pow(10. - 8.,2.) + 0.374 * pow(10. - 8.,3.)
    #    ref_value = a + b / R
    #    return ref_value + (x - 10.)
    #if (x >= 10.96):
    #    return 1000.
	  */
        
    return a + b / R_V;

}



