#ifndef _MAGNITUDES_H
#define _MAGNITUDES_H

#include <math.h>
#include <algorithm>
#include <gsl/gsl_sf_gamma.h>
#include "physical_constants.h"
#include "cosmology.h"
#include "galaxy.h"
#include "disruption.h"
#include "dust.h"


#define M_APPARENT_THRESHHOLD 20.0 // will need to distinguish how this is handled in various bands
#define R_PSF_ARCSEC 4.4 // 2 times r-band median PSF FWHM of 2.2 arcsec. will need to distinguish how this is handled in various bands
#define HOST_CONTRAST_CUT 0.5 // will need to distinguish how this is handled in various bands
#define RESOLUTION_FOR_NUKERGAMMA 0.04 // arsec. See Lauer et al 2007. Nick Stone's rate calculations were based on Nuker gamma as measured in this paper, so to convert n_sersic to nuker gamma we want to account for how they measured gamma

#define NU_GBAND 646243712006898.
#define NU_RBAND 489696925841228.

#define Z_MIN_CUT 0.003
#define Z_MAX_CUT 10.
#define Z_MIN_GENERAL 0.0001 // Just for setting numeric limits, but might need to think about this more 
#define Z_MAX_GENERAL 10. // Just for setting numeric limits, but might need to think about this more 



static double PlanckFunctionFrequency(double nu, double T)
{

  return 2. * H_PLANCK * pow(nu,3.)/( pow(C_LIGHT,2.) * (exp(H_PLANCK * nu /(K_BOLTZ * T)) - 1.) );

}

//includes K-correction
/* This function not needed except fo debugging?
double UnexctinctedBBGbandFlux(double z, double T, double L)
{
  double nu_corrected_g = NU_GBAND * (1. + z);

  double d_L = LuminosityDistance(z);

  double unobscured_Lnu =  PI * PlanckFunctionFrequency(nu_corrected_g, T) * L/(STEF_BOLTZ * pow(T,4.));

  //Now apply the "band stretch factor" 1 + z and convert to flux via luminosity distance

  return unobscured_Lnu * (1. + z)/(4. * PI * d_L * d_L);
  
}
*/

// Convert cgs F_nu to AB magnitude
// Assumes Fnu in cgs
static double mABFromFnu(double F_nu)
{
  if (F_nu <= 0) printf("ERROR: Fnu is negataive");

  return -2.5 * log10(F_nu) - 48.6;
}

// want to make this obsolete
static double Flare_m_r(double Lbol, double T, double z, double d_L)
{

  double nu_corrected_r = NU_RBAND * (1. + z);

  double unobscured_Lnu =  PI * PlanckFunctionFrequency(nu_corrected_r, T) * Lbol/(STEF_BOLTZ * pow(T,4.));

  //Now apply the "band stretch factor" 1 + z and convert to flux via luminosity distance

  return mABFromFnu(unobscured_Lnu * (1. + z)/(4. * PI * d_L * d_L));

}


// want to make this obsolete
static double Flare_m_g(double Lbol, double T, double z, double d_L)
{
  double nu_corrected_g = NU_GBAND * (1. + z);

  double unobscured_Lnu =  PI * PlanckFunctionFrequency(nu_corrected_g, T) * Lbol/(STEF_BOLTZ * pow(T,4.));

  //Now apply the "band stretch factor" 1 + z and convert to flux via luminosity distance

  return mABFromFnu(unobscured_Lnu * (1. + z)/(4. * PI * d_L * d_L));
}

static double Unobscured_Lnu(double nu_emit, double T, double Lbol)
{

  return PI * PlanckFunctionFrequency(nu_emit, T) * Lbol/(STEF_BOLTZ * pow(T,4.));
  
}

static double DustFluxFactorReduction(double nu_emit, double A_V, double R_V)
{

  double lambda_cm = C_LIGHT/nu_emit; // assumes nu_emit in Hz (rest frame)
  double lambda_micron = lambda_cm * 1.e4;
  double x = 1./lambda_micron;
  
  double cardelli_factor = Cardelli_Extinction(x,R_V);

  return pow(10., -1. * A_V * cardelli_factor/(2.5));

}

// can save time by precomputing the (1. + z)/(4. * PI * d_L * d_L)  for each galaxy and passing it heree
// can't precompute the dust flux factor reduction because A_V might be randomly generated
static double ExtinctedFluxObserved(double nu_emit, double cosmo_factor, double Lbol, double T, double A_V, double R_V)
{

  return Unobscured_Lnu(nu_emit, T,Lbol) * cosmo_factor * DustFluxFactorReduction(nu_emit,A_V,R_V);

}


// see reference mentioned at https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
static double get_approx_sersic_bn(double sersic_n)
{
  return 1.9992 * sersic_n - 0.3271;
}


// as found in https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
// will be in units of magnitudes / arcsec^2 assuming m_tot in magnitude and r_eff_arcsec in arsec
static double get_mu_eff(double sersic_n, double m_tot,double r_eff_arcsec)
{
  double bn = get_approx_sersic_bn(sersic_n);
    
  return m_tot + 5. * log10(r_eff_arcsec) + 2.5 * log10(2. * PI * sersic_n * exp(bn)/pow(bn,2. * sersic_n) * gsl_sf_gamma(2. * sersic_n)  );
}


// as found in https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
// I_e in mags / arcsec^2, r and r_e in arcsec
static double flux_enclosed_r_sersic(double r, double sersic_n,double r_e, double I_e)
{
  
  double bn = get_approx_sersic_bn(sersic_n);
      
  return I_e * r_e * r_e * 2. * PI * sersic_n * exp(bn)/pow(bn,2. * sersic_n) * gsl_sf_gamma_inc_P(2. *sersic_n, bn * (pow(r/r_e,1./sersic_n)))* gsl_sf_gamma(2. * sersic_n);

}

static double I_from_mu(double mu)
{    
  return pow(10., mu/(-2.5));
}

static double mu_from_I(double I)
{    
  return -2.5 * log10(I);
}

// maybe consider moving this to physical constants
static double arcsec_from_radians(double radians)
{    
  return radians * 180. * 3600./ PI;
}

// maybe consider moving this to physical constants
static double radian_from_arcsec(double r_arcsec)
{    
  return r_arcsec * PI / (180. * 3600.);
}

static double r_arcsec_from_kpc(double r_kpc, double z)
{
    
  double r_cm = r_kpc * 1000. * PARSEC; // convert to cm
  double d_A = AngularDiameterDistance(z); // in cm

  double r_radians = r_cm/d_A;
    
  return arcsec_from_radians(r_radians);
}

static double r_kpc_from_arcsec(double r_arcsec,double z)
{
  double r_rad = radian_from_arcsec(r_arcsec);
    
  double d_A = AngularDiameterDistance(z);
    
  double r_cm = r_rad * d_A;
    
  return r_cm / (1000. * PARSEC);
}

//double find_host_contrast_magnitude(double m_tot,double sersic_n,double re_kpc, double z)

static double find_host_contrast_magnitude(Galaxy gal)
{

  // REALLY NEED A GOOD WAY TO CHOOSE THE BAND HERE
  double m_tot = gal.Get_m_g();
  double sersic_n = gal.Get_sersic_n();
  double re_kpc = gal.Get_r50_kpc();
  double z = gal.Get_z();
  
  double re_arcsec = r_arcsec_from_kpc(re_kpc,z);
    
  double mu_e = get_mu_eff(sersic_n,m_tot,re_arcsec);
  double I_e = I_from_mu(mu_e);

  //magnitude (flux) enclosed in psf
  double m_psf = mu_from_I(flux_enclosed_r_sersic(R_PSF_ARCSEC,sersic_n,re_arcsec,I_e));
    
  return m_psf - 2.5 * log10(pow(10.,HOST_CONTRAST_CUT / 2.5) - 1.);
}


// Really need to find a way to put this in galaxy.h
static double Find_Nuker_Gammaprime_From_Sersic(double sersic_n, double re_kpc, double z)
{

  double bn = get_approx_sersic_bn(sersic_n);

  double rprime_kpc = r_kpc_from_arcsec(RESOLUTION_FOR_NUKERGAMMA,z);

  return bn/sersic_n * pow(rprime_kpc/re_kpc,1./sersic_n);

  // now use Graham 2003 eq 4

  
  
}

#endif
