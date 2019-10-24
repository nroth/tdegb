#ifndef _MAGNITUDES_H
#define _MAGNITUDES_H

#include <math.h>
#include <algorithm>
#include <gsl/gsl_sf_gamma.h>
#include "physical_constants.h"
#include "cosmology.h"
#include "galaxy.h"


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

// Given z and T, what is minimum Lbol that allows the k-corrected flux to exceed the value required to matach the survey apparent magnitude threshhold in the g band?
// Here it is easy enough to solve for L. More generally might need to solve an equation with Brent method the way you do for finding Zmax
static double LCriticalGband(double z, double T, double m_limit_contrast)
{

  double operating_limit = std::min(m_limit_contrast,M_APPARENT_THRESHHOLD);
  return pow(10., (operating_limit + 48.6)/-2.5) * STEF_BOLTZ * pow(T,4.) * 4. * PI * pow(LuminosityDistance(z),2.)/( (1. + z) * PI * PlanckFunctionFrequency(NU_GBAND * (1. + z), T));

}

static double LCriticalRband(double z, double T, double m_limit_contrast)
{

  double operating_limit = std::min(m_limit_contrast,M_APPARENT_THRESHHOLD);
  return pow(10., (operating_limit + 48.6)/-2.5) * STEF_BOLTZ * pow(T,4.) * 4. * PI * pow(LuminosityDistance(z),2.)/( (1. + z) * PI * PlanckFunctionFrequency(NU_RBAND * (1. + z), T));

}
  

// will need to specify which band(s), and distinguish how this is handled in various bands
static double KCorrection()
{
  return 0.;
}

// will need to distinguish how this is handled in various bands
static double AbsMagFromApparentMag(double apparent_mag, double z)
{

  return apparent_mag - DistanceModulus(z) - KCorrection();

}


static double ZmaxEqCondition(double this_z, double apparent_mag, double z_original)
{

  double abs_mag = AbsMagFromApparentMag(apparent_mag,z_original);

  return LuminosityDistance(this_z)/(10. * PARSEC) - pow(10., 1./5. * (M_APPARENT_THRESHHOLD - abs_mag - KCorrection()));
  
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double ZmaxBrentMethod(double apparent_mag, double z_original)
{

  double brent_solve_tolerance = 1.0e-2;
  double z_range_min = Z_MIN_GENERAL;
  double z_range_max = Z_MAX_GENERAL;

  int ITMAX = 100;
  double EPS = 3.0e-8;
  int iter;

  // Initial guesses
  double a=z_range_min;
  double b=z_range_max;
  double c=b;
  double d,e,min1,min2;
  double fa,fb = 0.;

  fa = ZmaxEqCondition(a, apparent_mag, z_original);
  fb = ZmaxEqCondition(b, apparent_mag, z_original);
  double fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    printf("Root must be bracketed in brent method");
  fc=fb;
  for (iter=1;iter<=ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*brent_solve_tolerance;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
         p=2.0*xm*s;
         q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);

    fb = ZmaxEqCondition(b, apparent_mag, z_original);
  }
  printf("Maximum number of iterations exceeded in zbrent");
  return 0.0;

}

static double FindZmax(double apparent_mag, double z_original)
{
  return ZmaxBrentMethod(apparent_mag,z_original);
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

//static double find_host_contrast_magnitude(double m_tot,double sersic_n,double re_kpc, double z)

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
