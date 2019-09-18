#ifndef MAGNITUDES_H
#define MAGNITUDES_H

#include <math.h>
#include "physical_constants.h"
#include "cosmology.h"


#define M_APPARENT_THRESHHOLD 22.0 // will need to distinguish how this is handled in various bands

#define Z_MIN_CUT 0.003
#define Z_MAX_CUT 10.
#define Z_MIN_GENERAL 0.0001 // Just for setting numeric limits, but might need to think about this more 
#define Z_MAX_GENERAL 10. // Just for setting numeric limits, but might need to think about this more 


#define NU_GBAND 646243712006898.
#define NU_RBAND 489696925841228.

double PlanckFunctionFrequency(double nu, double T)
{

  return 2. * H_PLANCK * pow(nu,3.)/( pow(C_LIGHT,2.) * (exp(H_PLANCK * nu /(K_BOLTZ * T)) - 1.) );

}

//includes K-correction
double UnexctinctedBBGbandFlux(double z, double T, double L)
{
  double nu_corrected_g = NU_GBAND * (1. + z);

  double d_L = LuminosityDistance(z);

  double unobscured_Lnu =  PI * PlanckFunctionFrequency(nu_corrected_g, T) * L/(STEF_BOLTZ * pow(T,4.));

  //Now apply the "band stretch factor" 1 + z and convert to flux via luminosity distance

  return unobscured_Lnu * (1. + z)/(4. * PI * d_L * d_L);
  
}

// Convert cgs F_nu to AB magnitude
// Assumes Fnu in cgs
double mABFromFnu(double F_nu)
{
  if (F_nu <= 0) printf("ERROR: Fnu is negataive");

  return -2.5 * log10(F_nu) - 48.6;
}

// Given z and T, what is minimum Lbol that allows the k-corrected flux to exceed the vaalue required to matach the survey apparent magnitude threshhold in the g band?
// Here it is easy enough to solve for L. More generally might need to solve an equation with Brent method the way you do for finding Zmax
double LCriticalGband(double z, double T)
{

  return pow(10., (M_APPARENT_THRESHHOLD + 48.6)/-2.5) * STEF_BOLTZ * pow(T,4.) * 4. * PI * pow(LuminosityDistance(z),2.)/( (1. + z) * PI * PlanckFunctionFrequency(NU_GBAND * (1. + z), T));

}
  

// will need to specify which band(s), and distinguish how this is handled in various bands
double KCorrection()
{
  return 0.;
}

// will need to distinguish how this is handled in various bands
double AbsMagFromApparentMag(double apparent_mag, double z)
{

  return apparent_mag - DistanceModulus(z) - KCorrection();

}


double ZmaxEqCondition(double this_z, double apparent_mag, double z_original)
{

  double abs_mag = AbsMagFromApparentMag(apparent_mag,z_original);

  return LuminosityDistance(this_z)/(10. * PARSEC) - pow(10., 1./5. * (M_APPARENT_THRESHHOLD - abs_mag - KCorrection()));
  
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
double ZmaxBrentMethod(double apparent_mag, double z_original)
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

double FindZmax(double apparent_mag, double z_original)
{
  return ZmaxBrentMethod(apparent_mag,z_original);
}

#endif
