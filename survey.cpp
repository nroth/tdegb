#include <math.h>
#include <algorithm>
#include "physical_constants.h"
#include "galaxy.h"
#include "survey.h"


Survey::Survey()
{

  nu_gband = 646243712006898.;
  nu_rband =  489696925841228.;

  m_r_threshhold = 20;
  m_g_threshhold = 20;

  host_contrast_cut = 0.5;
  
  r_psf_arcsec = 4.4; // 2 times r-band median PSF FWHM of 2.2 arcsec. will need to distinguish how this is handled in various bands

}

double Survey::Get_Nu_Gband()
{
  return nu_gband;
}

double Survey::Get_Nu_Rband()
{
  return nu_rband;
}

double Survey::Get_m_r_Threshhold()
{
  return m_r_threshhold;
}

double Survey::Get_m_g_Threshhold()
{
  return m_g_threshhold;
}

// Convert cgs F_nu to AB magnitude
// Assumes Fnu in cgs
double Survey::mAB_From_Fnu(double F_nu)
{
  if (F_nu <= 0) printf("ERROR: Fnu is negative\n");

  return -2.5 * log10(F_nu) - 48.6;
}

// find a way to make this passed galaxy constant?
double Survey::Find_Host_Contrast_Magnitude(Galaxy gal)
{

  // REALLY NEED A GOOD WAY TO CHOOSE THE BAND HERE
  double m_tot = gal.Get_m_g();
  double mu_e = gal.Get_Mu_Eff(m_tot);
  double I_e = I_From_Mu(mu_e);

  //magnitude (flux) enclosed in psf
  double m_psf = Mu_From_I(gal.Flux_Enclosed_R_Sersic(r_psf_arcsec,I_e));
    
  return m_psf - 2.5 * log10(pow(10.,host_contrast_cut / 2.5) - 1.);

}


//might want this as part of physical constants or standard definitions, not survey class
double Survey::I_From_Mu(double mu)
{    
  return pow(10., mu/(-2.5));
}

//might want this as part of physical constants or standard definitions, not survey class
double Survey::Mu_From_I(double I)
{    
  return -2.5 * log10(I);
}
