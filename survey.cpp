#include <math.h>
#include <stdio.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "survey.h"


Survey::Survey()
{

  // see e.g. http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Swift/UVOT.UVW2&&mode=browse&gname=Swift&gname2=UVOT#filter  
  nu_bands.push_back(C_LIGHT / (6436.92 * 1.e-8)); // ZTF_r
  nu_bands.push_back(C_LIGHT / (4804.79 * 1.e-8)); // ZTF_g
  nu_bands.push_back(C_LIGHT / (2688.46 * 1.e-8)); // Swift UVW1
  nu_bands.push_back(C_LIGHT / (2272.71 * 1.e-8)); // Swift UVM2
  nu_bands.push_back(C_LIGHT / (2140.26 * 1.e-8)); // Swift UVW2

  m_r_threshhold = 19;
  m_g_threshhold = 19;

  host_contrast_cut = 0.36385;

  g_psf_arcsec = 2.1; // g-band median PSF FWHM. Page 8 of Bellm et al 2019 https://iopscience.iop.org/article/10.1088/1538-3873/aaecbe/pdf
  r_psf_arcsec = 2.0; // r-band median PSF FWHM. Page 8 of Bellm et al 2019 https://iopscience.iop.org/article/10.1088/1538-3873/aaecbe/pdf
 

}

double Survey::Get_Band_Nu(int index)
{
  return nu_bands[index];
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
double Survey::Find_Host_Contrast_Magnitude(Galaxy gal, char band)
{

  double m_tot = gal.Get_m_r(); // default is r band
  if (band == 'g')
    {
      m_tot = gal.Get_m_g();
    }
  double mu_e = gal.Get_Mu_Eff(m_tot);
  double I_e = I_From_Mu(mu_e);

  double m_psf = 0.;
  //magnitude (flux) enclosed in psf
  if (band == 'g')
    {
      m_psf = Mu_From_I(gal.Flux_Enclosed_R_Sersic(g_psf_arcsec,I_e));
    }
  if (band == 'r')
    {
      m_psf = Mu_From_I(gal.Flux_Enclosed_R_Sersic(r_psf_arcsec,I_e));
    }
    
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
