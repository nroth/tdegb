/*
Copyright (c) 2020 Nathaniel Jacob Roth, licensed under GPL-3.0-or-later

This file is part of tdegb

tdegb is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <algorithm> // for std::max
#include <stdio.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_gamma.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "cosmology.h"


//***************************************************************
// Constructors
//***************************************************************


Galaxy::Galaxy(double* galaxy_info)
{

  mstar_max = 1.0;  // In the future these will depend on galaxy properties
  mstar_min = 0.08; // Only considering main sequence stars
  disruption_rate_normalization_combined =  pow(10.,-3.79);// for nuker gamma = 1 and
  //  disruption_rate_normalization_combined =  6.5e-5;// for nuker gamma = 1 and
  disruption_rate_powerlaw_mass = -0.223;
  disruption_rate_powerlaw_nuker = 0.852; // for galaxies such that black hole mass is below Hills mass for 1 solar mass star

  resolution_for_nuker_gamma = 0.04; // arsec. See Lauer et al 2007. Nick Stone's rate calculations were based on Nuker gamma as measured in this paper, so to convert n_sersic to nuker gamma we want to account for how they measured gamma
  resolution_for_central_sb_measurement = 0.04; // arcsec

  
  total_stellar_mass = pow(10.,galaxy_info[0]); // you passed it as log  
  mbh_sigma = pow(10.,galaxy_info[1]); // converting log to value
  mbh_bulge = pow(10.,galaxy_info[2]); // converting log to value
  z = galaxy_info[3];
  sersic_n = galaxy_info[4];
  r50_kpc = galaxy_info[5];
  m_g = galaxy_info[6];
  m_r = galaxy_info[7];
  ssfr = pow(10.,galaxy_info[8]);
  UminusR = galaxy_info[9];
  M_r = galaxy_info[10];


  sersic_bn = Get_Approx_Sersic_bn();
  
  nuker_gammaprime = Find_Nuker_Gammaprime_From_Sersic();

  mbh = mbh_sigma; // making this choice for now
  
  re_arcsec = R_Arcsec_From_Kpc(r50_kpc);

  Set_IMF_Normalization();

  median_R_V = 4.05; // Value recommended by Calzetti 2000. should depend on galaxy properties
  Set_Median_Extinction(); // values are hard-coded here, will probably want to make it more clear how to modif

}


// default constructor, included because the compiler needs it to be hre
Galaxy::Galaxy()
{
  total_stellar_mass = 0.;
}


double Galaxy::Get_Mbh() const
{
  return mbh;
}

double Galaxy::Get_Mbh_bulge() const
{
  return mbh_bulge;
}

double Galaxy::Get_z() const
{
  return z;
}


double Galaxy::Get_m_g() const
{
  return m_g;
}

double Galaxy::Get_m_r() const
{
  return m_r;
}

double Galaxy::Get_sersic_n() const
{
  return sersic_n;
}

double Galaxy::Get_nuker_gammaprime() const
{
  return nuker_gammaprime;
}

double Galaxy::Get_r50_kpc() const
{
  return r50_kpc;
}

double Galaxy::Get_ssfr() const
{
  return ssfr;
}

double Galaxy::Get_M_u() const
{
  return UminusR + M_r;
}

double Galaxy::Get_M_r() const
{
  return M_r;
}


double Galaxy::Get_Mstar_Min() const
{
  return mstar_min;
}

double Galaxy::Get_Mstar_Max() const
{
  return mstar_max;
}

double Galaxy::Get_Total_Stellar_Mass() const
{
  return total_stellar_mass;
}

double Galaxy::Get_imf_norm() const
{
  return imf_normalization;
}

double Galaxy::Get_median_R_V() const
{
  return median_R_V;
}

double Galaxy::Get_median_A_V() const
{
  return median_A_V;
}

double Galaxy::Get_sigma_A_V() const
{
  return sigma_A_V;
}

double Galaxy::Get_median_A_Ha() const
{
  return median_A_Ha;
}

double Galaxy::Get_sigma_A_Ha() const
{
  return sigma_A_Ha;
}

double Galaxy::Get_Disruption_Rate_Normalization_Combined() const
{
  return disruption_rate_normalization_combined;
}

double Galaxy::Get_Disruption_Rate_Powerlaw_Mass() const
{
  return disruption_rate_powerlaw_mass;
}

double Galaxy::Get_Disruption_Rate_Powerlaw_Nuker() const
{
  return disruption_rate_powerlaw_nuker;
}

void Galaxy::Set_z(double new_z)
{
  z = new_z;
}

bool Galaxy::Is_ETG_UminusR_Mstar(double u_minus_r, double log10_mstar)
{

  // For a given u_minus_r and mstar, does the galaxy fall on the ETG side of that parameter plane
  // The slope of this cut seems fairly consistent across the literature
  // The u-r intercept can be moved up and down depending on redshift range and how strict you want the cut to be
  // Remember, this is approximating a separation in morphology by using a related, but ultimately distinct, cut on color and magnitude

  return u_minus_r > -0.4 + 0.25 * log10_mstar;
  
}


// See equation 7 of https://arxiv.org/pdf/1212.0939.pdf
// the params has the overall normalization, so as not to recompute it each time
double Galaxy::Kroupa_IMF_for_integrating(double mstar, void * p)
{
  struct imf_params {double norm; double min_mass; double max_mass;}; // not sure I really need to pass in min and max values any more
  struct imf_params params = *(struct imf_params *)p;

  
  if (mstar < 0.08)
    return 1./params.norm * pow(mstar,-0.3);
  else if (mstar >= 0.08 and mstar < 0.5)
    return 1./params.norm * 0.08 * pow(mstar,-1.3);
  else 
    return 1./params.norm * 0.08 * pow(0.5,-1.3 + 2.3) * pow(mstar,-2.3);
}

double Galaxy::Kroupa_IMF_for_value(double mstar, double norm) const
{

  if (mstar < 0.08)
    return 1./norm * pow(mstar,-0.3);
  else if (mstar >= 0.08 and mstar < 0.5)
    return 1./norm * 0.08 * pow(mstar,-1.3);
  else 
    return 1./norm * 0.08 * pow(0.5,-1.3 + 2.3) * pow(mstar,-2.3);

}


void Galaxy::Set_IMF_Normalization()
{
  int grid_size = 128;
  double relative_error = 1.e-6;

  gsl_function F;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  F.function = Galaxy::Kroupa_IMF_for_integrating;

  struct imf_params {double temp_norm; double min_mass; double max_mass;};

  imf_params params = {1.,mstar_min, mstar_max};
  F.params = &params;

  double result, error;
  gsl_integration_qags(&F,mstar_min, mstar_max, 0, relative_error, grid_size, workspace, &result, &error);
  gsl_integration_workspace_free(workspace);
  
  imf_normalization = result;

}


// equation 5 of Garn & Best 2010. Median A_Halpha for star-forming galaxies 
double Galaxy::GarnBest_Median_SF_AHa(double logmstar)
{

  if (logmstar < 8.66) return 0.292;
  if (logmstar > 11.96) return 2.164;
  
  double x = logmstar - 10.;

  return 0.91 + 0.77 * x + 0.11 * x * x - 0.09 * x * x * x;

}

double Galaxy::Cardelli_Extinction(double x, double rv)
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
        
    return a + b / rv;

}


double Galaxy::Calzetti_Extinction(double x, double rv)
{

  // Calzetti 2000 gives a function for k(lambda)
  // To make this analogous to the Cardelli extinction function, here we return k(lambda) / R_V
  

  //    if (x < 5./11. || x > 25./3.)
  //    {
  //      printf("ERROR: Calzetti extinction not defined for x < 0.45 or x > 8.33...\n");
  //      return 0.;
  //    }

    if (x < 1.6)
      {
	
	return std::max((2.659 * (-1.857 + 1.040 * x) + rv)/rv,0.);
      }
    else
      {
	return (2.659 * (-2.156 + 1.509 * x - 0.198 * x *x + 0.011 * x * x * x) + rv)/rv;
      }

}



void Galaxy::Set_Median_Extinction()
{

  if (log10(ssfr) >= -11.3) // maybe make this more sophisticated?
    {
      median_A_Ha = GarnBest_Median_SF_AHa(log10(total_stellar_mass));
      median_A_V = median_A_Ha / Calzetti_Extinction(1./(0.65645),median_R_V);
      sigma_A_Ha = 0.28;
      sigma_A_V = sigma_A_Ha / Calzetti_Extinction(1./(0.65645),median_R_V); // sigma_A_V not used directly

    }
  else
    {
      median_A_V = 0.2;
      sigma_A_V = 0.06;
      sigma_A_Ha = sigma_A_V * Calzetti_Extinction(1./(0.65645),median_R_V); //sigma_A_Ha not used directly

    }
}


// maybe consider moving this to physical constants
double Galaxy::Arcsec_From_Radian(double radians) const
{    
  return radians * 180. * 3600./ PI;
}


// maybe consider moving this to physical constants
double Galaxy::Radian_From_Arcsec(double r_arcsec) const
{    
  return r_arcsec * PI / (180. * 3600.);
}



double Galaxy::R_Arcsec_From_Kpc(double r_kpc) const
{
    
  double r_cm = r_kpc * 1000. * PARSEC; // convert to cm
  double d_A = AngularDiameterDistance(z); // in cm

  double r_radians = r_cm/d_A;
    
  return Arcsec_From_Radian(r_radians);
}


double Galaxy::R_Kpc_From_Arcsec(double r_arcsec) const
{
  double r_rad = Radian_From_Arcsec(r_arcsec);
    
  double d_A = AngularDiameterDistance(z);
    
  double r_cm = r_rad * d_A;
    
  return r_cm / (1000. * PARSEC);
}


double Galaxy::Find_Nuker_Gammaprime_From_Sersic()
{

  double rprime_kpc = R_Kpc_From_Arcsec(resolution_for_nuker_gamma);

  return sersic_bn/sersic_n * pow(rprime_kpc/r50_kpc,1./sersic_n);

}

// see reference mentioned at https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
double Galaxy::Get_Approx_Sersic_bn()
{
  return 1.9992 * sersic_n - 0.3271;
}


// as found in https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
// will be in units of magnitudes / arcsec^2 assuming m_tot in magnitude and r_eff_arcsec in arsec
double Galaxy::Get_Mu_Eff(double m_tot) const
{

  return m_tot + 5. * log10(re_arcsec) + 2.5 * log10(2. * PI * sersic_n * exp(sersic_bn)/pow(sersic_bn,2. * sersic_n) * gsl_sf_gamma(2. * sersic_n)  );
}

// equation 6 https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
double Galaxy::Get_Mu_Central(double m_tot) const
{

  double mu_e = Get_Mu_Eff(m_tot);
  return mu_e + 2.5 * sersic_bn /log(10.) * (pow(resolution_for_central_sb_measurement/re_arcsec,1./sersic_n) - 1.);
    // from formula, you do in fact want this log to be natural log, which for c++ math.h is just log
}

// equation 6 https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
double Galaxy::Get_Mu_At_Rkpc_CosmoCorrected(double r_kpc, double m_tot) const
{

  double r_arcsec = R_Arcsec_From_Kpc(r_kpc);
  double mu_e = Get_Mu_Eff(m_tot);
  double cosmo_correction = -10. * (log10(1. + z));
  return mu_e + 2.5 * sersic_bn /log(10.) * (pow(r_arcsec/re_arcsec,1./sersic_n) - 1.) + cosmo_correction;
    // from formula, you do in fact want this log to be natural log, which for c++ math.h is just log
}

double Galaxy::Get_Luminosity_Distance() const
{
  return LuminosityDistance(z);
}


// as found in https://ned.ipac.caltech.edu/level5/March05/Graham/Graham2.html
// I_e has units intensity / arcsec^2, r has units arscsec
double Galaxy::Flux_Enclosed_R_Sersic(double r, double I_e)
{
  
  return I_e * re_arcsec * re_arcsec * 2. * PI * sersic_n * exp(sersic_bn)/pow(sersic_bn,2. * sersic_n) * gsl_sf_gamma_inc_P(2. *sersic_n, sersic_bn * (pow(r/re_arcsec,1./sersic_n)))* gsl_sf_gamma(2. * sersic_n);

}

