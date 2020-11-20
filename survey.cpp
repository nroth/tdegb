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
#include <stdio.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "survey.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>




Survey::Survey()
{

  // see e.g. http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Swift/UVOT.UVW2&&mode=browse&gname=Swift&gname2=UVOT#filter
  // later, when reading in, put this in a loop

  nu_bands.push_back(C_LIGHT / (6257.54 * 1.e-8)); // LSSt r
  nu_bands.push_back(C_LIGHT / (4840.83 * 1.e-8)); // LSST g
  nu_bands.push_back(C_LIGHT / (2688.46 * 1.e-8)); // Swift UVW1
  nu_bands.push_back(C_LIGHT / (2272.71 * 1.e-8)); // Swift UVM2
  nu_bands.push_back(C_LIGHT / (2140.26 * 1.e-8)); // Swift UVW2

  m_r_threshhold = 23.;
  m_g_threshhold = 23.;

  host_contrast_cut = 0.;

  g_psf_arcsec = 0.7; // g-band median PSF FWHM. Page 8 of Bellm et al 2019 https://iopscience.iop.org/article/10.1088/1538-3873/aaecbe/pdf
  r_psf_arcsec = 0.7; // r-band median PSF FWHM. Page 8 of Bellm et al 2019 https://iopscience.iop.org/article/10.1088/1538-3873/aaecbe/pdf

  ///// Stuff related to temperature fit /////

  num_fit_p = 2; // should read this in?
  mag_frac_error = 0.1; // 0.1 mag error corresponds to ~ 0.1 fractional error
  nu_ref = 1.e15; // Hz
  T_ref = 3.e4; // K
  R_ref = 3.e14; // cm
  lnu_ref = 4. * PI * PI * 2. * H_PLANCK /(C_LIGHT * C_LIGHT) * nu_ref * nu_ref * nu_ref * R_ref * R_ref;
  xtol = 1.e-8; // tolerance for small step sizes
  gtol = 1.e-8; // tolerance for gradient near cost function minimum
  ftol = 0.; // By setting this to zero, which is the tolerance in the actual residual value, you force the other tolerances to be used (I think)

  Type = gsl_multifit_nlinear_trust; // GSL only provides this one option for the type right now (although it could be multifit or multilarge).
  fdf_params = gsl_multifit_nlinear_default_parameters(); // you might want to tune some of these

  Nbands = nu_bands.size(); 

  // using malloc here so I could put the declarations in the header
  nuprimes = (double*) malloc(Nbands * sizeof(double));
  y = (double*) malloc(Nbands * sizeof(double));
  weights = (double*) malloc(Nbands * sizeof(double));
  
  // The struct is storing arrays, which really means storing pointers. When nuprime or y are updated, the struct d knows about the changes
  d.n = Nbands;
  d.nu_prime = nuprimes;
  d.l_prime = y;
  d.alpha = H_PLANCK * nu_ref/( K_BOLTZ * T_ref);
  
  // the starting guesses for the fit parameters, for every time the fit is run
  x_init =  (double*) malloc(num_fit_p * sizeof(double));
  for (int i = 0; i < num_fit_p; i++)
    {
      x_init[i] = 1.;
    }

  wts = gsl_vector_view_array(weights, Nbands);
  x = gsl_vector_view_array (x_init, num_fit_p); // these might seem redundant since I do this again later in the fit, but I find I need them 

  // define the function to be minimized //
  fdf.f = Bbfit_f;
  fdf.df = Bbfit_df;   // set to NULL for finite-difference Jacobian //
  fdf.fvv = NULL;     // not using geodesic acceleration //
  fdf.n = Nbands;
  fdf.p = num_fit_p;
  fdf.params = &d;

  /* allocate workspace with default parameters */
  w = gsl_multifit_nlinear_alloc (Type, &fdf_params, Nbands, num_fit_p);

}

Survey::~Survey()
{

  free(nuprimes);
  free(y);
  free(weights);
  free(x_init);
}

int Survey::Get_Nbands()
{
  return Nbands;
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

double Survey::Fnu_From_mAB(double mAB)
{
  return pow(10., (mAB + 48.6)/(-2.5));
}

// find a way to make this passed galaxy constant?
double Survey::Find_Host_Contrast_Magnitude(Galaxy gal, char band)
{

  double m_tot;
  if (band == 'g')
    {
      m_tot = gal.Get_m_g();
    }
  else if (band == 'r')
    {
      m_tot = gal.Get_m_r(); 
    }
  else
    {
      printf("ERROR, specify g or r band for specifying host contrat magniutde");
      exit(1);
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

//////// Functions relating to temperature fit ////////

// Model, before nondimensionalization: Lnu_e = 4 pi^2 R^2 B_nu(T) 
double Survey::Bb_nondim(double nu_prime, double R_prime, double T_prime, double alpha)
{
  return pow(nu_prime,3) * pow(R_prime,2) /(exp(alpha * nu_prime/T_prime) - 1.);
}

// x is vector of parameters, in this case [R_prime, T_prime]
int Survey::Bbfit_f (const gsl_vector * x, void *phot_data, gsl_vector * f)
{
  size_t n = ((struct phot_data *)phot_data)->n;
  double *nu_prime = ((struct phot_data *)phot_data)->nu_prime;
  double *l_prime = ((struct phot_data *)phot_data)->l_prime;
  double alpha = ((struct phot_data *)phot_data)->alpha;

  double R_prime = gsl_vector_get (x, 0);
  double T_prime = gsl_vector_get (x, 1);

  size_t i;

  for (i = 0; i < n; i++)
    {
      double Yi = Bb_nondim(nu_prime[i],R_prime,T_prime,alpha); 
      gsl_vector_set (f, i, Yi - l_prime[i]);
    }

  return GSL_SUCCESS;
}


int Survey::Bbfit_df (const gsl_vector * x, void *phot_data, gsl_matrix * J)
{
  size_t n = ((struct phot_data *)phot_data)->n;
  double *nu_prime = ((struct phot_data *)phot_data)->nu_prime;
  double alpha = ((struct phot_data *)phot_data)->alpha;

  double R_prime = gsl_vector_get (x, 0);
  double T_prime = gsl_vector_get (x, 1);

  size_t i;

  for (i = 0; i < n; i++)
    {

      double expf = exp(alpha*nu_prime[i]/T_prime);
      double denom = expf - 1.;
      gsl_matrix_set (J, i, 0, 2. * pow(nu_prime[i],3) * R_prime / denom);
      gsl_matrix_set (J, i, 1, alpha * pow(nu_prime[i],4) * pow(R_prime,2) * expf/(pow(T_prime * denom,2)));
    }

  return GSL_SUCCESS;
}

int Survey::Provide_Temperature_Fit_Data(double* nu_rests, double* lnus, gsl_rng *r)
{

  for (int i = 0; i < Nbands; i++)
    {

      double yi = lnus[i]/lnu_ref;
      double si = mag_frac_error * yi; 
      double dy = gsl_ran_gaussian(r, si);

      nuprimes[i] = nu_rests[i] / nu_ref;
      y[i] = yi + dy; 
      weights[i] = 1.0/(si * si); 

    }

  int status = Perform_Temperature_Fit(r);
  return status;

}

int Survey::Perform_Temperature_Fit(gsl_rng *r)
{

  
  int status, info;

  d.n = Nbands;
  for (int i = 0; i < Nbands; i++)
    {
      d.nu_prime[i] = nuprimes[i];
      d.l_prime[i] = y[i];
    }

  // initialize solver with starting point and weights
  wts = gsl_vector_view_array(weights, Nbands); // if weights gets updated, wts does too
  x = gsl_vector_view_array (x_init, num_fit_p);
  gsl_multifit_nlinear_winit (&x.vector, &wts.vector, &fdf, w); // The x vector that's used in the fit, which can be accssed as w->x, gets reverted to what was set in x_init here
  
  // solve the system with a maximum of 100 iterations 
  status = gsl_multifit_nlinear_driver(100, xtol, gtol, ftol, NULL, NULL, &info, w); // the NULL arguments avoid a callback function

  fitted_Rbb = gsl_vector_get(w->x,0) * R_ref;
  fitted_Tbb = gsl_vector_get(w->x,1) * T_ref;
  fitted_Lbol = 4. * PI * pow(fitted_Rbb,2) * STEF_BOLTZ * pow(fitted_Tbb,4);

  return status;
}

double Survey::Get_Tbb_Fit()
{
  return fitted_Tbb;
}

double Survey::Get_Rbb_Fit()
{
  return fitted_Rbb;
}

double Survey::Get_Lbol_Fit()
{
  return fitted_Lbol;
}

void Survey::Set_Tbb_Fit(double new_T)
{
  fitted_Tbb = new_T;
}

void Survey::Set_Rbb_Fit(double new_R)
{
  fitted_Rbb = new_R;
}

void Survey::Set_Lbol_Fit(double new_L)
{
  fitted_Lbol = new_L;
}
