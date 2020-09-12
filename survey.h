/*
MIT License

Copyright (c) 2020 Nathaniel Jacob Roth

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
*/

#ifndef SURVEY_H
#define SURVEY_H

#include "galaxy.h"
#include <vector>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_rng.h>

using std::vector;

enum survey_bands {ZTF_r, ZTF_g, UVW1, UVM2, UVW2};

class Survey 
{
 private:

  vector<double> nu_bands;

  double m_r_threshhold;
  double m_g_threshhold;

  double g_psf_arcsec;
  double r_psf_arcsec;

  double host_contrast_cut;


  //////// Stuff related to temperature fit ////////


  struct phot_data
  {
    size_t n; // number of bands
    double * nu_prime; // rest frame, normalized
    double * l_prime; // rest frame, normalized
    double alpha; // for non-dimensionalization
  };

  size_t num_fit_p;
  double mag_frac_error; 
    
  // just reference values to make fits nondimensional. Results shouldn't change if these are changed
  double nu_ref;
  double T_ref;
  double R_ref;
  double lnu_ref;
  
  double xtol; // tolerance for small step sizes                   
  double gtol; // tolerance for gradient near cost function minimum
  double ftol; // tolerance for residual

  double fitted_Tbb;
  double fitted_Rbb;
  double fitted_Lbol;

  const gsl_multifit_nlinear_type *Type; 
  gsl_multifit_nlinear_workspace *w;
  gsl_multifit_nlinear_fdf fdf;
  gsl_multifit_nlinear_parameters fdf_params;


  double* nuprimes;
  double* y;
  double* weights;
  double* x_init;
  
  int Nbands;
  struct phot_data d;
  gsl_vector_view x;
  gsl_vector_view wts;

  static double Bb_nondim(double, double, double,double); // needs to be static for GSL
  static int Bbfit_f (const gsl_vector *, void *, gsl_vector *); // needs to be static for GSL
  static int Bbfit_df (const gsl_vector *, void *, gsl_matrix *); // needs to be static for GSL

  int Perform_Temperature_Fit(gsl_rng *);
  

 public:

  Survey();
  ~Survey();
  double mAB_From_Fnu(double);
  double Fnu_From_mAB(double);
  double Find_Host_Contrast_Magnitude(Galaxy, char);

  int Get_Nbands();
  double Get_Band_Nu(int);
  double Get_m_r_Threshhold();
  double Get_m_g_Threshhold();

  double I_From_Mu(double); // might want to move these somewhere else
  double Mu_From_I(double); // might want to movee these somewhere else

  int Provide_Temperature_Fit_Data(double *, double *, gsl_rng *);

  double Get_Tbb_Fit();
  double Get_Rbb_Fit();
  double Get_Lbol_Fit();

  void Set_Tbb_Fit(double);
  void Set_Rbb_Fit(double);
  void Set_Lbol_Fit(double);

};

#endif
