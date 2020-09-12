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

#ifndef DISRUPTION_H
#define DISRUPTION_H

#include "galaxy.h"
#include <gsl/gsl_rng.h>

class Disruption {

 private:

  Galaxy* host_gal;
  double mbh; // could in principle randomly sample this from scatter in galaxy observables correlation

  double mstar;
  double beta;
  double T_opt;
  double A_V; // can vary along lines of sight so will vary between disruptions
  double R_V; // can vary along lines of sight so will vary between disruptions

  double beta_mean;
  double beta_sigma;
  double T_opt_min;
  double T_opt_max;
  double T_opt_log_powerlaw;

  double L_Edd;
  double peak_L; // beak Lbol for a given disruption
  double max_L; // maximum possible Lbol for any disruption with a given mstar in this galaxy

  double lf_log_powerlaw;
  double min_log_lbol;
  double max_radeff_streams;
  double max_edd_ratio;

  // randomly generated LOS obscuration? or put that in survey?

  double Max_Luminosity(double);
  double Peak_Mdot();
  double Main_Sequence_Radius();
  double Eddington_Luminosity();

  double Planck_Function_Frequency(double, double);

 public:

  Disruption();
  Disruption(Galaxy*);

  void Rejection_Sample_Mstar(gsl_rng *);
  void Determine_Max_L();
  void Sample_Peak_L(gsl_rng *);
  void Sample_Beta(gsl_rng *);
  void Sample_Topt(gsl_rng *);
  void Sample_A_V(gsl_rng *);
  void Sample_R_V(gsl_rng *);

  double Get_Mstar();
  double Get_Peak_L();
  double Get_Topt();
  double Get_beta();
  double Get_A_V();
  double Get_R_V();
  double Get_Hills_Mass();
  double Get_Max_L();
  double Get_Eddington_Luminosity();
  double Get_Min_Log_Lbol();

  void Set_Topt(double);
  void Set_Beta(double);
  void Set_Peak_L(double);
  void Set_A_V(double);
  void Set_R_V(double);
  
  double Unobscured_Lnu(double, double, double);
  double Dust_Flux_Factor_Reduction(double);
  double Extincted_Flux_Observed(double, double);

};

#endif
