#ifndef DISRUPTION_H
#define DISRUPTION_H

#include "galaxy.h"
#include <gsl/gsl_rng.h>

#define LF_LOG_POWERLAW 1.5

#define MIN_LOG_LBOL 43.0
#define MAX_RADEFF_STREAMS 0.1 // maximum radiative efficiency of stream collisions.

#define MAX_EDD_RATIO 2.

class Disruption {

 private:

  Galaxy host_gal;
  double mstar;
  double mbh; // could in principle randomly sample this from scatter in galaxy observables correlation
  double beta;
  double T_opt;
  double A_V; // can vary along lines of sight so will vary between disruptions
  double R_V; // can vary along lines of sight so will vary between disruptions

  double L_Edd;
  double peak_L; // beak Lbol for a given disruption
  double max_L; // maximum possible Lbol for any disruption in this host galaxy

  // randomly generated LOS obscuration? or put that in survey?

  double Max_Luminosity(double);
  double Peak_Mdot();
  double Main_Sequence_Radius();
  double Eddington_Luminosity();

  double Planck_Function_Frequency(double, double);
  double Unobscured_Lnu(double, double, double);
  double Cardelli_Extinction(double);
  double Dust_Flux_Factor_Reduction(double);

 public:

  Disruption(Galaxy);

  void Rejection_Sample_Mstar(gsl_rng *);
  void Determine_Max_L();
  void Sample_Peak_L(gsl_rng *);

  double Get_Mstar();
  double Get_Peak_L();
  double Get_Topt();
  double Get_beta();
  double Get_Hills_Mass();
  double Get_Max_L();
  

  double Extincted_Flux_Observed(double, double);

};

#endif
