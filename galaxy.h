#ifndef GALAXY_H
#define GALAXY_H

#include <math.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "individual_disruption.h"
#include "magnitudes.h"

using std::vector;

#define MSTAR_MAX 1.0 // maximum stellar mass, in solar units. Really, should consider making this a function of galaxy properties
#define MSTAR_MIN 0.08 // Only considering main sequence stars

// In the future these will depend on galaxy properties
#define RATE_NORMALIZATION_COMBINED 2.9e-5 
#define RATE_POWERLAW -0.404
#define RATE_POWERLAW_NUKER 0.705

// mstar in solar mass throughout


class Galaxy {

 private:

  double total_stellar_mass;
  double mbh_sigma;
  double mbh_bulge;
  double z;
  double sersic_n;
  double r50_kpc;
  double nuker_gammaprime;

  double mbh; // what you will take to be the actual mbh, not what would be inferred from a survey by correlations.

  // include star formation history information? Time since starburst, durtaion of starburst?

  double imf_normalization;

  // mbh spin

  Individual_Disruption disruptions;

 public:

  Galaxy(); //  including a pre-computed IMF normalization, and a temperature for the disrutpions. will want to clean this up
  
  Galaxy(double*, double, double, double); //  including a pre-computed IMF normalization, and a temperature for the disrutpions. will want to clean this up

  double Get_Mbh();

  static double Kroupa_IMF(double, void *); // static needed for gsl integration to work
  double Determine_IMF_Normalization();

  // these perform integrals using the static functions defined below
  double Total_Disruption_Rate();
  double Total_Disruption_Rate_Observed_Gband(double);
  double Total_Disruption_Rate_Observed_Rband(double);

  // helper functions. Made static for gsl integration. Need to have pointers to a galaxy instance to have access to data. Be careful here.
  static double Stellar_Rate_Integrand_GbandCut(double, void *);
  static double Stellar_Rate_Integrand_RbandCut(double, void *);
  static double Stellar_Disruption_Rate(double, Galaxy *);
  static double Stellar_Rate_Integrand(double, void *);
  static double Fraction_Observed (double, double, Galaxy*);


  
};


#endif
