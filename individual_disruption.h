#ifndef INDIVIDUAL_DISRUPTION_H
#define INDIVIDUAL_DISRUPTION_H

#include <math.h>
#include <vector>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"

using std::vector;


#define LF_LOG_POWERLAW 1.5

#define MIN_LOG_LBOL 43.0
#define MAX_RADEFF_STREAMS 0.1 // maximum radiative efficiency of stream collisions.

#define MAX_EDD_RATIO 2.

// mstar in solar mass throughout


class Individual_Disruption {

 private:

  //  Galaxy * host_galaxy;
  double mbh; // copied from host galaxy
  double beta;
  double temperature;
  //  double mstar; // randomly sample this, too? make a star class?
  // spin orbit inclination

 public:

  //  Individual_Disruption(double, double,double); // for now, passing in a temperataure and beta, will want to change that
  void Initialize(double, double, double);

  double Get_Temperature();

  // helper functions
  double Main_Sequence_Radius(double);
  double Hills_Mass(double);
  double Eddington_Luminosity();
  double Peak_Mdot(double); 
  double Max_Luminosity(double);

  
};

#endif
