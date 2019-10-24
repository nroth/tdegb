#ifndef DISRUPTION_H
#define DISRUPTION_H

#include "galaxy.h"

#define LF_LOG_POWERLAW 1.5

#define MIN_LOG_LBOL 43.0
#define MAX_RADEFF_STREAMS 0.1 // maximum radiative efficiency of stream collisions.

#define MAX_EDD_RATIO 2.

class Disruption {

 private:

  //  double mstar;
  double mbh; // could in principle randomly sample this from scatter in galaxy observables correlation
  double beta;
  double T;

  // randomly generated LOS obscuration? or put that in survey?

 public:

  Disruption(Galaxy);

  double Get_T();
  double Get_beta();
  
  double Max_Luminosity(double);
  
  //Helper functions. Make them private?
  double Peak_Mdot(double);
  double Main_Sequence_Radius(double);
  double Hills_Mass(double);
  double Eddington_Luminosity();
  


};

#endif
