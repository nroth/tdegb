#ifndef GALAXY_H
#define GALAXY_H

// In the future these will depend on galaxy properties
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
  double m_g;
  double m_r;
  double mbh_sigma;
  double mbh_bulge;
  double z;
  double sersic_n;
  double r50_kpc;
  double nuker_gammaprime;

  double mbh; // what you will take to be the actual mbh, not what would be inferred from a survey by correlations.

  // include star formation history information? Time since starburst, durtaion of starburst?
  // could be used to set the upper truncation of IMF to make the approximate present-day mass function

  double imf_normalization;

  // mbh spin

 public:

  Galaxy(double*); //  including a pre-computed IMF normalization, and a temperature for the disrutpions. will want to clean this up

  double Get_Mbh() const; // const member functions can be called by constant objects
  double Get_z() const;
  double Get_m_g() const;
  double Get_m_r() const;
  double Get_sersic_n() const;
  double Get_nuker_gammaprime() const;
  double Get_r50_kpc() const;
  double Get_imf_norm() const;

  static double Kroupa_IMF_for_integrating(double, void *); // static needed for gsl integration to work
  double Kroupa_IMF_for_value(double, double) const; // static needed for gsl integration to work
  void Set_IMF_Normalization();


  
};


#endif
