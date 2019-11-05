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

  double resolution_for_nuker_gamma; // arsec. Taken to be same for all galaxies. See Lauer et al 2007. Nick Stone's rate calculations were based on Nuker gamma as measured in this paper, so to convert n_sersic to nuker gamma we want to account for how they measured gamma

  double mbh; // what you will take to be the actual mbh, not what would be inferred from a survey by correlations.

  //derived quantities
  double re_arcsec;
  double sersic_bn;

  // include star formation history information? Time since starburst, durtaion of starburst?
  // could be used to set the upper truncation of IMF to make the approximate present-day mass function

  double imf_normalization;

  // mbh spin

  //helper functions
  double Get_Approx_Sersic_bn();
  double Find_Nuker_Gammaprime_From_Sersic();
  
  double R_Arcsec_From_Kpc(double);
  double R_Kpc_From_Arcsec(double);
  double Arcsec_From_Radian(double); // might want to move these somewhere else
  double Radian_From_Arcsec(double); // might want to move these somwhere else



public:
  
  Galaxy();  // default construtor, only here because it needs to be, filled with bogus values
  Galaxy(double*); //  including a pre-computed IMF normalization, and a temperature for the disrutpions. will want to clean this up

  double Get_Mbh() const; // const member functions can be called by constant objects
  double Get_z() const;
  double Get_m_g() const;
  double Get_m_r() const;
  double Get_sersic_n() const;
  double Get_nuker_gammaprime() const;
  double Get_r50_kpc() const;
  double Get_imf_norm() const;
  double Get_Luminosity_Distance() const;
  double Get_Mu_Eff(double) const;

  static double Kroupa_IMF_for_integrating(double, void *); // static needed for gsl integration to work
  double Kroupa_IMF_for_value(double, double) const; // static needed for gsl integration to work
  void Set_IMF_Normalization();

  double Flux_Enclosed_R_Sersic(double, double); // not totally sure this belongs here

  
};


#endif
