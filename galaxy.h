#ifndef GALAXY_H
#define GALAXY_H


// GSL ntuple requires a struct
struct galaxy_catalogue_data
{
  double attributes[17];
  double weight; // for determining volumetric disruption rate
};

// just for convenience
enum galaxy_catalogue_indices{ mstar_mendel_i, mbh_sigma_i,  mbh_bulge_i, z_i, sersic_n_i, r50_kpc_i, m_g_i, m_r_i,ssfr_i,UminusR_i,M_r_i,nuker_gamma_i,median_A_V_i,sigma_A_V_i };


class Galaxy {

 private:

  // entries contained in the catalogue data struct
  double total_stellar_mass; // in solar mass
  double m_g;
  double m_r;
  double mbh_sigma;
  double mbh_bulge;
  double z;
  double sersic_n;
  double r50_kpc;
  double ssfr;
  double UminusR;
  double M_r;
  // include star formation history information? Time since starburst, durtaion of starburst?
  // could be used to set the upper truncation of IMF to make the approximate present-day mass function

  double resolution_for_nuker_gamma; // arcsec. Taken to be same for all galaxies. See Lauer et al 2007. Nick Stone's rate calculations were based on Nuker gamma as measured in this paper, so to convert n_sersic to nuker gamma we want to account for how they measured gamma

  double resolution_for_central_sb_measurement; 
  
  double mbh; // what you will take to be the actual mbh, not what would be inferred from a survey by correlations.

  double mstar_max;
  double mstar_min;

  double disruption_rate_normalization_combined;
  double disruption_rate_powerlaw_mass;
  double disruption_rate_powerlaw_nuker;

  //derived quantities
  double re_arcsec;
  double sersic_bn;
  double nuker_gammaprime;
  double median_R_V;
  double median_A_Ha;
  double sigma_A_Ha;
  double median_A_V;
  double sigma_A_V;

  double imf_normalization;

  // mbh spin

  //helper functions
  double Get_Approx_Sersic_bn();
  double Find_Nuker_Gammaprime_From_Sersic();
  
  double R_Arcsec_From_Kpc(double) const;
  double R_Kpc_From_Arcsec(double) const;
  double Arcsec_From_Radian(double) const; // might want to move these somewhere else
  double Radian_From_Arcsec(double) const; // might want to move these somwhere else



public:
  
  Galaxy();  // default construtor, only here because it needs to be, filled with bogus values
  Galaxy(double*); //  including a pre-computed IMF normalization, and a temperature for the disrutpions. will want to clean this up

  double Get_Mbh() const; // const member functions can be called by constant objects
  double Get_Mbh_bulge() const; 
  double Get_z() const;
  double Get_m_g() const;
  double Get_m_r() const;
  double Get_M_u() const;
  double Get_M_r() const;
  double Get_sersic_n() const;
  double Get_nuker_gammaprime() const;
  double Get_r50_kpc() const;
  double Get_ssfr() const;
  double Get_Mstar_Min() const;
  double Get_Mstar_Max() const;
  double Get_Total_Stellar_Mass() const;
  double Get_imf_norm() const;
  double Get_median_R_V() const;
  double Get_median_A_Ha() const;
  double Get_sigma_A_Ha() const;
  double Get_median_A_V() const;
  double Get_sigma_A_V() const;
  double Get_Luminosity_Distance() const;
  double Get_Disruption_Rate_Normalization_Combined() const;
  double Get_Disruption_Rate_Powerlaw_Mass() const;
  double Get_Disruption_Rate_Powerlaw_Nuker() const;

  void Set_z(double);

  static bool Is_ETG_UminusR_Mstar(double, double);

  double Get_Mu_Eff(double) const;
  double Get_Mu_Central(double) const;
  double Get_Mu_At_Rkpc_CosmoCorrected(double, double) const;

  static double Kroupa_IMF_for_integrating(double, void *); // static needed for gsl integration to work
  double Kroupa_IMF_for_value(double, double) const; // static needed for gsl integration to work
  double GarnBest_Median_SF_AHa(double);

  double Calzetti_Extinction(double, double);
  double Cardelli_Extinction(double, double);  
  
  void Set_IMF_Normalization();
  void Set_Median_Extinction();

  double Flux_Enclosed_R_Sersic(double, double); // not totally sure this belongs here

  
};


#endif
