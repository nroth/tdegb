#ifndef SURVEY_H
#define SURVEY_H

#include "galaxy.h"


class Survey {

 private:

  double nu_gband;
  double nu_rband;

  double m_r_threshhold;
  double m_g_threshhold;

  double r_psf_arcsec;

  double host_contrast_cut;
  


 public:

  Survey();
  double mAB_From_Fnu(double);
  double Find_Host_Contrast_Magnitude(Galaxy, char);

  double Get_Nu_Gband();
  double Get_Nu_Rband();
  double Get_m_r_Threshhold();
  double Get_m_g_Threshhold();

  double I_From_Mu(double); // might want to move these somewhere else
  double Mu_From_I(double); // might want to movee these somewhere else

};

#endif
