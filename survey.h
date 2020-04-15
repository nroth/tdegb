#ifndef SURVEY_H
#define SURVEY_H

#include "galaxy.h"
#include <vector>

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
  


 public:

  Survey();
  double mAB_From_Fnu(double);
  double Find_Host_Contrast_Magnitude(Galaxy, char);

double Get_Band_Nu(int);
  double Get_m_r_Threshhold();
  double Get_m_g_Threshhold();

  double I_From_Mu(double); // might want to move these somewhere else
  double Mu_From_I(double); // might want to movee these somewhere else

};

#endif
