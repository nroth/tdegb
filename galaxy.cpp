#include <math.h>
#include <gsl/gsl_integration.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "magnitudes.h"


//***************************************************************
// Constructors
//***************************************************************


Galaxy::Galaxy(double* galaxy_info)
{
  total_stellar_mass = galaxy_info[0];  // log
  mbh_sigma = pow(10.,galaxy_info[1]); // converting log to value
  mbh_bulge = pow(10.,galaxy_info[2]); // converting log to value
  z = galaxy_info[3];
  sersic_n = galaxy_info[4];
  r50_kpc = galaxy_info[5];
  m_g = galaxy_info[6];
  m_r = galaxy_info[7];
  nuker_gammaprime = Find_Nuker_Gammaprime_From_Sersic(sersic_n, r50_kpc, z);

  mbh = mbh_sigma; // making this choice for now

  Set_IMF_Normalization();

}

double Galaxy::Get_Mbh() const
{
  return mbh;
}

double Galaxy::Get_z() const
{
  return z;
}


double Galaxy::Get_m_g() const
{
  return m_g;
}

double Galaxy::Get_m_r() const
{
  return m_r;
}

double Galaxy::Get_sersic_n() const
{
  return sersic_n;
}

double Galaxy::Get_nuker_gammaprime() const
{
  return nuker_gammaprime;
}

double Galaxy::Get_r50_kpc() const
{
  return r50_kpc;
}

double Galaxy::Get_imf_norm() const
{
  return imf_normalization;
}



// See equation 7 of https://arxiv.org/pdf/1212.0939.pdf
// the params has the overall normalization, so as not to recompute it each time
double Galaxy::Kroupa_IMF_for_integrating(double mstar, void * params)
{
  double norm = *(double *) params;
  
  double m1 = 0.08;
  double m2 = 0.5;
  double m3 = MSTAR_MAX;

  double k1 = pow(m1,-0.3 + 1.3);
  double k2 = k1 * pow(m2,-1.3 + 2.3);
  //  double k3 = k2 * pow(m3, -2.3 + 2.3);

  if (mstar >= m1 && mstar < m2)
    return 1./norm * k1 * pow(mstar,-1.3);
  else if (mstar >= m2 && mstar <= m3)
      return 1./norm * k2 * pow(mstar,-2.3);
  else
    return 0.;
  
}

double Galaxy::Kroupa_IMF_for_value(double mstar, double norm) const
{
  double m1 = 0.08;
  double m2 = 0.5;
  double m3 = MSTAR_MAX;

  double k1 = pow(m1,-0.3 + 1.3);
  double k2 = k1 * pow(m2,-1.3 + 2.3);
  //  double k3 = k2 * pow(m3, -2.3 + 2.3);

  if (mstar >= m1 && mstar < m2)
    return 1./norm * k1 * pow(mstar,-1.3);
  else if (mstar >= m2 && mstar <= m3)
      return 1./norm * k2 * pow(mstar,-2.3);
  else
    return 0.;
  
}


void Galaxy::Set_IMF_Normalization()
{
  int grid_size = 128;
  double relative_error = 1.e-6;

  gsl_function F;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  F.function = Galaxy::Kroupa_IMF_for_integrating;

  double temp_norm = 1.;
  F.params = &temp_norm;

  double result, error;
  gsl_integration_qags(&F,MSTAR_MIN, MSTAR_MAX, 0, relative_error, grid_size, workspace, &result, &error);
  gsl_integration_workspace_free(workspace);
  
  imf_normalization = result;

}





