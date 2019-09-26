#ifndef COSMOLOGY_H
#define COSMOLOGY_H

#include <gsl/gsl_integration.h>
#include <math.h>
#include "physical_constants.h"


// ASSUMES FLAT UNIVERSE THROUGHOUT

#define HUBBLE_PARAM  0.7 // i.e. little h
#define OMEGA_M  0.3
#define OMEGA_LAMBDA  0.7


//function needs to have this signature for params
double ComovingDistanceIntegrand(double z, void * params)
{

  //double param = *(double *) params;
  double H0_cgs = HUBBLE_PARAM * 100. * 1.e5 * 1./(1.e6 * PARSEC);
  return C_LIGHT / H0_cgs * 1./sqrt(OMEGA_M * pow(1. + z,3.) + OMEGA_LAMBDA);

}

double ComovingDistance(double z)
{

  //    double ComovingDistanceIntegrand(double, void *); // function needs to have this signature for gsl

  int grid_size = 256;
  double relative_error = 1.e-7;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = &ComovingDistanceIntegrand;

  //  double example_parameter = 1.;
  //  F.params = &example_parameter;

  double result, error; 
  
  gsl_integration_qags(&F, 0, z, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return result;

}

double ComovingVolumeIntegrand(double z, void * params)
{

  double d_c = ComovingDistance(z);

  return d_c * d_c * ComovingDistanceIntegrand(z,params);

}

// return in Mpc^3
double ComovingVolume(double z_min, double z_max, double omega)
{

  int grid_size = 256;
  double relative_error = 1.e-7;
  gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(grid_size);
  gsl_function F;
  F.function = &ComovingVolumeIntegrand;

  //  double example_parameter = 1.;
  //  F.params = &example_parameter;

  double result, error; 
  
  gsl_integration_qags(&F, z_min, z_max, 0, relative_error, grid_size, workspace, &result, &error);

  gsl_integration_workspace_free(workspace);

  return omega * result / pow(1.e6 * PARSEC,3.);


}

double LuminosityDistance(double z)
{

  return (1. + z) * ComovingDistance(z);
  
}

double AngularDiameterDistance(double z)
{

  return 1./(1. + z) * ComovingDistance(z);
  
}

double DistanceModulus(double z)
{

  double ld = LuminosityDistance(z);
  return 5. * log10(ld/(10. * PARSEC));
  
}


#endif
