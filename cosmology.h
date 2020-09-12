/*
MIT License

Copyright (c) 2020 Nathaniel Jacob Roth

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE
*/

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
static double ComovingDistanceIntegrand(double z, void * params)
{

  //double param = *(double *) params;
  double H0_cgs = HUBBLE_PARAM * 100. * 1.e5 * 1./(1.e6 * PARSEC);
  return C_LIGHT / H0_cgs * 1./sqrt(OMEGA_M * pow(1. + z,3.) + OMEGA_LAMBDA);

}

static double ComovingDistance(double z)
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

static double ComovingVolumeIntegrand(double z, void * params)
{

  double d_c = ComovingDistance(z);

  return d_c * d_c * ComovingDistanceIntegrand(z,params);

}

// return in Mpc^3
static double ComovingVolume(double z_min, double z_max, double omega)
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

static double LuminosityDistance(double z)
{

  return (1. + z) * ComovingDistance(z);
  
}

static double AngularDiameterDistance(double z)
{

  return 1./(1. + z) * ComovingDistance(z);
  
}

static double DistanceModulus(double z)
{

  double ld = LuminosityDistance(z);
  return 5. * log10(ld/(10. * PARSEC));
  
}


#endif
