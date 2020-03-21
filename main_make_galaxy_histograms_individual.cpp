#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_ntuple.h>
#include "galaxy.h"
#include "histogram1d_ntuple.h"
#include "histogram2d_ntuple.h"


using std::vector;
using std::string;

int main(int argc, char **argv)
{


  clock_t begin;
  clock_t end;
  float elapsed_secs;

  begin = clock();

  vector<int> num_bins(2);
  vector<int> icols(2);
  vector<double> spec(2);
  int ibin;
  vector<vector<double>> bin_specs;
  string base_name;
  string v_name("default"); 
  string h_name("default");

  
  string ntuple_filename = "gal_catalogue_ntuple_combined.dat";
  

  base_name = "gals";


  // create 1d histogram
  begin = clock();
  spec[0] = 0.;
  spec[1] = 1.;
  h_name = "z";
  icols[0] = z_i;
  num_bins[0] = 30;
  Histogram1dNtuple this_hist1d(num_bins[0],spec,base_name, h_name,icols[0],ntuple_filename);
  this_hist1d.Print_Histogram_1D(1);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);



  begin = clock();
  v_name = "UminusR"; // for now, if you're going to be computing a function, put it here
  h_name = "mstar_mendel";
  num_bins[0] = 27;
  spec[0] = 0.75;
  spec[1] = 3.45;
  bin_specs.push_back(spec);
  num_bins[1] = 25;
  spec[0] = 8.;
  spec[1] = 12.;
  bin_specs.push_back(spec);
  icols[0] = -1; // function
  icols[1] = mstar_mendel_i;
  ibin = 0;
  Histogram2dNtuple this_hist2d(num_bins,bin_specs,base_name, v_name,h_name,icols,ibin, ntuple_filename);
  bin_specs.clear();
  this_hist2d.Print_Histogram_2D_With_Header(1);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);


  

}

