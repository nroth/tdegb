#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <ctime>
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
  


  // create 1d histogram
  begin = clock();
  base_name = "gals";
  spec[0] = 0.;
  spec[1] = 1.;
  h_name = "z";
  icols[0] = z_i;
  num_bins[0] = 25;
  Histogram1dNtuple<galaxy_catalogue_data> this_hist1d(num_bins[0],spec,base_name, h_name,icols[0],ntuple_filename);
  this_hist1d.Print_Histogram_1D(0);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  begin = clock();
  base_name = "vol_disrupt";
  v_name = "z"; // for now, if you're going to be computing a function, put it here
  num_bins[0] = 20;
  spec[0] = 0.;
  spec[1] = 1.;
  icols[0] = z_i; // function
  bin_specs.push_back(spec);
  h_name = "log_mbh_sigma";
  num_bins[1] = 40;
  spec[0] = 3.;
  spec[1] = 9.;
  icols[1] = mbh_sigma_i;
  bin_specs.push_back(spec);
  ibin = 0;
  Histogram2dNtuple<galaxy_catalogue_data> this_hist2d(num_bins,bin_specs,base_name, v_name,h_name,icols,ibin, ntuple_filename);
  bin_specs.clear();
  this_hist2d.Print_Histogram_2D_With_Header(1);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  

}

