#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <string>
#include <ctime>
#include "galaxy.h"
#include "histogram1d_ntuple.h"
#include "histogram2d_ntuple.h"
#include "integration_helper_functions.h" // for now just to get flare data struct


using std::vector;
using std::string;


void make_1d_bins(vector<string> identifiers, int num_bins, vector<double> bin_specs, int bin_index)
{

  string data_path = identifiers[0];
  string save_path = identifiers[1];
  string base_name = identifiers[2];
  string axis_name = identifiers[3];

  string ntuple_filename = data_path;

  if (base_name == "gals" || base_name == "vol_disrupt")
    {

      ntuple_filename.append("gal_catalogue_ntuple_combined.dat");
      Histogram1dNtuple<galaxy_catalogue_data> hist1d(num_bins,bin_specs,base_name, axis_name,bin_index,ntuple_filename);

      if (base_name == "vol_disrupt")
	hist1d.Print_Histogram_1D(1, save_path);
      if (base_name == "gals")
	hist1d.Print_Histogram_1D(0, save_path);

    }
  else if (base_name == "flares" || base_name == "flares_unweighted")
    {
      ntuple_filename.append("flare_ntuple_combined.dat");
      Histogram1dNtuple<flare_data> hist1d(num_bins,bin_specs,base_name, axis_name,bin_index,ntuple_filename);

      if (base_name == "flares")
	hist1d.Print_Histogram_1D(1, save_path);
      if (base_name == "flares_unweighted")
	hist1d.Print_Histogram_1D(0, save_path);
    }
  else
    {
      printf("Error: histogram type (the type should be one of: gals, vol_disrupt, flares, or flares_unweighted)\n");
      exit(1);
    }

}


void make_2d_bins(vector<string> identifiers, vector<int> num_bins, vector<vector<double>> bin_specs_2d, vector<int> bin_indices)
{

  string data_path = identifiers[0];
  string save_path = identifiers[1];
  string base_name = identifiers[2];
  string axis1_name = identifiers[3];
  string axis2_name = identifiers[4];

  string ntuple_filename = data_path;

  if (base_name == "gals" || base_name == "vol_disrupt")
    {

      ntuple_filename.append("gal_catalogue_ntuple_combined.dat");
      Histogram2dNtuple<galaxy_catalogue_data> hist2d(num_bins,bin_specs_2d,base_name,axis1_name,axis2_name,bin_indices,ntuple_filename);

      if (base_name == "vol_disrupt")
	hist2d.Print_Histogram_2D_With_Header(1, save_path);
      if (base_name == "gals")
	hist2d.Print_Histogram_2D_With_Header(0, save_path);

    }
  else if (base_name == "flares" || base_name == "flares_unweighted")
    {
      ntuple_filename.append("flare_ntuple_combined.dat");
      Histogram2dNtuple<flare_data> hist2d(num_bins,bin_specs_2d,base_name, axis1_name,  axis2_name, bin_indices,ntuple_filename);

      if (base_name == "flares")
	hist2d.Print_Histogram_2D_With_Header(1, save_path);
      if (base_name == "flares_unweighted")
	hist2d.Print_Histogram_2D_With_Header(0, save_path);
    }
  else
    {
      printf("Error: histogram type (the type should be one of: gals, vol_disrupt, flares, or flares_unweighted)\n");
      exit(1);
    }

}


int main(int argc, char **argv)
{

  clock_t begin;
  clock_t end;
  float elapsed_secs;

  vector<int> num_bins(2);
  vector<int> icols(2);
  vector<double> bin_specs(2);
  vector<vector<double>> bin_specs_2d;
  vector<string> identifiers(4); // careful, though, different for 1d vs 2d
  string data_path;
  string save_path;


  // make 1d histogram
  data_path = "results/third_production_runs/v5/";
  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "A_V";
  num_bins[0] = 40;
  bin_specs[0] = 0.;
  bin_specs[1] = 2.;
  icols[0] = 18; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);



  // make 2d histogram
  identifiers.resize(5);
  data_path = "results/third_production_runs/v7/";
  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "UminusR";
  num_bins[0] = 27;
  bin_specs[0] = 0.75;
  bin_specs[1] = 3.45;
  icols[0] = -1; // function
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "mstar_mendel";
  num_bins[1] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[1] = mstar_mendel_i;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

}

