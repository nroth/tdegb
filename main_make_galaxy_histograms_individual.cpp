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

  else if (base_name == "vol_flares")
    {
      ntuple_filename.append("vol_flare_ntuple_combined.dat");
      Histogram1dNtuple<vol_flare_data> hist1d(num_bins,bin_specs,base_name, axis_name,bin_index,ntuple_filename);

     	hist1d.Print_Histogram_1D(1, save_path);
     
    }

  else
    {
      printf("Error: histogram type (the type should be one of: gals, vol_disrupt, flares, or flares_unweighted)\n");
      exit(1);
    }

}

void make_Lg_bins(vector<string> identifiers, int num_bins, vector<double> bin_specs, int bin_index)
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
	hist1d.Print_Lg_Histogram_1D(1, save_path);
      if (base_name == "gals")
	hist1d.Print_Lg_Histogram_1D(0, save_path);

    }


    else if (base_name == "flares" || base_name == "flares_unweighted")
    {
      ntuple_filename.append("flare_ntuple_combined.dat");
      Histogram1dNtuple<flare_data> hist1d(num_bins,bin_specs,base_name, axis_name,bin_index,ntuple_filename);

      if (base_name == "flares")
	hist1d.Print_Lg_Histogram_1D(1, save_path);
      if (base_name == "flares_unweighted")
	hist1d.Print_Lg_Histogram_1D(0, save_path);
    }
  
    else if (base_name == "vol_flares")
    {
      ntuple_filename.append("vol_flare_ntuple_combined.dat");
      Histogram1dNtuple<vol_flare_data> hist1d(num_bins,bin_specs,base_name, axis_name,bin_index,ntuple_filename);

	hist1d.Print_volLg_Histogram_1D(1, save_path);
    }

}



void make_2d_bins(vector<string> identifiers, vector<size_t> num_bins, vector<vector<double>> bin_specs_2d, vector<int> bin_indices)
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


  vector<size_t> num_bins(2);
  vector<int> icols(2);
  vector<double> bin_specs(2);
  vector<vector<double>> bin_specs_2d;
  vector<string> identifiers(4); // careful, though, different for 1d vs 2d
  string data_path;
  string save_path;



  data_path = "./results/fifth_production_local/nuker_gamma_rate/fiducial/";
  save_path = "./results/fifth_production_local/nuker_gamma_rate/fiducial/";
  // make 1d histogram
  identifiers[2] = "flares";
  identifiers[3] = "log_Rbb_fit";
  //  identifiers[3] = "gammaprime";
  num_bins[0] = 40;
  bin_specs[0] = 12.;
  bin_specs[1] = 16.;
  icols[0] = 26; 
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  
}

