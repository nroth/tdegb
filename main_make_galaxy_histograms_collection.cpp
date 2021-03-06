/*

Copyright (c) 2020 Nathaniel Jacob Roth, licensed under GPL-3.0-or-later

This file is part of tdegb

tdegb is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

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


  // MODIFY FOR YOUR NEEDS
  data_path = "./results/fifth_production_local/nuker_gamma_rate/fiducial/";
  save_path = "./results/fifth_production_local/nuker_gamma_rate/fiducial/";



    // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Lg";
  num_bins[0] = 75;
  bin_specs[0] = 40.5;
  bin_specs[1] = 45.;
  icols[0] = 0; // Lg is special
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_Lg_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

      // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "A_V";
  num_bins[0] = 75;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.75;
  icols[0] = 23;
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

        // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "A_V_median";
  num_bins[0] = 75;
  bin_specs[0] = 0.;
  bin_specs[1] = 2.5;
  icols[0] = 12;
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  /*

  //  identifiers[2] = "vol_flares";
  //identifiers[3] = "log_volLg";
  //num_bins[0] = 40;
  //bin_specs[0] = 41.;
  //bin_specs[1] = 45.;
  //icols[0] = 0; // Lg is special
  //identifiers[0] = data_path;
  //identifiers[1] = save_path;
  //begin = clock();
  //make_Lg_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  //end = clock();
  //elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  //printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  */

  // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "mstar_mendel";
  num_bins[0] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[0] = 0; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Lopt_fit";
  num_bins[0] = 25;
  bin_specs[0] = 41;
  bin_specs[1] = 45;
  icols[0] = 25; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

    

  // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Lopt";
  num_bins[0] = 25;
  bin_specs[0] = 41;
  bin_specs[1] = 45;
  icols[0] = 19; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

    // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Edd_ratio_fit";
  num_bins[0] = 28;
  bin_specs[0] = -4.;
  bin_specs[1] = 0.;
  icols[0] = 27; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  
  // make 1d histogram
  //  save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Edd_ratio";
  num_bins[0] = 25;
  bin_specs[0] = -4.;
  bin_specs[1] = 0.;
  icols[0] = 27; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

    // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "Topt_fit";
  num_bins[0] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 5.e4;
  icols[0] = 24; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  
  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Topt_fit";
  num_bins[0] = 25;
  bin_specs[0] = log10(1.e3);
  bin_specs[1] = log10(5.e4);
  icols[0] = 30; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  
    // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "Topt";
  num_bins[0] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 5.e4;
  icols[0] = 20; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

      // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Topt";
  num_bins[0] = 25;
  bin_specs[0] = log10(1.e3);
  bin_specs[1] = log10(5.e4);
  icols[0] = 29; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 25;
  bin_specs[0] = 3.;
  bin_specs[1] = 9.;
  icols[0] = 1; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


    // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 25;
  bin_specs[0] = 3.;
  bin_specs[1] = 9.;
  icols[0] = 1; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  

  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Rbb_fit";
  num_bins[0] = 25;
  bin_specs[0] = 12.;
  bin_specs[1] = 16.;
  icols[0] = 26; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

    // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "mstar";
  num_bins[0] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.;
  icols[0] = 17; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_ssfr";
  num_bins[0] = 25;
  bin_specs[0] = -13.;
  bin_specs[1] = -9.;
  icols[0] = 8; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

    // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "z";
  num_bins[0] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[0] = 3; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "Mu_100pc";
  num_bins[0] = 25;
  bin_specs[0] = 15.;
  bin_specs[1] = 25;
  icols[0] = 16; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);

  
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Lopt";
  num_bins[0] = 50;
  bin_specs[0] = 42.;
  bin_specs[1] = 45.;
  icols[0] = 19; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Lopt_fit";
  num_bins[1] = 50;
  bin_specs[0] = 41.;
  bin_specs[1] = 45.;
  icols[1] = 25;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again





  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_mstar_mendel";
  num_bins[0] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[0] = 0; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_ssfr";
  num_bins[1] = 25;
  bin_specs[0] = -13.;
  bin_specs[1] = -9.;
  icols[1] = 8;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  
    // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "log_mstar_mendel";
  num_bins[0] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[0] = 0; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_ssfr";
  num_bins[1] = 25;
  bin_specs[0] = -13.;
  bin_specs[1] = -9.;
  icols[1] = 8;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again


      // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  
  //save_path = data_path;
  identifiers[2] = "vol_disrupt";
  identifiers[3] = "log_mstar_mendel";
  num_bins[0] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[0] = 0; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_ssfr";
  num_bins[1] = 25;
  bin_specs[0] = -13.;
  bin_specs[1] = -9.;
  icols[1] = 8;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again
  

  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "log_mstar_mendel";
  num_bins[0] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[0] = 0; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "z";
  num_bins[1] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[1] = 3;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

    // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  
  //save_path = data_path;
  identifiers[2] = "vol_disrupt";
  identifiers[3] = "log_mstar_mendel";
  num_bins[0] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[0] = 0; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "z";
  num_bins[1] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[1] = 3;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  

  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_mstar_mendel";
  num_bins[0] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[0] = 0; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "z";
  num_bins[1] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[1] = 3;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again



      // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Topt_fit";
  num_bins[0] = 25;
  bin_specs[0] = 4.;
  bin_specs[1] = 4.8;
  icols[0] = 30; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Lopt_fit";
  num_bins[1] = 25;
  bin_specs[0] = 41.;
  bin_specs[1] = 45.;
  icols[1] = 25;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

        // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Topt";
  num_bins[0] = 25;
  bin_specs[0] = 4.0;
  bin_specs[1] = 4.8;
  icols[0] = 29; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Lopt";
  num_bins[1] = 25;
  bin_specs[0] = 41.;
  bin_specs[1] = 45.;
  icols[1] = 19;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again


  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Rbb_fit";
  num_bins[0] = 25;
  bin_specs[0] = 14.;
  bin_specs[1] = 15.5;
  icols[0] = 26; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Lopt_fit";
  num_bins[1] = 25;
  bin_specs[0] = 41.;
  bin_specs[1] = 45.;
  icols[1] = 25;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

        // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Rbb_fit";
  num_bins[0] = 25;
  bin_specs[0] = 14.;
  bin_specs[1] = 15.5;
  icols[0] = 26; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Lopt";
  num_bins[1] = 25;
  bin_specs[0] = 41.;
  bin_specs[1] = 45.;
  icols[1] = 19;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  
  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Rbb_fit";
  num_bins[0] = 25;
  bin_specs[0] = 14.;
  bin_specs[1] = 15.5;
  icols[0] = 26; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Topt_fit";
  num_bins[1] = 25;
  bin_specs[0] = 4.;
  bin_specs[1] = 4.8;
  icols[1] = 30;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_Rbb_fit";
  num_bins[0] = 25;
  bin_specs[0] = 14.;
  bin_specs[1] = 15.5;
  icols[0] = 26; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Topt";
  num_bins[1] = 25;
  bin_specs[0] = 4.;
  bin_specs[1] = 4.8;
  icols[1] = 29;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again


  
      // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "UminusR";
  num_bins[0] = 27;
  bin_specs[0] = 0.75;
  bin_specs[1] = 3.;
  icols[0] = 9; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_mstar_mendel";
  num_bins[1] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[1] = 0;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

    
  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "vol_disrupt";
  identifiers[3] = "UminusR";
  num_bins[0] = 27;
  bin_specs[0] = 0.75;
  bin_specs[1] = 3.;
  icols[0] = 9; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_mstar_mendel";
  num_bins[1] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[1] = 0;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again
  


  // make 2d histogram
    if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "UminusR";
  num_bins[0] = 27;
  bin_specs[0] = 0.75;
  bin_specs[1] = 3.0;
  icols[0] = 9; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_mstar_mendel";
  num_bins[1] = 25;
  bin_specs[0] = 8.;
  bin_specs[1] = 12.;
  icols[1] = 0;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  
  
  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 25;
  bin_specs[0] = 3.;
  bin_specs[1] = 10.;
  icols[0] = 1; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "z";
  num_bins[1] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[1] = 3;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

    // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "vol_disrupt";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 25;
  bin_specs[0] = 3.;
  bin_specs[1] = 10.;
  icols[0] = 1; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "z";
  num_bins[1] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[1] = 3;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  

  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 25;
  bin_specs[0] = 3.;
  bin_specs[1] = 10.;
  icols[0] = 1; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "z";
  num_bins[1] = 25;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[1] = 3;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again


  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "Topt";
  num_bins[0] = 50;
  bin_specs[0] = 1.e4;
  bin_specs[1] = 3.e4;
  icols[0] = 20; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "Topt_fit";
  num_bins[1] = 50;
  bin_specs[0] = 5.e3;
  bin_specs[1] = 3.e4;
  icols[1] = 24;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  
  // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 25;
  bin_specs[0] = 3.;
  bin_specs[1] = 9.;
  icols[0] = 1; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_ssfr";
  num_bins[1] = 25;
  bin_specs[0] = -13.;
  bin_specs[1] = -9.;
  icols[1] = 8;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  
    // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 25;
  bin_specs[0] = 3.;
  bin_specs[1] = 9.;
  icols[0] = 1; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_ssfr";
  num_bins[1] = 25;
  bin_specs[0] = -13.;
  bin_specs[1] = -9.;
  icols[1] = 8;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again
  

      // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "log_mbh_sigma";
  num_bins[0] = 30;
  bin_specs[0] = 4.;
  bin_specs[1] = 8.5;
  icols[0] = 1; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "log_Lopt_fit";
  num_bins[1] = 30;
  bin_specs[0] = 42.;
  bin_specs[1] = 45.5;
  icols[1] = 25;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

  
      // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "z";
  num_bins[0] = 30;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[0] = 3; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "Topt_fit";
  num_bins[1] = 30;
  bin_specs[0] = 5.e3;
  bin_specs[1] = 5.e4;
  icols[1] = 24;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again


    
      // make 2d histogram
  if (identifiers.size() != 5) identifiers.resize(5);
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "z";
  num_bins[0] = 30;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[0] = 3; 
  bin_specs_2d.push_back(bin_specs);
  identifiers[4] = "Mu_100pc_cut";
  num_bins[1] = 30;
  bin_specs[0] = 12.;
  bin_specs[1] = 17.;
  icols[1] = 16;
  bin_specs_2d.push_back(bin_specs);
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_2d_bins(identifiers,num_bins,bin_specs_2d,icols);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
  bin_specs_2d.clear(); // had put this here in case you wanted to use it again

    



  // CAREFUL ABOUT ZCUT HERE
  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "Mu_100pc";
  num_bins[0] = 25;
  bin_specs[0] = 15.;
  bin_specs[1] = 25;
  icols[0] = 16; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  

  // CAREFUL ABOUT ZCUT HERE
  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "Mu_100pc";
  num_bins[0] = 25;
  bin_specs[0] = 15.;
  bin_specs[1] = 25.;
  icols[0] = 16; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  
  identifiers[2] = "flares";
  identifiers[3] = "z_coarse";
  num_bins[0] = 8;
  bin_specs[0] = 0.;
  bin_specs[1] = 1.0;
  icols[0] = 3; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);



  identifiers[2] = "flares";
  identifiers[3] = "log_mstar_coarse";
  num_bins[0] = 14;
  bin_specs[0] = 8.25;
  bin_specs[1] = 11.75;
  icols[0] = 0; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  identifiers[2] = "flares";
  identifiers[3] = "log_Lopt_coarse";
  num_bins[0] = 8;
  bin_specs[0] = 41.8;
  bin_specs[1] = 45.;
  icols[0] = 25; 
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);



  identifiers[2] = "flares";
  identifiers[3] = "log_Topt_fit_coarse";
  num_bins[0] = 7;
  bin_specs[0] = 4. - 4./30.;
  bin_specs[1] = 4.8;
  icols[0] = 30; 
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  identifiers[2] = "flares";
  identifiers[3] = "log_Ropt_fit_coarse";
  num_bins[0] = 12;
  bin_specs[0] =13.7 ;
  bin_specs[1] = 16.1;
  icols[0] = 26; 
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);


  identifiers[2] = "flares";
  identifiers[3] = "log_Lg_coarse";
  num_bins[0] = 8;
  bin_specs[0] = 41.00 ;
  bin_specs[1] = 45.00;
  icols[0] = 0; // Lg is special
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_Lg_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);



  // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "gals";
  identifiers[3] = "m_r";
  num_bins[0] = 150;
  bin_specs[0] = 10.;
  bin_specs[1] = 25.;
  icols[0] = 7; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);



    // make 1d histogram
  //save_path = data_path;
  identifiers[2] = "flares";
  identifiers[3] = "m_r";
  num_bins[0] = 150;
  bin_specs[0] = 10.;
  bin_specs[1] = 25.;
  icols[0] = 7; // use an enum, avoid specifying this again? Also end of input
  identifiers[0] = data_path;
  identifiers[1] = save_path;
  begin = clock();
  make_1d_bins(identifiers,num_bins[0],bin_specs,icols[0]);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  
}






