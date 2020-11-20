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
#include <math.h>
#include <vector>
#include <string>
#include <ctime>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_ntuple.h>
#include "histogram1d_ntuple.h"
#include "histogram2d_ntuple.h"

using std::vector;
using std::string;

int main(int argc, char **argv)
{

  MPI_Init( &argc, &argv );
  int my_rank,n_procs;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);

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

  if (my_rank == 0)
    {
      // create 1d histogram
      begin = clock();
      spec[0] = 8.;
      spec[1] = 12.;
      h_name = "mstar_mendel";
      icols[0] = 0;
      num_bins[0] = 24;
      Histogram1dNtuple this_hist1d(num_bins[0],spec,base_name, h_name,icols[0],ntuple_filename);
      this_hist1d.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);


      // create 1d histogram
      begin = clock();
      spec[0] = 0.;
      spec[1] = 1.;
      h_name = "z";
      icols[0] = 3;
      num_bins[0] = 20;
      this_hist1d.Init(num_bins[0],spec,base_name,h_name,icols[0],ntuple_filename);
      this_hist1d.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);

      
      // create 1d histogram
      begin = clock();
      spec[0] = 3.;
      spec[1] = 9.;
      h_name = "log_mbh_sigma";
      icols[0] = 1;
      num_bins[0] = 40;
      this_hist1d.Init(num_bins[0],spec,base_name, h_name,icols[0],ntuple_filename);
      this_hist1d.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);


    }

  if (my_rank == 1)
    {


      // create 2d histogram
      begin = clock();
      v_name = "z";
      h_name = "log_mbh_sigma";
      num_bins[0] = 20;
      spec[0] = 0;
      spec[1] = 1.;
      bin_specs.push_back(spec);
      num_bins[1] = 40;
      spec[0] = 3.;
      spec[1] = 9.;
      bin_specs.push_back(spec);
      icols[0] = 3;
      icols[1] = 1;
      ibin = 0;
      Histogram2dNtuple this_hist2d(num_bins,bin_specs,base_name, v_name,h_name,icols,ibin, ntuple_filename);
      bin_specs.clear();
      //this_hist2d.Print_Histogram_2D_With_Header();
      this_hist2d.Print_Weighted_Histogram_2D_With_Header();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram on rank %d\n",elapsed_secs, my_rank);


    }

  
  if (my_rank == 2)
    {


      // create 2d histogram
      begin = clock();
      v_name = "mstar_mendel";
      h_name = "UminusR";
      num_bins[0] = 25;
      spec[0] = 8.;
      spec[1] = 12.;
      bin_specs.push_back(spec);
      num_bins[1] = 27;
      spec[0] = 0.75;
      spec[1] = 3.45;
      bin_specs.push_back(spec);
      icols[0] = 0;
      icols[1] = 11;
      ibin = 0;
      Histogram2dNtuple this_hist2d(num_bins,bin_specs,base_name, v_name,h_name,icols,ibin, ntuple_filename);
      bin_specs.clear();
      //this_hist2d.Print_Histogram_2D_With_Header();
      this_hist2d.Print_Weighted_Histogram_2D_With_Header();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram on rank %d\n",elapsed_secs, my_rank);

    }

      
  MPI_Finalize();

}

