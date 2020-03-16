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
      num_bins[0] = 12;
      Histogram1dNtuple this_hist1d(num_bins[0],spec,base_name, h_name,icols[0],ntuple_filename);
      this_hist1d.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);

      /*

      // create 1d histogram
      begin = clock();
      spec[0] = 0.;
      spec[1] = 0.2;
      h_name = "z";
      icols[0] = 1;
      num_bins[0] = 40;
      this_hist1d.Init(num_bins[0],spec,base_name,h_name,icols[0],combined_ntuple_filename);
      this_hist1d.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);

      
      // create 1d histogram
      begin = clock();
      spec[0] = 5.;
      spec[1] = 8.;
      h_name = "log_mbh";
      icols[0] = 2;
      num_bins[0] = 30;
      this_hist1d.Init(num_bins[0],spec,base_name, h_name,icols[0],combined_ntuple_filename);
      this_hist1d.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);
      */

    }

  if (my_rank == 1)
    {

      /*
      // create 2d histogram
      begin = clock();
      v_name = "m_r";
      h_name = "z";
      num_bins[0] = 32;
      spec[0] = 10.;
      spec[1] = 26.;
      bin_specs.push_back(spec);
      num_bins[1] = 20;
      spec[0] = 0.;
      spec[1] = 0.2;
      bin_specs.push_back(spec);
      icols[0] = 0;
      icols[1] = 1;
      ibin = 0;
      Histogram2dNtuple this_hist2d(num_bins,bin_specs,base_name, v_name,h_name,icols,ibin, combined_ntuple_filename);
      bin_specs.clear();
      this_hist2d.Print_Histogram_2D_With_Header();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram on rank %d\n",elapsed_secs, my_rank);
      */     

    }

  
  if (my_rank == 2)
    {

      /*
      // create 2d histogram
      begin = clock();
      v_name = "m_r";
      h_name = "log_mbh";
      num_bins[0] = 32;
      spec[0] = 10.;
      spec[1] = 26.;
      bin_specs.push_back(spec);
      num_bins[1] = 31;
      spec[0] = 5.;
      spec[1] = 8.;
      bin_specs.push_back(spec);
      icols[0] = 0;
      icols[1] = 2;
      ibin = 0;
      Histogram2dNtuple this_hist2d(num_bins,bin_specs,base_name, v_name,h_name,icols,ibin, combined_ntuple_filename);
      bin_specs.clear();
      this_hist2d.Print_Histogram_2D_With_Header();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram on rank %d\n",elapsed_secs, my_rank);
      */
    }

      
  MPI_Finalize();

}

