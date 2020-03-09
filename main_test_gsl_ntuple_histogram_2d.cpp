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

  vector<int> num_bins(2);
  vector<int> icols(2);
  vector<double> spec(2);
  int ibin;
  vector<vector<double>> bin_specs;
  string base_name;
  string v_name("default"); 
  string h_name("default");

  // maybe clean this up?
  string filename_prefix = "gal_ntuple_";
  string extension = ".dat";
  char ntuple_filename_array[35];
  sprintf(ntuple_filename_array, "%s%d%s", filename_prefix.c_str(),my_rank,extension.c_str());


  double mean_m_r = 10. + 0.5 * (26. - 10.);
  double sigma_m_r = 0.2 * (26. - 10.);

  double mean_z = 0. + 0.5 * (0.2 - 0.);
  double sigma_z = 0.2 * (0.2 - 0.);

  double mean_log_mbh = 5. + 0.5 * (8. - 5.);
  double sigma_log_mbh = 0.2 * (8. - 5.);

  // setup random num generator
  gsl_rng *rangen;
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  //gsl_rng_default_seed = (unsigned int)time(NULL);
  //gsl_rng_default_seed = 2;
  gsl_rng_default_seed = my_rank;
  //gsl_rng_default_seed = my_rank + (unsigned int)time(NULL);
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);


  int num_samples = 1000000;

  clock_t begin;
  clock_t end;
  
  struct data gal_row;
  
  gsl_ntuple *gal_ntuple  = gsl_ntuple_create(ntuple_filename_array, &gal_row, sizeof (gal_row));


  begin = clock();
  for (int t = 0; t < num_samples; t++)
    {

      gal_row.attributes[0] = mean_m_r + gsl_ran_gaussian(rangen, sigma_m_r);
      gal_row.attributes[1] = mean_z + gsl_ran_gaussian(rangen, sigma_z);
      gal_row.attributes[2] = mean_log_mbh + gsl_ran_gaussian(rangen, sigma_log_mbh);

      gsl_ntuple_write(gal_ntuple);
    }

  
  gsl_ntuple_close (gal_ntuple);

  end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to generate %d samples on rank %d\n",elapsed_secs, num_samples, my_rank);
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much


  char combined_ntuple_filename[35];
  string filename = "gal_ntuple_combined.dat";
  strcpy(combined_ntuple_filename, filename.c_str());

  

  // combine the ntuples
  if ( my_rank == 0)
    {

      begin = clock();

      struct data combined_gal_row;

      char working_ntuple_filename[35];
      filename_prefix = "gal_ntuple_";
      extension = ".dat";
      gsl_ntuple *working_ntuple;

      gsl_ntuple *combined_ntuple  = gsl_ntuple_create(combined_ntuple_filename, &combined_gal_row, sizeof (combined_gal_row));
	
      // loop over files
      for (int i =0; i < n_procs; i++)
	{
	  
	sprintf(working_ntuple_filename, "%s%d%s", filename_prefix.c_str(),i,extension.c_str());
	working_ntuple = gsl_ntuple_open(working_ntuple_filename, &gal_row, sizeof (gal_row));
      
	for (int t = 0; t < num_samples; t++)
	  {

	    gsl_ntuple_read(working_ntuple);
	    memcpy(combined_gal_row.attributes, gal_row.attributes, sizeof(combined_gal_row.attributes));
	    gsl_ntuple_write(combined_ntuple);
	  
	  }
      
	gsl_ntuple_close (working_ntuple);
	
	}

      gsl_ntuple_close (combined_ntuple);

      //delete files with "remove?" http://www.cplusplus.com/reference/cstdio/remove/
      
      end = clock();
      float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to copy events into single file\n",elapsed_secs);

    }

  
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much

  
  // can this be split up among processes?

  base_name = "gals";

  if (my_rank == 0)
    {
      // create 1d histogram
      begin = clock();
      spec[0] = 10.;
      spec[1] = 26.;
      h_name = "m_r";
      icols[0] = 0;
      num_bins[0] = 50;
      Histogram1dNtuple this_hist1d(num_bins[0],spec,base_name, h_name,icols[0],combined_ntuple_filename);
      this_hist1d.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);

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

    }

  if (my_rank == 1)
    {

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
     

    }

  
  if (my_rank == 2)
    {

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

    }

      
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much

  gsl_rng_free(rangen);

  MPI_Finalize();

}

