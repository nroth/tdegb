#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "physical_constants.h"
#include "histogramNd.h"
#include "cdf.h"
#include "galaxy.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "integration_helper_functions.h"

using std::vector;
using std::string;

//#define SURVEY_SQ_DEG 14555.

int main(int argc, char **argv)
{

  //  double Omega = SURVEY_SQ_DEG / (41252.96);

  vector<vector<double> > bin_specs;
  vector<string> dimension_names;
  string base_name;
  vector<double> spec(3);

  // Define the axes for a histogram of the host galaxy properties
  
  string  name = "mstar";
  spec[0] = 8.;   // start, stop, delta
  spec[1] = 12.;
  spec[2] = 0.16;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  
  name = "g_minus_r";
  spec[0] = 0.2;   // start, stop, delta
  spec[1] = 2.;
  spec[2] = 0.08;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  
  name = "mbh_bulge";
  spec[0] = 3.;   // start, stop, delta
  spec[1] = 9.;
  spec[2] = 0.3;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  name = "mbh_sigma";
  spec[0] = 3.;   // start, stop, delta
  spec[1] = 9.;
  spec[2] = 0.3;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  name = "z";
  spec[0] = 0.;   // start, stop, delta
  spec[1] = 1.02;
  spec[2] = 0.05;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  /*
  name = "sersic";
  spec[0] = 0.;   // start, stop, delta
  spec[1] = 5.;
  spec[2] = 0.25;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);
  */

  base_name = "gals";
  HistogramNd hist_gals(bin_specs, dimension_names, base_name);

  int hist_gals_dimension = hist_gals.Get_Dimension();
  int n_bins_total = hist_gals.Get_Num_Bins_Total();

  //  printf("Galaxy histogram dimension is %d\n", hist_gals_dimension);
  //  printf("Total number of bins in galaxy histogram is %d\n", n_bins_total);

  // Define the axes for a histogram of the volumetric disruption properties. For now, give it the same axes as for the host galaxies, but later will need to give it other axes like flare luminosity, flare temperature, maybe redshift

  HistogramNd hist_vol_disrupt;
  hist_vol_disrupt = hist_gals;
  base_name = "vol_disrupt";
  hist_vol_disrupt.Set_Base_Name(base_name);

  HistogramNd hist_detected_disrupt;
  hist_detected_disrupt = hist_gals;
  base_name = "detected_disrupt";
  hist_detected_disrupt.Set_Base_Name(base_name);

  // And maybe later you'll make another histogram for the *observed* flares

  // Start reading in the catalogue

  clock_t begin = clock();

  int num_galaxies = 6101944; // should think more about how to make this flexible
  
  // need to define these everywhere but only want to allocate memory to them on rank 0, hence the use of malloc later for rank 0
  double* z_vector_big;
  double* m_g_vector_big; 
  double* m_r_vector_big; 
  double* mbh_bulge_vector_big;
  double* mbh_sigma_vector_big;
  double* mass_vector_big;
  double* sersic_n_vector_big;
  double* r50_kpc_vector_big;
  
  MPI_Init( &argc, &argv );

  int my_rank,n_procs;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);

  int base_gals_per_proc = num_galaxies / n_procs; // integer division is automatically floored
  int remainder = num_galaxies % n_procs;

  int* displacements = new int[n_procs];
  int* send_counts = new int[n_procs];

  int offset = 0;
  int remainder_count = 0;

  
  for (int i = 0; i < n_procs; ++i)

    {
      if (remainder_count < remainder)
	{
	  send_counts[i] = base_gals_per_proc + 1;
	  remainder_count++;
	}
      else
	{
	  send_counts[i] = base_gals_per_proc;
	}
      displacements[i] = offset;
      offset += send_counts[i];
    }

  int my_num_gals = send_counts[my_rank];

  if (my_rank == 0)
    {
      printf("num_galaxies %d, n_procs %d, base_gals_per_proc %d, remainder %d\n",num_galaxies,n_procs, base_gals_per_proc, remainder);

      string catalogue_filename = "/Users/nathanielroth/Dropbox/research/TDE/host_galaxies/sjoert_catalogue/van_velzen_2018_catalogue.h5";

      z_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      m_g_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      m_r_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      mbh_bulge_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      mbh_sigma_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      mass_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      sersic_n_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      r50_kpc_vector_big = (double*) malloc(num_galaxies * sizeof(double));

      hid_t file_id = H5Fopen(catalogue_filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
      
      H5LTread_dataset_double(file_id,"z",z_vector_big);
      H5LTread_dataset_double(file_id,"m_g",m_g_vector_big);
      H5LTread_dataset_double(file_id,"m_r",m_r_vector_big);
      H5LTread_dataset_double(file_id,"mbh_bulge",mbh_bulge_vector_big);
      H5LTread_dataset_double(file_id,"mbh_sigma",mbh_sigma_vector_big);
      H5LTread_dataset_double(file_id,"mass",mass_vector_big);
      H5LTread_dataset_double(file_id,"sersic_n",sersic_n_vector_big);
      H5LTread_dataset_double(file_id,"r50_kpc",r50_kpc_vector_big);

      H5Fclose(file_id);

    }


  double* z = new double[my_num_gals]; // for some reason, the hdf5 read will only work if you do the manual memory allocation this way (and when using new operator, you should elete later)
  double* m_g = new double[my_num_gals];
  double* m_r = new double[my_num_gals];
  double* mbh_bulge = new double[my_num_gals];
  double* mbh_sigma = new double[my_num_gals];
  double* mass = new double[my_num_gals];
  double* sersic_n = new double[my_num_gals];
  double* r50_kpc = new double[my_num_gals];

  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much

  MPI_Scatterv(z_vector_big, send_counts, displacements, MPI_DOUBLE, z, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(m_g_vector_big, send_counts, displacements, MPI_DOUBLE, m_g, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(m_r_vector_big, send_counts, displacements, MPI_DOUBLE, m_r, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mbh_bulge_vector_big, send_counts, displacements, MPI_DOUBLE, mbh_bulge, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mbh_sigma_vector_big, send_counts, displacements, MPI_DOUBLE, mbh_sigma, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mass_vector_big, send_counts, displacements, MPI_DOUBLE, mass, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(sersic_n_vector_big, send_counts, displacements, MPI_DOUBLE, sersic_n, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(r50_kpc_vector_big, send_counts, displacements, MPI_DOUBLE, r50_kpc, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
      free(z_vector_big);
    }

  delete[] send_counts;
  delete[] displacements;

  vector<double> catalogue_data(hist_gals_dimension);

  for (int i = 0; i < my_num_gals; i++)
    {

      catalogue_data[0] = log10(mass[i]); // stored as value, we want log
      catalogue_data[1] = m_g[i] - m_r[i];
      catalogue_data[2] = mbh_bulge[i];// stored as log
      catalogue_data[3] = mbh_sigma[i]; // stored as log
      catalogue_data[4] = z[i];

      hist_gals.Count(catalogue_data);

      double galaxy_info[8] = { mass[i], mbh_sigma[i], mbh_bulge[i], z[i], sersic_n[i], r50_kpc[i], m_g[i], m_r[i]};


      Galaxy this_galaxy(galaxy_info);

      double total_disruption_rate_volumetric = Total_Disruption_Rate(this_galaxy);

      hist_vol_disrupt.Count(catalogue_data,total_disruption_rate_volumetric); // the volumetric disruption rate histogram is just like the host galaxy histogram, but weighted by per-galaxy disruption rate. The z is needed to convert from galaxy time frame to observer rest framee


      double T = 3.e4;
      double beta = 1.;

      double total_rate_obs_rband = Total_Disruption_Rate_Observed_Rband(this_galaxy,T,beta);

      hist_detected_disrupt.Count(catalogue_data,total_rate_obs_rband);

    }


  // we're done with the input now, only need the histograms
  delete[] z;
  delete[] m_g;
  delete[] m_r;
  delete[] mbh_bulge;
  delete[] mbh_sigma;
  delete[] mass;
  delete[] sersic_n;
  delete[] r50_kpc;

  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much

  // now reduce

  // make arrays for communication
  vector<double> src_mpi_vec(hist_gals.Get_All_Counts());
  vector<double> dst_mpi_vec(n_bins_total,0.);
  double* src_mpi = &src_mpi_vec[0];
  double* dst_mpi = new double[n_bins_total]; // remember to delete later
  for (int i =0; i < n_bins_total; ++i) dst_mpi[i] = 0.;

  MPI_Allreduce(src_mpi,dst_mpi,n_bins_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  dst_mpi_vec.assign(dst_mpi, dst_mpi + n_bins_total);
  hist_gals.Set_All_Counts(dst_mpi_vec);

  // reset
  src_mpi_vec = hist_vol_disrupt.Get_All_Counts();
  for (int i =0; i < n_bins_total; ++i) dst_mpi[i] = 0.;
  src_mpi = &src_mpi_vec[0]; // not needed?

  MPI_Allreduce(src_mpi,dst_mpi,n_bins_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  dst_mpi_vec.assign(dst_mpi, dst_mpi + n_bins_total);
  hist_vol_disrupt.Set_All_Counts(dst_mpi_vec);

  // reset
  src_mpi_vec = hist_detected_disrupt.Get_All_Counts();
  for (int i =0; i < n_bins_total; ++i) dst_mpi[i] = 0.;
  src_mpi = &src_mpi_vec[0]; // not needed?

  MPI_Allreduce(src_mpi,dst_mpi,n_bins_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  dst_mpi_vec.assign(dst_mpi, dst_mpi + n_bins_total);
  hist_detected_disrupt.Set_All_Counts(dst_mpi_vec);

  delete[] dst_mpi;

  clock_t end = clock();

  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to bin %d entries\n",elapsed_secs, num_galaxies);

  vector<int> kept_axes(2);
  int ka[2];

  if (my_rank == 0)
    {
  
  // Print out some histograms of galaxy properties
  HistogramNd hist_projected_gals;
  

    
  ka[0] = 0;
  ka[1] = 1;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 2;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 2;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);
  
  ka[0] = 3;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

    }
  /*

  ka[0] = 0;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 3;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  ka[0] = 4;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
  hist_projected_gals.Print_Histogram_2D(0,1);

  */

  // Print out some histograms of volumetric disruption properties


  if (my_rank == 1)
    {
  HistogramNd hist_projected_vol_disrupt;
  
  ka[0] = 0;
  ka[1] = 1;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 2;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 2;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 3;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);
    }

  /*

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 3;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 4;
  ka[1] = 5;
  kept_axes.assign(ka, ka + 2);

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);


  */

  /////////

  // Print out some histograms of detectable disruption

  if (my_rank == 2)
    {
  HistogramNd hist_projected_detected_disrupt;
  
  ka[0] = 0;
  ka[1] = 1;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 2;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 2;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 3;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 0;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 1;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 2;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

  ka[0] = 3;
  ka[1] = 4;
  kept_axes.assign(ka, ka + 2);

  hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_detected_disrupt.Print_Histogram_2D(0,1);

    }

    MPI_Finalize();

}
  
