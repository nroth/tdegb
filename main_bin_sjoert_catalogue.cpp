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
#include <ctime>
#include <mpi.h>
#include <gsl/gsl_rng.h>
#include "histogramNd.h"
#include "galaxy.h"
#include "survey.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include "integration_helper_functions.h"

using std::vector;
using std::string;

//#define SURVEY_SQ_DEG 14555.

int main(int argc, char **argv)
{

  clock_t begin = clock();

  //  double Omega = SURVEY_SQ_DEG / (41252.96);

  vector<vector<double> > bin_specs;
  vector<string> dimension_names;
  string base_name;
  vector<double> spec(3);

  // Define the axes for a histogram of the host galaxy properties
  
  string  name = "mstar_mendel";
  spec[0] = 8.;   // start, stop, delta
  spec[1] = 12.;
  spec[2] = 0.16;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  
  name = "u_minus_r";
  spec[0] = 0.75;   // start, stop, delta
  spec[1] = 3.5;
  spec[2] = 0.1;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  /*  
  name = "mbh_bulge";
  spec[0] = 3.;   // start, stop, delta
  spec[1] = 9.;
  spec[2] = 0.3;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);
  */

  name = "mbh_sigma";
  spec[0] = 3.;   // start, stop, delta
  spec[1] = 9.;
  spec[2] = 0.15;

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
  spec[1] = 6.;
  spec[2] = 0.25;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);
  */

  name = "sfr";
  spec[0] = -3.;   // start, stop, delta
  spec[1] = 3.;
  spec[2] = 0.25;

  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  base_name = "gals";
  HistogramNd hist_gals(bin_specs, dimension_names, base_name);

  int hist_gals_dimension = hist_gals.Get_Dimension();
  int n_bins_total = hist_gals.Get_Num_Bins_Total();

 
  ///// DISRUPTION HISTOGRAMS

  HistogramNd hist_vol_disrupt;
  hist_vol_disrupt = hist_gals;
  base_name = "vol_disrupt";
  hist_vol_disrupt.Set_Base_Name(base_name);

  HistogramNd hist_detected_disrupt;
  hist_detected_disrupt = hist_gals;
  base_name = "detected_disrupt";
  hist_detected_disrupt.Set_Base_Name(base_name);

  ///////// FLARE HISTOGRAMS
  bin_specs.clear();
  dimension_names.clear();
  
  name = "m_r";
  spec[0] = 12.;   // start, stop, delta
  spec[1] = 20.;
  spec[2] = 0.25;
  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  name = "g_minus_r";
  spec[0] = -1.;   // start, stop, delta
  spec[1] = 0.5;
  spec[2] = 0.05;
  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  name = "z";
  spec[0] = 0.;   // start, stop, delta
  spec[1] = 1.02;
  spec[2] = 0.05;
  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  name = "Lbol_peak";
  spec[0] = 43.;
  spec[1] = 47.;
  spec[2] = 0.1;
  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  name = "mbh"; // might be redundant, but might as well
  spec[0] = 3.;   // start, stop, delta
  spec[1] = 9.;
  spec[2] = 0.3;
  bin_specs.push_back(spec);
  dimension_names.push_back(name);

  base_name = "flares_detected";
  HistogramNd hist_detected_flares(bin_specs, dimension_names, base_name);

  ////////// Start reading in the catalogue


  int num_galaxies = 6092542; // should think more about how to make this flexible
  
  // need to define these everywhere but only want to allocate memory to them on rank 0, hence the use of malloc later for rank 0
  double* z_vector_big;
  double* m_g_vector_big; 
  double* m_r_vector_big;
  double* M_u_vector_big;
  double* M_r_vector_big; 
  double* mbh_bulge_vector_big;
  double* mbh_sigma_vector_big;
  double* mass_mendel_vector_big;
  double* sersic_n_vector_big;
  double* r50_kpc_vector_big;
  double* ssfr_vector_big;
  
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

      //      string catalogue_filename = "/Users/nathanielroth/Dropbox/research/TDE/host_galaxies/sjoert_catalogue/van_velzen_2018_catalogue.h5";
      //      string catalogue_filename = "/lustre/nroth1/TDE/host_galaxy_binning/van_velzen_2018_catalogue.h5";
      string catalogue_filename = "/lustre/nroth1/TDE/host_galaxy_binning/van_velzen_Nov2019_catalogue.h5";
      

      z_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      m_g_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      m_r_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      M_u_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      M_r_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      mbh_bulge_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      mbh_sigma_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      mass_mendel_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      sersic_n_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      r50_kpc_vector_big = (double*) malloc(num_galaxies * sizeof(double));
      ssfr_vector_big = (double*) malloc(num_galaxies * sizeof(double));

      hid_t file_id = H5Fopen(catalogue_filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
      
      H5LTread_dataset_double(file_id,"z",z_vector_big);
      H5LTread_dataset_double(file_id,"m_g",m_g_vector_big);
      H5LTread_dataset_double(file_id,"m_r",m_r_vector_big);
      H5LTread_dataset_double(file_id,"M_u",M_u_vector_big);
      H5LTread_dataset_double(file_id,"M_r",M_r_vector_big);
      H5LTread_dataset_double(file_id,"mbh_bulge",mbh_bulge_vector_big);
      H5LTread_dataset_double(file_id,"mbh_sigma",mbh_sigma_vector_big);
      H5LTread_dataset_double(file_id,"mass_mendel",mass_mendel_vector_big);
      H5LTread_dataset_double(file_id,"sersic_n",sersic_n_vector_big);
      H5LTread_dataset_double(file_id,"r50_kpc",r50_kpc_vector_big);
      H5LTread_dataset_double(file_id,"ssfr",ssfr_vector_big);

      H5Fclose(file_id);

    }


  double* z = new double[my_num_gals]; // for some reason, the hdf5 read will only work if you do the manual memory allocation this way (and when using new operator, you should delete later)
  double* m_g = new double[my_num_gals];
  double* m_r = new double[my_num_gals];
  double* M_u = new double[my_num_gals];
  double* M_r = new double[my_num_gals];
  double* mbh_bulge = new double[my_num_gals];
  double* mbh_sigma = new double[my_num_gals];
  double* mass_mendel = new double[my_num_gals];
  double* sersic_n = new double[my_num_gals];
  double* r50_kpc = new double[my_num_gals];
  double* ssfr = new double[my_num_gals];

  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much

  MPI_Scatterv(z_vector_big, send_counts, displacements, MPI_DOUBLE, z, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(m_g_vector_big, send_counts, displacements, MPI_DOUBLE, m_g, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(m_r_vector_big, send_counts, displacements, MPI_DOUBLE, m_r, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(M_u_vector_big, send_counts, displacements, MPI_DOUBLE, M_u, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(M_r_vector_big, send_counts, displacements, MPI_DOUBLE, M_r, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mbh_bulge_vector_big, send_counts, displacements, MPI_DOUBLE, mbh_bulge, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mbh_sigma_vector_big, send_counts, displacements, MPI_DOUBLE, mbh_sigma, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mass_mendel_vector_big, send_counts, displacements, MPI_DOUBLE, mass_mendel, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(sersic_n_vector_big, send_counts, displacements, MPI_DOUBLE, sersic_n, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(r50_kpc_vector_big, send_counts, displacements, MPI_DOUBLE, r50_kpc, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(ssfr_vector_big, send_counts, displacements, MPI_DOUBLE, ssfr, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
      free(z_vector_big);
      free(m_g_vector_big); 
      free(m_r_vector_big);
      free(M_u_vector_big); 
      free(M_r_vector_big); 
      free(mbh_bulge_vector_big);
      free(mbh_sigma_vector_big);
      free(mass_mendel_vector_big);
      free(sersic_n_vector_big);
      free(r50_kpc_vector_big);
      free(ssfr_vector_big);
    }

  delete[] send_counts;
  delete[] displacements;

  vector<double> catalogue_data(hist_gals_dimension);

  //Set up random number stuff
  gsl_rng *rangen;
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL);
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);

  Survey surv;

  clock_t end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to prepare everything on rank %d\n",elapsed_secs,my_rank );

  begin = clock();

  /////// MAIN LOOP
  for (int i = 0; i < my_num_gals; i++)
    {

      catalogue_data[0] = log10(mass_mendel[i]); // stored as value, we want log
      catalogue_data[1] = M_u[i] - M_r[i];
      catalogue_data[2] = mbh_sigma[i]; // stored as log
      catalogue_data[3] = z[i];
      //      catalogue_data[4] = sersic_n[i];
      catalogue_data[4] = log10(ssfr[i] * mass_mendel[i]);
					

      hist_gals.Count(catalogue_data);

      double galaxy_info[9] = { mass_mendel[i], mbh_sigma[i], mbh_bulge[i], z[i], sersic_n[i], r50_kpc[i], m_g[i], m_r[i],ssfr[i]};

      Galaxy this_galaxy(galaxy_info);

      double vol_rate_weight;
      double detected_rate_weight;

      Sample_Disruption_Parameters(rangen,surv,this_galaxy,vol_rate_weight, detected_rate_weight,hist_detected_flares);

      hist_vol_disrupt.Count(catalogue_data, vol_rate_weight);
      hist_detected_disrupt.Count(catalogue_data, detected_rate_weight);


    }

  gsl_rng_free(rangen); // we don't need more random numbers
  
  // we've finished with the input now, only need the histograms
  delete[] z;
  delete[] m_g;
  delete[] m_r;
  delete[] M_u;
  delete[] M_r;
  delete[] mbh_bulge;
  delete[] mbh_sigma;
  delete[] mass_mendel;
  delete[] sersic_n;
  delete[] r50_kpc;
  delete[] ssfr;

  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to finish the flares on rank %d\n",elapsed_secs,my_rank );

  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much
  
  // now reduce
  begin = clock();

  // make arrays for communication
  vector<double> src_mpi_vec(hist_gals.Get_All_Counts());
  vector<double> dst_mpi_vec(n_bins_total,0.);
  double* src_mpi = &src_mpi_vec[0];
  double* dst_mpi = new double[n_bins_total]; // remember to delete later
  for (int i =0; i < n_bins_total; ++i) dst_mpi[i] = 0.;

  MPI_Allreduce(src_mpi,dst_mpi,n_bins_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much
  dst_mpi_vec.assign(dst_mpi, dst_mpi + n_bins_total);
  hist_gals.Set_All_Counts(dst_mpi_vec);

  // reset
  src_mpi_vec = hist_vol_disrupt.Get_All_Counts();
  for (int i =0; i < n_bins_total; ++i) dst_mpi[i] = 0.;
  src_mpi = &src_mpi_vec[0]; // not needed?

  MPI_Allreduce(src_mpi,dst_mpi,n_bins_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much
  dst_mpi_vec.assign(dst_mpi, dst_mpi + n_bins_total);
  hist_vol_disrupt.Set_All_Counts(dst_mpi_vec);

  // reset
  src_mpi_vec = hist_detected_disrupt.Get_All_Counts();
  for (int i =0; i < n_bins_total; ++i) dst_mpi[i] = 0.;
  src_mpi = &src_mpi_vec[0]; // not needed?

  MPI_Allreduce(src_mpi,dst_mpi,n_bins_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much
  dst_mpi_vec.assign(dst_mpi, dst_mpi + n_bins_total);
  hist_detected_disrupt.Set_All_Counts(dst_mpi_vec);

  delete[] dst_mpi;

  // now do it for flares
  src_mpi_vec = hist_detected_flares.Get_All_Counts();
  dst_mpi_vec.clear(); // not needed?

  n_bins_total = hist_detected_flares.Get_Num_Bins_Total();
  dst_mpi_vec.resize(n_bins_total);

  double* src_mpi_flares = &src_mpi_vec[0];
  double* dst_mpi_flares = new double[n_bins_total]; // remember to delete later
  for (int i =0; i < n_bins_total; ++i) dst_mpi_flares[i] = 0.;

  MPI_Allreduce(src_mpi_flares,dst_mpi_flares,n_bins_total,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much
  dst_mpi_vec.assign(dst_mpi_flares, dst_mpi_flares + n_bins_total);
  hist_detected_flares.Set_All_Counts(dst_mpi_vec);

  delete[] dst_mpi_flares;
  
  end = clock();
  if (my_rank == 0)
    {
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to reduce everything\n",elapsed_secs);
    }

  vector<int> kept_axes(2);
  int ka[2];


  int assigned_ranks = 6;
  int rank_counter = 0;
  for (int i = 1; i < hist_gals.Get_Dimension(); i++)
    for (int j = 0; j < i; j++)
      {

	if (rank_counter % assigned_ranks == my_rank)
	  {
	    printf("Rank %d outputting galaxy info set %d %d\n", my_rank,i,j);
	    HistogramNd hist_projected_gals;
	    ka[0] = j;
	    ka[1] = i;
	    kept_axes.assign(ka, ka + 2);
	    hist_projected_gals = hist_gals.Create_Projected_Histogram(kept_axes);
	    hist_projected_gals.Print_Histogram_2D(0,1);
	  }
	rank_counter++;
	    
      }

  
  assigned_ranks = 6;
  int rank_displacement = 6;
  rank_counter = 0;
  for (int i = 1; i < hist_vol_disrupt.Get_Dimension(); i++)
	for (int j = 0; j < i; j++)
	  {

	    if (rank_counter % assigned_ranks + rank_displacement == my_rank)
	      {
		printf("Rank %d outputting vol disrupt set %d %d\n", my_rank,i,j);
		  HistogramNd hist_projected_vol_disrupt;
	    ka[0] = j;
	    ka[1] = i;
	    kept_axes.assign(ka, ka + 2);
	    hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
	    hist_projected_vol_disrupt.Print_Histogram_2D(0,1);
	      }

	    rank_counter++;
	  }

  assigned_ranks = 6;
  rank_displacement = 12;
  rank_counter = 0;

  for (int i = 1; i < hist_detected_disrupt.Get_Dimension(); i++)
    for (int j = 0; j < i; j++)
      {
	
	if (rank_counter % assigned_ranks + rank_displacement == my_rank)
	  {
	    printf("Rank %d outputting detected disrupt set %d %d\n", my_rank,i,j);
	    HistogramNd hist_projected_detected_disrupt;
	ka[0] = j;
	ka[1] = i;
	kept_axes.assign(ka, ka + 2);
	hist_projected_detected_disrupt = hist_detected_disrupt.Create_Projected_Histogram(kept_axes);
	hist_projected_detected_disrupt.Print_Histogram_2D(0,1);
	  }

	rank_counter++;
      }


  assigned_ranks = 1;
  rank_displacement = 18;
  rank_counter = 0;
  for (int i = 1; i < hist_detected_flares.Get_Dimension(); i++)
    for (int j = 0; j < i; j++)
      {
	
	if (rank_counter % assigned_ranks + rank_displacement == my_rank)
	  {
	    printf("Rank %d outputting flare set %d %d\n", my_rank,i,j);
	      HistogramNd hist_projected_flares;
	ka[0] = j;
	ka[1] = i;
	kept_axes.assign(ka, ka + 2);
	hist_projected_flares = hist_detected_flares.Create_Projected_Histogram(kept_axes);
	hist_projected_flares.Print_Histogram_2D(0,1);
	  }

	rank_counter++;
      }






  MPI_Finalize();

}
  
