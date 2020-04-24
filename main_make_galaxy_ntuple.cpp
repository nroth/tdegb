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
#include <hdf5.h>
#include <hdf5_hl.h>
#include "galaxy.h"
#include "histogram1d_ntuple.h"
#include "histogram2d_ntuple.h"
#include "integration_helper_functions.h"

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
  struct galaxy_catalogue_data gal_row;
  struct flare_data flare_row;

  string full_filename_string;
  int name_length;

  string filename_prefix_gals = "gal_catalogue_ntuple_";
  string extension_gals = ".dat";
  full_filename_string = filename_prefix_gals + std::to_string(my_rank);
  full_filename_string += extension_gals;
  name_length = full_filename_string.length();
  char* gal_ntuple_filename_array = new char[name_length];
  strcpy(gal_ntuple_filename_array, full_filename_string.c_str());
  gsl_ntuple *gal_ntuple  = gsl_ntuple_create(gal_ntuple_filename_array, &gal_row, sizeof (gal_row));
  delete [] gal_ntuple_filename_array;

  string filename_prefix_flares = "flare_ntuple_";
  string extension_flares = ".dat";
  full_filename_string = filename_prefix_flares + std::to_string(my_rank);
  full_filename_string += extension_flares;
  name_length = full_filename_string.length();
  char* flare_ntuple_filename_array = new char[name_length];
  strcpy(flare_ntuple_filename_array, full_filename_string.c_str());
  gsl_ntuple *flare_ntuple  = gsl_ntuple_create(flare_ntuple_filename_array, &flare_row, sizeof (flare_row));
  delete [] flare_ntuple_filename_array;

  // setup random num generator
  gsl_rng *rangen;
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  //gsl_rng_default_seed = (unsigned int)time(NULL);
  //gsl_rng_default_seed = 2;
  //gsl_rng_default_seed = my_rank;
  gsl_rng_default_seed = my_rank + (unsigned int)time(NULL);
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);


  string catalogue_filename = "/Users/nathanielroth/Dropbox/research/TDE/host_galaxies/sjoert_catalogue/van_velzen_Nov2019_catalogue.h5";
  hid_t catalogue_id = H5Fopen(catalogue_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  hsize_t num_galaxies;
  H5LTget_dataset_info(catalogue_id,"/z",&num_galaxies,NULL,NULL);
  //  int num_galaxies = catalogue_length;
  printf("number of galaxies is %llu\n", num_galaxies);
  H5Fclose (catalogue_id);

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
  
  
  int base_gals_per_proc = num_galaxies / n_procs; // integer division is automatically floored
  int remainder = num_galaxies % n_procs;

  int* displacements = new int[n_procs];
  int* num_gals_per_proc = new int[n_procs];

  int offset = 0;
  int remainder_count = 0;

  
  for (int i = 0; i < n_procs; ++i)

    {
      if (remainder_count < remainder)
	{
	  num_gals_per_proc[i] = base_gals_per_proc + 1;
	  remainder_count++;
	}
      else
	{
	  num_gals_per_proc[i] = base_gals_per_proc;
	}
      displacements[i] = offset;
      offset += num_gals_per_proc[i];
    }


  int my_num_gals = num_gals_per_proc[my_rank];

  

  if (my_rank == 0)
    {
      printf("num_galaxies %llu, n_procs %d, base_gals_per_proc %d, remainder %d\n",num_galaxies,n_procs, base_gals_per_proc, remainder);

      //read this in


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

  MPI_Scatterv(z_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, z, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(m_g_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, m_g, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(m_r_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, m_r, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(M_u_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, M_u, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(M_r_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, M_r, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mbh_bulge_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, mbh_bulge, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mbh_sigma_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, mbh_sigma, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(mass_mendel_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, mass_mendel, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(sersic_n_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, sersic_n, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(r50_kpc_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, r50_kpc, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Scatterv(ssfr_vector_big, num_gals_per_proc, displacements, MPI_DOUBLE, ssfr, my_num_gals, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

  delete[] displacements;

  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to prepare everything on rank %d\n",elapsed_secs,my_rank );

  Survey* surv = new Survey; // using a pointer here so as to be sure that constructor and destructor are only called once, outside the main integration function.
  // Before, when I just declared a survey and then passed it to Sample_Disruption_Parameters, I think funny things were happening with copy constructors, leading
  // to memory allocation errors.

  printf("starting to read catalogue on rank %d\n",my_rank);
  begin = clock();

  /////// MAIN LOOP
  for (int i = 0; i < my_num_gals; i++)
    {

      double galaxy_info[11] = { log10(mass_mendel[i]), mbh_sigma[i], mbh_bulge[i], z[i], sersic_n[i], r50_kpc[i], m_g[i], m_r[i],ssfr[i],M_u[i]- M_r[i],M_r[i]};

      Galaxy this_galaxy(galaxy_info);
      
      double vol_rate_weight;

      Sample_Disruption_Parameters(rangen,surv,this_galaxy,vol_rate_weight,flare_ntuple,&flare_row);
      
      for (int j = 0; j < 11; j++)
	{
	  gal_row.attributes[j] = galaxy_info[j];
	}
      gal_row.attributes[11] = this_galaxy.Get_nuker_gammaprime();
      gal_row.attributes[12] = this_galaxy.Get_median_A_V();
      gal_row.attributes[13] = this_galaxy.Get_sigma_A_V();
      gal_row.attributes[14] = this_galaxy.Get_Mu_Eff(this_galaxy.Get_m_g());
      gal_row.attributes[15] = this_galaxy.Get_Mu_Central(this_galaxy.Get_m_g());
      gal_row.attributes[16] = this_galaxy.Get_Mu_At_Rkpc_CosmoCorrected(0.1,this_galaxy.Get_m_g());
      gal_row.weight = vol_rate_weight;

      gsl_ntuple_write(gal_ntuple);

      if (i % 100000 == 0) printf("gal %d vol rate %e\n",i,vol_rate_weight);

    }

  gsl_ntuple_close (gal_ntuple);
  gsl_ntuple_close (flare_ntuple);

  gsl_rng_free(rangen); // we don't need any more random numbers

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
  printf("#\n# It took %f seconds to read %d galaxies from catalogue and write ntuple portion on rank %d\n",elapsed_secs,my_num_gals,my_rank );

  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much
  
  char combined_gal_ntuple_filename[35];
  string filename = "gal_catalogue_ntuple_combined.dat";
  strcpy(combined_gal_ntuple_filename, filename.c_str());

  char combined_flare_ntuple_filename[35];
  filename = "flare_ntuple_combined.dat";
  strcpy(combined_flare_ntuple_filename, filename.c_str());


  // combine the galaxy ntuples
  if ( my_rank == 0)
    {

      begin = clock();

      struct galaxy_catalogue_data combined_gal_row;

      char working_ntuple_filename[35];

      gsl_ntuple *working_ntuple;

      gsl_ntuple *combined_ntuple  = gsl_ntuple_create(combined_gal_ntuple_filename, &combined_gal_row, sizeof (combined_gal_row));
	
      // loop over files
      for (int i =0; i < n_procs; i++)
	{

	  //using the filename_prefix and extension that you defined earlier when you originally wrote the files
	sprintf(working_ntuple_filename, "%s%d%s", filename_prefix_gals.c_str(),i,extension_gals.c_str());
	working_ntuple = gsl_ntuple_open(working_ntuple_filename, &gal_row, sizeof (gal_row));
      
	while(gsl_ntuple_read(working_ntuple) != GSL_EOF)
	  {

	    memcpy(combined_gal_row.attributes, gal_row.attributes, sizeof(combined_gal_row.attributes));
	    combined_gal_row.weight = gal_row.weight;
	    
	    gsl_ntuple_write(combined_ntuple);
	  }


	gsl_ntuple_close (working_ntuple);
	
	}

      gsl_ntuple_close (combined_ntuple);

      //delete files with "remove?" http://www.cplusplus.com/reference/cstdio/remove/
      
      end = clock();
      float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to copy galaxy entries into single file\n",elapsed_secs);

    }


    // combine the flare ntuples
  if ( my_rank == 1)
    {

      begin = clock();

      struct flare_data combined_flare_row;

      char working_ntuple_filename[35];

      gsl_ntuple *working_ntuple;

      gsl_ntuple *combined_ntuple  = gsl_ntuple_create(combined_flare_ntuple_filename, &combined_flare_row, sizeof (combined_flare_row));
	
      // loop over files
      for (int i =0; i < n_procs; i++)
	{

	  //using the filename_prefix and extension that you defined earlier when you originally wrote the files
	sprintf(working_ntuple_filename, "%s%d%s", filename_prefix_flares.c_str(),i,extension_flares.c_str());
	working_ntuple = gsl_ntuple_open(working_ntuple_filename, &flare_row, sizeof (flare_row));


	while(gsl_ntuple_read(working_ntuple) != GSL_EOF)
	  {

	    memcpy(combined_flare_row.attributes, flare_row.attributes, sizeof(combined_flare_row.attributes));
	    combined_flare_row.weight = flare_row.weight;
	    
	    gsl_ntuple_write(combined_ntuple);
	  }


	gsl_ntuple_close (working_ntuple);
	
	}

      gsl_ntuple_close (combined_ntuple);

      //delete files with "remove?" http://www.cplusplus.com/reference/cstdio/remove/
      
      end = clock();
      float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to copy flare events into single file\n",elapsed_secs);

    }


  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much

  delete[] num_gals_per_proc;
  MPI_Finalize();

}
  
