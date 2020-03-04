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

using std::vector;
using std::string;

struct bin_params_2d
{
  gsl_histogram * hist1;
  gsl_histogram * hist2;
  int icol1; // which column in the ntuple this is tracking
  int icol2; // which column in the ntuple this is tracking
  int ibin;  // which (row?, column?) you're currently working on
};

int main(int argc, char **argv)
{

  MPI_Init( &argc, &argv );
  int my_rank,n_procs;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);

  vector<double> bin_specs;
  string base_name;
  string axis_name;
  int icol;
  int num_bins;

  //  void Print_Hist1d(char * , string, bin_params_1d);
  void Print_Hist2d_With_Header(char * , string, string, bin_params_2d);

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

  string v_name("default");
  string h_name("default");
  
  //  bin_params_1d binfo1d;
  bin_params_2d binfo2d;


  if (my_rank == 0)
    {
      // create 1d histogram
      begin = clock();
      bin_specs.push_back(10.);
      bin_specs.push_back(26.);
      base_name = "gals";
      axis_name = "m_r";
      icol = 0;
      num_bins = 50;
      Histogram1dNtuple hist_gals_m_r(num_bins,bin_specs,base_name, axis_name,icol,combined_ntuple_filename);
      bin_specs.clear();
      hist_gals_m_r.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);

      // create 1d histogram
      begin = clock();
      bin_specs.push_back(0.);
      bin_specs.push_back(0.2);
      base_name = "gals";
      axis_name = "z";
      icol = 1;
      num_bins = 40;
      Histogram1dNtuple hist_gals_z(num_bins,bin_specs,base_name, axis_name,icol,combined_ntuple_filename);
      bin_specs.clear();
      hist_gals_z.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);

      
      // create 1d histogram
      begin = clock();
      bin_specs.push_back(5.);
      bin_specs.push_back(8.);
      base_name = "gals";
      axis_name = "log_mbh";
      icol = 2;
      num_bins = 30;
      Histogram1dNtuple hist_gals_mbh(num_bins,bin_specs,base_name, axis_name,icol,combined_ntuple_filename);
      bin_specs.clear();
      hist_gals_mbh.Print_Histogram_1D();
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram on rank %d\n",elapsed_secs, my_rank);

    }

  if (my_rank == 1)
    {
      /*
      // create 2d histogram
      v_name = "m_r";
      h_name = "z";
      binfo2d.hist1 = hist_gals_m_r; 
      binfo2d.hist2 = hist_gals_z;
      binfo2d.icol1 = 0;
      binfo2d.icol2 = 1;
      binfo2d.ibin = 0;
      begin = clock();
      Print_Hist2d_With_Header(combined_ntuple_filename, v_name, h_name, binfo2d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram on rank %d\n",elapsed_secs, my_rank);
      */

      // create 2d histogram
      /*
      v_name = "z";
      h_name = "log_mbh";
      binfo2d.hist1 = hist_gals_z; 
      binfo2d.hist2 = hist_gals_mbh;
      binfo2d.icol1 = 1;
      binfo2d.icol2 = 2;
      binfo2d.ibin = 0;
      begin = clock();
      Print_Hist2d_With_Header(combined_ntuple_filename, v_name, h_name, binfo2d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram on rank %d\n",elapsed_secs, my_rank);
      */
    }

  /*
  if (my_rank == 2)
    {
      // create 2d histogram
      v_name = "m_r";
      h_name = "log_mbh";
      binfo2d.hist1 = hist_gals_m_r; 
      binfo2d.hist2 = hist_gals_mbh;
      binfo2d.icol1 = 0;
      binfo2d.icol2 = 2;
      binfo2d.ibin = 0;
      begin = clock();
      Print_Hist2d_With_Header(combined_ntuple_filename, v_name, h_name,binfo2d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram on rank %d\n",elapsed_secs,my_rank);
}
  */
      
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt much

  gsl_rng_free(rangen);

  MPI_Finalize();

}

/*
int sel_func_1d (void *this_data, void *p)
{
  return 1;
}
*/


int sel_func_2d (void *this_data, void *p)
{

  bin_params_2d binfo = *(struct bin_params_2d *)p;

  gsl_histogram * horizontal_hist = binfo.hist2;

  struct data * data_pointer = (struct data *) this_data;
  double this_horizontal_data = data_pointer->attributes[binfo.icol2];

  // handle extreme values
  if (this_horizontal_data >= gsl_histogram_max(horizontal_hist))
    {
      this_horizontal_data = (1. - 1.e-15) * gsl_histogram_max(horizontal_hist);
    }
  if (this_horizontal_data < gsl_histogram_min(horizontal_hist))
    {
      this_horizontal_data = gsl_histogram_min(horizontal_hist);
    }

  return ( this_horizontal_data >= horizontal_hist->range[binfo.ibin] && this_horizontal_data < horizontal_hist->range[binfo.ibin + 1]);
  //  return ( this_horizontal_data >= 0.1);
  // return true;
}


/*
double val_func_1d (void *this_data, void *p)
{
  bin_params_1d binfo = *(struct bin_params_1d *)p;

  int icol = binfo.icol;
  
  struct data * data_pointer = (struct data *) this_data;
  double this_col_value;

  this_col_value  = data_pointer->attributes[icol];

  // handle extreme values
  if (this_col_value < gsl_histogram_min(binfo.hist))
    {
      this_col_value = gsl_histogram_min(binfo.hist);
    }

  if (this_col_value >= gsl_histogram_max(binfo.hist))
    {
      this_col_value = (1. - 1.e-15) * gsl_histogram_max(binfo.hist);
    }


  return this_col_value;
}
*/

double val_func_2d (void *this_data, void *p)
{

  struct data * data_pointer = (struct data *) this_data;
  double this_col_data; //, v_data;

  bin_params_2d binfo = *(struct bin_params_2d *)p;

  int icol = binfo.icol1;
  this_col_data = data_pointer->attributes[icol];

  if (this_col_data < gsl_histogram_min(binfo.hist1))
    {
      this_col_data = gsl_histogram_min(binfo.hist1);
    }

  if (this_col_data >= gsl_histogram_max(binfo.hist1))
    {
      this_col_data = (1. - 1.e-15) * gsl_histogram_max(binfo.hist1);
    }


  return this_col_data;
}


//void Print_Hist2d_With_Header(string ntuple_filename, string v_name, string h_name, bin_params_2d binfo)
void Print_Hist2d_With_Header(char *ntuple_filename, string v_name, string h_name, bin_params_2d binfo)
{

  //write header
  string outfilename = "gsl_hist_2d_" + v_name + "_" + h_name + ".hist";
  FILE * outfile = fopen(outfilename.c_str(),"w");

  //  int n = ntuple_filename.length(); // doing this so that gsl_ntuple_create doesn't complain
  //  char ntuple_filename_array[n + 1]; // '' ''
  //  strcpy(ntuple_filename_array, ntuple_filename.c_str());  // '' ''

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;

  S.function = &sel_func_2d;
  V.function = &val_func_2d;

  struct data gal_row;
  
  // write header
  for (int iy = 0; iy < binfo.hist2->n + 1; iy++)
    {    
      fprintf(outfile,"%g ",binfo.hist2->range[iy]);
    }
  fprintf(outfile,"\n");
  for (int ix = 0; ix < binfo.hist1->n + 1; ix++)
    {    
      fprintf(outfile,"%g ",binfo.hist1->range[ix]);
    }
  fprintf(outfile,"\n");
  fclose(outfile);


  bin_params_2d bin_params_2_V = {binfo.hist1,binfo.hist2, binfo.icol1,binfo.icol2, 0}; // ibin of zero not used
  bin_params_2d bin_params_2_S = {binfo.hist1,binfo.hist2, binfo.icol1,binfo.icol2, 0}; // ibin  of 0 is just initialization

  V.params = &bin_params_2_V;
  
  for (int ibin = 0; ibin < (int)binfo.hist2->n; ibin++)
    {

      bin_params_2_S.ibin = ibin;
      S.params = &bin_params_2_S;

      //      gsl_ntuple *gal_ntuple = gsl_ntuple_open (ntuple_filename_array, &gal_row, sizeof (gal_row));
      gsl_ntuple *gal_ntuple = gsl_ntuple_open (ntuple_filename, &gal_row, sizeof (gal_row));
  
      gsl_histogram_reset(binfo.hist1);
      gsl_ntuple_project(binfo.hist1,gal_ntuple,&V,&S);

      outfile = fopen(outfilename.c_str(),"a"); 
      gsl_histogram_fprintf(outfile, binfo.hist1, "%-16.6e", "%-16.6e");
      fclose(outfile);

      gsl_ntuple_close(gal_ntuple);
      
    }

}

/*
//void Print_Hist1d(string ntuple_filename, string h_name, bin_params_1d binfo)
void Print_Hist1d(char *ntuple_filename, string h_name, bin_params_1d binfo)
{

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;

  S.function = &sel_func_1d;
  V.function = &val_func_1d;

  V.params = &binfo;

  //  int n = ntuple_filename.length(); // doing this so that gsl_ntuple_create doesn't complain
  //  char ntuple_filename_array[n + 1]; // '' ''
  //  strcpy(ntuple_filename_array, ntuple_filename.c_str());  // '' ''

  struct data gal_row;
  //  gsl_ntuple *gal_ntuple = gsl_ntuple_open (ntuple_filename_array, &gal_row, sizeof (gal_row));
  gsl_ntuple *gal_ntuple = gsl_ntuple_open (ntuple_filename, &gal_row, sizeof (gal_row));  
  gsl_ntuple_project(binfo.hist,gal_ntuple,&V,&S);
  gsl_ntuple_close(gal_ntuple);

  string outfilename = "gsl_hist_1d_" + h_name + ".hist";
  FILE * outfile = fopen(outfilename.c_str(),"w");
  gsl_histogram_fprintf(outfile, binfo.hist, "%-16.6e", "%-16.6e");
  fclose(outfile);

}
*/
