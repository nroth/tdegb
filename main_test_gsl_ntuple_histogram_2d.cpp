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

using std::vector;
using std::string;

struct ntuple_data
  {
    double attributes[3];
  };

struct bin_params_1d
{

  gsl_histogram * hist;
  int icol; // which column in the ntuple this is tracking
};


struct bin_params_2d
{
  gsl_histogram * hist1;
  gsl_histogram * hist2;
  int icol1; // which column in the ntuple this is tracking
  int icol2; // which column in the ntuple this is tracking
  int ibin;  // which (row?, column?) you're currently working on
};

//gsl_ntuple_select_fn S; ? 
//gsl_ntuple_value_fn V; ?


int main(int argc, char **argv)
{

  MPI_Init( &argc, &argv );
  int my_rank,n_procs;
  MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &n_procs);

  void Print_Hist1d(char * , string, bin_params_1d);
  void Print_Hist2d_With_Header(char * , string, string, bin_params_2d);


  int sel_func_2d (void *ntuple_data, void *params);
  double val_func_2d (void *ntuple_data, void *params);

  int sel_func_1d (void *ntuple_data, void *params);
  double val_func_1d (void *ntuple_data, void *params);


  // maybe clean this up?
  string filename_prefix = "gal_ntuple_";
  string extension = ".dat";
  char ntuple_filename_array[35];
  sprintf(ntuple_filename_array, "%s%d%s", filename_prefix.c_str(),my_rank,extension.c_str());

  double min_m_r = 10;
  double max_m_r = 26;
  int num_bins_m_r = 50;

  double mean_m_r = min_m_r + 0.5 * (max_m_r - min_m_r);
  double sigma_m_r = 0.2 * (max_m_r - min_m_r);

  gsl_histogram * hist_gals_m_r = gsl_histogram_alloc(num_bins_m_r);
  gsl_histogram_set_ranges_uniform(hist_gals_m_r, min_m_r, max_m_r);
    
  double min_z = 0;
  double max_z = 0.2;
  int num_bins_z = 40;

  double mean_z = min_z + 0.5 * (max_z - min_z);
  double sigma_z = 0.2 * (max_z - min_z);

  gsl_histogram * hist_gals_z = gsl_histogram_alloc(num_bins_z);
  gsl_histogram_set_ranges_uniform(hist_gals_z, min_z, max_z);

  double min_log_mbh = 5.;
  double max_log_mbh = 8.;
  int num_bins_mbh = 30;

  double mean_log_mbh = min_log_mbh + 0.5 * (max_log_mbh - min_log_mbh);
  double sigma_log_mbh = 0.2 * (max_log_mbh - min_log_mbh);

  gsl_histogram * hist_gals_mbh = gsl_histogram_alloc(num_bins_mbh);
  gsl_histogram_set_ranges_uniform(hist_gals_mbh, min_log_mbh, max_log_mbh);
  

  // setup random num generator
  gsl_rng *rangen;
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL);
  //gsl_rng_default_seed = 2;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);


  int num_samples = 1000000;

  clock_t begin;
  clock_t end;
  
  struct ntuple_data gal_row;
  
  gsl_ntuple *gal_ntuple  = gsl_ntuple_create(ntuple_filename_array, &gal_row, sizeof (gal_row));


  begin = clock();
  for (int t = 0; t < num_samples; t++)
    {

      gal_row.attributes[0] = mean_m_r + gsl_ran_gaussian(rangen, sigma_m_r);
      gal_row.attributes[1] = mean_z + gsl_ran_gaussian(rangen, sigma_z);
      gal_row.attributes[2] = mean_log_mbh + gsl_ran_gaussian(rangen, sigma_log_mbh);;

      gsl_ntuple_write(gal_ntuple);
    }

  
  gsl_ntuple_close (gal_ntuple);

  end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to generate %d samples on rank %d\n",elapsed_secs, num_samples, my_rank);
  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt muchw


  if ( my_rank == 0)
    {

      string filename = "gal_ntuple_combined.dat";
      char combined_ntuple_filename[35];
      strcpy(combined_ntuple_filename, filename.c_str());  

      // combine
      //      for (int i = 0 ; i < n_procs; i++) 

      std::ifstream if_a("gal_ntuple_0.dat", std::ios_base::binary);
      std::ifstream if_b("gal_ntuple_1.dat", std::ios_base::binary);
      std::ofstream of_c(combined_ntuple_filename, std::ios_base::binary);

      of_c << if_a.rdbuf() << if_b.rdbuf();


      string v_name("default");
      string h_name("default");


      // create 1d histogram
      begin = clock();
      h_name = "m_r";
      bin_params_1d binfo1d = {hist_gals_m_r, 0};
      Print_Hist1d(combined_ntuple_filename,h_name,binfo1d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

      // create 1d histogram
      begin = clock();
      h_name = "z";
      binfo1d.hist = hist_gals_z;
      binfo1d.icol = 1;
      Print_Hist1d(combined_ntuple_filename,h_name,binfo1d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

      // create 1d histogram
      begin = clock();
      h_name = "log_mbh";
      binfo1d.hist = hist_gals_mbh;
      binfo1d.icol = 2;
      Print_Hist1d(combined_ntuple_filename,h_name,binfo1d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

      // create 2d histogram
      v_name = "m_r";
      h_name = "z";
      bin_params_2d binfo2d = {hist_gals_m_r, hist_gals_z,0,1,0};
      begin = clock();
      Print_Hist2d_With_Header(combined_ntuple_filename, v_name, h_name, binfo2d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);

      // create 2d histogram
      v_name = "m_r";
      h_name = "log_mbh";
      binfo2d = {hist_gals_m_r,hist_gals_mbh,0,2,0};
      begin = clock();
      Print_Hist2d_With_Header(combined_ntuple_filename, v_name, h_name,binfo2d);
      end = clock();
      elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
      printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);
      
    }

  MPI_Barrier(MPI_COMM_WORLD); // might not be necessary, but doesn't hurt muchw
  
  gsl_histogram_free(hist_gals_m_r);
  gsl_histogram_free(hist_gals_z);
  gsl_histogram_free(hist_gals_mbh);
  
  gsl_rng_free(rangen);

  MPI_Finalize();

}

int sel_func_1d (void *this_data, void *p)
{
  return 1;
}


int sel_func_2d (void *this_data, void *p)
{

  bin_params_2d binfo = *(struct bin_params_2d *)p;

  gsl_histogram * horizontal_hist = binfo.hist2;

  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double this_horizontal_data = data->attributes[binfo.icol2];

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


double val_func_1d (void *this_data, void *p)
{
  bin_params_1d binfo = *(struct bin_params_1d *)p;

  int icol = binfo.icol;
  
  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double this_col_value;

  this_col_value  = data->attributes[icol];

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


double val_func_2d (void *this_data, void *p)
{

  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double this_col_data; //, v_data;

  bin_params_2d binfo = *(struct bin_params_2d *)p;

  int icol = binfo.icol1;
  this_col_data = data->attributes[icol];

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

  struct ntuple_data gal_row;
  
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

  struct ntuple_data gal_row;
  //  gsl_ntuple *gal_ntuple = gsl_ntuple_open (ntuple_filename_array, &gal_row, sizeof (gal_row));
  gsl_ntuple *gal_ntuple = gsl_ntuple_open (ntuple_filename, &gal_row, sizeof (gal_row));  
  gsl_ntuple_project(binfo.hist,gal_ntuple,&V,&S);
  gsl_ntuple_close(gal_ntuple);

  string outfilename = "gsl_hist_1d_" + h_name + ".hist";
  FILE * outfile = fopen(outfilename.c_str(),"w");
  gsl_histogram_fprintf(outfile, binfo.hist, "%-16.6e", "%-16.6e");
  fclose(outfile);

}



  
