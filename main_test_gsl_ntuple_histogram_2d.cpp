#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_ntuple.h>

using std::vector;
using std::string;

struct ntuple_data
  {

    // m_r, z, mbh, mstar
    double attributes[3];

    /*
    double m_r;
    double z;
    double mbh;
    double mstar;
    */
    
  };


int main(int argc, char **argv)
{

  int Bin_Index_Hist2d(int, int,size_t,size_t); // not currently used?
  void Print_Hist2d_With_Header(string, gsl_histogram *, gsl_histogram *, string, string, int, int);
  void Print_Hist1d(gsl_histogram *, string);

  int sel_func_2d (void *ntuple_data, void *params);
  double val_func_2d (void *ntuple_data, void *params);

  int sel_func_1d (void *ntuple_data, void *params);
  double val_func_1d (void *ntuple_data, void *params);


  string gal_ntuple_filename = "test_ntuple2d.dat";
  int n = gal_ntuple_filename.length(); // doing this so that gsl_ntuple_create doesn't complain
  char filename_array[n + 1]; // '' ''
  strcpy(filename_array, gal_ntuple_filename.c_str());  // '' ''

  struct ntuple_data gal_row;

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;


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
  //  gsl_rng_default_seed = (unsigned int)time(NULL);
  gsl_rng_default_seed = 2;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);


  int num_samples = 1000000;

  clock_t begin;
  clock_t end;
  

  gsl_ntuple *gal_ntuple  = gsl_ntuple_create(filename_array, &gal_row, sizeof (gal_row));


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
  printf("#\n# It took %f seconds to generate %d samples\n",elapsed_secs, num_samples);

  struct bin_param_1d {gsl_histogram * hist; int icol;};


  string v_name("default");
  string h_name("default");

  // create 1d histogram
  begin = clock();

  S.function = &sel_func_1d;
  V.function = &val_func_1d;

  bin_param_1d bin_params = {hist_gals_m_r, 0};
  V.params = &bin_params;

  gal_ntuple = gsl_ntuple_open (filename_array, &gal_row, sizeof (gal_row));
  gsl_ntuple_project(hist_gals_m_r,gal_ntuple,&V,&S);
  gsl_ntuple_close(gal_ntuple);

  h_name = "m_r";
  string outfilename = "gsl_hist_1d_" + h_name + ".hist";
  Print_Hist1d(hist_gals_m_r, outfilename);

  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  // create 1d histogram
  begin = clock();

  //  S.function = &sel_func_1d;
  //  V.function = &val_func_1d;

  bin_params.hist = hist_gals_z;
  bin_params.icol = 1;
  V.params = &bin_params;

  gal_ntuple = gsl_ntuple_open(filename_array, &gal_row, sizeof (gal_row));
  gsl_ntuple_project(hist_gals_z,gal_ntuple,&V,&S);
  gsl_ntuple_close(gal_ntuple);

  h_name = "z";
  outfilename = "gsl_hist_1d_" + h_name + ".hist";

  Print_Hist1d(hist_gals_z, outfilename);

  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  // create 1d histogram
  begin = clock();

  //  S.function = &sel_func_1d;
  //  V.function = &val_func_1d;

  bin_params.hist = hist_gals_mbh;
  bin_params.icol = 2;
  V.params = &bin_params;

  gal_ntuple = gsl_ntuple_open(filename_array, &gal_row, sizeof (gal_row));
  gsl_ntuple_project(hist_gals_mbh,gal_ntuple,&V,&S);
  gsl_ntuple_close(gal_ntuple);

  h_name = "log_mbh";
  outfilename = "gsl_hist_1d_" + h_name + ".hist";

  Print_Hist1d(hist_gals_mbh, outfilename);

  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  // create 2d histogram
  v_name = "m_r";
  h_name = "z";
  int v_icol = 0; // do this better
  int h_icol = 1;
  begin = clock();
  Print_Hist2d_With_Header(gal_ntuple_filename, hist_gals_m_r, hist_gals_z, v_name, h_name, v_icol,h_icol);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);

  // create 2d histogram
  v_name = "m_r";
  h_name = "log_mbh";
  v_icol = 0; // do this better
  h_icol = 2;
  begin = clock();
  Print_Hist2d_With_Header(gal_ntuple_filename, hist_gals_m_r, hist_gals_mbh, v_name, h_name,v_icol,h_icol);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 2d histogram\n",elapsed_secs);

  gsl_histogram_free(hist_gals_m_r);
  gsl_histogram_free(hist_gals_z);
  gsl_histogram_free(hist_gals_mbh);
  
  gsl_rng_free(rangen);


}

int sel_func_1d (void *this_data, void *p)
{
  return 1;
}


int sel_func_2d (void *this_data, void *p)
{

  struct bin_struct { gsl_histogram * hist; gsl_histogram * hist_h; int icol; int icol_h; int ibin_h;};
  bin_struct binfo = *(struct bin_struct *)p;

  gsl_histogram * horizontal_hist = binfo.hist_h;

  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double this_horizontal_data = data->attributes[binfo.icol_h];

  // handle extreme values
  if (this_horizontal_data >= gsl_histogram_max(horizontal_hist))
    {
      this_horizontal_data = (1. - 1.e-15) * gsl_histogram_max(horizontal_hist);
    }
  if (this_horizontal_data < gsl_histogram_min(horizontal_hist))
    {
      this_horizontal_data = gsl_histogram_min(horizontal_hist);
    }

  return ( this_horizontal_data >= horizontal_hist->range[binfo.ibin_h] && this_horizontal_data < horizontal_hist->range[binfo.ibin_h + 1]);
  //  return ( this_horizontal_data >= 0.1);
  // return true;
}


double val_func_1d (void *this_data, void *p)
{
  struct bin_struct { gsl_histogram * hist; int icol; };
  bin_struct binfo = *(struct bin_struct *)p;

  int icol = binfo.icol;
  
  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double this_col_value;

  this_col_value  = data->attributes[icol];

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
  //  (void)(params); /* avoid unused parameter warning */

  struct bin_struct { gsl_histogram * hist; gsl_histogram * hist_h; int icol; int icol_v; int ibin_v;};
  
  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double this_col_data; //, v_data;

  bin_struct ph = *(struct bin_struct *)p;

  int icol = ph.icol;
  this_col_data = data->attributes[icol];
  //  v_data = data->attributes[ibin_v];
  //  mbh = data->mbh;
  //  mstar = data->mstar;

  // replace z with generic name
  if (this_col_data < gsl_histogram_min(ph.hist))
    {
      this_col_data = gsl_histogram_min(ph.hist);
    }

  if (this_col_data >= gsl_histogram_max(ph.hist))
    {
      this_col_data = (1. - 1.e-15) * gsl_histogram_max(ph.hist);
    }


  return this_col_data;
}



int Bin_Index_Hist2d(int i, int j,size_t nx,size_t ny)
{
  return i * ny + j;

}

void Print_Hist2d_With_Header(string ntuple_filename, gsl_histogram * hist1, gsl_histogram * hist2, string v_name, string h_name, int i1, int i2)
{

  //write header
  string outfilename = "gsl_hist_2d_" + v_name + "_" + h_name + ".hist";
  FILE * outfile = fopen(outfilename.c_str(),"w");

  int n = ntuple_filename.length(); // doing this so that gsl_ntuple_create doesn't complain
  char ntuple_filename_array[n + 1]; // '' ''
  strcpy(ntuple_filename_array, ntuple_filename.c_str());  // '' ''

  gsl_ntuple_select_fn S;
  gsl_ntuple_value_fn V;

  struct bin_param_2d {gsl_histogram * hist; gsl_histogram * hist_h; int icol; int icol_v; int ibin;};

  S.function = &sel_func_2d;
  V.function = &val_func_2d;

  struct ntuple_data gal_row;
  
  // write header
  for (int iy = 0; iy < hist2->n + 1; iy++)
    {    
      fprintf(outfile,"%g ",hist2->range[iy]);
    }
  fprintf(outfile,"\n");
  for (int ix = 0; ix < hist1->n + 1; ix++)
    {    
      fprintf(outfile,"%g ",hist1->range[ix]);
    }
  fprintf(outfile,"\n");
  fclose(outfile);


  //clean this up, pass these in?
  bin_param_2d bin_params_2_V = {hist1,hist2, i1,i2, 0};
  bin_param_2d bin_params_2_S = {hist1,hist2, i1,i2, 0}; // just starting definition

  V.params = &bin_params_2_V;
  
  for (int i2 = 0; i2 < (int)hist2->n; i2++)
    {

      bin_params_2_S.ibin = i2;
      S.params = &bin_params_2_S;

      gsl_ntuple *gal_ntuple = gsl_ntuple_open (ntuple_filename_array, &gal_row, sizeof (gal_row));
  
      gsl_histogram_reset(hist1);
      gsl_ntuple_project(hist1,gal_ntuple,&V,&S);

      outfile = fopen(outfilename.c_str(),"a"); 
      gsl_histogram_fprintf(outfile, hist1, "%-16.6e", "%-16.6e");
      fclose(outfile);

      gsl_ntuple_close(gal_ntuple);
      
    }


}

void Print_Hist1d(gsl_histogram * hist1d, string name)
{

  FILE * outfile = fopen(name.c_str(),"w"); 

  gsl_histogram_fprintf(outfile, hist1d, "%-16.6e", "%-16.6e");

  fclose(outfile);

}



  
