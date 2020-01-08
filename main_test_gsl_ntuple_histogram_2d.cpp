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
    double attributes[2];

    /*
    double m_r;
    double z;
    double mbh;
    double mstar;
    */
    
  };


int main(int argc, char **argv)
{

  int Bin_Index_Hist2d(int, int,size_t,size_t);
  void Handle_Extreme_Values_Hist2d(gsl_histogram2d *, double, double);
  void Print_Hist2d_With_Header(gsl_histogram2d *, string);
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
  int num_bins_m_r = 20;

  double mean_m_r = min_m_r + 0.5 * (max_m_r - min_m_r);
  double sigma_m_r = 0.2 * (max_m_r - min_m_r);



  
  double min_z = 0;
  double max_z = 0.2;
  int num_bins_z = 25;

  double mean_z = min_z + 0.5 * (max_z - min_z);
  double sigma_z = 0.2 * (max_z - min_z);

  gsl_histogram * hist_gals_m_r = gsl_histogram_alloc(num_bins_m_r);
  gsl_histogram * hist_gals_z = gsl_histogram_alloc(num_bins_z);
  
  gsl_histogram_set_ranges_uniform(hist_gals_m_r, min_m_r, max_m_r);
  gsl_histogram_set_ranges_uniform(hist_gals_z, min_z, max_z);


  // setup random num generator
  gsl_rng *rangen;
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  //  gsl_rng_default_seed = (unsigned int)time(NULL);
  gsl_rng_default_seed = 0;
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);


  int num_samples = 10000000;

  clock_t begin;
  

  gsl_ntuple *gal_ntuple  = gsl_ntuple_create(filename_array, &gal_row, sizeof (gal_row));



  begin = clock();
  for (int t = 0; t < num_samples; t++)
    {

      gal_row.attributes[0] = mean_m_r + gsl_ran_gaussian(rangen, sigma_m_r);
      gal_row.attributes[1] = mean_z + gsl_ran_gaussian(rangen, sigma_z);

      //gal_row.mbh = 6.;
      //      gal_row.mstar = 8.;

      gsl_ntuple_write(gal_ntuple);
    }

  
  gsl_ntuple_close (gal_ntuple);


  
  clock_t end = clock();
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to generate %d samples\n",elapsed_secs, num_samples);


  // create 1d histogram
  begin = clock();

  S.function = &sel_func_1d;
  V.function = &val_func_1d;

  struct bin_param {gsl_histogram * hist; int icol;};

  bin_param bin_params = {hist_gals_m_r, 0};
  V.params = &bin_params;

  gal_ntuple = gsl_ntuple_open (filename_array, &gal_row, sizeof (gal_row));
  gsl_ntuple_project(hist_gals_m_r,gal_ntuple,&V,&S);
  gsl_ntuple_close(gal_ntuple);

  string x_name("m_r");
  string outfilename = "gsl_hist_1d_" + x_name + ".hist";
  Print_Hist1d(hist_gals_m_r, outfilename);

  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to 1d histogram\n",elapsed_secs);

  // create 2d histogram
  begin = clock();

  x_name = "m_r";
  string y_name = "z";

  //write header
  outfilename = "gsl_hist_2d_" + x_name + "_" + y_name + ".hist";
  FILE * outfile = fopen(outfilename.c_str(),"w");

  // header
    for (int ix = 0; ix < hist_gals_z->n + 1; ix++)
    {    
      fprintf(outfile,"%g ",hist_gals_z->range[ix]);
    }
  fprintf(outfile,"\n");
  for (int iy = 0; iy < hist_gals_m_r->n + 1; iy++)
    {    
      fprintf(outfile,"%g ",hist_gals_m_r->range[iy]);
    }
  fprintf(outfile,"\n");
  fclose(outfile);


  struct hist_param {gsl_histogram * hist;};
  //bin param already defined
  
  hist_param h_params = {hist_gals_m_r};
  V.params = &h_params;

  S.function = &sel_func_2d;
  V.function = &val_func_2d;



  for (int iz = 0; iz < num_bins_z; iz++)
    {

      bin_param b_params = {hist_gals_z, iz};
      S.params = &b_params;

      gal_ntuple = gsl_ntuple_open (filename_array, &gal_row, sizeof (gal_row));
  
      gsl_histogram_reset(hist_gals_m_r);
      gsl_ntuple_project(hist_gals_m_r,gal_ntuple,&V,&S);
      
      // will also need to supply info about which ntuple variable is the vertical bin
      Print_Hist1d(hist_gals_m_r,outfilename);

      gsl_ntuple_close(gal_ntuple);
      
    }


  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to project to histograms\n",elapsed_secs);


  begin = clock();
  gal_ntuple = gsl_ntuple_open (filename_array, &gal_row, sizeof (gal_row));
  gsl_ntuple_close(gal_ntuple);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds just to open and close the n-tuple\n",elapsed_secs);


  begin = clock();
  gal_ntuple = gsl_ntuple_open (filename_array, &gal_row, sizeof (gal_row));
  gsl_histogram_reset(hist_gals_m_r);
  gsl_ntuple_project(hist_gals_m_r,gal_ntuple,&V,&S);
  gsl_ntuple_close(gal_ntuple);
  end = clock();
  elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to open, close, and do a single histogram projectione\n",elapsed_secs);


  //string outfilename = "gsl_hist_2d_" + x_name + "_" + y_name + ".hist";

  //Print_Hist1d(hist_gals_z);
  
  gsl_histogram_free(hist_gals_m_r);
  gsl_histogram_free(hist_gals_z);
  
  gsl_rng_free(rangen);




}

int sel_func_1d (void *this_data, void *p)
{

  /*
  struct bin_struct { int column_index; };
  bin_struct binfo = *(struct bin_struct *)p;

  int column_index = binfo.column_index;
  */

  /*
  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double m_r, z, mbh,mstar;
  //  scale = *(double *) params;

  m_r = data->m_r;
  z = data->z;
  mbh = data -> mbh;
  mstar = data -> mstar;
  */

  return 1;

}


// will also need to pass info about which ntuple variable is the vertical bin
int sel_func_2d (void *this_data, void *p)
{

  struct bin_struct { gsl_histogram * hist; int ibin; };
  bin_struct binfo = *(struct bin_struct *)p;

  gsl_histogram * vertical_hist = binfo.hist;
  int vertical_bin_index = binfo.ibin;
  

  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double m_r, z, mbh,mstar;
  //  scale = *(double *) params;

  m_r = data->attributes[0];
  z = data->attributes[1];
  //  mbh = data -> mbh;
  //  mstar = data -> mstar;


  // need to replace z with generic vertical variable
  if (z >= gsl_histogram_max(vertical_hist))
    {
      z = (1. - 1.e-15) * gsl_histogram_max(vertical_hist);
    }
  if (z < gsl_histogram_min(vertical_hist))
    {
      z = gsl_histogram_min(vertical_hist);
    }

  return ( z >= vertical_hist->range[vertical_bin_index] &&  z < vertical_hist->range[vertical_bin_index + 1]);
  //  return ( z >= 0.1);
  // return true;
}

double val_func_2d (void *this_data, void *p)
{
  //  (void)(params); /* avoid unused parameter warning */

  struct hist_struct { gsl_histogram * hist; };

  
  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double m_r, z,mbh,mstar;

  hist_struct ph = *(struct hist_struct *)p;

  m_r = data->attributes[0];
  z = data->attributes[1];
  //  mbh = data->mbh;
  //  mstar = data->mstar;

  // replace z with generic name
  if (m_r < gsl_histogram_min(ph.hist))
    {
      m_r = gsl_histogram_min(ph.hist);
    }

  if (m_r >= gsl_histogram_max(ph.hist))
    {
      m_r = (1. - 1.e-15) * gsl_histogram_max(ph.hist);
    }


  return m_r;
}

double val_func_1d (void *this_data, void *p)
{
  //  (void)(params); /* avoid unused parameter warning */

  struct bin_struct { gsl_histogram * hist; int column_index; };
  bin_struct binfo = *(struct bin_struct *)p;

  int column_index = binfo.column_index;

  
  struct ntuple_data * data = (struct ntuple_data *) this_data;
  double this_col_value;

  this_col_value  = data->attributes[column_index];
  //  z = data->attributes[1];
  //  mbh = data->mbh;
  //  mstar = data->mstar;

  // replace z with generic name
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



int Bin_Index_Hist2d(int i, int j,size_t nx,size_t ny)
{
  return i * ny + j;

}

void Handle_Extreme_Values_Hist2d(gsl_histogram2d * hist2d, double x_value, double y_value)
{

  size_t i = 0,j = 0;
  size_t dummy = 0;
  size_t nx = hist2d->nx;
  size_t ny = hist2d->ny;

  if (x_value < gsl_histogram2d_xmin(hist2d))
    {
      i = 0;
      //left_errors_x++;
    }
  else if (x_value < gsl_histogram2d_xmax(hist2d))
    {
      gsl_histogram2d_find(hist2d,x_value, gsl_histogram2d_ymin(hist2d),&i,&dummy);
    }
  if (x_value >= gsl_histogram2d_xmax(hist2d))
    {
      i = nx - 1;
      //right_errors_x++;
    }
  else if (x_value >= gsl_histogram2d_xmin(hist2d))
    {
      gsl_histogram2d_find(hist2d,x_value, gsl_histogram2d_ymin(hist2d),&i,&dummy);
    }
	  
  if (y_value < gsl_histogram2d_ymin(hist2d))
    {
      j = 0;
      //left_errors_z++;
    }
  else if (y_value < gsl_histogram2d_ymax(hist2d))
    {
      gsl_histogram2d_find(hist2d,gsl_histogram2d_xmin(hist2d), y_value,&dummy,&j);
    }
  if (y_value >= gsl_histogram2d_ymax(hist2d))
    {
      j = ny - 1;
      //right_errors_z++;
    }
  else if (y_value >= gsl_histogram2d_ymin(hist2d))
    {
      gsl_histogram2d_find(hist2d,gsl_histogram2d_xmin(hist2d), y_value,&dummy,&j);
    }

  hist2d->bin[Bin_Index_Hist2d(i,j,nx,ny)]++;

}

void Print_Hist2d_With_Header(gsl_histogram2d * hist2d, string name)
{

  
  FILE * outfile = fopen(name.c_str(),"w");
  
  //  fprintf(outfile,"%s %d\n",x_name.c_str(), num_bins_m_r);
  //  fprintf(outfile,"%s %d\n",y_name.c_str(), num_bins_z);

  for (int ix = 0; ix < hist2d->nx + 1; ix++)
    {    
      fprintf(outfile,"%g ",hist2d->xrange[ix]);
    }
  fprintf(outfile,"\n");
  for (int iy = 0; iy < hist2d->ny + 1; iy++)
    {    
      fprintf(outfile,"%g ",hist2d->yrange[iy]);
    }
  fprintf(outfile,"\n");

  gsl_histogram2d_fprintf(outfile, hist2d, "%-16.6e", "%-16.6e");

  fclose(outfile);


}

void Print_Hist1d(gsl_histogram * hist1d, string name)
{

  FILE * outfile = fopen(name.c_str(),"a");

  gsl_histogram_fprintf(outfile, hist1d, "%-16.6e", "%-16.6e");

  fclose(outfile);

}



  
