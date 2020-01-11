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

using std::vector;
using std::string;

int main(int argc, char **argv)
{

  int Bin_Index_Hist2d(int, int,size_t,size_t);
  void Handle_Extreme_Values_Hist2d(gsl_histogram2d *, double, double);
  void Print_Hist2d_With_Header(gsl_histogram2d *, string);
  
  double min_m_r = 12;
  double max_m_r = 28;
  int num_bins_m_r = 80;

  double mean_m_r = min_m_r + 0.5 * (max_m_r - min_m_r);
  double sigma_m_r = 0.2 * (max_m_r - min_m_r);

  string x_name("m_r");
  string y_name("z");
  
  double min_z = 0;
  double max_z = 0.2;
  int num_bins_z = 100;

  double mean_z = min_z + 0.5 * (max_z - min_z);
  double sigma_z = 0.2 * (max_z - min_z);

  gsl_histogram2d * hist_gals = gsl_histogram2d_alloc(num_bins_m_r,num_bins_z);
  gsl_histogram2d_set_ranges_uniform(hist_gals, min_m_r, max_m_r,min_z,max_z);

  string outfilename = "gsl_hist_2d_" + x_name + "_" + y_name + ".hist";


  // setup random num generator
  gsl_rng *rangen;
  const gsl_rng_type * TypeR;
  gsl_rng_env_setup();
  gsl_rng_default_seed = (unsigned int)time(NULL);
  TypeR = gsl_rng_default;
  rangen = gsl_rng_alloc (TypeR);


  int num_samples = 10000000;



  double m_r_value = 0.;
  double z_value = 0.;
  /*
  int left_errors_m_r = 0;
  int right_errors_m_r = 0;
  int left_errors_z = 0;
  int right_errors_z = 0;
  int general_errors = 0;
  */
  clock_t begin = clock();
  for (int t = 0; t < num_samples; t++)
    {

      m_r_value = mean_m_r + gsl_ran_gaussian(rangen, sigma_m_r);
      z_value = mean_z + gsl_ran_gaussian(rangen, sigma_z);
								   
      int increment_status = gsl_histogram2d_increment(hist_gals,m_r_value,z_value);

      if (increment_status == GSL_EDOM)
	{
	  Handle_Extreme_Values_Hist2d(hist_gals,m_r_value,z_value);
	}

    }

  clock_t end = clock();
    
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to generate %d samples\n",elapsed_secs, num_samples);

  /*
  printf("# of general errors: %d\n",general_errors);
  printf("# of left errors m_r: %d\n",left_errors_m_r);
  printf("# of right errors m_r: %d\n",right_errors_m_r);
  printf("# of left errors z: %d\n",left_errors_z);
  printf("# of right errors z: %d\n",right_errors_z);
  */

  Print_Hist2d_With_Header(hist_gals, outfilename);

  gsl_histogram2d_free(hist_gals);

  gsl_rng_free(rangen);


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



  
