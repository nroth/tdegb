#include <fstream>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "physical_constants.h"
#include "histogramNd.h"
#include "cdf.h"
#include "individual_disruption.h"
#include "magnitudes.h"
#include "hdf5.h"
#include "hdf5_hl.h"


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

  printf("Galaxy histogram dimension is %d\n", hist_gals_dimension);
  printf("Total number of bins in galaxy histogram is %d\n", hist_gals.Get_Num_Bins_Total());

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

  string catalogue_filename = "/Users/nathanielroth/Dropbox/research/TDE/host_galaxies/sjoert_catalogue/van_velzen_2018_catalogue.h5";

  hid_t file_id = H5Fopen(catalogue_filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);

  int num_galaxies = 6101944; // should think more about how to make this flexible

  double* z = new double[num_galaxies]; // for some reason, the hdf5 read will only work if you do the manual memory allocation this way (and you must delete later)
  double* m_g = new double[num_galaxies];
  double* m_r = new double[num_galaxies];
  double* mbh_bulge = new double[num_galaxies];
  double* mbh_sigma = new double[num_galaxies];
  double* mass = new double[num_galaxies];
  double* sersic_n = new double[num_galaxies];
  double* r50_kpc = new double[num_galaxies];

  H5LTread_dataset_double(file_id,"z",z);
  H5LTread_dataset_double(file_id,"m_g",m_g);
  H5LTread_dataset_double(file_id,"m_r",m_r);
  H5LTread_dataset_double(file_id,"mbh_bulge",mbh_bulge);
  H5LTread_dataset_double(file_id,"mbh_sigma",mbh_sigma);
  H5LTread_dataset_double(file_id,"mass",mass);
  H5LTread_dataset_double(file_id,"sersic_n",sersic_n);
  H5LTread_dataset_double(file_id,"r50_kpc",r50_kpc);

  vector<double> catalogue_data(hist_gals_dimension);

  Initialize_IMF(); // for computing volumetric disruption rates

  for (int i = 0; i < num_galaxies; i++)
    {

      catalogue_data[0] = log10(mass[i]);
      catalogue_data[1] = m_g[i] - m_r[i];
      catalogue_data[2] = mbh_bulge[i];
      catalogue_data[3] = mbh_sigma[i];
      catalogue_data[4] = z[i];

      hist_gals.Count(catalogue_data);

      double nuker_gamma = find_nuker_gammaprime_from_sersic(sersic_n[i],r50_kpc[i],z[i]);

      hist_vol_disrupt.Count(catalogue_data,Total_Disruption_Rate(pow(10.,mbh_sigma[i]),z[i],nuker_gamma)); // the volumetric disruption rate histogram is just like the host galaxy histogram, but weighted by per-galaxy disruption rate. The z is needed to convert from galaxy time frame to observer rest framee
      double beta = 1.;
      double T = 3.e4;

      double m_limit_contrast = find_host_contrast_magnitude(m_g[i],sersic_n[i],r50_kpc[i],z[i]);
      
      hist_detected_disrupt.Count(catalogue_data,Total_Disruption_Rate_Observed_Rband(pow(10.,mbh_sigma[i]),beta,z[i],T,23.,nuker_gamma));
      //hist_detected_disrupt.Count(catalogue_data,Total_Disruption_Rate_Observed_Gband(pow(10.,mbh_sigma),beta,z,T,23.));
      

      /*
      if (counter % 10000 == 0)
      	{
	  printf("at galaxy %d\n",counter);
	}
	    

      counter++;
      */

    }

  H5Fclose(file_id);
  delete[] z;
  delete[] m_g;
  delete[] m_r;
  delete[] mbh_bulge;
  delete[] mbh_sigma;
  delete[] mass;
  delete[] sersic_n;
  delete[] r50_kpc;

  clock_t end = clock();
    
  float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
  printf("#\n# It took %f seconds to bin %d entries\n",elapsed_secs, num_galaxies);


  // Print out some histograms of galaxy properties
  HistogramNd hist_projected_gals;
  
  vector<int> kept_axes(2);
  int ka[2];
    
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
  
