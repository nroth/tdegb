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

  string galaxy_catalogue_filename = "/Users/nathanielroth/Dropbox/research/TDE/host_galaxies/sjoert_catalogue/sjoert_catalogue_ascii.dat";

  float z, ra, dec, mass, b300, b1000, ssfr, bt, r50_kpc, sersic_n, sigma, sigma_SDSS,sigma_SDSS_err, mbh_sigma, mbh_bulge, m_r, Mabs_r, m_g, M_g;

  int num_galaxies = 0;
  
  clock_t begin = clock();

  FILE * pFile;
  pFile = fopen(galaxy_catalogue_filename.c_str(),"r");

  vector<double> catalogue_data(hist_gals_dimension);

  Initialize_IMF(); // for computing volumetric disruption rates

  int counter = 0;
  while(fscanf(pFile,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &z, &ra, &dec, &mass, &b300, &b1000, &ssfr, &bt, &r50_kpc, &sersic_n, &sigma, &sigma_SDSS,&sigma_SDSS_err, &mbh_sigma, &mbh_bulge, &m_r, &Mabs_r, &m_g, &M_g) != EOF)
    {


      //      if (counter > 1007270 && counter < 1007291)
      //	{
      
      catalogue_data[0] = log10(mass);
      catalogue_data[1] = m_g - m_r;
      catalogue_data[2] = mbh_bulge;
      catalogue_data[3] = mbh_sigma;
      catalogue_data[4] = z;

      hist_gals.Count(catalogue_data);

      hist_vol_disrupt.Count(catalogue_data,Total_Disruption_Rate(pow(10.,mbh_bulge),z)); // the volumetric disruption rate histogram is just like the host galaxy histogram, but weighted by per-galaxy disruption rate. The z is needed to convert from galaxy time frame to observer rest framee
      double beta = 1.;
      double T = 3.e4;
      //      hist_detected_disrupt.Count(catalogue_data,Total_Disruption_Rate_Observed_Gband(pow(10.,mbh_bulge),beta,z,T));
      
      num_galaxies++;


      //      if (counter % 1 == 0)
      //	{
      //	  //printf("at galaxy %d\n",num_galaxies);
      //	  printf("at galaxy %d\n",counter);
      //	}
	    

      //      	}
      counter++;



    }

  fclose(pFile);

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

  hist_projected_vol_disrupt = hist_vol_disrupt.Create_Projected_Histogram(kept_axes);
  hist_projected_vol_disrupt.Print_Histogram_2D(0,1);


  /////////

  // Print out some histograms of detectable disruption
  /*

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


  */  

}
  
