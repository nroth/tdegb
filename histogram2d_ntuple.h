#ifndef HISTOGRAM2DNTUPLE_H
#define HISTOGRAM2DNTUPLE_H 

#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_ntuple.h>
#include "ntuple_data.h"


using std::string;
using std::vector;


template <class data_struct>
class Histogram2dNtuple {
 
private:

  gsl_histogram * hist1;
  gsl_histogram * hist2; 
  string base_name;
  string v_axis_name; 
  string h_axis_name; 
  int icol1;
  int icol2;
  int ibin;

  string ntuple_filename_string;

  int sel_func_2d (void *);
  double val_func_2d (void *);
  

public:

  // constructors
  Histogram2dNtuple();
  Histogram2dNtuple(vector<size_t>, vector<vector<double>> , string, string, string, vector<int>, string);
  
  // Initialize
  void Init(vector<size_t>, vector<vector<double>>, string, string, string, vector<int>, string);
  void Print_Histogram_2D_With_Header(bool, string);

  ~Histogram2dNtuple();
  
 
};


//***************************************************************
// Constructors
//***************************************************************

template <class data_struct>
Histogram2dNtuple<data_struct>::Histogram2dNtuple()
{
  base_name = "default";
  v_axis_name = "default";
  h_axis_name = "default";
  icol1 = 0;
  icol2 = 0;
  ibin = 0;
}

template <class data_struct>
Histogram2dNtuple<data_struct>::Histogram2dNtuple(vector<size_t> num_bins, vector<vector<double>> bin_specs, string b_name, string v_name, string h_name, vector<int> indices, string ntuple_filename)
{
  Init(num_bins,bin_specs, b_name,v_name,h_name, indices, ntuple_filename);
}

//***************************************************************
// Destructor
//***************************************************************

template <class data_struct>
Histogram2dNtuple<data_struct>::~Histogram2dNtuple()
{
  gsl_histogram_free(hist1);
  gsl_histogram_free(hist2);
}

//***************************************************************
// Initialization
//***************************************************************

template <class data_struct>
void Histogram2dNtuple<data_struct>::Init (vector<size_t> num_bins, vector<vector<double>> bin_specs, string b_name, string v_name, string h_name, vector<int> indices, string ntuple_filename)
{

  hist1 = gsl_histogram_alloc(num_bins[0]);
  hist2 = gsl_histogram_alloc(num_bins[1]);
  gsl_histogram_set_ranges_uniform(hist1, bin_specs[0][0], bin_specs[0][1]);
  gsl_histogram_set_ranges_uniform(hist2, bin_specs[1][0], bin_specs[1][1]);
  icol1 = indices[0];
  icol2 = indices[1];
  
  base_name = b_name;
  v_axis_name = v_name;
  h_axis_name = h_name;

  ntuple_filename_string = ntuple_filename;

}

//***************************************************************
// Functions for specifying how histograms are constructed (with GSL)
// Selection function used to make cuts
// Value function to determine the quantity that is binned
//***************************************************************

template <class data_struct>
int Histogram2dNtuple<data_struct>::sel_func_2d (void *this_data)
{

  data_struct * data_pointer = (data_struct *) this_data;
  double this_horizontal_data = data_pointer->attributes[icol2];
  //  double this_horizontal_data = log10(data_pointer->attributes[icol2]);
  //  double this_horizontal_data = data_pointer->attributes[23] - data_pointer->attributes[19];


  // handle extreme values
  if (this_horizontal_data >= gsl_histogram_max(hist2))
    {
      if ( gsl_histogram_max(hist2) >=  0.)
	this_horizontal_data = (1. - 1.e-15) * gsl_histogram_max(hist2);
      else
	this_horizontal_data = (1. + 1.e-15) * gsl_histogram_max(hist2);
    }
  else 
    {
      if (this_horizontal_data < gsl_histogram_min(hist2))
	{
	  this_horizontal_data = gsl_histogram_min(hist2);
	}
    }
      

  
  return ( this_horizontal_data >= hist2->range[ibin] && this_horizontal_data < hist2->range[ibin + 1]);
  //  return ( this_horizontal_data >= hist2->range[ibin] && this_horizontal_data < hist2->range[ibin + 1] && data_pointer->attributes[z_i] < 0.4);

}

template <class data_struct>
double Histogram2dNtuple<data_struct>::val_func_2d (void *this_data)
{
  
  data_struct * data_pointer = (data_struct *) this_data;
  double this_col_data;
  
  this_col_data = data_pointer->attributes[icol1];
  //  this_col_data = log10(data_pointer->attributes[icol1]);
  //this_col_data = data_pointer->attributes[23] - data_pointer->attributes[19];


  if (this_col_data < gsl_histogram_min(hist1))
    {
      this_col_data = gsl_histogram_min(hist1);
    }

  else
    {
      if (this_col_data >= gsl_histogram_max(hist1))
	{
	  if ( gsl_histogram_max(hist1) >= 0.)
	    this_col_data = (1. - 1.e-15) * gsl_histogram_max(hist1);
	  else
	    this_col_data = (1. + 1.e-15) * gsl_histogram_max(hist1);
	}
    }


  return this_col_data;

}


//***************************************************************
// n-tuple projection to histogram, and output
//***************************************************************
template <class data_struct>
void Histogram2dNtuple<data_struct>::Print_Histogram_2D_With_Header(bool weighted, string save_path)
{

  string outfilename = save_path + '/' + base_name + "_" + v_axis_name + "_" + h_axis_name + "_2d.hist";

  Histogram2dNtuple* ptr2S = this;
  auto ptrS = [=](void *ntuple_data)->int{return ptr2S->sel_func_2d(ntuple_data);};
  gsl_ntuple_select_fn_pp<decltype(ptrS)> FpS(ptrS);     
  gsl_ntuple_select_fn *S = static_cast<gsl_ntuple_select_fn*>(&FpS);

  Histogram2dNtuple* ptr2V = this;
  auto ptrV = [=](void *ntuple_data)->double{return ptr2V->val_func_2d(ntuple_data);};
  gsl_ntuple_value_fn_pp<decltype(ptrV)> FpV(ptrV);     
  gsl_ntuple_value_fn *V = static_cast<gsl_ntuple_value_fn*>(&FpV);

  data_struct data_row;

  int name_length = ntuple_filename_string.length();
  char* ntuple_filename_array = new char[name_length + 1];
  strcpy(ntuple_filename_array, ntuple_filename_string.c_str());

  // write header
  FILE * outfile = fopen(outfilename.c_str(),"w");
  for (size_t iy = 0; iy < hist2->n + 1; iy++)
    {    
      fprintf(outfile,"%g ",hist2->range[iy]);
    }
  fprintf(outfile,"\n");

  for (size_t ix = 0; ix < hist1->n + 1; ix++)
    {    
      fprintf(outfile,"%g ",hist1->range[ix]);
    }
  fprintf(outfile,"\n");
  fclose(outfile);

  outfile = fopen(outfilename.c_str(),"a");

    for (size_t i = 0; i < hist2->n; i++)
    {
      ibin = i;

      gsl_ntuple *this_ntuple = gsl_ntuple_open (ntuple_filename_array, &data_row, sizeof (data_row));

  
      gsl_histogram_reset(hist1);
      if (weighted)
	{
	  gsl_ntuple_project_weighted<data_struct>(hist1,this_ntuple,V,S);
	}
      else
	{
	  gsl_ntuple_project(hist1,this_ntuple,V,S);
	}
      gsl_histogram_fprintf(outfile, hist1, "%-16.6e", "%-16.6e");

      gsl_ntuple_close(this_ntuple);
      
    }

    delete [] ntuple_filename_array;
    fclose(outfile);
}



#endif
