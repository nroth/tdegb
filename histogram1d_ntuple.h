#ifndef HISTOGRAM1DNTUPLE_H
#define HISTOGRAM1DNTUPLE_H 

#include <stdio.h>
#include <string>
#include <vector>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_ntuple.h>
#include "ntuple_data.h"

using std::string;
using std::vector;

template <class data_struct>
class Histogram1dNtuple {
 
private:

  gsl_histogram * hist; 
  string base_name;
  string axis_name; 
  int icol;

  string ntuple_filename_string;

  int sel_func_1d (void *);
  double val_func_1d (void *);


public:

  // constructors
  Histogram1dNtuple();
  Histogram1dNtuple(int, vector<double> , string, string, int, string);
  
  // Initialize
  void Init(int, vector<double>, string, string, int, string);
  void Print_Histogram_1D(bool, string);

  ~Histogram1dNtuple();
  
 
};



//***************************************************************
// Constructors
//***************************************************************

template <class data_struct>
Histogram1dNtuple<data_struct>::Histogram1dNtuple()
{
  base_name = "default";
  axis_name = "default";
  icol = 0;
}

template <class data_struct>
Histogram1dNtuple<data_struct>::Histogram1dNtuple(int num_bins, vector<double> bin_specs, string b_name, string a_name, int index, string ntuple_filename)
{
  Init(num_bins,bin_specs, b_name,a_name, index, ntuple_filename);
}

//***************************************************************
// Destructor
//***************************************************************

template <class data_struct>
Histogram1dNtuple<data_struct>::~Histogram1dNtuple()
{
  gsl_histogram_free(hist);
}

//***************************************************************
// Initialization
//***************************************************************

template <class data_struct>
void Histogram1dNtuple<data_struct>::Init (int num_bins, vector<double> bin_specs, string b_name, string a_name, int index, string ntuple_filename)
{

  hist = gsl_histogram_alloc(num_bins);
  gsl_histogram_set_ranges_uniform(hist, bin_specs[0], bin_specs[1]);
  icol = index;
  
  base_name = b_name;
  axis_name = a_name;

  ntuple_filename_string = ntuple_filename;

}


//***************************************************************
// Functions for specifying how histograms are constructed (with GSL)
// Selection function used to make cuts
// Value function to determine the quantity that is binned
//***************************************************************

template <class data_struct>
int Histogram1dNtuple<data_struct>::sel_func_1d (void *this_data)
{
  return 1;
}


template <class data_struct>  
double Histogram1dNtuple<data_struct>::val_func_1d (void *this_data)
{

  data_struct * data_pointer = (data_struct *) this_data;
  double this_col_value  = data_pointer->attributes[icol];

  // handle extreme values
  if (this_col_value < gsl_histogram_min(hist))
    {
      this_col_value = gsl_histogram_min(hist);
    }

  if (this_col_value >= gsl_histogram_max(hist))
    {
      this_col_value = (1. - 1.e-15) * gsl_histogram_max(hist);
    }


  return this_col_value;
}



//***************************************************************
// n-tuple projection to histogram, and output
//***************************************************************
template <class data_struct>
void Histogram1dNtuple<data_struct>::Print_Histogram_1D(bool weighted, string save_path)
{

  Histogram1dNtuple* ptr2S = this;
  auto ptrS = [=](void *ntuple_data)->int{return ptr2S->sel_func_1d(ntuple_data);};
  gsl_ntuple_select_fn_pp<decltype(ptrS)> FpS(ptrS);     
  gsl_ntuple_select_fn *S = static_cast<gsl_ntuple_select_fn*>(&FpS);

  Histogram1dNtuple* ptr2V = this;
  auto ptrV = [=](void *ntuple_data)->double{return ptr2V->val_func_1d(ntuple_data);};
  gsl_ntuple_value_fn_pp<decltype(ptrV)> FpV(ptrV);     
  gsl_ntuple_value_fn *V = static_cast<gsl_ntuple_value_fn*>(&FpV);

  int name_length = ntuple_filename_string.length();
  char* ntuple_filename_array = new char[name_length];
  strcpy(ntuple_filename_array, ntuple_filename_string.c_str());
  
  data_struct data_row;
  gsl_ntuple *this_ntuple = gsl_ntuple_open(ntuple_filename_array, &data_row, sizeof (data_row));
  delete [] ntuple_filename_array;

  if (weighted)
    {
      gsl_ntuple_project_weighted<data_struct>(hist,this_ntuple,V,S);
    }
  else
    {
      gsl_ntuple_project(hist,this_ntuple,V,S);
    }
  gsl_ntuple_close(this_ntuple);

  string outfilename = save_path + '/' + base_name + "_" + axis_name + "_1d.hist";
  FILE * outfile = fopen(outfilename.c_str(),"w");
  gsl_histogram_fprintf(outfile, hist, "%-16.6e", "%-16.6e");
  fclose(outfile);

}

#endif
