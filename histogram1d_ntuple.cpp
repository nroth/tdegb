#include <stdio.h>
#include <math.h>
#include "histogram1d_ntuple.h"

using std::vector;
using std::string;



//***************************************************************
// Constructors
//***************************************************************

Histogram1dNtuple::Histogram1dNtuple()
{
  base_name = "default";
  axis_name = "default";
  icol = 0;
}

Histogram1dNtuple::Histogram1dNtuple(int num_bins, vector<double> bin_specs, string b_name, string a_name, int index, string ntuple_filename)
{
  Init(num_bins,bin_specs, b_name,a_name, index, ntuple_filename);
}

//***************************************************************
// Destructor
//***************************************************************

Histogram1dNtuple::~Histogram1dNtuple()
{
  gsl_histogram_free(hist);
}

//***************************************************************
// Initialization
//***************************************************************
void Histogram1dNtuple::Init (int num_bins, vector<double> bin_specs, string b_name, string a_name, int index, string ntuple_filename)
{

  hist = gsl_histogram_alloc(num_bins);
  gsl_histogram_set_ranges_uniform(hist, bin_specs[0], bin_specs[1]);
  icol = index;
  
  base_name = b_name;
  axis_name = a_name;

  //  gsl_function_pp Fp( std::bind(&Class::member_function, &(*this),  std::placeholders::_1) );
  //  gsl_function *F = static_cast<gsl_function*>(&Fp);

//  S.function = &Histogram1dNtuple::sel_func_1d;
//  V.function = &Histogram1dNtuple::val_func_1d;

  strcpy(ntuple_filename_array, ntuple_filename.c_str());

}

//***************************************************************
// GSL helper functions
//***************************************************************

int Histogram1dNtuple::sel_func_1d (void *this_data)
{
  return 1;
}


  
double Histogram1dNtuple::val_func_1d (void *this_data)
{
  
  struct data * data_pointer = (struct data *) this_data;

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
void Histogram1dNtuple::Print_Histogram_1D()
{

  Histogram1dNtuple* ptr2S = this;
  auto ptrS = [=](void *ntuple_data)->int{return ptr2S->sel_func_1d(ntuple_data);};
  gsl_ntuple_select_fn_pp<decltype(ptrS)> FpS(ptrS);     
  gsl_ntuple_select_fn *S = static_cast<gsl_ntuple_select_fn*>(&FpS);

  Histogram1dNtuple* ptr2V = this;
  auto ptrV = [=](void *ntuple_data)->double{return ptr2V->val_func_1d(ntuple_data);};
  gsl_ntuple_value_fn_pp<decltype(ptrV)> FpV(ptrV);     
  gsl_ntuple_value_fn *V = static_cast<gsl_ntuple_value_fn*>(&FpV);

  struct data data_row;
  gsl_ntuple *this_ntuple = gsl_ntuple_open(ntuple_filename_array, &data_row, sizeof (data_row));  
  //gsl_ntuple_project(hist,this_ntuple,&V,&S);
  gsl_ntuple_project(hist,this_ntuple,V,S);
  gsl_ntuple_close(this_ntuple);

  string outfilename = base_name + "_" + axis_name + "_1d.hist";
  FILE * outfile = fopen(outfilename.c_str(),"w");
  gsl_histogram_fprintf(outfile, hist, "%-16.6e", "%-16.6e");
  fclose(outfile);

}

