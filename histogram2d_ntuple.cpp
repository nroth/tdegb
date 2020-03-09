#include <stdio.h>
#include <math.h>
#include "histogram2d_ntuple.h"

using std::vector;
using std::string;



//***************************************************************
// Constructors
//***************************************************************

Histogram2dNtuple::Histogram2dNtuple()
{
  base_name = "default";
  v_axis_name = "default";
  h_axis_name = "default";
  icol1 = 0;
  icol2 = 0;
  ibin = 0;
}

Histogram2dNtuple::Histogram2dNtuple(vector<int> num_bins, vector<vector<double>> bin_specs, string b_name, string v_name, string h_name, vector<int> indices, int ibin, string ntuple_filename)
{
  Init(num_bins,bin_specs, b_name,v_name,h_name, indices, ibin, ntuple_filename);
}

//***************************************************************
// Destructor
//***************************************************************

Histogram2dNtuple::~Histogram2dNtuple()
{
  gsl_histogram_free(hist1);
  gsl_histogram_free(hist2);
}

//***************************************************************
// Initialization
//***************************************************************
void Histogram2dNtuple::Init (vector<int> num_bins, vector<vector<double>> bin_specs, string b_name, string v_name, string h_name, vector<int> indices, int index, string ntuple_filename)
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

  strcpy(ntuple_filename_array, ntuple_filename.c_str());

}

//***************************************************************
// GSL helper functions
//***************************************************************

int Histogram2dNtuple::sel_func_2d (void *this_data)
{
  
  struct data * data_pointer = (struct data *) this_data;
  double this_horizontal_data = data_pointer->attributes[icol2];

  // handle extreme values
  if (this_horizontal_data >= gsl_histogram_max(hist2))
    {
      this_horizontal_data = (1. - 1.e-15) * gsl_histogram_max(hist2);
    }
  if (this_horizontal_data < gsl_histogram_min(hist2))
    {
      this_horizontal_data = gsl_histogram_min(hist2);
    }

  return ( this_horizontal_data >= hist2->range[ibin] && this_horizontal_data < hist2->range[ibin + 1]);
  //  return ( this_horizontal_data >= 0.1);
  // return true;

}

double Histogram2dNtuple::val_func_2d (void *this_data)
{
  
  struct data * data_pointer = (struct data *) this_data;
  double this_col_data; //, v_data;

  this_col_data = data_pointer->attributes[icol1];

  if (this_col_data < gsl_histogram_min(hist1))
    {
      this_col_data = gsl_histogram_min(hist1);
    }

  if (this_col_data >= gsl_histogram_max(hist1))
    {
      this_col_data = (1. - 1.e-15) * gsl_histogram_max(hist1);
    }

  return this_col_data;

}



//***************************************************************
// n-tuple projection to histogram, and output
//***************************************************************
void Histogram2dNtuple::Print_Histogram_2D_With_Header()
{

  string outfilename = "gsl_hist_2d_" + v_axis_name + "_" + h_axis_name + ".hist";


  Histogram2dNtuple* ptr2S = this;
  auto ptrS = [=](void *ntuple_data)->int{return ptr2S->sel_func_2d(ntuple_data);};
  gsl_ntuple_select_fn_pp<decltype(ptrS)> FpS(ptrS);     
  gsl_ntuple_select_fn *S = static_cast<gsl_ntuple_select_fn*>(&FpS);

  Histogram2dNtuple* ptr2V = this;
  auto ptrV = [=](void *ntuple_data)->double{return ptr2V->val_func_2d(ntuple_data);};
  gsl_ntuple_value_fn_pp<decltype(ptrV)> FpV(ptrV);     
  gsl_ntuple_value_fn *V = static_cast<gsl_ntuple_value_fn*>(&FpV);

  struct data data_row;

  // write header
  FILE * outfile = fopen(outfilename.c_str(),"w");
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

  outfile = fopen(outfilename.c_str(),"a");

    for (int i = 0; i < hist2->n; i++)
    {
      ibin = i;

      gsl_ntuple *this_ntuple = gsl_ntuple_open (ntuple_filename_array, &data_row, sizeof (data_row));
  
      gsl_histogram_reset(hist1);
      gsl_ntuple_project(hist1,this_ntuple,V,S);

      outfile = fopen(outfilename.c_str(),"a"); 
      gsl_histogram_fprintf(outfile, hist1, "%-16.6e", "%-16.6e");

      gsl_ntuple_close(this_ntuple);
      
    }

  fclose(outfile);

}

