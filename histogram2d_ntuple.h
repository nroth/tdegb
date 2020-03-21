#ifndef HISTOGRAM2DNTUPLE_H
#define HISTOGRAM2DNTUPLE_H 

#include <string>
#include <vector>
#include <gsl/gsl_histogram.h>


using std::string;
using std::vector;



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

  char ntuple_filename_array[35];

  int sel_func_2d (void *);
  double val_func_2d (void *);
  

public:

  // constructors
  Histogram2dNtuple();
  Histogram2dNtuple(vector<int>, vector<vector<double>> , string, string, string, vector<int>, int, string);
  
  // Initialize
  void Init(vector<int>, vector<vector<double>>, string, string, string, vector<int>, int, string);
  void Print_Histogram_2D_With_Header(bool);

  ~Histogram2dNtuple();
  
 
};


#endif
