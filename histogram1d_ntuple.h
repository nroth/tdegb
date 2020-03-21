#ifndef HISTOGRAM1DNTUPLE_H
#define HISTOGRAM1DNTUPLE_H 

#include <string>
#include <vector>
#include <gsl/gsl_histogram.h>

using std::string;
using std::vector;


class Histogram1dNtuple {
 
private:

  gsl_histogram * hist; 
  string base_name;
  string axis_name; 
  int icol;

  char ntuple_filename_array[35];

  int sel_func_1d (void *);
  double val_func_1d (void *);


public:

  // constructors
  Histogram1dNtuple();
  Histogram1dNtuple(int, vector<double> , string, string, int, string);
  
  // Initialize
  void Init(int, vector<double>, string, string, int, string);
  void Print_Histogram_1D(bool);

  ~Histogram1dNtuple();
  
 
};



#endif
