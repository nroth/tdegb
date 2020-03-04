#ifndef HISTOGRAM1DNTUPLE_H
#define HISTOGRAM1DNTUPLE_H 

#include <string>
#include <vector>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_ntuple.h>

using std::string;
using std::vector;


struct data
  {
    double attributes[3];
  };


class Histogram1dNtuple {
 
private:

  gsl_histogram * hist; //= gsl_histogram_alloc(num_bins_m_r);
  string base_name;
  string axis_name; // histogram name
  int icol;

  char ntuple_filename_array[35];

  // GSL ntuple requires a struct
    struct data 
  {
    double attributes[3];
  };

  //  gsl_ntuple_select_fn S;
  //  gsl_ntuple_value_fn V; 

  int sel_func_1d (void *);
  double val_func_1d (void *);


public:

  // constructors
  Histogram1dNtuple();
  Histogram1dNtuple(int, vector<double> , string, string, int, string);
  
  // Initialize
  void Init(int, vector<double>, string, string, int, string);
  void Print_Histogram_1D();

  ~Histogram1dNtuple();
  
 
};



// This is required so that you can use a class member function as a GSL histogram selection function
// See https://stackoverflow.com/questions/13074756/how-to-avoid-static-member-function-when-using-gsl-with-c/18181494#18181494
template<typename F>
class gsl_ntuple_select_fn_pp : public gsl_ntuple_select_fn {
 public:
  gsl_ntuple_select_fn_pp(const F& func) : _func(func) {
    function = &gsl_ntuple_select_fn_pp::invoke;
    params=this;
  }
  private:
  const F& _func;
  static int invoke(void *ntuple_data, void *params) {
    return static_cast<gsl_ntuple_select_fn_pp*>(params)->_func(ntuple_data);
  }
};

// This is required so that you can use a class member function as a GSL histogram value function
template<typename F>
class gsl_ntuple_value_fn_pp : public gsl_ntuple_value_fn {
 public:
  gsl_ntuple_value_fn_pp(const F& func) : _func(func) {
    function = &gsl_ntuple_value_fn_pp::invoke;
    params=this;
  }
  private:
  const F& _func;
  static double invoke(void *ntuple_data, void *params) {
    return static_cast<gsl_ntuple_value_fn_pp*>(params)->_func(ntuple_data);
  }
};


#endif
