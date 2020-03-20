#ifndef NTUPLE_DATA_H
#define NTUPLE_DATA_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_ntuple.h>

// GSL ntuple requires a struct
struct data
  {
    double attributes[12];
    double weight; 
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



int gsl_ntuple_project_weighted (gsl_histogram *, gsl_ntuple *,gsl_ntuple_value_fn *, 
				 gsl_ntuple_select_fn*);


#endif
