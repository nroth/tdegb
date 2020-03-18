#ifndef NTUPLE_DATA_H
#define NTUPLE_DATA_H

// GSL ntuple requires a struct
struct data
  {
    double attributes[12];
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
