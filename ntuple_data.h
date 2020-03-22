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

// Similar to above, this time for value function
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


// this is just a small modification from the original ntuple.c
// Uses gsl_histogram_accumulate instead of gsl_histogram_increment
// and is templated for different data structures that might compose the ntuple
#define EVAL(f,x) ((*((f)->function))(x,(f)->params))

template <class data_struct>
int gsl_ntuple_project_weighted (gsl_histogram * h, gsl_ntuple * ntuple,
                     gsl_ntuple_value_fn * value_func, 
                     gsl_ntuple_select_fn * select_func)
 {
   size_t nread;

   do
     {
       nread = fread (ntuple->ntuple_data, ntuple->size,
                      1, ntuple->file);
  
       if (nread == 0 && feof(ntuple->file))
         {
           break ;
         }
       
       if (nread != 1) 
         {
           GSL_ERROR ("failed to read ntuple for projection", GSL_EFAILED);
         }
  
       if (EVAL(select_func, ntuple->ntuple_data))
         {
	   data_struct * data_pointer = (data_struct *) ntuple->ntuple_data;
	   
	   gsl_histogram_accumulate (h, EVAL(value_func, ntuple->ntuple_data),data_pointer->weight);
         }
     }
   while (1);
  
   return GSL_SUCCESS;
 }



#endif
