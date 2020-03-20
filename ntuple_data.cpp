#include "ntuple_data.h"

// this was in the original ntuple.c
#define EVAL(f,x) ((*((f)->function))(x,(f)->params))

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
	   struct data * data_pointer = (data *) ntuple->ntuple_data;
	   
	   gsl_histogram_accumulate (h, EVAL(value_func, ntuple->ntuple_data),data_pointer->weight);
         }
     }
   while (1);
  
   return GSL_SUCCESS;
 }

