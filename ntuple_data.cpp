#include "ntuple_data.h"
#include "galaxy.h" // to use the galaxy catalogue data struct

// this was in the original ntuple.c
#define EVAL(f,x) ((*((f)->function))(x,(f)->params))

// How do you easily specify here what type of data struct to use? Templating?
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
	   struct galaxy_catalogue_data * data_pointer = (galaxy_catalogue_data *) ntuple->ntuple_data;
	   
	   gsl_histogram_accumulate (h, EVAL(value_func, ntuple->ntuple_data),data_pointer->weight);
         }
     }
   while (1);
  
   return GSL_SUCCESS;
 }

