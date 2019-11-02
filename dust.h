#ifndef _DUST_H
#define _DUST_H

#include <math.h>
#include <stdio.h>
#include "physical_constants.h"


static double Cardelli_Extinction(double x, double R_v)
{

  double a = 0.;
  double b = 0.;

  if (x < 0.3 || x > 10.)
    {
      printf("ERROR: Cardelli extinction not defined for x < 0.3 or x > 10.\n");
      return 0.;
    }
  if (x >= 0.3 && x < 1.1)
    {
      a = 0.574 * pow(x,1.61);
      b = -0.527 * pow(x,1.61);
    }
  if (x >= 1.1 and x < 3.3)
    {
      double y = x - 1.82;
      a = 1. + 0.17699 * y - 0.50447 * pow(y,2.) -0.02427 * pow(y,3.) + 0.72085 * pow(y,4.) + 0.01979 * pow(y,5.) - 0.77530 * pow(y,6.) + 0.32999 * pow(y,7.);
      b = 1.41338 * y + 2.28305 * pow(y,2.) + 1.07233 * pow(y,3.) - 5.38434 * pow(y,4.) - 0.62251 * pow(y,5.) + 5.3026 * pow(y,6) - 2.09002 * pow(y,7.);
    }
    if (x >= 3.3 and x < 8.)
      {
	double Fa = 0.;
	double Fb = 0.;
	if (x >= 5.9)
	  {
            Fa = -0.04473 * pow(x - 5.9,2.) - 0.009779 *pow(x - 5.9,3.);
            Fb = 0.213 * pow(x - 5.9,2.) + 0.1207 *pow(x - 5.9,3.);
	  }
        a = 1.752 - 0.316 * x - 0.104/(pow(x - 4.67,2.) + 0.341) + Fa;
        b = -3.09 + 1.825 * x + 1.206/(pow(x - 4.62,2.) + 0.263) + Fb;
      }
    if (x >= 8. && x <= 10.)
      {
        a = -1.073 - 0.628 * (x - 8.) + 0.137 * pow(x - 8.,2.) - 0.070 * pow(x - 8.,3.);
        b = 13.67 + 4.257 * (x - 8.) - 0.42 * pow(x - 8.,2.) + 0.374 * pow(x - 8.,3.);
      }
	  /*
    #if (x > 10. and x < 10.96):
    #    #print "ERROR, Cardelli extinction not defined"
    #    a = -1.073 - 0.628 * (10. - 8.) + 0.137 * pow(10. - 8.,2.) - 0.070 * pow(10. - 8.,3.)
    #    b = 13.67 + 4.257 * (10. - 8.) - 0.42 * pow(10. - 8.,2.) + 0.374 * pow(10. - 8.,3.)
    #    ref_value = a + b / R
    #    return ref_value + (x - 10.)
    #if (x >= 10.96):
    #    return 1000.
	  */
        
    return a + b / R_v;

}

#endif
