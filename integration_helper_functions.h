#include <math.h>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_ntuple.h>
#include "physical_constants.h"
#include "galaxy.h"
#include "disruption.h"
#include "histogramNd.h"
#include "survey.h"

// move this to where appropriate. 
struct flare_data
{
  double attributes[28]; // 17 for galaxy, 11 for flare
  double weight; // for determining volumetric disruption rate
};


void Sample_Disruption_Parameters(gsl_rng *rangen, Survey* surv, Galaxy gal, double& vol_rate_accumulator, gsl_ntuple *flare_ntuple, flare_data *flare_row )
{

  double m_r_limit_contrast = surv->Find_Host_Contrast_Magnitude(gal,'r');
  double m_g_limit_contrast = surv->Find_Host_Contrast_Magnitude(gal,'g'); 

  // for all flares in this galaxy
  double operating_m_r_limit = std::min(m_r_limit_contrast,surv->Get_m_r_Threshhold());
  double operating_m_g_limit = std::min(m_g_limit_contrast,surv->Get_m_g_Threshhold());

  double mbh = gal.Get_Mbh();
  double z = gal.Get_z();

  double cosmo_factor = (1. + z)/(4. * PI * pow(gal.Get_Luminosity_Distance(),2.));

  const int Nbands = surv->Get_Nbands();

  double * nu_emit = new double[Nbands];
  double * lnu_obsc_rest = new double[Nbands];

  // these are used to see if the event passes the survey flux limit
  // call in the enum, rather than using it globally?
  double nu_g_emit = (1. + z) * surv->Get_Band_Nu(ZTF_g);
  double nu_r_emit = (1. + z) * surv->Get_Band_Nu(ZTF_r);

  // save these for later for the temperature fit
  nu_emit[0] = nu_r_emit;
  nu_emit[1] = nu_g_emit;
  nu_emit[2] = (1. + z) * surv->Get_Band_Nu(UVW1);
  nu_emit[3] = (1. + z) * surv->Get_Band_Nu(UVM2);
  nu_emit[4] = (1. + z) * surv->Get_Band_Nu(UVW2);

  int num_trials = 100;
  vol_rate_accumulator = 0.;

  //  double rate_normalization = 1./( (double) num_trials) * gal.Get_Disruption_Rate_Normalization_Combined() * pow(gal.Get_nuker_gammaprime()/1.0,gal.Get_Disruption_Rate_Powerlaw_Nuker()) * 1./(1. + z);

  double rate_normalization = 1./( (double) num_trials) * gal.Get_Disruption_Rate_Normalization_Combined() * pow(mbh/1.e8,gal.Get_Disruption_Rate_Powerlaw_Mass()) * 1./(1. + z);

  Galaxy* gal_pointer = &gal;
  Disruption disrupt(gal_pointer); // default values filled in now
  //  Disruption disrupt(gal); // default values filled in now

  for (int i = 0; i < num_trials; ++i)
    {
      // sample mstar
      disrupt.Rejection_Sample_Mstar(rangen);
      double mhills = disrupt.Get_Hills_Mass();

      if (mbh > mhills) continue;
      else
	{
	  // will then have to do the flare observability criteria for each flare when doing detected disrupt. For volumteric disrupt, accept all of these

	  vol_rate_accumulator += 1.;

	  disrupt.Sample_Beta(rangen);
	  
	  // need to have sampled mstar, beta (and mbh) already 
	  disrupt.Determine_Max_L();
      	  double max_L = disrupt.Get_Max_L();  // will also be used to convert the sampled x to a phsyical L
	  if (log10(max_L) < disrupt.Get_Min_Log_Lbol()) continue;

	  disrupt.Sample_Peak_L(rangen);
	  disrupt.Sample_Topt(rangen);
	  disrupt.Sample_R_V(rangen);
	  disrupt.Sample_A_V(rangen);
	    
	  double r_mag_observed = surv->mAB_From_Fnu(disrupt.Extincted_Flux_Observed(nu_r_emit,cosmo_factor));
	  double g_mag_observed = surv->mAB_From_Fnu(disrupt.Extincted_Flux_Observed(nu_g_emit,cosmo_factor));

	  if (r_mag_observed < operating_m_r_limit && g_mag_observed < operating_m_g_limit && (g_mag_observed - r_mag_observed) < 0.)
	    {

	      /////// temperature fit ///////
	      // collect the values from disruption::Unobscured_Lnu in an array
	      for (int nu_index = 0; nu_index < Nbands; nu_index++)
		{
		  lnu_obsc_rest[nu_index] = disrupt.Unobscured_Lnu(nu_emit[nu_index],disrupt.Get_Topt(),disrupt.Get_Peak_L()) * disrupt.Dust_Flux_Factor_Reduction(nu_emit[nu_index]); // consider storing the flux factor reduction (and maybe enven the unobscured lnu) from earlier so you don't need to repeat that here.
		}

	      int status = surv->Provide_Temperature_Fit_Data(nu_emit,lnu_obsc_rest, rangen); // and perform the fit
	      if (status != GSL_SUCCESS) printf("ERROR: Temperature fit failed!\n");

	      //could use your enum(s)
	      flare_row->attributes[0] = log10(gal.Get_Total_Stellar_Mass());
	      flare_row->attributes[1] = log10(gal.Get_Mbh());
	      flare_row->attributes[2] = log10(gal.Get_Mbh_bulge());
	      flare_row->attributes[3] = gal.Get_z();
	      flare_row->attributes[4] = gal.Get_sersic_n();
	      flare_row->attributes[5] = gal.Get_r50_kpc();
	      flare_row->attributes[6] = gal.Get_m_g();
	      flare_row->attributes[7] = gal.Get_m_r();
	      flare_row->attributes[8] = log10(gal.Get_ssfr());
	      flare_row->attributes[9] = gal.Get_M_u() - gal.Get_M_r();
	      flare_row->attributes[10] = gal.Get_M_r();
	      flare_row->attributes[11] = gal.Get_nuker_gammaprime();
	      flare_row->attributes[12] = gal.Get_median_A_V();
	      flare_row->attributes[13] = gal.Get_sigma_A_V();
	      flare_row->attributes[14] = gal.Get_Mu_Eff(gal.Get_m_g());
	      flare_row->attributes[15] = gal.Get_Mu_Central(gal.Get_m_g());
	      flare_row->attributes[16] = gal.Get_Mu_At_Rkpc_CosmoCorrected(0.1,gal.Get_m_g());

	      flare_row->attributes[17] = disrupt.Get_Mstar();
	      flare_row->attributes[18] = disrupt.Get_beta();
	      flare_row->attributes[19] = log10(disrupt.Get_Peak_L());
	      flare_row->attributes[20] = disrupt.Get_Topt();
	      flare_row->attributes[21] = g_mag_observed;
	      flare_row->attributes[22] = g_mag_observed - r_mag_observed;
	      flare_row->attributes[23] = disrupt.Get_A_V();
	      flare_row->attributes[24] = surv->Get_Tbb_Fit();
	      flare_row->attributes[25] = log10(surv->Get_Lbol_Fit());
	      flare_row->attributes[26] = log10(surv->Get_Rbb_Fit()); // in the future, improve  capability to bin functions of columns so you don't need to store extra columns
	      flare_row->attributes[27] = log10(disrupt.Get_Peak_L() / disrupt.Get_Eddington_Luminosity());
	      			    
	      flare_row->weight = rate_normalization;

	      gsl_ntuple_write(flare_ntuple);

	    }	    

	}

    }

    // this can be done at the end to save time
  // normalize by number of trials
  vol_rate_accumulator *=  rate_normalization;

  delete[] nu_emit;
  delete[] lnu_obsc_rest;

}

