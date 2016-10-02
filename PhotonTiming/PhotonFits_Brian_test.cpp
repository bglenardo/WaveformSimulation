///////////////////////////////////////////////////////////////////////////////
//
//    2015-07-13  Implemented timing constraints on photon times to be within
//                the boundaries of the waveform. Did it work? (Brian L.)
//
//    2015-08-04  Tried to clean up the code a bit. Reorganized the fitting 
//                routines so that they're only written once. Also, changed
//                name of Stats() function to FindPeaks(), which is more
//                accurate. The numbering of all the results should also be
//                fixed, but I can leave that to a later date. (Brian L.)
//   
//
//   2015-08-05   Note: Currently parameter stats[12] doesn't actually refer
//                to anything, and instead stats[13] is the location of
//                peak 4, (see line 59). (Rose)
//
//   2015-09-01   Fixed the parameter constraints in the FCN_S1 function.
//                Added a FCN_BL to allow baseline fits.
//
//   2015-09-17   Added FindSortedPeaks() which just sorts the peaks from 
//                smallest to largest.
//
//   2015-10-21   Set initializations so that multi-photon spikes have the 
//                additinal photons set *before* the peak, rather than after.
//                Due to the long tail, this is what we expect.
//
//   2015-10-25   Set area constraints to be on fitted values rather than
//                actual values. This includes changing it so that weight
//                applied to the Bayes factor is due to the integral of 
//                the DPH PDF from fit-error to fit+error, divided by the error.
//                This is an approximation of the read Bayes factor.
//
//   2015-10-25   Changed fit to run mnseek and THEN Migrad. Also, reduced 
//                the step size of the fitting routines.
//
//////////////////////////////////////////////////////////////////////////////


#include <vector>
#include <cmath>
#include <iostream>
#include <TMinuit.h>
#include "TF1.h"
#include "TMath.h"
#include "SPEFunctions_Rose_test.h"
#include "PhotonFits.h"
#include "PhotonAreaDistributions.h"
#include "SortIndices.h"

#define PI 3.14159265358979323   

#define _DEBUG_ 0

using namespace std;

////////////////////////////////////////////////////////////////////////////////
// Global variables 
////////////////////////////////////////////////////////////////////////////////
vector<double>::iterator pulse_start;
vector<double>::iterator pulse_end;
TMinuit *minuit = new TMinuit(10);
double negLogLikelihood = 0;
double gerror = 1;
double minimizer_small = 1e-10;


//////////////////////////////////////////////////////////////////////////////
//
//  FindSortedPeaks()
//
//////////////////////////////////////////////////////////////////////////////

void FindSortedPeaks(double * stats) {
  // calculates basic stats quickly
  // [0] max
  // [1] max2
  // [2] max3 
  // [3] max4
  // [4] max5 
  // [5] location of peak 1 
  // [6] location of peak 2
  // [7] location of peak 3
  // [8] location of peak 4
  // [9] location of peak 5
  // [10] num_spikes
  // [11] area
  stats[0] = stats[1] = stats[2] = stats[3] = 0.;
  stats[4] = stats[5] = stats[6] = stats[7] = 0.;
  stats[8] = stats[9] = stats[10] = 0.;
  stats[11] = 0.;

  double maximums[50] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  double max_positions[50] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
                          0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  int num_spikes=0;
  int pulse_ind=0; 

  int counter = 0;
 
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++) {

    if( it > pulse_start && it < pulse_end-1 ){
      vector<double>::iterator previous = it, next = it;
      previous--;
      next++;
      if(*previous < *it && *next <= *it && *it > 0.1) { 
        num_spikes += 1;
        maximums[counter] = *it;
   //     cerr << "Maximum found at: " << *it << endl;
        max_positions[counter] = (double) pulse_ind;
        counter++;
      } // end if(*previous..
    } // end if( it >...
    stats[11] += *it;
    pulse_ind++;
  } // end for (vector<double>::iterator...


  int indices[5] = {0,1,2,3,4};
  if(num_spikes < 5)
    std::sort(indices,indices+num_spikes,sort_indices(maximums));
  if(num_spikes >= 5)
    std::sort(indices,indices+5,sort_indices(maximums));
  stats[0] = maximums[indices[0]];
  stats[1] = maximums[indices[1]];
  stats[2] = maximums[indices[2]];
  stats[3] = maximums[indices[3]];
  stats[4] = maximums[indices[4]];

  //cerr << "Maximums: " << endl;
  //cerr << stats[0] << "\t" << stats[1] << "\t" << stats[2] << "\t" 
  //     << stats[3] << "\t" << stats[4] << endl;

  stats[5] = max_positions[indices[0]];
  stats[6] = max_positions[indices[1]];
  stats[7] = max_positions[indices[2]];
  stats[8] = max_positions[indices[3]];
  stats[9] = max_positions[indices[4]];

  //cerr << stats[5] << "\t" << stats[6] << "\t" << stats[7] << "\t"
 //      << stats[8] << "\t" << stats[9] << endl;
 
  stats[10] = num_spikes;

  
} // end function


////////////////////////////////////////////////////////////////////////////////
// TMinuit Functions
////////////////////////////////////////////////////////////////////////////////
void FCN1(int &npar, double *grad, double &fval, double *par, int flag) {
  // This function calculates the negative log likelihood. It is passed to a minimizer 
  // (TMinuit) to maximize the likelihood. 
  // npar = number of parameters (defined in TMinuit)
  // grad = first derivitives of the function (optional). Not used. 
  // fval = Chisquare value 
  // par  = parameters. 
  // flag = TMinuit flag corresponding to its options 

  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  double x[1];
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++) {
    x[0] = (double) i;
    theory = singleSPE(x, par);
    error = sqrt(gerror*gerror)+minimizer_small;
    fval += (theory - *it)*(theory - *it)/(2*error*error);// + TMath::Log(2*PI*error*error)/2;
    i++;
  }
//  cerr << "fval = " << fval << endl;
  if (!fval) fval = 1e100;
  if (par[1] > (double) i || par[1] < 0.) fval = 1e100;
  negLogLikelihood = fval;
}

////////////////////////////////////////////////////////////////////////////////
void FCN2(int &npar, double *grad, double &fval, double *par, int flag) {
  // This function calculates the negative log likelihood. It is passed to a minimizer 
  // (TMinuit) to maximize the likelihood. 
  // npar = number of parameters (defined in TMinuit)
  // grad = first derivitives of the function (optional). Not used. 
  // fval = Chisquare value 
  // par  = parameters. The first parameter (par[0]) is set to a constant and is
  //        used to choose which function to fit to, ExpFunc or GausFunc.
  // flag = TMinuit flag corresponding to its options 
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  double x[1];
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++) {
    x[0] = (double) i;
    theory = doubleSPE(x, par);
    error = sqrt(gerror*gerror)+minimizer_small;
    
    fval += (theory - *it)*(theory - *it)/(2*error*error); // + TMath::Log(2*PI*error*error)/2
    i++;
  }
  if (!fval) fval = 1e100;
  if ((par[1] > (double) i || par[1] < 0.) &&
      (par[3] > (double) i || par[3] < 0.)) 
         fval = 1e100;
  negLogLikelihood = fval;
}

////////////////////////////////////////////////////////////////////////////////
void FCN3(int &npar, double *grad, double &fval, double *par, int flag) {
  // This function calculates the negative log likelihood. It is passed to a minimizer 
  // (TMinuit) to maximize the likelihood. 
  // npar = number of parameters (defined in TMinuit)
  // grad = first derivitives of the function (optional). Not used. 
  // fval = Chisquare value 
  // par  = parameters. The first parameter (par[0]) is set to a constant and is
  //        used to choose which function to fit to, ExpFunc or GausFunc.
  // flag = TMinuit flag corresponding to its options 
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  double x[1];
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++) {
    x[0] = (double) i;
    theory = tripleSPE(x, par);
    error = sqrt(gerror*gerror)+minimizer_small;
    
    fval += (theory - *it)*(theory - *it)/(2*error*error); // + TMath::Log(2*PI*error*error)/2
    i++;
  }
  if (!fval) fval = 1e100;
  if ((par[1] > (double) i || par[1] < 0.) &&
      (par[3] > (double) i || par[3] < 0.) &&
      (par[5] > (double) i || par[5] < 0.)) 
         fval = 1e100;
  negLogLikelihood = fval;
}

//////////////////////////////////////////////////////////////////////////////////
void FCN4(int &npar, double *grad, double &fval, double *par, int flag) {
  // see comments for FCN3
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  double x[1];
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++){
        x[0] = (double) i;
        theory = quadSPE(x, par); 
        error = sqrt(gerror*gerror)+minimizer_small;

        fval += (theory - *it)*(theory - *it)/(2*error*error); // + TMath::Log(2*PI*error*error)/2
        i++;
  }
  if (!fval) fval = 1e100;
  if ((par[1] > (double) i || par[1] < 0.) &&
      (par[3] > (double) i || par[3] < 0.) &&
      (par[5] > (double) i || par[5] < 0.) && 
      (par[7] > (double) i || par[7] < 0.)) 
         fval = 1e100;
  negLogLikelihood = fval;
}
/////////////////////////////////////////////////////////////////////////////////////
void FCN5(int &npar, double *grad, double &fval, double *par, int flag) {
  //see comments for FCN3
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  double x[1];
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++){
        x[0] = (double) i;
        theory = quintSPE(x, par); 
        error = sqrt(gerror*gerror)+minimizer_small;

        fval += (theory - *it)*(theory - *it)/(2*error*error); // + TMath::Log(2*PI*error*error)/2
        i++;
  }
  if (!fval) fval = 1e100;
  if ((par[1] > (double) i || par[1] < 0.) &&
      (par[3] > (double) i || par[3] < 0.) &&
      (par[5] > (double) i || par[5] < 0.) && 
      (par[7] > (double) i || par[7] < 0.) && 
      (par[9] > (double) i || par[9] < 0.)) 
         fval = 1e100;
  negLogLikelihood = fval;
}
/////////////////////////////////////////////////////////////////////////////////////
void FCN_S1(int &npar, double *grad, double &fval, double *par, int flag) {
  //see comments for FCN3
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  double x[1];
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++){
        x[0] = (double) i;
        theory = S1Template(x, par); 
        error = sqrt(gerror*gerror)+minimizer_small;

        fval += (theory - *it)*(theory - *it)/(2*error*error); // + TMath::Log(2*PI*error*error)/2
        i++;
  }
  if (!fval) fval = 1e100;
  if ((par[1] > (double) i || par[1] < 0.))
         fval = 1e100;
  negLogLikelihood = fval;
}
/////////////////////////////////////////////////////////////////////////////////////
void FCN_BL(int &npar, double *grad, double &fval, double *par, int flag) {
  //see comments for FCN3
  
  fval = 0;
  double theory = 0, error;
  size_t i = 0;
  double x[1];
  for (vector<double>::iterator it=pulse_start; it != pulse_end; it++){
        x[0] = (double) i;
        theory = Baseline(x, par); 
        error = sqrt(gerror*gerror)+minimizer_small;

        fval += (theory - *it)*(theory - *it)/(2*error*error); // + TMath::Log(2*PI*error*error)/2
        i++;
  }
  if (!fval) fval = 1e100;
  negLogLikelihood = fval;
}




////////////////////////////////////////////////////////////////////////////////
// Fitting 
////////////////////////////////////////////////////////////////////////////////
void Fit(vector<double> wave, double trace_start, double perror, double *results,double *results_err) {
  // results [0] = number of photons 
  // results [1] = 
  // results [2] = exp_fit_tau_rise_samples 
 
  // Set values to 0
  double temp_results[40], likelihoods[10], occam_factors[10];// bayes_factors[10];  
  
  results[0] = results[1] = results[2] = results[3] = results[4] = 1.e7;
  results[5] = results[6] = results[7] = results[8] = results[9] = 1.e7;
  results[10] = results[11] = 1.e7;
  
  // Set Globals
  gerror = perror;
  pulse_start = wave.begin();
  pulse_end = wave.end();
  double time_err = 0.01;
  double norm_err = 0.02;
  double ph_height_lim = 0.02;
  double ph_height_max = 7.;
  double default_height_init = 0.2;
  double SPEfitNorm = 1.;

  minuit->SetMaxIterations(200);
  minuit->SetPrintLevel(-1); // (-1) = quiet mode 
  
  // get guesses
  double stats[12];
  FindSortedPeaks(stats);
  double max = stats[0];
  double max2 = stats[1];
  double max3 = stats[2]; 
  double max4 = stats[3];
  double max5 = stats[4]; 

  double max_loc = stats[5]+0.8;
  double max_loc2 = stats[6]+0.8;
  double max_loc3 = stats[7]+0.8;
  double max_loc4 = stats[8]+0.8;
  double max_loc5 = stats[9]+0.8;

  double num_spikes = stats[10];
  double area = stats[11];

  // If there are no spikes, return 0 photons.
  if( num_spikes < 0.5 ) {
     results[0] = 0.;
     return;
  }  

  double error[40] = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 
                       0., 0., 0., 0., 0., 0., 0., 0., 0., 0. }; 

  //cerr << "Num spikes: " << num_spikes << endl;
  //cerr << "Maximums: " << endl;
  //cerr << stats[0] << "\t" << stats[1] << "\t" << stats[2] << "\t" 
  //     << stats[3] << "\t" << stats[4] << endl;
  //cerr << stats[5] << "\t" << stats[6] << "\t" << stats[7] << "\t"
  //     << stats[8] << "\t" << stats[9] << endl;
  
//  double first_guess[2] = {stats[0],stats[5]+0.8};
//  for(int a=0; a <(int) end; a++){
//     double x[1] = {(double) a};
//     cerr << a << "\t" <<  wave[0][a] << "\t" <<  singleSPE(x,first_guess) << endl;
//  }



  //-------- Perform the SPE fit ------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------


  if( num_spikes > 1.5 || num_spikes < 0.5 )
     occam_factors[0] = -1.; 

  else { 
  if( _DEBUG_ == 1 ) cerr << "Performing the SPE fit..." << endl;

    // Perform the SPE fit
    minuit->SetFCN(&FCN1);
    minuit->DefineParameter(0, "Normalization", max, norm_err, ph_height_lim, ph_height_max);
    minuit->DefineParameter(1, "t0", max_loc, time_err, 0., (double) wave.size() - 3.);
//    minuit->Migrad();
    minuit->mnseek();
    minuit->GetParameter(0, temp_results[1], error[1]);
    minuit->GetParameter(1, temp_results[2], error[2]);

    minuit->DefineParameter(0, "Normalization", temp_results[1], norm_err, ph_height_lim, ph_height_max);
    minuit->DefineParameter(1, "t0", temp_results[2], time_err, 0., (double) wave.size() - 3.);
    minuit->Migrad();
    minuit->GetParameter(0, temp_results[1], error[1]);
    minuit->GetParameter(1, temp_results[2], error[2]);

    likelihoods[0] = TMath::Exp( -negLogLikelihood );
  //  if( _DEBUG_ == 1 ) cerr << "Best Fit: n1 = " << temp_results[1] << "\t t1 = " << temp_results[2] << endl; 

//    occam_factors[0] = likelihoods[0]*error[1]*SPE_area_PDF(area,channel); 
    occam_factors[0] = likelihoods[0] *
                           SPE_area_CDF_int((temp_results[1]-error[1])*SPEfitNorm,(temp_results[1]+error[1])*SPEfitNorm) / (error[1]*SPEfitNorm); 

    if( _DEBUG_ == 1) cerr << "Occam_factor = " << occam_factors[0] << 
                              "\t Likelihood: " << likelihoods[0] <<
                              "\t error[1] = " << error[1] <<
                              "\t SPE_area_PDF() = " << SPE_area_PDF(area) << endl;
    
  }


  //-------- Perform the DPE fit ------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

  if( num_spikes > 2.5 || num_spikes < 0.5)
     occam_factors[1] = -1.;

  else {

   // Perform the DPE fit
   if( _DEBUG_ == 1 ) cerr << "Performing the DPE fit" << endl;
   minuit->SetFCN(&FCN2);

   if( num_spikes > 1.5 ) {
   minuit->DefineParameter(0,"Norm 1",max, norm_err, ph_height_lim, ph_height_max);
   minuit->DefineParameter(1,"t1",max_loc,time_err,0.,(double) wave.size() - 3.);
   minuit->DefineParameter(2,"Norm 2",max2,norm_err, ph_height_lim,ph_height_max);
   minuit->DefineParameter(3,"t2",max_loc2,time_err,0.,(double) wave.size()- 3.);
   } else {
   double height_init;
   if( max/1.7 <= ph_height_lim )
     height_init = default_height_init;
   else
     height_init = max/1.7;
   minuit->DefineParameter(0,"Norm 1",height_init, norm_err, ph_height_lim, ph_height_max);
   minuit->DefineParameter(1,"t1",max_loc-0.805,time_err,0.,(double) wave.size() - 3.);
   minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max); 
   minuit->DefineParameter(3,"t2",max_loc+0.305,time_err,0.,(double) wave.size()- 3.);
   }

   minuit->mnseek();
   minuit->GetParameter(0, temp_results[3], error[3]);
   minuit->GetParameter(1, temp_results[4], error[4]);
   minuit->GetParameter(2, temp_results[5], error[5]);
   minuit->GetParameter(3, temp_results[6], error[6]);
   // Now run MIGRAD
   minuit->DefineParameter(0,"Norm 1",temp_results[3], norm_err, ph_height_lim, ph_height_max);
   minuit->DefineParameter(1,"t1",temp_results[4],time_err,0.,(double) wave.size()- 3.);
   minuit->DefineParameter(2,"Norm 2",temp_results[5],norm_err, ph_height_lim,ph_height_max);
   minuit->DefineParameter(3,"t2",temp_results[6],time_err,0.,(double) wave.size() -3.);
   minuit->Migrad();
   minuit->GetParameter(0, temp_results[3], error[3]);
   minuit->GetParameter(1, temp_results[4], error[4]);
   minuit->GetParameter(2, temp_results[5], error[5]);
   minuit->GetParameter(3, temp_results[6], error[6]);
  
   likelihoods[1] = TMath::Exp(-negLogLikelihood);

//   occam_factors[1] = likelihoods[1]*error[3]*error[5]*DPE_area_PDF(area,channel);
   occam_factors[1] = likelihoods[1] *
                           (SPE_area_CDF_int((temp_results[3]-error[3])*SPEfitNorm,(temp_results[3]+error[3])*SPEfitNorm) / (error[3]*SPEfitNorm) +
                           SPE_area_CDF_int((temp_results[5]-error[5])*SPEfitNorm,(temp_results[5]+error[5])*SPEfitNorm) / (error[5]*SPEfitNorm)); 
   //occam_factors[1] = likelihoods[1]*DPE_area_PDF(area,channel);

    if( _DEBUG_ == 1) cerr << "Occam_factor = " << occam_factors[1] << 
                              "\t Likelihood: " << likelihoods[1] <<
                              "\t error[3] = " << error[3] <<
                              "\t error[5] = " << error[5] <<
                              "\t DPE_area_PDF() = " << DPE_area_PDF(area) << endl;
   }

  //-------- Perform the TPE fit ------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

     if( num_spikes > 3.5 || num_spikes < 0.5 )
        occam_factors[2] = -1.;

     else {
     // Perform the TPE fit
     if( _DEBUG_ == 1 ) cerr << "Performing the TPE fit" << endl;
     minuit->SetFCN(&FCN3);

     if( num_spikes < 1.5 && num_spikes > 0.5 ) {
       double height_init;
       if(max/2.4 <= ph_height_lim )
         height_init = default_height_init;
       else
         height_init = max/2.4;
     minuit->DefineParameter(0,"Norm 1",height_init, norm_err, ph_height_lim, ph_height_max);
     minuit->DefineParameter(1,"t1",max_loc-1.515,time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(3,"t2",max_loc-0.505,time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(4,"Norm 3",height_init,norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(5,"t3",max_loc+0.505,time_err,0.,(double) wave.size()-3.);
     } 
     else if ( stats[8] < 2.5 && stats[8] > 1.5) {
       double height_init;
       if( max/1.7 <= ph_height_lim)
         height_init = default_height_init;
       else
         height_init = max/1.7;
     minuit->DefineParameter(0,"Norm 1",height_init, norm_err, ph_height_lim, ph_height_max);
     minuit->DefineParameter(1,"t1",max_loc-0.805,time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(3,"t2",max_loc+0.305,time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(4,"Norm 3",max2,norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(5,"t3",max_loc2,time_err,0.,(double) wave.size()-3.);
     } else if ( stats[8] > 2.5) {
     minuit->DefineParameter(0,"Norm 1",max, norm_err, ph_height_lim, ph_height_max);
     minuit->DefineParameter(1,"t1",max_loc,time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(2,"Norm 2",max2,norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(3,"t2",max_loc2,time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(4,"Norm 3",max3,norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(5,"t3",max_loc3,time_err,0.,(double) wave.size()-3.);
    } 
     minuit->mnseek();
     minuit->GetParameter(0, temp_results[7], error[7]);
     minuit->GetParameter(1, temp_results[8], error[8]);
     minuit->GetParameter(2, temp_results[9], error[9]);
     minuit->GetParameter(3, temp_results[10], error[10]);
     minuit->GetParameter(4, temp_results[11], error[11]);
     minuit->GetParameter(5, temp_results[12], error[12]);

     minuit->DefineParameter(0,"Norm 1",temp_results[7], norm_err, ph_height_lim, ph_height_max);
     minuit->DefineParameter(1,"t1",temp_results[8],time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(2,"Norm 2",temp_results[9],norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(3,"t2",temp_results[10],time_err,0.,(double) wave.size()-3.);
     minuit->DefineParameter(4,"Norm 3",temp_results[11],norm_err, ph_height_lim,ph_height_max);
     minuit->DefineParameter(5,"t3",temp_results[12],time_err,0.,(double) wave.size()-3.);
     minuit->Migrad();
     minuit->GetParameter(0, temp_results[7], error[7]);
     minuit->GetParameter(1, temp_results[8], error[8]);
     minuit->GetParameter(2, temp_results[9], error[9]);
     minuit->GetParameter(3, temp_results[10], error[10]);
     minuit->GetParameter(4, temp_results[11], error[11]);
     minuit->GetParameter(5, temp_results[12], error[12]);

     likelihoods[2] = TMath::Exp(-negLogLikelihood);
    occam_factors[2] = likelihoods[2] *
                          (SPE_area_CDF_int((temp_results[7]-error[7])*SPEfitNorm,(temp_results[7]+error[7])*SPEfitNorm) / (error[7]*SPEfitNorm) + 
                           SPE_area_CDF_int((temp_results[9]-error[9])*SPEfitNorm,(temp_results[9]+error[9])*SPEfitNorm) / (error[9]*SPEfitNorm) +
                           SPE_area_CDF_int((temp_results[11]-error[11])*SPEfitNorm,(temp_results[11]+error[11])*SPEfitNorm) / (error[11]*SPEfitNorm)); 
   // occam_factors[2] = likelihoods[2]*error[7]*error[9]*error[11]*TPE_area_PDF(area,channel);
    //occam_factors[2] = likelihoods[2]*TPE_area_PDF(area,channel);

    if( _DEBUG_ == 1) cerr << "Occam_factor = " << occam_factors[2] << 
                              "\t Likelihood: " << likelihoods[2] <<
                              "\t error[7] = " << error[7] <<
                              "\t error[9] = " << error[9] <<
                              "\t error[11] = " << error[11] <<
                              "\t TPE_area_PDF() = " << TPE_area_PDF(area) << endl;
    }


  //-------- Perform the QPE fit ------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

    if( num_spikes > 4.5 || num_spikes < 0.5 )
      occam_factors[3] = -1.;
   
    else {
     if( _DEBUG_ == 1 ) cerr << "Performing the QPE fit" << endl;
    minuit->SetFCN(&FCN4);

    if( num_spikes > 0.5 && num_spikes < 1.5 ) {
         double height_init;
         if( max/3.2 <= ph_height_lim )
           height_init = default_height_init;
         else
           height_init = max/3.2;
       minuit->DefineParameter(0,"Norm 1", height_init, norm_err,ph_height_lim,ph_height_max);
       minuit->DefineParameter(1,"t1",max_loc-2.515,time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
       minuit->DefineParameter(3,"t2",max_loc-1.51,time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(4,"Norm 3",height_init, norm_err, ph_height_lim,ph_height_max);
       minuit->DefineParameter(5,"t3",max_loc-0.505,time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(6,"Norm 4",height_init, norm_err, ph_height_lim, ph_height_max);
       minuit->DefineParameter(7,"t4",max_loc+0.505,time_err,0.,(double) wave.size()-3.);
    } else if ( num_spikes > 1.5 && num_spikes < 2.5 ) {
       if( max > 2 * max2 ) {
            double height_init;
            if( max/2.7 <= ph_height_lim )
              height_init = default_height_init;
            else 
              height_init = max/2.7;
          minuit->DefineParameter(0,"Norm 1", height_init, norm_err,ph_height_lim,ph_height_max);
          minuit->DefineParameter(1,"t1",max_loc-1.515,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(2,"Norm 2", height_init,norm_err, ph_height_lim,ph_height_max);
          minuit->DefineParameter(3,"t2",max_loc-0.505,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(4,"Norm 3", height_init, norm_err, ph_height_lim,ph_height_max);
          minuit->DefineParameter(5,"t3",max_loc+0.505,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(6,"Norm 4",max2, norm_err, ph_height_lim, ph_height_max);
          minuit->DefineParameter(7,"t4",max_loc2,time_err,0.,(double) wave.size()-3.);
       } else {
            double height_init[2];
            if( max/1.7 <= ph_height_lim )
              height_init[0] = default_height_init;
            else
              height_init[0] = max/1.7;
            if( max2/1.7 <= ph_height_lim )
              height_init[1] = default_height_init;
            else
              height_init[1] = max2/1.7;
          minuit->DefineParameter(0,"Norm 1", height_init[0], norm_err,ph_height_lim,ph_height_max);
          minuit->DefineParameter(1,"t1",max_loc-0.805,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(2,"Norm 2", height_init[0],norm_err, ph_height_lim,ph_height_max);
          minuit->DefineParameter(3,"t2",max_loc+0.305,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(4,"Norm 3", height_init[1], norm_err, ph_height_lim,ph_height_max);
          minuit->DefineParameter(5,"t3",max_loc2-0.805,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(6,"Norm 4", height_init[1], norm_err, ph_height_lim, ph_height_max);
          minuit->DefineParameter(7,"t4",max_loc2+0.305,time_err,0.,(double) wave.size()-3.);
       }
    } else if( num_spikes > 2.5 && num_spikes < 3.5 ) {
            double height_init;
            if( max/1.7 <= ph_height_lim )
              height_init = default_height_init;
            else
              height_init = max/1.7;
          minuit->DefineParameter(0,"Norm 1", height_init, norm_err,ph_height_lim,ph_height_max);
          minuit->DefineParameter(1,"t1",max_loc-0.805,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
          minuit->DefineParameter(3,"t2",max_loc+0.305,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(4,"Norm 3",max2, norm_err, ph_height_lim,ph_height_max);
          minuit->DefineParameter(5,"t3",max_loc2,time_err,0.,(double) wave.size()-3.);
          minuit->DefineParameter(6,"Norm 4",max3, norm_err, ph_height_lim, ph_height_max);
          minuit->DefineParameter(7,"t4",max_loc3,time_err,0.,(double) wave.size()-3.);
   } else if( num_spikes > 3.5 ) {
       minuit->DefineParameter(0,"Norm 1", max, norm_err,ph_height_lim,ph_height_max);
       minuit->DefineParameter(1,"t1",max_loc,time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(2,"Norm 2",max2,norm_err, ph_height_lim,ph_height_max);
       minuit->DefineParameter(3,"t2",max_loc2,time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(4,"Norm 3",max3, norm_err, ph_height_lim,ph_height_max);
       minuit->DefineParameter(5,"t3",max_loc3,time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(6,"Norm 4",max4, norm_err, ph_height_lim, ph_height_max);
       minuit->DefineParameter(7,"t4",max_loc4,time_err,0.,(double) wave.size()-3.);
   }

    minuit->mnseek();
    minuit->GetParameter(0, temp_results[13], error[13]);
    minuit->GetParameter(1, temp_results[14], error[14]);
    minuit->GetParameter(2, temp_results[15], error[15]);
    minuit->GetParameter(3, temp_results[16], error[16]);
    minuit->GetParameter(4, temp_results[17], error[17]);
    minuit->GetParameter(5, temp_results[18], error[18]);
    minuit->GetParameter(6, temp_results[19], error[19]);
    minuit->GetParameter(7, temp_results[20], error[20]);

       minuit->DefineParameter(0,"Norm 1", temp_results[13], norm_err,ph_height_lim,ph_height_max);
       minuit->DefineParameter(1,"t1",temp_results[14],time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(2,"Norm 2",temp_results[15],norm_err, ph_height_lim,ph_height_max);
       minuit->DefineParameter(3,"t2",temp_results[16],time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(4,"Norm 3",temp_results[17], norm_err, ph_height_lim,ph_height_max);
       minuit->DefineParameter(5,"t3",temp_results[18],time_err,0.,(double) wave.size()-3.);
       minuit->DefineParameter(6,"Norm 4",temp_results[19], norm_err, ph_height_lim, ph_height_max);
       minuit->DefineParameter(7,"t4",temp_results[20],time_err,0.,(double) wave.size()-3.);
    minuit->Migrad();
    minuit->GetParameter(0, temp_results[13], error[13]);
    minuit->GetParameter(1, temp_results[14], error[14]);
    minuit->GetParameter(2, temp_results[15], error[15]);
    minuit->GetParameter(3, temp_results[16], error[16]);
    minuit->GetParameter(4, temp_results[17], error[17]);
    minuit->GetParameter(5, temp_results[18], error[18]);
    minuit->GetParameter(6, temp_results[19], error[19]);
    minuit->GetParameter(7, temp_results[20], error[20]);
    

    likelihoods[3] = TMath::Exp(-negLogLikelihood);
//    occam_factors[3] = likelihoods[3]*error[13]*error[15]*error[17]*error[19]*QPE_area_PDF(area,channel);
    occam_factors[3] = likelihoods[3] *
                          (SPE_area_CDF_int((temp_results[13]-error[13])*SPEfitNorm,(temp_results[13]+error[13])*SPEfitNorm) / (error[13]*SPEfitNorm) + 
                           SPE_area_CDF_int((temp_results[15]-error[15])*SPEfitNorm,(temp_results[15]+error[15])*SPEfitNorm) / (error[15]*SPEfitNorm) +
                           SPE_area_CDF_int((temp_results[17]-error[17])*SPEfitNorm,(temp_results[17]+error[17])*SPEfitNorm) / (error[17]*SPEfitNorm) +
                           SPE_area_CDF_int((temp_results[19]-error[19])*SPEfitNorm,(temp_results[19]+error[19])*SPEfitNorm) / (error[19]*SPEfitNorm)); 
    //occam_factors[3] = likelihoods[3]*QPE_area_PDF(area,channel);

    if( _DEBUG_ == 1) cerr << "Occam_factor = " << occam_factors[3] << 
                              "\t Likelihood: " << likelihoods[3] <<
                              "\t error[13] = " << error[13] <<
                              "\t error[15] = " << error[15] <<
                              "\t error[17] = " << error[17] <<
                             "\t error[19] = " << error[19] <<
                              "\t QPE_area_PDF() = " << QPE_area_PDF(area) << endl;
    }



  //-------- Perform the QIPE fit ------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

   //Perform the QIPE fit
    
    if( num_spikes < 0.5 )
      occam_factors[4] = -1.;
   
    else {
     if( _DEBUG_ == 1 ) cerr << "Performing the QIPE fit" << endl;
   minuit->SetFCN(&FCN5);

   if( num_spikes > 0.5 && num_spikes < 1.5 ) {
        double height_init;
        if( max/4. <= ph_height_lim )
           height_init = default_height_init;
        else
           height_init = max/4.;
      minuit->DefineParameter(0,"Norm 1",height_init, norm_err, ph_height_lim, ph_height_max);
      minuit->DefineParameter(1,"t1",max_loc-3.52,time_err,0.,(double) wave.size()-3.);
      minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
      minuit->DefineParameter(3,"t2",max_loc-2.515,time_err,0.,(double) wave.size()-3.);
      minuit->DefineParameter(4,"Norm 3",height_init,norm_err, ph_height_lim,ph_height_max);
      minuit->DefineParameter(5,"t3",max_loc-1.510,time_err,0.,(double) wave.size()-3.);
      minuit->DefineParameter(6,"Norm 4",height_init, norm_err, ph_height_lim, ph_height_max);
      minuit->DefineParameter(7,"t4",max_loc-0.505,time_err,0.,(double) wave.size()-3.);
      minuit->DefineParameter(8,"Norm 5",height_init,norm_err, ph_height_lim,ph_height_max);
      minuit->DefineParameter(9,"t5",max_loc+0.505,time_err,0.,(double) wave.size()-3.);
   } else if (num_spikes > 1.5 && num_spikes < 2.5) {
      if (max < 3*max2) {
           double height_init[2];
           if( max/2.7 <= ph_height_lim )
             height_init[0] = default_height_init;
           else
             height_init[0] = max/2.7;
           if( max2/1.7 <= ph_height_lim )
             height_init[1] = default_height_init;
           else 
             height_init[1] = max2/1.7;
         minuit->DefineParameter(0,"Norm 1",height_init[0], norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(1,"t1",max_loc-1.51,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(2,"Norm 2",height_init[0],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(3,"t2",max_loc-0.505,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(4,"Norm 3",height_init[0],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(5,"t3",max_loc+0.505,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(6,"Norm 4",height_init[1], norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(7,"t4",max_loc2-0.805,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(8,"Norm 5",height_init[1],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(9,"t5",max_loc2+0.305,time_err,0.,(double) wave.size()-3.);
      } else {
           double height_init;
           if( max/3.2 <= ph_height_lim )
             height_init = default_height_init;
           else
             height_init = max/3.2;
         minuit->DefineParameter(0,"Norm 1",height_init, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(1,"t1",max_loc-2.52,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(3,"t2",max_loc-1.51,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(4,"Norm 3",height_init,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(5,"t3",max_loc-0.505,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(6,"Norm 4",height_init, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(7,"t4",max_loc+0.505,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(8,"Norm 5",max2,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(9,"t5",max_loc2,time_err,0.,(double) wave.size()-3.);
      }
  } else if( num_spikes > 2.5 && num_spikes < 3.5 ) {
      if(max > 2.*max2) {
           double height_init;
           if( max/2.4 <= ph_height_lim )
             height_init = default_height_init;
           else
             height_init = max/2.4;
         minuit->DefineParameter(0,"Norm 1",height_init, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(1,"t1",max_loc-1.51,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(3,"t2",max_loc-0.505,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(4,"Norm 3",height_init,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(5,"t3",max_loc+0.505,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(6,"Norm 4",max2, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(7,"t4",max_loc2,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(8,"Norm 5",max3,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(9,"t5",max_loc3,time_err,0.,(double) wave.size()-3.);
      } else {
           double height_init[2];
           if( max/1.7 <= ph_height_lim )
             height_init[0] = default_height_init;
           else
             height_init[0] = max/1.7;
           if( max2/1.7 <= ph_height_lim )
             height_init[1] = default_height_init;
           else 
             height_init[1] = max2/1.7;
         minuit->DefineParameter(0,"Norm 1",height_init[0], norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(1,"t1",max_loc-0.805,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(2,"Norm 2",height_init[0],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(3,"t2",max_loc+0.305,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(4,"Norm 3",height_init[1],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(5,"t3",max_loc2-0.805,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(6,"Norm 4",height_init[1], norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(7,"t4",max_loc2+0.305,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(8,"Norm 5",max3,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(9,"t5",max_loc3,time_err,0.,(double) wave.size()-3.);
      }
  } else if(num_spikes > 3.5 && num_spikes < 4.5) {
           double height_init;
           if( max/1.7 <= ph_height_lim )
             height_init = default_height_init;
           else
             height_init = max/1.7;
         minuit->DefineParameter(0,"Norm 1",height_init, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(1,"t1",max_loc-0.805,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(2,"Norm 2",height_init,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(3,"t2",max_loc+0.305,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(4,"Norm 3",max2,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(5,"t3",max_loc2,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(6,"Norm 4",max3, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(7,"t4",max_loc3,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(8,"Norm 5",max4,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(9,"t5",max_loc4,time_err,0.,(double) wave.size()-3.);
  } else if(num_spikes > 4.5) { 
         minuit->DefineParameter(0,"Norm 1",max, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(1,"t1",max_loc,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(2,"Norm 2",max2,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(3,"t2",max_loc2,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(4,"Norm 3",max3,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(5,"t3",max_loc3,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(6,"Norm 4",max4, norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(7,"t4",max_loc4,time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(8,"Norm 5",max5,norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(9,"t5",max_loc5,time_err,0.,(double) wave.size()-3.);
  }
   minuit->mnseek();
   minuit->GetParameter(0, temp_results[21], error[21]);
   minuit->GetParameter(1, temp_results[22], error[22]);
   minuit->GetParameter(2, temp_results[23], error[23]);
   minuit->GetParameter(3, temp_results[24], error[24]);
   minuit->GetParameter(4, temp_results[25], error[25]);
   minuit->GetParameter(5, temp_results[26], error[26]);
   minuit->GetParameter(6, temp_results[27], error[27]);
   minuit->GetParameter(7, temp_results[28], error[28]);
   minuit->GetParameter(8, temp_results[29], error[29]);
   minuit->GetParameter(9, temp_results[30], error[30]);

         minuit->DefineParameter(0,"Norm 1",temp_results[21], norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(1,"t1",temp_results[22],time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(2,"Norm 2",temp_results[23],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(3,"t2",temp_results[24],time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(4,"Norm 3",temp_results[25],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(5,"t3",temp_results[26],time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(6,"Norm 4",temp_results[27], norm_err, ph_height_lim, ph_height_max);
         minuit->DefineParameter(7,"t4",temp_results[28],time_err,0.,(double) wave.size()-3.);
         minuit->DefineParameter(8,"Norm 5",temp_results[29],norm_err, ph_height_lim,ph_height_max);
         minuit->DefineParameter(9,"t5",temp_results[30],time_err,0.,(double) wave.size()-3.);
   minuit->Migrad();
   minuit->GetParameter(0, temp_results[21], error[21]);
   minuit->GetParameter(1, temp_results[22], error[22]);
   minuit->GetParameter(2, temp_results[23], error[23]);
   minuit->GetParameter(3, temp_results[24], error[24]);
   minuit->GetParameter(4, temp_results[25], error[25]);
   minuit->GetParameter(5, temp_results[26], error[26]);
   minuit->GetParameter(6, temp_results[27], error[27]);
   minuit->GetParameter(7, temp_results[28], error[28]);
   minuit->GetParameter(8, temp_results[29], error[29]);
   minuit->GetParameter(9, temp_results[30], error[30]);
   
   likelihoods[4] = TMath::Exp(-negLogLikelihood);
   //occam_factors[4] = likelihoods[4]*error[21]*error[23]*error[25]*error[27]*error[29]*QIPE_area_PDF(area,channel);
   occam_factors[4] = likelihoods[4] *
                          (SPE_area_CDF_int((temp_results[21]-error[21])*SPEfitNorm,(temp_results[21]+error[21])*SPEfitNorm) / (error[21]*SPEfitNorm) + 
                           SPE_area_CDF_int((temp_results[23]-error[23])*SPEfitNorm,(temp_results[23]+error[23])*SPEfitNorm) / (error[23]*SPEfitNorm) +
                           SPE_area_CDF_int((temp_results[25]-error[25])*SPEfitNorm,(temp_results[25]+error[25])*SPEfitNorm) / (error[25]*SPEfitNorm) +
                           SPE_area_CDF_int((temp_results[27]-error[27])*SPEfitNorm,(temp_results[27]+error[27])*SPEfitNorm) / (error[27]*SPEfitNorm) + 
                           SPE_area_CDF_int((temp_results[29]-error[29])*SPEfitNorm,(temp_results[29]+error[29])*SPEfitNorm) / (error[29]*SPEfitNorm)); 
   //occam_factors[4] = likelihoods[4]*QIPE_area_PDF(area,channel);
    if( _DEBUG_ == 1) cerr << "Occam_factor = " << occam_factors[4] << 
                              "\t Likelihood: " << likelihoods[4] <<
                              "\t error[31] = " << error[31] <<
                              "\t QIPE_area_PDF() = " << QIPE_area_PDF(area) << endl;



  }

  
  //-------- Perform the Baseline fit ------------------------------------------------------------------
  //-----------------------------------------------------------------------------------------------

  if( _DEBUG_ == 1 ) cerr << "Performing the baseline fit..." << endl;

  if( num_spikes > 0.5 )
     occam_factors[5] = -1.; 

  else { 

    // Perform the SPE fit
    minuit->SetFCN(&FCN_BL);
    minuit->DefineParameter(0, "Normalization", 0., 0.01, -10., 10.);
    if( _DEBUG_ == 1 ) cerr << "Initial values: n1 = " << max << "\t t1 = " << stats[9] << endl; 
//    minuit->Migrad();
    minuit->mnseek();
    minuit->GetParameter(0, temp_results[31], error[31]);
    likelihoods[5] = TMath::Exp( -negLogLikelihood );
    if( _DEBUG_ == 1 ) cerr << "Best Fit: norm = " << temp_results[31] << endl;

    //occam_factors[5] = likelihoods[5]*error[31]*Baseline_area_PDF(area); 
    occam_factors[5] = likelihoods[5]*Baseline_area_PDF(area); 

    if( _DEBUG_ == 1) cerr << "Occam_factor = " << occam_factors[5] << 
                              "\t Likelihood: " << likelihoods[5] <<
                              "\t error[31] = " << error[31] << endl;
    
  }






   // NEED TO MAKE THE DECISION OF WHICH TO USE  

     // This block implements the too-close cuts, the minimum height cuts, and the 
     // waveform boundary cuts
     double maximum = -1.;// occam_factors[5];
     bool photons_separated[5] = {false,false,false,false,false};
     //bool in_bounds[5] = {false,false,false,false,false};
     bool in_bounds[5] = {true, true, true, true, true};
     bool tall_enough[5] = {false,false,false,false,false};
    // bool photons_separated[5] = {true, true, true, true, true};

     for (int k=0; k<5; k++) {
       if( k == 0 && num_spikes > 0.5 && num_spikes < 1.5){
          photons_separated[k] = true;
          if( temp_results[1] > 0.05 )
              tall_enough[0] = true;
       }
       if( k == 1 && num_spikes > 0.5 && num_spikes < 2.5){
         if( TMath::Abs(temp_results[4]-temp_results[6]) > 1.)
             photons_separated[1] = true;
         if( temp_results[3] > 0.05 && temp_results[5] > 0.05 )
             tall_enough[1] = true;
       }
       if( k == 2 && num_spikes > 0.5 && num_spikes < 3.5){
         if( TMath::Abs(temp_results[8]-temp_results[10]) > 1. &&
             TMath::Abs(temp_results[10]-temp_results[12]) > 1. &&
             TMath::Abs(temp_results[8]-temp_results[12]) > 1)
              photons_separated[k] = true;
         if( temp_results[7] > 0.05 && temp_results[9] > 0.05 && 
                   temp_results[11] > 0.05 )
              tall_enough[2] = true;
       }
       if( k == 3 && num_spikes > 0.5 && num_spikes < 4.5){  
         if( TMath::Abs(temp_results[14]-temp_results[16]) > 1. &&
             TMath::Abs(temp_results[14]-temp_results[18]) > 1. &&
             TMath::Abs(temp_results[14]-temp_results[20]) > 1. && 
             TMath::Abs(temp_results[16]-temp_results[18]) > 1. && 
             TMath::Abs(temp_results[16]-temp_results[20]) > 1. && 
             TMath::Abs(temp_results[18]-temp_results[20]) > 1.)
          photons_separated[k] = true; 
         if( temp_results[13] > 0.05 && temp_results[15] > 0.05 &&
             temp_results[17] > 0.05 && temp_results[19] > 0.05)
          tall_enough[3] = true;
       }
       if( k == 4 && num_spikes > 0.5 && num_spikes < 5.5) {
         if( TMath::Abs(temp_results[22]-temp_results[24]) > 1. &&
             TMath::Abs(temp_results[22]-temp_results[26]) > 1. &&
             TMath::Abs(temp_results[22]-temp_results[28]) > 1. && 
             TMath::Abs(temp_results[22]-temp_results[30]) > 1. && 
             TMath::Abs(temp_results[24]-temp_results[26]) > 1. && 
             TMath::Abs(temp_results[24]-temp_results[28]) > 1. && 
             TMath::Abs(temp_results[24]-temp_results[30]) > 1. && 
             TMath::Abs(temp_results[26]-temp_results[28]) > 1. && 
             TMath::Abs(temp_results[26]-temp_results[30]) > 1. && 
             TMath::Abs(temp_results[28]-temp_results[30]) > 1.)
          photons_separated[k] = true; 
         if( temp_results[21] > 0.05 && temp_results[23] > 0.05 && 
             temp_results[25] > 0.05 && temp_results[27] > 0.05 &&
             temp_results[29] > 0.05 )
          tall_enough[k] = true;
      }

        if (maximum < occam_factors[k] && photons_separated[k] && tall_enough[k] && in_bounds[k])
           maximum = occam_factors[k];
     }


     if( maximum < 0. || occam_factors[5]==maximum ) {
        if( _DEBUG_ == 1) cerr << "No fit - wave below threshold or it's best fit by baseline."
                               << endl;
        results[0] = 0.;
        return;
     } 
     else if( occam_factors[0]==maximum ){
     results[0] = 1.;
     results[1] = temp_results[1]*SPEfitNorm;
     results[2] = temp_results[2] + trace_start;
     results_err[2] = error[2];
     }
     else if ( occam_factors[1]==maximum ) {
     results[0] = 2.;
     results[1] = temp_results[3]*SPEfitNorm;
     results[2] = temp_results[4] + trace_start;
     results[3] = temp_results[5]*SPEfitNorm;
     results[4] = temp_results[6] + trace_start;     
     results_err[2] = error[4];
     results_err[4] = error[6];
     } 
     else if (occam_factors[2]==maximum)  {
     results[0] = 3.;
     results[1] = temp_results[7]*SPEfitNorm;
     results[2] = temp_results[8] + trace_start;
     results[3] = temp_results[9]*SPEfitNorm;
     results[4] = temp_results[10] + trace_start;
     results[5] = temp_results[11]*SPEfitNorm;
     results[6] = temp_results[12] + trace_start;   
     results_err[2] = error[8];
     results_err[4] = error[10];
     results_err[6] = error[12]; 
     } 
     else if (occam_factors[3]==maximum ) {
     results[0] = 4.;
     results[1] = temp_results[13]*SPEfitNorm;
     results[2] = temp_results[14] + trace_start;
     results[3] = temp_results[15]*SPEfitNorm;
     results[4] = temp_results[16] + trace_start;
     results[5] = temp_results[17]*SPEfitNorm;
     results[6] = temp_results[18] + trace_start;
     results[7] = temp_results[19]*SPEfitNorm;
     results[8] = temp_results[20] + trace_start;
     results_err[2] = error[14];
     results_err[4] = error[16];
     results_err[6] = error[18];
     results_err[8] = error[20];
     }  
     else {
     results[0] = 5.;
     results[1] = temp_results[21]*SPEfitNorm;
     results[2] = temp_results[22] + trace_start;
     results[3] = temp_results[23]*SPEfitNorm;
     results[4] = temp_results[24] + trace_start;
     results[5] = temp_results[25]*SPEfitNorm;
     results[6] = temp_results[26] + trace_start;
     results[7] = temp_results[27]*SPEfitNorm;
     results[8] = temp_results[28] + trace_start;
     results[9] = temp_results[29]*SPEfitNorm;
     results[10] = temp_results[30] + trace_start;
     results_err[2] = error[22];
     results_err[4] = error[24];
     results_err[6] = error[26];
     results_err[8] = error[28];
     results_err[10] = error[30];
     }

     if( _DEBUG_ ==1 ) {
        cerr << "Photons found: " << results[0] << endl;

       if( results[0] > 3.5 && results[0] < 5.5) {
          int kk = 0;
          for (vector<double>::iterator it=pulse_start; it != pulse_end; it++){
             cerr << max << "\t " <<   kk << "\t " << *it << endl;
             kk++;
          }
       } 
     }


return;

}


void FitFullS1(vector<double> *wave, double perror, double *results,double *results_err) {


    double stats[20], temp_results[2], error[2], likelihoods[2];
    FindSortedPeaks(stats);
    
    minuit->SetFCN(&FCN_S1);
    minuit->DefineParameter(0, "Normalization", stats[0], 0.01, 0.001, 1000.);
    minuit->DefineParameter(1, "t0", stats[9], 0.01, 0., (double) wave->size()-3.);
    if( _DEBUG_ == 1 ) cerr << "Initial values: n1 = " << stats[0] << "\t t1 = " << stats[9] << endl; 
   // minuit->Migrad();
    minuit->mnseek();
    minuit->GetParameter(0, temp_results[0], error[0]);
    minuit->GetParameter(1, temp_results[1], error[1]);
    likelihoods[0] = TMath::Exp( -negLogLikelihood );


    results[0] = temp_results[0];
    results[1] = temp_results[1]; 

    results_err[0] = error[0];
    results_err[1] = error[1];

  return;
}




