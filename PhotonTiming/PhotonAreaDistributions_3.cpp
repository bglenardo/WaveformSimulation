///////////////////////////////////////////////////////////////////////////////////////
//
//   PhotonAreaDistributions.cpp
//
//   2015-07-14 - Initial submission. (Rose)
//
//   2015-07-31 - Added "delete" commands for the TF1's that we're creating, so they
//                don't hang around in memory. (Brian L.)
//   2015-09-02 - Added a Baseline_area_PDF() to constrain baseline fits.
//
///////////////////////////////////////////////////////////////////////////////////////



#include <cmath>
#include <iostream>
#include <vector>
#include "SPEparameters.h"
#include "TF1.h"
#include "TMath.h"

using namespace std;


// Global variables
double spe_mu = 0.8327873877969729;
double spe_sig = 0.35;
double sqrt2 = TMath::Sqrt(2);
double pi = TMath::Pi();
double sqrt2pi = TMath::Sqrt(2*pi);

/////////////////////////////////////////////////////////////////////////////////
//
//  Baseline_area_PDF (Rose)
//
/////////////////////////////////////////////////////////////////////////////////

double Baseline_area_PDF(double peak_area_phe) {

   TF1 * baselineFit = new TF1("baselineFit","gaus(0) + gaus(3)",-10.,100.);
   baselineFit->SetParameters(1.3437,-0.002051,0.1721,0.577,0.03223,0.2903);
   double y = baselineFit->Eval(peak_area_phe);
   delete baselineFit;
   return y;

}





/////////////////////////////////////////////////////////////////////////////////
//
//  SPE_area_PDF (Rose)
//
/////////////////////////////////////////////////////////////////////////////////
double SPE_area_PDF(double peak_area_phe) {

  double y = 0.8/spe_sig/sqrt2pi*TMath::Exp( -TMath::Power( (peak_area_phe-spe_mu)/spe_sig, 2.)*0.5 ) + 
      0.2/spe_sig/sqrt2pi/sqrt2*TMath::Exp( -TMath::Power( (peak_area_phe-2*spe_mu)/(sqrt2 * spe_sig), 2. )*0.5 );	

  return y;
	
}


/////////////////////////////////////////////////////////////////////////////////
//
//  SPE_area_CDF (Brian)
//
/////////////////////////////////////////////////////////////////////////////////

double SPE_area_CDF(double x) {

   double y = 
      0.8/spe_sig/sqrt2pi /(2*sqrt2/sqrt2pi) * sqrt2 * spe_sig *
      (TMath::Erf( (x - spe_mu)/spe_sig ) + TMath::Erf( spe_mu/spe_sig ) ) + 
      
      0.2/spe_sig/sqrt2pi/sqrt2 /(2*sqrt2/sqrt2pi) * 2. * spe_sig *
      (TMath::Erf( (x - 2*spe_mu)/(sqrt2*spe_sig) ) + TMath::Erf( 2*spe_mu/spe_sig/sqrt2 ) );

   return y;
}


/////////////////////////////////////////////////////////////////////////////////
//
//  SPE_area_PDF_int (Brian)
//
/////////////////////////////////////////////////////////////////////////////////

double SPE_area_CDF_int(double a, double b) {
 
        return SPE_area_CDF(b) - SPE_area_CDF(a);
}


/////////////////////////////////////////////////////////////////////////////////
//
//  DPE_area_PDF (Rose)
//
/////////////////////////////////////////////////////////////////////////////////
double DPE_area_PDF(double peak_area_phe) {

  double y =  
  0.8*0.8 /spe_sig/sqrt2pi/sqrt2 * 
        TMath::Exp( -TMath::Power( (peak_area_phe-2*spe_mu)/(sqrt2*spe_sig), 2.)*0.5 ) + 
  2 * 0.8 * 0.2 /spe_sig/sqrt2pi/TMath::Sqrt(3) *
        TMath::Exp( - TMath::Power( (peak_area_phe-3*spe_mu)/(TMath::Sqrt(3)*spe_sig), 2.)*0.5 ) + 
  0.2*0.2 /spe_sig/sqrt2pi/2.  *
        TMath::Exp( - TMath::Power( (peak_area_phe - 4*spe_mu)/(4*spe_sig), 2.)*0.5);

  return y;

}

/////////////////////////////////////////////////////////////////////////////////
//
//  TPE_area_PDF (Rose)
//
/////////////////////////////////////////////////////////////////////////////////
double TPE_area_PDF(double peak_area_phe) {

  double y = 
  0.8*0.8*0.8 /spe_sig/sqrt2pi/TMath::Sqrt(3) * 
        TMath::Exp( -TMath::Power( (peak_area_phe-3*spe_mu)/(TMath::Sqrt(3)*spe_sig), 2.)*0.5 ) + 
  3 * 0.8 * 0.8 * 0.2 /spe_sig/sqrt2pi/TMath::Sqrt(4) *
        TMath::Exp( -TMath::Power( (peak_area_phe-4*spe_mu)/(TMath::Sqrt(4)*spe_sig), 2.)*0.5 ) + 
  3 * 0.8 * 0.2 * 0.2 / spe_sig/sqrt2pi/TMath::Sqrt(5) *
        TMath::Exp( -TMath::Power( (peak_area_phe-5*spe_mu)/(TMath::Sqrt(5)*spe_sig), 2.)*0.5 ) +
  0.2*0.2*0.2 / spe_sig/sqrt2pi/TMath::Sqrt(6) *
        TMath::Exp( -TMath::Power( (peak_area_phe-6*spe_mu)/(TMath::Sqrt(6)*spe_sig), 2.)*0.5 );
 
  return y; 


}


/////////////////////////////////////////////////////////////////////////////////
//
//  QPE_area_PDF (Rose)
//
/////////////////////////////////////////////////////////////////////////////////
double QPE_area_PDF(double peak_area_phe, int channel) {

  double y = 
  0.8*0.8*0.8*0.8 / spe_sig/sqrt2pi/TMath::Sqrt(4.) *
     TMath::Exp( -TMath::Power( (peak_area_phe-4*spe_mu)/(TMath::Sqrt(4)*spe_sig), 2.)*0.5 ) +
  4 * 0.8 * 0.8 * 0.8 * 0.2/spe_sig/sqrt2pi/TMath::Sqrt(5.) *
     TMath::Exp( -TMath::Power( (peak_area_phe-5*spe_mu)/(TMath::Sqrt(5)*spe_sig), 2.)*0.5 ) +
  6 * 0.8 * 0.8 * 0.2 * 0.2/spe_sig/sqrt2pi/TMath::Sqrt(6.) *
     TMath::Exp( -TMath::Power( (peak_area_phe-6*spe_mu)/(TMath::Sqrt(6)*spe_sig), 2.)*0.5 ) +
  4 * 0.8 * 0.2 * 0.2 * 0.2/spe_sig/sqrt2pi/TMath::Sqrt(7.) *
     TMath::Exp( -TMath::Power( (peak_area_phe-7*spe_mu)/(TMath::Sqrt(7)*spe_sig), 2.)*0.5 ) +
  0.2*0.2*0.2*0.2 / spe_sig/sqrt2pi/TMath::Sqrt(8.) *
     TMath::Exp( -TMath::Power( (peak_area_phe-8*spe_mu)/(TMath::Sqrt(8)*spe_sig), 2.)*0.5 );
  
  return y;
}   


/////////////////////////////////////////////////////////////////////////////////
//
//  QIPE_area_PDF (Rose)
//
/////////////////////////////////////////////////////////////////////////////////
double QIPE_area_PDF(double peak_area_phe, int channel) {

  double y = 
  0.8*0.8*0.8*0.8*0.8 / spe_sig/sqrt2pi/TMath::Sqrt(5.) *
     TMath::Exp( -TMath::Power( (peak_area_phe-5*spe_mu)/(TMath::Sqrt(5)*spe_sig), 2.)*0.5 ) + 
  5 * 0.8 * 0.8 * 0.8 * 0.8 * 0.2 /spe_sig/sqrt2pi/TMath::Sqrt(6.)*
     TMath::Exp( -TMath::Power( (peak_area_phe-6*spe_mu)/(TMath::Sqrt(6)*spe_sig), 2.)*0.5 ) + 
  10 * 0.8 * 0.8 * 0.8 * 0.2 * 0.2 /spe_sig/sqrt2pi/TMath::Sqrt(7.)*
     TMath::Exp( -TMath::Power( (peak_area_phe-7*spe_mu)/(TMath::Sqrt(7)*spe_sig), 2.)*0.5 ) + 
  10 * 0.8 * 0.8 * 0.2 * 0.2 * 0.2 / spe_sig/sqrt2pi/TMath::Sqrt(8.)*
     TMath::Exp( -TMath::Power( (peak_area_phe-8*spe_mu)/(TMath::Sqrt(8)*spe_sig), 2.)*0.5 ) + 
  5 * 0.8 * 0.2 * 0.2 * 0.2 * 0.2 / spe_sig/sqrt2pi/TMath::Sqrt(9.)*
     TMath::Exp( -TMath::Power( (peak_area_phe-9*spe_mu)/(TMath::Sqrt(9)*spe_sig), 2.)*0.5 ) +
  0.2*0.2*0.2*0.2*0.2 / spe_sig/sqrt2pi/TMath::Sqrt(10.)*
     TMath::Exp( -TMath::Power( (peak_area_phe-10*spe_mu)/(TMath::Sqrt(10)*spe_sig), 2.)*0.5 );

  return y;

}

