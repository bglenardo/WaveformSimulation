////////////////////////////////////////////////////////////////////////////////////////////
//
//  SPEFunctions_Rose_test.cpp 
//
//     This file contains the parameterized function for SPE response in LUX, along with
//     the sum of doubles and triples and quadruples and quintuples.
//  
//     2015-05-28 - Created. (Brian L.)
//     2015-07-22 - Modified to include functions of up to 5 photons. (Rose)
//     2015-07-30 - The height threshold doesn't seem to be doing much. I've raised it to
//                  0.095 just to see if that helps at all. I also changed the
//                  large-negative-value to -1.e100. (Brian L.) 
//     2015-08-18 - Added the summed S1 function to allow fits to get start time.
//     2015-09-01 - Added the Baseline function to fit traces and allow for 0-ph cases.
//
////////////////////////////////////////////////////////////////////////////////////////////

#include "TMath.h"

// Global shape parameters;
double first_width = 7.876/10.;
double second_width = 12.57/10.;
double first_offset = 6.644/10.;
double second_offset = 6.36/10.;
double first_exp_norm = 411./523.;
double second_exp_norm = (80.418/523.);
double first_exp_dec = 0.090*10.;
double second_exp_dec = 0.02847*10.;
double SPEfitNorm = 3.04;


//double first_width = 7.176/10;
//double second_width = 1.08;
//double first_offset = 6.644/10.;
//double second_offset = 6.36/10.;
//double first_exp_norm = 340./523.;
//double second_exp_norm = (85.418/523.);
//double first_exp_dec = 1.19;
//double second_exp_dec = 0.19;
//double SPEfitNorm = 2.982674469;
//------------------------------------------------------------------------------------------
//   S1Template  
//------------------------------------------------------------------------------------------

double S1Template(double * x, double * par) {



  double l = 0.28;
  double s = 1.3;

  double y = par[0] * l/(2.*0.1544) * TMath::Exp(l/2*(2*par[1] + l*s*s - 2*x[0])) * 
                      TMath::Erfc( (par[1] + l*s*s - x[0]) / (TMath::Sqrt(2)*s) );

  return y;

}


//------------------------------------------------------------------------------------------
//   Baseline  
//------------------------------------------------------------------------------------------

double Baseline(double * x, double * par) {

  double y = par[0]; 
  return y;

}



//------------------------------------------------------------------------------------------
//    singleSPE - fit function for single photons
//------------------------------------------------------------------------------------------

double singleSPE(double * x, double * par) {

  // par[0] is normalization
  // par[1] is t0

  double y;
  
  if ( x[0] < (par[1] - 0.7) )
     y = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y = par[0]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[1])) +
         par[0]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[1]));
  }

 
//  if (par[0] < 0.085)
//   y = -1.e100;
  
  return y/3.0479881347475857;

}


//------------------------------------------------------------------------------------------
//    doubleSPE - fit function for double SPE
//------------------------------------------------------------------------------------------

double doubleSPE(double * x, double * par) {

  
 double y1;
 double y2;
 
  // First photon 
  if ( x[0] < (par[1] - 0.7) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y1 = par[0]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[1])) +
         par[0]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[1]));
  }
 
//  if (par[0] < 0.085)
//    y1 = -1.e100;

  
  // Second photon
  if ( x[0] < (par[3] - 0.7) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[3] + 1.) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[3] + 1.) ) {
     y2 = par[2]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[3])) +
         par[2]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[3]));
  }

//  if (par[2] < 0.085)
//    y2 = -1.e100;

  // If photons are closer than a single sample,
  // return garbage. 
//  if ( TMath::Abs(par[3] - par[1]) < 1. )
//    y2 = -1.e100;

  return (y1 + y2)/3.0479881347475857;

}


//------------------------------------------------------------------------------------------
//    tripleSPE - fit function for triple SPE
//------------------------------------------------------------------------------------------


double tripleSPE(double * x, double * par) {

  
 double y1;
 double y2;
 double y3;


  // First photon  
  if ( x[0] < (par[1] - 0.7) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y1 = par[0]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[1])) +
         par[0]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[1]));
  }
 
  //if (par[0] < 0.085)
  //  y1 = -1.e100;
  

  // Second Photon
  if ( x[0] < (par[3] - 0.7) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[3] + 1.) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[3] + 1.) ) {
     y2 = par[2]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[3])) +
         par[2]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[3]));
  }
 
 // if (par[2] < 0.085)
 //   y2 = -1.e100;


  // Third Photon
  if ( x[0] < (par[5] - 0.7) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[5] + 1.) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[5] + 1.) ) {
     y3 = par[4]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[5])) +
         par[4]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[5]));
  }
 
  // if (par[4] < 0.085)
   // y3 = -1.e100;

//  if ( TMath::Abs(par[3]-par[1]) < 1. ||
//       TMath::Abs(par[5]-par[1]) < 1. ||
//       TMath::Abs(par[5]-par[3]) < 1. )
//    y3 = -1.e100;

  return (y1 + y2 + y3)/3.0479881347475857;

}

//----------------------------------------------------------------
//    quadSPE - fit function for four PHE
//----------------------------------------------------------------

double quadSPE(double * x, double * par) {

double y1;
double y2;
double y3;
double y4;

//First Photon
if ( x[0] < (par[1] - 0.7) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y1 = par[0]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[1])) +
         par[0]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[1]));
  }

//  if (par[0] < 0.085)
//    y1 = -1.e100;

//Second Photon
if ( x[0] < (par[3] - 0.7) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[3] + 1.) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[3] + 1.) ) {
     y2 = par[2]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[3])) +
         par[2]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[3]));
  }

  //if (par[2] < 0.085)
 //   y2 = -1.e100;

//Third Photon
if ( x[0] < (par[5] - 0.7) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[5] + 1.) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[5] + 1.) ) {
     y3 = par[4]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[5])) +
         par[4]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[5]));
  }

//  if (par[4] < 0.085)
//    y3 = -1.e100;

//Fourth Photon
if ( x[0] < (par[7] - 0.7) )
     y4 = par[6] * TMath::Exp( - TMath::Power( ( x[0]-(par[7]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[7] + 1.) )
     y4 = par[6] * TMath::Exp( - TMath::Power( ( x[0]-(par[7]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[7] + 1.) ) {
     y4 = par[6]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[7])) +
         par[6]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[7]));
  }

// if (par[6] < 0.085)
//    y4 = -1.e100;

//  if (TMath::Abs(par[3]-par[1]) < 1. ||
//      TMath::Abs(par[5]-par[1]) < 1. ||
//      TMath::Abs(par[7]-par[1]) < 1. ||
//      TMath::Abs(par[7]-par[5]) < 1. ||
//      TMath::Abs(par[7]-par[3]) < 1. ||
//      TMath::Abs(par[5]-par[3]) <1.)
//    y4 = -1.e100;

  return (y1 + y2 + y3 + y4)/3.0479881347475857;
}

//-------------------------------------------------------------------
//       quintSPE - fit function for 5 PHE
//-------------------------------------------------------------------


double quintSPE(double * x, double * par) {

double y1;
double y2;
double y3;
double y4;
double y5;

//First photon
if ( x[0] < (par[1] - 0.7) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y1 = par[0]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[1])) +
         par[0]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[1]));
  }

//  if (par[0] < 0.085)
//    y1 = -1.e100;

//Second photon
 if ( x[0] < (par[3] - 0.7) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[3] + 1.) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[3] + 1.) ) {
     y2 = par[2]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[3])) +
         par[2]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[3]));
  }

//  if (par[2] < 0.085)
//    y2 = -1.e100;

//Third photon
if ( x[0] < (par[5] - 0.7) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[5] + 1.) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[5] + 1.) ) {
     y3 = par[4]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[5])) +
         par[4]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[5]));
  }

//  if (par[4] < 0.085)
//    y3 = -1.e100;

//Fourth photon
if ( x[0] < (par[7] - 0.7) )
     y4 = par[6] * TMath::Exp( - TMath::Power( ( x[0]-(par[7]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[7] + 1.) )
     y4 = par[6] * TMath::Exp( - TMath::Power( ( x[0]-(par[7]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[7] + 1.) ) {
     y4 = par[6]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[7])) +
         par[6]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[7]));
  }

//  if (par[6] < 0.085)
//    y4 = -1.e100;
    
//Fifth photon
if ( x[0] < (par[9] - 0.7) )
     y5 = par[8] * TMath::Exp( - TMath::Power( ( x[0]-(par[9]-first_offset) )/first_width, 2.)/2.);
  else if ( x[0] < (par[9] + 1.) )
     y5 = par[8] * TMath::Exp( - TMath::Power( ( x[0]-(par[9]-second_offset) )/second_width, 2.)/2.);
  else if ( x[0] >= (par[9] + 1.) ) {
     y5 = par[8]*first_exp_norm*TMath::Exp( - first_exp_dec*(x[0] - par[9])) +
         par[8]*second_exp_norm*TMath::Exp( - second_exp_dec*(x[0] - par[9]));
  }

//  if (par[8] < 0.085)
//    y5 = -1.e100;

//  if (TMath::Abs(par[3]-par[1]) < 1. ||
//      TMath::Abs(par[5]-par[1]) < 1. ||
//      TMath::Abs(par[7]-par[1]) < 1. ||
//      TMath::Abs(par[9]-par[1]) < 1. ||
//      TMath::Abs(par[9]-par[7]) < 1. ||
//      TMath::Abs(par[9]-par[5]) < 1. ||
//      TMath::Abs(par[9]-par[3]) < 1. ||
//      TMath::Abs(par[7]-par[5]) < 1. ||
//      TMath::Abs(par[7]-par[3]) < 1. ||
//      TMath::Abs(par[5]-par[3]) < 1. )
//   y5 = -1.e100;

  return (y1 + y2 + y3 + y4 + y5)/3.0479881347475857;

}

