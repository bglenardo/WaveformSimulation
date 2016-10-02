////////////////////////////////////////////////////////////////////////////////////////////
//
//  SPEFunctions.cpp 
//
//     This file contains the parameterized funciton for SPE response in LUX, along with
//     the sum of doubles and triples and quadruples.
//  
//   2015-05-28 - Created. (Brian L.)
//
////////////////////////////////////////////////////////////////////////////////////////////

#include "TMath.h"
#include "SPEFunctions.hh"

//------------------------------------------------------------------------------------------
//    singleSPE - fit function for single photons
//------------------------------------------------------------------------------------------

double singleSPE(double * x, double * par) {

  // par[0] is normalization
  // par[1] is t0

  double y;
  
  if ( x[0] < (par[1] - 0.7) )
     y = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-6.644/10) )/(7.876/10.), 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-6.36/10) )/(12.57/10.), 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y = par[0]*(411./523.)*TMath::Exp( - 0.090*10.*(x[0] - par[1])) +
         par[0]*(80.418/523.)*TMath::Exp( - 0.02847*10.*(x[0] - par[1]));
  }

 
  if (par[0] < 0.)
    y = -1.e6;
  
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
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-6.644/10) )/(7.876/10.), 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-6.36/10) )/(12.57/10.), 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y1 = par[0]*(411./523.)*TMath::Exp( - 0.090*10.*(x[0] - par[1])) +
         par[0]*(80.418/523.)*TMath::Exp( - 0.02847*10.*(x[0] - par[1]));
  }
 
  if (par[0] < 0.)
    y1 = -1.e6;

  
  // Second photon
  if ( x[0] < (par[3] - 0.7) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-6.644/10) )/(7.876/10.), 2.)/2.);
  else if ( x[0] < (par[3] + 1.) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-6.36/10) )/(12.57/10.), 2.)/2.);
  else if ( x[0] >= (par[3] + 1.) ) {
     y2 = par[2]*(411./523.)*TMath::Exp( - 0.090*10.*(x[0] - par[3])) +
         par[2]*(80.418/523.)*TMath::Exp( - 0.02847*10.*(x[0] - par[3]));
  }

  if (par[2] < 0.)
    y2 = -1.e6;


  return (y1 + y2)/3.0479881347475857;

}


//------------------------------------------------------------------------------------------
//    tripleSPE - fit function for double SPE
//------------------------------------------------------------------------------------------


double tripleSPE(double * x, double * par) {

  
 double y1;
 double y2;
 double y3;


  // First photon  
  if ( x[0] < (par[1] - 0.7) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-6.644/10) )/(7.876/10.), 2.)/2.);
  else if ( x[0] < (par[1] + 1.) )
     y1 = par[0] * TMath::Exp( - TMath::Power( ( x[0]-(par[1]-6.36/10) )/(12.57/10.), 2.)/2.);
  else if ( x[0] >= (par[1] + 1.) ) {
     y1 = par[0]*(411./523.)*TMath::Exp( - 0.090*10.*(x[0] - par[1])) +
         par[0]*(80.418/523.)*TMath::Exp( - 0.02847*10.*(x[0] - par[1]));
  }
 
  if (par[0] < 0.)
    y1 = -1.e6;
  

  // Second Photon
  if ( x[0] < (par[3] - 0.7) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-6.644/10) )/(7.876/10.), 2.)/2.);
  else if ( x[0] < (par[3] + 1.) )
     y2 = par[2] * TMath::Exp( - TMath::Power( ( x[0]-(par[3]-6.36/10) )/(12.57/10.), 2.)/2.);
  else if ( x[0] >= (par[3] + 1.) ) {
     y2 = par[2]*(411./523.)*TMath::Exp( - 0.090*10.*(x[0] - par[3])) +
         par[2]*(80.418/523.)*TMath::Exp( - 0.02847*10.*(x[0] - par[3]));
  }
 
  if (par[2] < 0.)
    y2 = -1.e6;


  // Third Photon
  if ( x[0] < (par[5] - 0.7) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-6.644/10) )/(7.876/10.), 2.)/2.);
  else if ( x[0] < (par[5] + 1.) )
     y3 = par[4] * TMath::Exp( - TMath::Power( ( x[0]-(par[5]-6.36/10) )/(12.57/10.), 2.)/2.);
  else if ( x[0] >= (par[5] + 1.) ) {
     y3 = par[4]*(411./523.)*TMath::Exp( - 0.090*10.*(x[0] - par[5])) +
         par[4]*(80.418/523.)*TMath::Exp( - 0.02847*10.*(x[0] - par[5]));
  }
 
  if (par[4] < 0.)
    y3 = -1.e6;


  return (y1 + y2 + y3)/3.0479881347475857;

}


