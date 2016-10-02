////////////////////////////////////////////////////////////////////////////////////////////
//
//  SPEFunctions_Rose_test.h
//
//     This file contains the parameterized funciton for SPE response in LUX, along with
//     the sum of doubles and triples and quadruples.
//  
//   2015-05-29 - Created. (Brian L.)
//   2015-07-28 - Modified (Rose)
//   2015-08-18 - Added S1Template function (Brian L.)
//   2015-09-01 - Added Baseline function (Brian L.)
//
////////////////////////////////////////////////////////////////////////////////////////////

double S1Template(double * x, double * par);
double Baseline(  double * x, double * par);
double singleSPE( double * x, double * par);
double doubleSPE( double * x, double * par);
double tripleSPE( double * x, double * par);
double quadSPE(   double * x, double * par);
double quintSPE(  double * x, double * par);
