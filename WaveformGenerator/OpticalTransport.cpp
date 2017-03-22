#include "TMath.h"


double OpticalTransport(double *x, double*p) {

   double y;

   if( x[0] >= 0.) 
       y = p[0]/p[1] * TMath::Exp(-x[0]/p[1]) + (1-p[0])/p[2] * TMath::Exp(-x[0]/p[2]);
   else
       y = 0.;

   return y;
  
}
