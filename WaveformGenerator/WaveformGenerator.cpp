#include "TMath.h"
#include "TRandom3.h"
#include <algorithm>
#include "SPEFunctions.hh"
#include "WaveformGenerator.hh" 

/************************************************************************
     WaveformGenerator (constructor)
************************************************************************/

WaveformGenerator::WaveformGenerator() {

   t1 = 0.31;
   t3 = 2.4;
   tr = 1.8;
   sig = 0.55;
   singlet_fraction = 0.4;
 
   r.SetSeed(0);

   photons_in_array = 0;
   photons_in_ch = 0;
      
}



/************************************************************************
     GenPhotonsInCh
************************************************************************/

void WaveformGenerator::GenPhotonsInCh() {

  double p = 1./61.;
  photons_in_ch = r.Binomial( photons_in_array, p );

  return;
}


/************************************************************************
     GenPhotonArrivalTimes
************************************************************************/

void WaveformGenerator::GenPhotonArrivalTimes() {

   if( photons_in_ch == 0 ) return;   

   sorted_times.clear();

   double temp_time = 0.;

   for(int i=0; i<photons_in_ch; i++){
     temp_time = 0.;     
     temp_time += r.Gaus( 0., sig ) + r.Exp( tr );
     if( r.Uniform() < (1. - singlet_fraction) )
        temp_time += r.Exp( t3 );
     else
        temp_time += r.Exp( t1 );
     sorted_times.push_back(temp_time);
   }
   
   std::sort( sorted_times.begin(), sorted_times.end() );

   return;
}

/************************************************************************
     GenRandomPhdArea
************************************************************************/

void WaveformGenerator::GenRandomPhdArea() {

  if( photons_in_ch == 0 ) return;   

  areas.clear();

  double mu = 0.8327873877969729;
  double spe_sig = 0.35;

  double temp_area = 0.;

  for(int i=0; i<photons_in_ch; i++) {
    if(r.Uniform() < 0.2)
      temp_area = r.Gaus( 2 * mu, TMath::Sqrt(2) * spe_sig );
    else
      temp_area = r.Gaus( mu, spe_sig);

    areas.push_back(temp_area);
  }

}


/************************************************************************
     GenerateWaveform
************************************************************************/

void WaveformGenerator::GenerateWaveform() {

   double start = -20.;
   double end = 30.;
   wave_vec.clear();
   waveform.Set(0);

   if(photons_in_ch > 0) {
      start = sorted_times[0] - 20. + r.Uniform();
      end = sorted_times[photons_in_ch-1] + 30.;
   }

   trace_start = start;

   double x[1];
   double par[2];
   double y = 0;
   double pt_idx = 0;


   for(double xval=start; xval<end; xval += 1.){
      y = 0.;
      for(int ii=0; ii<photons_in_ch; ii++) {
         par[0] = areas[ii]; par[1] = sorted_times[ii];
         x[0] = xval;
         y += singleSPE( x, par );
      }
     y += r.Gaus(0.,0.01); 
     waveform.SetPoint(pt_idx, xval, y);
     wave_vec.push_back(y);
     pt_idx++;

  }

}


/************************************************************************
     GenerateNewWaveform
************************************************************************/

void WaveformGenerator::GenerateNewWaveform() {

   GenPhotonsInCh();
   GenPhotonArrivalTimes();
   GenRandomPhdArea();
   GenerateWaveform();

}
