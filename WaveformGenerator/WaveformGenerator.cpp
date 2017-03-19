#include "TMath.h"
#include "TRandom3.h"
#include <algorithm>
#include "SPEFunctions.hh"
#include "WaveformGenerator.hh" 
#include "NoiseSpectrum.hh"
#include <iostream>
#include "OpticalTransport.hh"
#include "TF1.h"

/************************************************************************
     WaveformGenerator (constructor)
************************************************************************/

WaveformGenerator::WaveformGenerator() {

   t1 = 0.31;
   t3 = 2.7;
   ta = 1.118;
   tb = 2.0;
   A = 1.;
   optical = TF1("optical",OpticalTransport,0.,200.,3);
   optical.SetParameters(A,ta,tb);
   optical.SetNpx(400);
   sig = 0.75;
   singlet_fraction = 0.4;
 
   r.SetSeed(0);

   photons_in_array = 0;
   photons_in_ch = 0;
   aft_t25_samples = -100.;
   aft_t05_samples = -100.;

      
}


/************************************************************************
     SetTa 
************************************************************************/
void WaveformGenerator::SetTa( double ta_ ) {
    ta = ta_;
    optical.SetParameter(1,ta_);
}

/************************************************************************
     SetTb 
************************************************************************/
void WaveformGenerator::SetTb( double tb_ ) {
    tb = tb_;
    optical.SetParameter(2,tb_);
}

/************************************************************************
     SetA
************************************************************************/
void WaveformGenerator::SetA( double A_ ) {
    A = A_;
    optical.SetParameter(0,A_);
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
     temp_time += r.Gaus( 0., sig );

     if( r.Uniform() < A && A <= 1.)
        temp_time += r.Exp( ta );
     else if ( r.Uniform() >= A && A <= 1. )
        temp_time += r.Exp( tb );
     else if ( A > 1. ) 
        temp_time += optical.GetRandom();  

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
   trace_start = -100.;
   aft_t25_samples = -100.;
   aft_t05_samples = -100.;

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
//     y += r.Gaus(0.,0.01); 
     wave_vec.push_back(y);
     pt_idx++;

  }
  GenerateBaseline();
  for(int i=0; i < (int) wave_vec.size(); i++) {
//     wave_vec[i] += r.Gaus(0.,0.01); // baseline_vec[i];
     wave_vec[i] += baseline_vec[i];
//     std::cout << baseline_vec[i] << std::endl;
     waveform.SetPoint(i, i+start, wave_vec[i]);
//     std::cout << baseline_vec[i] << std::endl;
  }
  GenPeakArea();
  GenAftTimes();

}

/************************************************************************
     GenerateBaseline
************************************************************************/
void WaveformGenerator::GenerateBaseline() {
   baseline_vec.clear();

   int N = (sizeof(NoiseBinCenter)/sizeof(*NoiseBinCenter)); // number of bins of log histogram. 
   int num_samples = (int) wave_vec.size();
   if (num_samples == 0) { 
         std::cerr << "Trying to generate baseline_vec before waveform! Abort!" << std::endl;
         return;
   }
   double randomPhase = 0.;// r.Uniform() * 2 * TMath::Pi();

   double phe_scale = 22.0/1000.; // (4.6 phe/mV) / (1000 mV/V)


   for(int i=0; i<N; i++) {

      randomPhase = r.Uniform() * 2 * TMath::Pi();
//      if( NoiseBinCenter[i] > 3.e7 ) break;

      for(int s=0; s<num_samples; s++) {
          if(i == 0) {
              baseline_vec.push_back( TMath::Sqrt( 2*NoiseBinWidth[i] ) * NoiseBinHeight[i] / phe_scale *
                                  TMath::Sin( 2 * TMath::Pi() * NoiseBinCenter[i] * s * 1e-8 + randomPhase));
          } else {
              baseline_vec[s] += TMath::Sqrt( 2*NoiseBinWidth[i] ) * NoiseBinHeight[i] / phe_scale *
                             TMath::Sin( 2 * TMath::Pi() * NoiseBinCenter[i] * s * 1e-8 + randomPhase);
          }
      }  

   }

   for(int i=0; i<num_samples; i++) {
      baseline.SetPoint(i,(double) i, baseline_vec[i]);
   }
   
}



/************************************************************************
     GenerateWaveform
************************************************************************/

void WaveformGenerator::GenPeakArea() {
   peak_area = 0.;
   if( wave_vec.size() == 0) {
      peak_area = 0.;
      return;
   }

   for( std::vector<double>::size_type i=0; i<wave_vec.size(); i++) {
     peak_area += wave_vec[i];
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

/************************************************************************
     GenAftTimes
************************************************************************/

void WaveformGenerator::GenAftTimes() {
   
   double sum = 0.;
   bool found_aft_t05 = false, 
        found_aft_t25 = false;

   for(int i=0; i< (int)wave_vec.size(); i++) {
      sum += wave_vec[i];
      if( sum > 0.05*peak_area && !found_aft_t05) {
        aft_t05_samples = (double) i + trace_start;
        found_aft_t05 = true;
      }
      if( sum > 0.25*peak_area && !found_aft_t25) {
        aft_t25_samples = (double) i + trace_start;
        found_aft_t25 = true;
      }

   }
}


