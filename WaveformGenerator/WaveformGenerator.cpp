#include "TMath.h"
#include "TRandom3.h"
#include <algorithm>
#include "SPEFunctions.hh"
#include "WaveformGenerator.hh" 
#include "NoiseSpectrum.hh"
#include <iostream>
#include "OpticalTransport.hh"
#include "TF1.h"
#include "corrections/PMTInfo.h"
#include "corrections/DAQ_Corrections.h"
#include "corrections/DAQ_Errors_04.h"

/************************************************************************
     WaveformGenerator (constructor)
************************************************************************/

WaveformGenerator::WaveformGenerator() {

   t1 = 0.31;
   t3 = 2.7;
   ta = 1.118;
   tb = 2.0;
   A = 0.088;
   B = 1.;
   optical = TF1("optical",OpticalTransport,0.,200.,3);
   optical.SetParameters(B,ta,tb);
   optical.SetNpx(400);
   sig = 0.175;
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
}

/************************************************************************
     SetB
************************************************************************/
void WaveformGenerator::SetB( double B_ ) {
    B = B_;
    optical.SetParameter(0,B_);
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
   uncorrected_sorted_times.clear();

   double temp_time = 0.;

   for(int i=0; i<photons_in_ch; i++){
     temp_time = 0.;     
     temp_time += r.Gaus( 0., sig );

     if( r.Uniform() > A ) {
       if( r.Uniform() < B && B <= 1.)
          temp_time += r.Exp( ta );
       else if ( r.Uniform() >= B && B <= 1. )
          temp_time += r.Exp( tb );
       else if ( B > 1. ) 
          temp_time += optical.GetRandom();  
     }

     if( r.Uniform() < (1. - singlet_fraction) )
        temp_time += r.Exp( t3 );
     else
        temp_time += r.Exp( t1 );
     temp_ch = GenChannelNum();

     sorted_times.push_back(temp_time + daq_errors[temp_ch]/10.);

     uncorrected_sorted_times.push_back(temp_time + GenUncorrectionTime( temp_ch ) );

   }
   
   std::sort( sorted_times.begin(), sorted_times.end() );
   std::sort( uncorrected_sorted_times.begin(), uncorrected_sorted_times.end() );

   return;
}


/************************************************************************
     GenUncorrectionTime
************************************************************************/
double WaveformGenerator::GenUncorrectionTime( int ch ) {

     double z_pmt;
     if( ch < 60 || ch == 120 )
         z_pmt = 59.;
     else
         z_pmt = -1.; 
     double z = drift_time * 1.5 / 10.; // converts us to cm
     // The following block uniformly distributes events in x/y
     double t, u, ra;
     t = 2 * TMath::Pi() * r.Uniform();
     u = (r.Uniform() + r.Uniform())*20.;
     if( u > 20. ) ra = 40. - u; 
     else  ra = u;
     double x = ra * TMath::Cos(t);
     double y = ra * TMath::Sin(t);

     // Compute direct-path time (which is subtracted off in real data)
     double path_time = TMath::Sqrt( (0. - PMT[ch][0])*(0. - PMT[ch][0]) +
                         (0. - PMT[ch][1])*(0. - PMT[ch][1]) +
                         (46. - z_pmt)*( 46. -z_pmt) ) / 19.2 / 10.;
 
//    printf("Correction: %f\t Ch: %d\n",path_time + daq_correction[ch],ch);
    return path_time + daq_correction[ch];

}


/************************************************************************
     GenChannelNum
************************************************************************/
int WaveformGenerator::GenChannelNum() {

     // Coefficients for probability for photons to arrive in the top
     // PMT array. Fit by Evan Pease to tritium data.
     topArrayProbabilityCoeffs[0] = 1.10677683e-07;
     topArrayProbabilityCoeffs[1] = -7.23082313e-04; 
     topArrayProbabilityCoeffs[2] = 3.36092810e-01; 

     topProb = topArrayProbabilityCoeffs[0]*drift_time*drift_time +
               topArrayProbabilityCoeffs[1]*drift_time +
               topArrayProbabilityCoeffs[0];

     int ch = 106;     
     while( daq_correction[ch] < -1000 ) {
       if( r.Uniform() < 1. - topProb ) {
         ch = TMath::Floor( r.Uniform() * 61 ) + 60;
         if( ch == 120 ) ch = 121;
       } else {
         ch = TMath::Floor( r.Uniform() * 61 );
         if( ch == 60 ) ch = 120;
       }
     }
     return ch;
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
   uncorrected_wave_vec.clear();
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
   double y_un = 0;
   double pt_idx = 0;


   for(double xval=start; xval<end; xval += 1.){
      y = 0.;
      y_un = 0.;
      for(int ii=0; ii<photons_in_ch; ii++) {
         par[0] = areas[ii]; par[1] = sorted_times[ii];
         x[0] = xval;
         y += singleSPE( x, par );
         par[0] = areas[ii]; par[1] = uncorrected_sorted_times[ii];
         y_un += singleSPE( x, par );
      }
     wave_vec.push_back(y);
     uncorrected_wave_vec.push_back(y_un);
     pt_idx++;

  }
  GenerateBaseline();
  for(int i=0; i < (int) wave_vec.size(); i++) {
//     wave_vec[i] += r.Gaus(0.,0.01); // baseline_vec[i];
     wave_vec[i] += baseline_vec[i];
     uncorrected_wave_vec[i] += baseline_vec[i];
//     std::cout << baseline_vec[i] << std::endl;
     waveform.SetPoint(i, i+start, wave_vec[i]);
     uncorrected_waveform.SetPoint(i, i+start, uncorrected_wave_vec[i]);
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

   for( std::vector<double>::size_type i=0; i<uncorrected_wave_vec.size(); i++) {
     peak_area += uncorrected_wave_vec[i];
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

   for(int i=0; i< (int)uncorrected_wave_vec.size(); i++) {
      sum += uncorrected_wave_vec[i];
      if( sum < 0.05*peak_area /*&& !found_aft_t05*/ ) {
        aft_t05_samples = (double) i + trace_start;
        found_aft_t05 = true;
      }
      //printf("%f\t%f\n",uncorrected_wave_vec[i],aft_t05_samples);
      if( sum > 0.25*peak_area && !found_aft_t25) {
        aft_t25_samples = (double) i + trace_start;
        found_aft_t25 = true;
      }

   }
}


