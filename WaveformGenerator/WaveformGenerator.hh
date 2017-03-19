/***********************************************************************************

   WaveformGenerator.hh
    
     Class definition for WaveformGenerator

***********************************************************************************/

#ifndef _WaveformGenerator_
#define _WaveformGenerator_

#include "TMath.h"
#include "TRandom3.h"
#include <vector>
#include "TGraph.h"
#include "TF1.h"
#include "OpticalTransport.hh"

class WaveformGenerator {

   private:
      int photons_in_array;
      int photons_in_ch;
      std::vector<double> sorted_times;
      std::vector<double> areas;     
      double peak_area;
      double aft_t05_samples;
      double aft_t25_samples;
 
      double t1;  // Singlet time constant
      double t3;  // Triplet time constant
      double ta;  // Reflection time constant 1
      double tb;  // Reflection time constant 2
      double A;  // Reflection time weightin
      double sig;  // Normal smearing
      double singlet_fraction;

      TRandom3 r;

      TGraph waveform;
      TGraph baseline;
      std::vector<double> wave_vec;
      std::vector<double> baseline_vec;
      double trace_start;
      TF1 optical;
     
   public:
      WaveformGenerator();

      void GenPhotonsInCh();
      void GenRandomPhdArea();
      void GenPhotonArrivalTimes();
      void GenerateBaseline();
      void GenerateWaveform();  
      void GenerateNewWaveform();
      void GenPeakArea();
      void GenAftTimes();

      double GetT1() { return t1; }
      double GetT3() { return t3; }
      double GetTa() { return ta; }
      double GetTb() { return tb; }
      double GetA() { return A; }
      double GetSig() { return sig; }
      int GetPhotonsInArray() { return photons_in_array; }
      int GetPhotonsInCh()    { return photons_in_ch; }
      std::vector<double> GetSortedTimes() { return sorted_times; }
      std::vector<double> GetAreas() { return areas; }
      double GetPeakArea() { return peak_area; }
      TGraph GetWaveform() { return waveform; }
      TGraph GetBaseline() { return baseline; }
      std::vector<double> GetWaveVec() { return wave_vec; } 
      std::vector<double> GetBaselineVec() { return baseline_vec; } 
      double GetTraceStart() { return trace_start; }
      double GetAftT25Samples() { return aft_t25_samples; }
      double GetAftT05Samples() { return aft_t05_samples; }

      void SetT1( double t1_set ) { t1 = t1_set; }
      void SetT3( double t3_set ) { t3 = t3_set; }
      void SetTa( double ta_set );
      void SetTb( double tb_set );
      void SetA( double A_set );
      void SetSig( double sig_set ) { sig = sig_set; }
      void SetSingletFraction( double singlet_fraction_set ) { singlet_fraction = singlet_fraction_set; }

      void SetPhotonsInArray( int photons_in_array_set ) { photons_in_array = photons_in_array_set; }
      void SetPhotonsInCh( int photons_in_ch_set ) { photons_in_ch = photons_in_ch_set; }
      void SetSortedTimes( std::vector<double> sorted_times_set ) { sorted_times = sorted_times_set; }
      
};

#endif
