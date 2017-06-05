#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include <vector>

using namespace std;

int main(int argc, char * argv[]) {

  double offset;// = atof(argv[1]);

  bool tritium = false;
  bool DD = true;  // True if you want to simulate DD
  int num_events_per_bin = 10000; // Number of events to simulate

  // Create a rootfile to store the histograms
  char fname[40];
  sprintf(fname,"nr_band_4-12-2017.root");
  TFile * file = new TFile( fname, "RECREATE");
  TGraph waveform;


  // Create the histograms

  vector<double> psd_params;

  double psd_mean[100],
         psd_lower[100],
         psd_upper[100];

  TGraph * g_lower = new TGraph();
  TGraph * g_upper = new TGraph();
  TGraph * g_median = new TGraph();



  TRandom3 rando;
  rando.SetSeed(0);

  // Create the waveform generator object with the right
  // parameters for either DD or tritium
  WaveformGenerator w;
  //double offset;
  if( tritium ) {
    w.SetSingletFraction(1. - 1./(1 + 0.044));
    w.SetT1(2.68858e-01);
    w.SetT3(2.58778e+00);
    w.SetSig(3.82367e-01);
    w.SetA(1.07);
    w.SetTa(1.12);
    w.SetTb(0.235);
    offset = 0.05;
  }
  
  if( DD ) {
    w.SetSingletFraction(1. - 1./(1 + 0.265));
    w.SetT1(2.68858e-01);
    w.SetT3(2.38718e+00);
    w.SetSig(3.82367e-01);
    w.SetA(1.07);
    w.SetTa(1.12);
    w.SetTb(0.235);
    offset = 0.15;
  }

  double x;

  // Simulated variables I'll need to reconstruct PSD parameter
  vector<double> sorted_times;
  vector<double> areas;
  double T05;


  int n_ph, index;  
  double prompt_sum = 0.;
  double total_sum = 0.;

  // S1 bin loop
  for(int bin=1; bin<69; bin++) {
    printf("Running S1 = %d\n",bin*3);
    psd_params.clear();

    // Event loop
    for(int i=0; i<num_events_per_bin; i++) {
      prompt_sum = 0.;
      total_sum = 0.;
      n_ph = bin*3;
      if( i % 1000 == 0 )
         std::cout << "\tRunning event " << i << " of " << num_events_per_bin <<  std::endl;
      
      w.SetPhotonsInCh( n_ph );
      w.GenPhotonArrivalTimes();
      w.GenRandomPhdArea();
      w.GenerateWaveform();
  
      sorted_times = w.GetSortedTimes();
      areas = w.GetAreas();    
      T05 = w.GetAftT05Samples();
  
      for( int j=0; j<n_ph; j++){
        double t_ph = sorted_times[j]-offset;
        if( (t_ph >= T05-1.4) && (t_ph <= T05+13.) )
          total_sum += areas[j]; 
        if( (t_ph >= T05-1.) && (t_ph <= T05+3.6) )
          prompt_sum += areas[j]; 
      }
      psd_params.push_back( prompt_sum/total_sum );
  
    }
    std::sort( psd_params.begin(), psd_params.end() );
    double low_ind = TMath::Floor( num_events_per_bin * 0.16 );
    double med_ind = TMath::Floor( num_events_per_bin * 0.5 );
    double hi_ind = TMath::Floor( num_events_per_bin * 0.84 );
    g_lower->SetPoint(bin,(double) bin*3., psd_params[low_ind]);  
    g_median->SetPoint(bin,(double) bin*3., psd_params[med_ind]);  
    g_upper->SetPoint(bin,(double) bin*3., psd_params[hi_ind]);

  }




  printf("\nPSD distribution:\n\n");
  printf("PF\tCounts\n");
  for(int i=0; i<100; i++) {
    double x,ylo,yme,yhi;
    g_lower->GetPoint(i,x,ylo);
    g_median->GetPoint(i,x,yme);
    g_upper->GetPoint(i,x,yhi);
    printf("%f\t%f\t%f\t%f\n",i*3.,ylo,yme,yhi);
  }
  
  g_lower->SetName("g_lower");
  g_median->SetName("g_median");
  g_upper->SetName("g_upper");
  g_lower->Write();
  g_median->Write();
  g_upper->Write();
  file->Close();

  return 0;

}
