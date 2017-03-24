#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include <vector>

using namespace std;

int main(int argc, char * argv[]) {

  bool tritium = false;
  bool DD = true;  // True if you want to simulate DD
  int num_events = 50000; // Number of events to simulate
  int range_min = 40; // Min of S1 size range in phd
  int range_max = 50; // Max of S1 size range in phd


  // Create a rootfile to store the histograms
  TFile * file = new TFile("dd_PSD_40-50_with_uncorrection_1_5ns.root","RECREATE");
  TGraph waveform;


  // Create the histograms
  TH1F * h_PSD = new TH1F("h_PSD","h_PSD",50,0.,1.);
  TH1F * h_pulseShape = new TH1F("h_pulseShape","h_pulseShape",100,-20.,100.);
  TH1F * h_pulseShape25 = new TH1F("h_pulseShape25","h_pulseShape25",200,-20.,100.);




  TRandom3 rando;
  rando.SetSeed(0);

  // Create the waveform generator object with the right
  // parameters for either DD or tritium
  WaveformGenerator w;
  double offset;
  if( tritium ) {
    w.SetSingletFraction(1. - 1./(1 + 0.0908));
    w.SetT1(0.362);
    w.SetT3(2.55);
    w.SetSig(0.47);
    w.SetA(1.07);
    w.SetTa(1.12);
    w.SetTb(0.235);
    offset = 0.1;
  }
  
  if( DD ) {
    w.SetSingletFraction(1. - 1./(1 + 0.352));
    w.SetT1(0.358);
    w.SetT3(2.40);
    w.SetSig(0.47);
    w.SetA(1.07);
    w.SetTa(1.12);
    w.SetTb(0.235);
    offset = 0.15;
  }

  double x;

  // Simulated variables I'll need to reconstruct PSD parameter
  vector<double> sorted_times;
  vector<double> areas;
  double T05, T25;


  int n_ph, index;  
  double prompt_sum = 0.;
  double total_sum = 0.;


  // Event loop
  for(int i=0; i<num_events; i++) {

    prompt_sum = 0.;
    total_sum = 0.;

    if( i % 5000 == 0 )
       std::cout << "Running event " << i << " of " << num_events <<  std::endl;
    
    n_ph = (int) TMath::Floor( rando.Uniform() * (range_max-range_min) ) + range_min;
   
    w.SetPhotonsInCh( n_ph );
    w.GenPhotonArrivalTimes();
    w.GenRandomPhdArea();
    w.GenerateWaveform();

    sorted_times = w.GetSortedTimes();
    areas = w.GetAreas();    
    T05 = w.GetAftT05Samples();
    T25 = w.GetAftT25Samples();

    for( int j=0; j<n_ph; j++){
      double t_ph = sorted_times[j]-offset;
      if( (t_ph >= T05-1.4) && (t_ph <= T05+13.) )
        total_sum += areas[j]; 
      if( (t_ph >= T05-1.) && (t_ph <= T05+3.6) )
        prompt_sum += areas[j]; 

      h_pulseShape->Fill( (t_ph - T05)*10. );
      h_pulseShape25->Fill( (t_ph - T25)*10. );
    }
    h_PSD->Fill(prompt_sum/total_sum);

  }

  printf("\nPSD distribution:\n\n");
  printf("PF\tCounts\n");
  for(int i=0; i<h_PSD->GetNbinsX(); i++) {
    printf("%f\t%f\n",h_PSD->GetBinLowEdge(i),h_PSD->GetBinContent(i));
  }
  
  printf("\nPulse shape distribution:\n\n");
  printf("time(ns)\tcounts\n");
  for(int i=0; i<h_PSD->GetNbinsX(); i++) {
    printf("%f\t%f\n",h_PSD->GetBinLowEdge(i),h_PSD->GetBinContent(i));
  }


  file->Write();
  file->Close();

  return 0;

}
