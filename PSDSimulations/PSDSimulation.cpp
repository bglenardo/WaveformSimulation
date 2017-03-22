#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include <vector>

using namespace std;

int main(int argc, char * argv[]) {

  TFile * file = new TFile("dd_PSD_40-50.root","RECREATE");
  TGraph waveform;

  TH1F * h_PSD = new TH1F("h_PSD","h_PSD",50,0.,1.);
  TH1F * h_pulseShape = new TH1F("h_pulseShape","h_pulseShape",200,-20.,100.);
  TH1F * h_pulseShape25 = new TH1F("h_pulseShape25","h_pulseShape25",200,-20.,100.);

  TRandom3 rando;
  rando.SetSeed(0);
  WaveformGenerator w;
    w.SetSingletFraction(1. - 1./(1 + 0.347));
    w.SetT1(0.362);
    w.SetT3(2.40);
    w.SetSig(0.47);
    w.SetA(1.07);
    w.SetTa(1.12);
    w.SetTb(0.235);
  double x;

  // Simulated variables I'll need to reconstruct PSD parameter
  vector<double> sorted_times;
  vector<double> areas;
  double T05, T25;


  int n_ph, index;  
  double prompt_sum = 0.;
  double total_sum = 0.;

  for(int i=0; i<100000; i++) {

    prompt_sum = 0.;
    total_sum = 0.;

    if( i % 10000 == 0 )
       std::cout << "Running event " << i << std::endl;
    
    n_ph = (int) TMath::Floor( rando.Uniform() * 10. ) + 40;
   
    w.SetPhotonsInCh( n_ph );
    w.GenPhotonArrivalTimes();
    w.GenRandomPhdArea();
    w.GenerateWaveform();

    sorted_times = w.GetSortedTimes();
    areas = w.GetAreas();    
    T05 = w.GetAftT05Samples() + 0.8;
    T25 = w.GetAftT25Samples() + 0.8;

    for( int j=0; j<n_ph; j++){
      if( (sorted_times[j] >= T05-1.4) && (sorted_times[j] <= T05+13.) )
        total_sum += areas[j]; 
      if( (sorted_times[j] >= T05-1.) && (sorted_times[j] <= T05+3.6) )
        prompt_sum += areas[j]; 

      h_pulseShape->Fill( (sorted_times[j] - T05)*10. );
      h_pulseShape25->Fill( (sorted_times[j] - T25)*10. );
    }
    h_PSD->Fill(prompt_sum/total_sum);

  }

  for(int i=0; i<h_PSD->GetNbinsX(); i++) {
    printf("%f\t%f\n",h_PSD->GetBinLowEdge(i),h_PSD->GetBinContent(i));
  }

  file->Write();
  file->Close();

  return 0;

}
