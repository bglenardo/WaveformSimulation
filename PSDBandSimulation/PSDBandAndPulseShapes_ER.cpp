#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include <vector>
#include "nr_band_sim_5-9-2017.h"

using namespace std;

int main(int argc, char * argv[]) {

  double offset;// = atof(argv[1]);

  bool tritium = true;
  bool DD = false;  // True if you want to simulate DD
  int num_events_per_bin = 100000; // Number of events to simulate

  // Create a rootfile to store the histograms
  char fname[40];
  sprintf(fname,"er_band_dd_5-9-2017.root");
  TFile * file = new TFile( fname, "RECREATE");
  TGraph waveform;

  TH1F * h_PSD_40 = new TH1F("h_PSD_40","h_PSD_40",100,0.,1.);
  TH1F * h_pulseShape_40 = new TH1F("h_pulseShape_40","h_pulseShape_40",300,-10.,50.);

  // Create the histograms

  vector<double> psd_params;

  double psd_mean[100],
         psd_lower[100],
         psd_upper[100];

  TGraph * g_lower = new TGraph();
  TGraph * g_upper = new TGraph();
  TGraph * g_median = new TGraph();
  TGraph * g_leakage = new TGraph();
  double above_nr_median, total;


  TRandom3 rando;
  rando.SetSeed(0);

  // Create the waveform generator object with the right
  // parameters for either DD or tritium
  WaveformGenerator w;
  //double offset;
  if( tritium ) {
    w.SetSingletFraction(1. - 1./(1 + 0.042));
    w.SetT1(3.27441e-01);
    w.SetT3(2.58851e+00);
    w.SetSig(0.359220);
//    w.SetA(1.07);
//    w.SetTa(1.12);
//    w.SetTb(0.236);
    w.SetA(1.062);
    w.SetB(0.0574);
    w.SetTa(1.12);
    w.SetTb(0.27);
    offset = 0.10;
  }
  
  if( DD ) {
    w.SetSingletFraction(1. - 1./(1 + 0.269));
    w.SetT1(3.27441e-01);
    w.SetT3(2.58851e+00);
    w.SetSig(3.84173e-01);
    w.SetA(1.062);
    w.SetB(0.0574);
    w.SetTa(1.12);
    w.SetTb(0.27);
    offset = 0.12;
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
  for(int bin=1; bin<41; bin++) {
    above_nr_median = 0.; total = 0.;
    printf("Running S1 = %d\n",bin*5-2);
    psd_params.clear();

    // Event loop
    for(int i=0; i<num_events_per_bin; i++) {
      prompt_sum = 0.;
      total_sum = 0.;
      n_ph = bin*5-2;
      if( i % 10000 == 0 )
         std::cout << "\tRunning event " << i << " of " << num_events_per_bin <<  std::endl;
      
      w.SetPhotonsInCh( n_ph );
      w.GenPhotonArrivalTimes();
      w.GenRandomPhdArea();
      w.GenerateWaveform();
  
      sorted_times = w.GetSortedTimes();
      areas = w.GetAreas();    
      T05 = w.GetAftT05Samples();
  
      for( int j=0; j<n_ph; j++){
        double t_ph = sorted_times[j]-offset-T05;
        if( (t_ph >= -1.4) && (t_ph <= 13.) )
          total_sum += areas[j]; 
        if( (t_ph >= -1.) && (t_ph <= 3.4) )
          prompt_sum += areas[j];
        if( n_ph > 40. && n_ph < 50. )
          h_pulseShape_40->Fill( t_ph ); 
      }
      psd_params.push_back( prompt_sum/total_sum );
      if( n_ph > 40. && n_ph < 50. )
        h_PSD_40->Fill( prompt_sum/total_sum );
      if( prompt_sum/total_sum > nr_band_sim_array[bin-1][2] ) 
        above_nr_median = above_nr_median + 1.;
      total = total + 1.;
  
    }
    std::sort( psd_params.begin(), psd_params.end() );
    double low_ind = TMath::Floor( num_events_per_bin * 0.16 );
    double med_ind = TMath::Floor( num_events_per_bin * 0.5 );
    double hi_ind = TMath::Floor( num_events_per_bin * 0.84 );
    g_lower->SetPoint(bin-1,(double) bin*5.-2., psd_params[low_ind]);  
    g_median->SetPoint(bin-1,(double) bin*5.-2., psd_params[med_ind]);  
    g_upper->SetPoint(bin-1,(double) bin*5.-2., psd_params[hi_ind]);
    g_leakage->SetPoint(bin-1,(double) bin*5.-2., above_nr_median/total);

  }




  printf("\nPSD distribution:\n\n");
  printf("PF\tCounts\n");
  for(int i=1; i<41; i++) {
    double x,ylo,yme,yhi;
    g_lower->GetPoint(i,x,ylo);
    g_median->GetPoint(i,x,yme);
    g_upper->GetPoint(i,x,yhi);
    printf("%f\t%f\t%f\t%f\n",i*5.-2.,ylo,yme,yhi);
  }
  
  g_lower->SetName("g_lower");
  g_median->SetName("g_median");
  g_upper->SetName("g_upper");
  g_leakage->SetName("g_leakage");
  g_lower->Write();
  g_median->Write();
  g_upper->Write();
  g_leakage->Write();
  file->Write();
  file->Close();

  return 0;

}
