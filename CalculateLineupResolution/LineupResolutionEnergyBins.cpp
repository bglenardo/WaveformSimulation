#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"

double tritium[6][2] = {
          {20.2,3.93},
          {28.72,4.46},
          {37.17,5.28},
          {45.44,5.44},
          {54.32,6.74},
          {64.84,7.76}};

double dd[6][2] = {
          {26.79,3.88},
          {37.4,4.04},
          {49.13,4.51},
          {60.72,4.51},
          {71.97,4.99},
          {84.44,6.4}};



int main(int argc, char * argv[]) {

  TFile * file = new TFile("energy_bin_sigmas_offsets_it02_transit_time_mod_noise_18.root","RECREATE");
  TGraph waveform;

  TH1F * tt_h_aft_25[6];
  TH1F * tt_h_aft_05[6];
  TH1F * dd_h_aft_25[6];
  TH1F * dd_h_aft_05[6];
  char hist_name[100];

  for(int hh=0; hh<6; hh++){
     sprintf(hist_name,"tt_h_aft_25_%02d",hh);
     tt_h_aft_25[hh] = new TH1F(hist_name,hist_name,400,-10.,20.);
     sprintf(hist_name,"tt_h_aft_05_%02d",hh);
     tt_h_aft_05[hh] = new TH1F(hist_name,hist_name,400,-10.,20.); 
     sprintf(hist_name,"dd_h_aft_25_%02d",hh);
     dd_h_aft_25[hh] = new TH1F(hist_name,hist_name,400,-10.,20.);
     sprintf(hist_name,"dd_h_aft_05_%02d",hh);
     dd_h_aft_05[hh] = new TH1F(hist_name,hist_name,400,-10.,20.); 
  }

  TRandom3 rando;
  rando.SetSeed(0);
  WaveformGenerator w;
    w.SetSingletFraction(1. - 1./(1 + 1.*0.363/2.545));
    w.SetT1(0.363);
    w.SetT3(2.545);
    w.SetTa(1.118);
    w.SetTb(0.2358);
    w.SetA(1.070);
    w.SetSig(0.17);
  int n_ph;  

  for(int k=4; k<5; k++) { 
    for(int i=0; i<10000; i++) {
      if( i % 1000 == 0 )
         std::cout << "Running event " << i << std::endl;
      
      n_ph = (int) TMath::Floor( rando.Gaus(tritium[k][0],tritium[k][1])+0.5 );
      
      w.SetPhotonsInCh( n_ph );
      w.GenPhotonArrivalTimes();
      w.GenRandomPhdArea();
      w.GenerateWaveform();
   
      tt_h_aft_25[k]->Fill( w.GetAftT25Samples() );
      tt_h_aft_05[k]->Fill( w.GetAftT05Samples() );
  
    }
  }

    w.SetSingletFraction(1. - 1./(1 + 2.14*0.363/2.40));
    w.SetT1(0.363);
    w.SetT3(2.40);
    w.SetSig(0.345);
  


  for(int k=0; k<0; k++) { 
    for(int i=0; i<10000; i++) {
      if( i % 100 == 0 )
         std::cout << "Running event " << i << std::endl;
      
      n_ph = (int) TMath::Floor( rando.Gaus(dd[k][0],dd[k][1])+0.5 );
      
      w.SetPhotonsInCh( n_ph );
      w.GenPhotonArrivalTimes();
      w.GenRandomPhdArea();
      w.GenerateWaveform();
   
      dd_h_aft_25[k]->Fill( w.GetAftT25Samples() );
      dd_h_aft_05[k]->Fill( w.GetAftT05Samples() );
  
    }
  }

  file->Write();
  file->Close();

  return 0;

}
