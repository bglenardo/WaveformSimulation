#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"

int main(int argc, char * argv[]) {

  TFile * file = new TFile("outfile.root","RECREATE");
  TGraph waveform;

  TH1F * h_aft_25[60];
  TH1F * h_aft_05[60];
  char hist_name[100];

  for(int hh=0; hh<60; hh++){
     sprintf(hist_name,"h_aft_25_%02d",hh);
     h_aft_25[hh] = new TH1F(hist_name,hist_name,400,-10.,20.);
     sprintf(hist_name,"h_aft_05_%02d",hh);
     h_aft_05[hh] = new TH1F(hist_name,hist_name,400,-10.,20.); 
  }

  TRandom3 rando;
  rando.SetSeed(0);
  WaveformGenerator w;
  double x;
  int n_ph, index;  

  for(int i=0; i<500000; i++) {
    if( i % 10000 == 0 )
       std::cout << "Running event " << i << std::endl;
    
    n_ph = (int) TMath::Floor( rando.Uniform() * 300 );
    index = (int) TMath::Floor(n_ph/5.);

    w.SetPhotonsInCh( n_ph );
    w.GenPhotonArrivalTimes();
    w.GenRandomPhdArea();
    w.GenerateWaveform();
 
    h_aft_25[index]->Fill( w.GetAftT25Samples() );
    h_aft_05[index]->Fill( w.GetAftT05Samples() );

  }

  file->Write();
  file->Close();

  return 0;

}
