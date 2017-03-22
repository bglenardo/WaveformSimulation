#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include <cstdlib>



int main(int argc, char * argv[]) {
  
  if (argc != 2) {
      std::cout << "Incorrect number of arguments! " << std::endl;
      std::cout << "\nUsage:" << std::endl;
      std::cout << "\tLineupResolutionSingleBin <bin>" << std::endl;
      return 1;
  }

  double bin = atof(argv[1]);

  char outfile_name[50];
  sprintf(outfile_name,"outfile_%2.2f_dd.root",bin);

  TFile * file = new TFile(outfile_name,"RECREATE");
  TGraph waveform;

  TH1F * h_aft_25 = new TH1F("h_aft_25","h_aft_25",400,-10.,20.);
  TH1F * h_aft_05 = new TH1F("h_aft_05","h_aft_05",400,-10.,20.);

  TRandom3 rando;
  rando.SetSeed(0);
  WaveformGenerator w;
    w.SetSingletFraction(1. - 1./(1 + 2.20*0.389/2.46));
    w.SetT1(0.389);
    w.SetT3(2.46);
    w.SetSig(0.365);
    w.SetTr(1.08);
  double x;
  int n_ph, index;  

  for(int i=0; i<10000; i++) {
    if( i % 100 == 0 )
       std::cout << "Running event " << i << std::endl;
    
    n_ph = (int) TMath::Floor( rando.Uniform()*10. + bin );
   

    w.SetPhotonsInCh( n_ph );
    w.GenPhotonArrivalTimes();
    w.GenRandomPhdArea();
    w.GenerateWaveform();
 
    h_aft_25->Fill( w.GetAftT25Samples() );
    h_aft_05->Fill( w.GetAftT05Samples() );

  }

  file->Write();
  file->Close();

  return 0;

}
