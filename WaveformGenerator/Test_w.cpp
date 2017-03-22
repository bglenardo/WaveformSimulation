#include "WaveformGenerator.hh"
#include "SPEFunctions.hh"
#include "TH1F.h"
#include "TFile.h"
#include <vector>
#include "TRandom3.h"

using namespace std;

int main() {

  WaveformGenerator w;
  TRandom3 rand;
  rand.SetSeed(0);
  double offset;
  offset = rand.Uniform();

  TFile * outfile = new TFile("outfile_no_sig.root","RECREATE"); 
  TH1F * hist = new TH1F("hist","hist",400,-20.,40.);

  w.SetSig(0.);
  w.SetSingletFraction(0.);
  w.SetPhotonsInCh(2000000);
  w.GenPhotonArrivalTimes();
  
  vector<double> vec = w.GetSortedTimes();

  for(int i=0; i<2000000; i++) {
     if( i % 10000 == 0 ) printf("Event: %d\n",i);
     if( i % 100 == 0 ) { offset = 0.; } // rand.Uniform(); }//printf("offset = %f\n",offset); }
     hist->Fill(vec[i] + offset);
  }

  outfile->Write();
  outfile->Close();

  return 0;

}  
