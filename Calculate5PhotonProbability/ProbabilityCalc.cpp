#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TFile.h"
#include "TH1.h"
#include <iostream>
#include "TRandom3.h"
#include "TMath.h"
#include <vector>

using namespace std;

int main(int argc, char * argv[]) {

  int N = 2000000;

  TFile * file = new TFile("outfile.root","RECREATE");
  TH1F * h_nph = new TH1F("h_nph","h_nph",20,0.,20.);


  vector<double> sorted_times;

  TRandom3 rando;
  rando.SetSeed(0);
  WaveformGenerator w;
  double x;
  int n_ph, index;  

  int num_without_5_independent_hits = 0;

  for(int i=0; i<N; i++) {
    if( i % 1000000 == 0 )
       cout << "Running event " << i << endl;
    
    w.SetPhotonsInArray( 100 );
    w.GenPhotonsInCh();
    w.GenPhotonArrivalTimes();
    sorted_times = w.GetSortedTimes();

    h_nph->Fill( (double) w.GetPhotonsInCh() );

    if (sorted_times.size() > 4) {
      for(int i=1; i<sorted_times.size(); i++) {
        if( TMath::Abs((sorted_times[i] - sorted_times[i-1])) < 1 ) {
           num_without_5_independent_hits++;
           break;           
        }   
      }
    } else {
        num_without_5_independent_hits++;
    }

  }


  cout << "The probability of seeing 5 independent photons (none overlapping within 1 sample) is: " << endl;
  cout << "\t" << (double) (N - num_without_5_independent_hits)/N << endl;

  file->Write();
  file->Close();

  return 0;

}
