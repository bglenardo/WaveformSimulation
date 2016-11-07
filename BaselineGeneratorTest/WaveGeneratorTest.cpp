#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TGraph.h"
#include "TFile.h"
#include <vector>

void Test() {

  WaveformGenerator w;
  w.SetPhotonsInArray(90);
  w.GenPhotonsInCh();
  w.GenPhotonArrivalTimes();
  w.GenRandomPhdArea();
  w.GenerateWaveform();
//
//
//
//
//  w.GenerateNewWaveform();
//  std::vector<double> times = w.GetSortedTimes();
//  std::vector<double> areas = w.GetAreas();
//  for(int i=0; i < times.size(); i++){
//    printf("%f \t %f\n",times[i], areas[i]);
//  }
//  
  TGraph g = w.GetWaveform();
  TGraph b = w.GetBaseline();

  g.SetName("g");
  g.SetTitle("g");

  b.SetName("b");
  b.SetTitle("b");

  TFile * outfile = new TFile("WaveformGeneratorTest_30Mhz_cutoff.root","RECREATE");

  g.Write();
  b.Write();
  outfile->Close();

}

int main() {

  Test();

  return 1;
}

