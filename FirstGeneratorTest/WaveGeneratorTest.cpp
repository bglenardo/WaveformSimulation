#include "WaveformGenerator.hh"
#include "TGraph.h"
#include "TFile.h"
#include <vector>

void Test() {

  WaveformGenerator w;
  w.SetPhotonsInArray(120);
  w.GenerateNewWaveform();
  std::vector<double> times = w.GetSortedTimes();
  std::vector<double> areas = w.GetAreas();
  for(int i=0; i < times.size(); i++){
    printf("%f \t %f\n",times[i], areas[i]);
  }
  
  TGraph g = w.GetWaveform();


  g.SetName("g");
  g.SetTitle("g");

  TFile * outfile = new TFile("WaveformGeneratorTest.root","RECREATE");

  g.Write();
  outfile->Close();


}

int main() {

  Test();

  return 1;
}

