#include "../WaveformGenerator/WaveformGenerator.hh"
#include "TGraph.h"
#include "TFile.h"
#include <vector>

void Test() {

  WaveformGenerator w;
  w.SetPhotonsInArray(1200);
  w.GenerateNewWaveform();
  std::vector<double> times = w.GetSortedTimes();
  std::vector<double> times2 = w.GetUncorrectedSortedTimes();
  std::vector<double> areas = w.GetAreas();
  for(int i=0; i < times.size(); i++){
    printf("%f \t %f \t %f\n",times[i], times2[i], areas[i]);
  }
  
  TGraph g = w.GetWaveform();
  TGraph g2 = w.GetUncorrectedWaveform();

  g.SetName("g");
  g.SetTitle("g");

  g2.SetName("g_un");

  TFile * outfile = new TFile("WaveformGeneratorTest.root","RECREATE");

  g.Write();
  g2.Write();
  outfile->Close();


}

int main() {

  Test();

  return 1;
}

