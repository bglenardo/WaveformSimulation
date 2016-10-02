#include "../PhotonTiming/PhotonFits.h"
#include "../WaveformGenerator/WaveformGenerator.hh"
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include "../PhotonTiming/SPEFunctions_Rose_test.h"


int main() {

   double results[20];
   double results_err[20];
   double threshold;


   TGraph g_waveform;
   std::vector<double> v_waveform;
   std::vector<double> times;

   WaveformGenerator w;
   w.SetPhotonsInArray(120);
   w.GenerateNewWaveform();
   g_waveform = w.GetWaveform();
   v_waveform = w.GetWaveVec();
   Fit( v_waveform, w.GetTraceStart(), 0.01, results, results_err);

   int n_phot = TMath::Floor(results[0]+0.5);

   printf("True num photons = %d, Fitted num photons = %d\n",w.GetPhotonsInCh(),n_phot);
   TGraph * g_fittedWaveform = new TGraph();

   printf("True photon times: ");
   times = w.GetSortedTimes();
   for(int i=0; i<w.GetPhotonsInCh(); i++) {
      printf("\t%f",times[i]);
   }
   printf("\n");

   printf("Fitted times: ");
   for(int i=0; i<n_phot; i++) {
     int ind = i*2 + 2;
     printf("\t%f",results[ind]);
   }
   printf("\n");

   double x[1], y, par[2];

   for(int i=0; i<5000; i++){

      x[0] = (double) i/100. - 20.;
      y = 0.;
      for(int j=0; j<n_phot; j++) {
          int ind = j*2 + 1;
          par[0] = results[ind];
          par[1] = results[ind+1];
          y += singleSPE( x, par );             
      }
      g_fittedWaveform->SetPoint(i,x[0],y);
   }

   

   TFile * outfile = new TFile("FitTest.root","RECREATE");
   g_fittedWaveform->SetName("fit");
   g_fittedWaveform->SetLineColor(1);
   g_fittedWaveform->Write();

   g_waveform.SetName("data");
   g_waveform.SetMarkerStyle(21);
   g_waveform.SetLineColor(4);
   g_waveform.SetMarkerColor(4);
   g_waveform.SetMarkerSize(0.5);
   g_waveform.Write();

   outfile->Close();

   return 1;
}
