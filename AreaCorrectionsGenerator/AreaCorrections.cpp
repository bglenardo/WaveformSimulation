#include "../PhotonTiming/PhotonFits.h"
#include "../WaveformGenerator/WaveformGenerator.hh"
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include "TH1.h"
#include "../PhotonTiming/SPEFunctions_Rose_test.h"
#include <ctime>

std::vector<int> AssignTrueToFitPhotons( std::vector<double> true_times, std::vector<double> fit_times);

int main() {

   time_t start_time;
   time_t end_time;
   time(&start_time);
   bool print_assignments = false;
   TFile * outfile = new TFile("AreaCorrections.root","RECREATE");
   
   double results[20];
   double results_err[20];
   double threshold;
   int n_phot,
       n_fits_near_phot,
       n_other_near_phot,
       fit_area_index;

   std::vector<double> v_waveform;
   std::vector<double> true_times;
   std::vector<double> fit_times;
   std::vector<double> fit_areas;
   std::vector<int> fit_num; 

   

   TH1D * fit_area_bins[100];
   char hist_name[100];
   for(int hh=0; hh<100; hh++){
      sprintf(hist_name,"h_fit_area_%02d",hh);
      fit_area_bins[hh] = new TH1D(hist_name,hist_name,200,0.,20.);
   }

      WaveformGenerator w;
      w.SetPhotonsInArray(100);

   for(int ii=0; ii<1000; ii++) {
      if( ii % 100 == 0 ) printf("Running %d...\n",ii);
      true_times.clear();
      fit_times.clear();
      fit_num.clear();
      fit_areas.clear();

      w.GenerateNewWaveform();
      v_waveform = w.GetWaveVec();
      true_times = w.GetSortedTimes();

      Fit( v_waveform, w.GetTraceStart(), 0.01, results, results_err);
      n_phot = TMath::Floor(results[0]+0.5);
      n_fits_near_phot = 0.; n_other_near_phot = 0.;
      for(int ph=0; ph<n_phot; ph++) {
            int ind = ph*2 + 1; 
            fit_areas.push_back(results[ind]);
            fit_times.push_back(results[ind+1]);  
      }
      if( n_phot > 0 )
        fit_num = AssignTrueToFitPhotons( true_times, fit_times );
      if( print_assignments ) {
         printf("True times: ");
         for(int t=0; t<true_times.size(); t++) printf("\t%4.4f",true_times[t]);
         printf("\nFit times: ");
         for(int f=0; f<fit_times.size(); f++) printf("\t%4.4f",fit_times[f]);
         printf("\nFit areas: ");
         for(int f=0; f<fit_areas.size(); f++) printf("\t%4.4f",fit_areas[f]);
         printf("\nAssignments: ");
         for(int f=0; f<fit_num.size(); f++) printf("\t%d",fit_num[f]); 
         printf("\n");
     }

     if( fit_areas.size() > 0 ) {
       for(int ph=0; ph<n_phot; ph++){
          fit_area_index = TMath::Floor( fit_areas[ph]/0.2 );
          fit_area_bins[fit_area_index]->Fill( fit_num[ph] );
       }
     }
         
     
   
   }




   outfile->Write();
   outfile->Close();

   time(&end_time);
   double seconds = difftime(end_time,start_time);
   printf("Run time: %4.4f minutes.\n",seconds/60.);


   return 1;
}


std::vector<int> AssignTrueToFitPhotons( std::vector<double> true_times, std::vector<double> fit_times) {

  double dt, min_dt;
  int min_dt_ind = 0;
  std::vector<int> assignments;
  for(int t=0; t<fit_times.size(); t++){
    assignments.push_back(0.);
  }


  for(int t=0; t<true_times.size(); t++) {
     min_dt = 1000.;
     for(int f=0; f<fit_times.size(); f++) {
        dt = TMath::Abs( fit_times[f] - true_times[t] );
        if( dt < min_dt ) {
           min_dt = dt;
           min_dt_ind = f;
        }
     }
     assignments[min_dt_ind] += 1.;
  }

  return assignments;
}





