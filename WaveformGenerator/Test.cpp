#include "../PhotonTiming/PhotonFits.h"
#include "../WaveformGenerator/WaveformGenerator.hh"
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include "TH1.h"
#include "../PhotonTiming/SPEFunctions_Rose_test.h"
#include <ctime>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>

std::vector<int> AssignTrueToFitPhotons( std::vector<double> true_times, 
                                         std::vector<double> fit_times,
                                         std::vector<double> fit_areas);
int GetNumSpikes(TGraph graph);

int main(int argc, char * argv[]) {

   int photons_in_array, run_num;
   if(argc > 1) {
     std::stringstream arg_in( argv[1] );
     if(!(arg_in >> photons_in_array)){
        std::cerr << "Invalid number of photons...\nExiting...." << std::endl;
        return 1;
     }
     std::stringstream run_arg_in( argv[2] );
     if(!(run_arg_in >> run_num)){
        std::cerr << "Invalid run number...\nExiting..." << std::endl;
        return 1;
     }
   }
   else
      photons_in_array = 150;
   
   std::cout << "\nPhotons in array: "  << photons_in_array << "\n" << std::endl;

   time_t start_time;
   time_t end_time;
   time(&start_time);
   bool print_assignments = true; //false;
   char outfile_name[100];
   sprintf(outfile_name,"graph.root");
   TFile * outfile = new TFile(outfile_name,"RECREATE");
   
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

   
      WaveformGenerator w;
      w.SetPhotonsInArray( photons_in_array );
      w.SetT3(3.5);

   for(int ii=0; ii<1; ii++) {
//     for(int ii=0; ii<100; ii++) {
      if( ii % 100== 0 ) std::cout << "Running " << ii << "..." << std::endl;
      true_times.clear();
      fit_times.clear();
      fit_num.clear();
      fit_areas.clear();

      std::vector<double> my_times;
      my_times.push_back(-0.75);
      my_times.push_back(1.6);
      my_times.push_back(6.5);

      w.SetPhotonsInCh(3);
      w.SetSortedTimes(my_times);
      w.GenRandomPhdArea();
      w.GenerateWaveform();
      v_waveform = w.GetWaveVec();
      true_times = w.GetSortedTimes();

      Fit( v_waveform, w.GetTraceStart(), 0.025, results, results_err);
      n_phot = TMath::Floor(results[0]+0.5);
      n_fits_near_phot = 0.; n_other_near_phot = 0.;
      for(int ph=0; ph<n_phot; ph++) {
            int ind = ph*2 + 1;
            if(results[ind] > 0.15) { 
              fit_areas.push_back(results[ind]);
              fit_times.push_back(results[ind+1]);  
            }
      }
      if( print_assignments ) {
         printf("True times: ");
         for(int t=0; t<true_times.size(); t++) printf("\t%4.4f",true_times[t]);
         printf("\nFit times: ");
         for(int f=0; f<fit_times.size(); f++) printf("\t%4.4f",fit_times[f]);
         printf("\nFit areas: ");
         for(int f=0; f<fit_areas.size(); f++) printf("\t%4.4f",fit_areas[f]);
     }

   }

   printf("\n");

   TGraph graph = w.GetWaveform();
   graph.SetName("graph");

   printf("Num spikes: %d\n",GetNumSpikes(graph));

    graph.Write();

   outfile->Close();

   time(&end_time);
   double seconds = difftime(end_time,start_time);
   std::cout << "Run time: " << seconds/60. << " minutes." << std::endl;


   return 1;
}


std::vector<int> AssignTrueToFitPhotons( std::vector<double> true_times, 
                                         std::vector<double> fit_times, 
                                         std::vector<double> fit_areas) {

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
        if( dt < min_dt  && fit_areas[f] > 0.16) {
           min_dt = dt;
           min_dt_ind = f;
        }
     }
     assignments[min_dt_ind] += 1.;
  }

  return assignments;
}


int GetNumSpikes( TGraph graph ) {

  int num_spikes = 0;
  double x,y,yprev,ynext;
  for(int i=1; i<graph.GetN()-1; i++) {

     graph.GetPoint(i-1,x,yprev);
     graph.GetPoint(i,x,y);
     graph.GetPoint(i+1,x,ynext);
   
     if( y > yprev && y >= ynext && y >= 0.08 )
        num_spikes++; 

  }

  return num_spikes;

}


