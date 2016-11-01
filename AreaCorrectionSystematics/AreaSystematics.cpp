#include "../PhotonTiming/PhotonFits.h"
#include "../WaveformGenerator/WaveformGenerator.hh"
#include <vector>
#include "TGraph.h"
#include "TFile.h"
#include "TH1.h"
#include "TMath.h"
#include "TTree.h"
#include "../PhotonTiming/SPEFunctions_Rose_test.h"
#include <ctime>
#include <sstream>
#include <string>
#include <iostream>
#include "TRandom3.h"
#include "Area_Corrections_NewMC.h"

std::vector<int> AssignTrueToFitPhotons( std::vector<double> true_times, std::vector<double> fit_times);
double SampleFromPhdDist();

int main(int argc, char * argv[]) {

   if(argc > 2) {
        std::cerr << "Too many inputs...\nExiting...." << std::endl;
        return 1;
   } else if (argc < 2) {
        std::cerr << "Please specify run #...\nExiting..." << std::endl;
        return 1;
   }

   int Nevents=5;

   // Measure runtime
   time_t start_time;
   time_t end_time;
   time(&start_time);
   bool print_assignments = false;
   
   // Generate output file
   char outfile_name[100];
   sprintf(outfile_name,"./AreaCorSystematics_10-120_%d.root",atoi(argv[1]));
   TFile * outfile = new TFile(outfile_name,"RECREATE");
   std::cout << "Output file: "  << outfile_name << std::endl;
   
   // Create TTree to store important variables
   double pulse_area,
          photon_times[300],
          best_weights[300],
          fit_area[300]; 
   int photon_count,
       true_photon_count;
   double pulse_area_mean, pulse_area_width;

   TTree * outtree = new TTree("data","data");
   outtree->Branch("pulse_area",&pulse_area,"pulse_area/D");
   outtree->Branch("photon_count",&photon_count,"photon_count/I");
   outtree->Branch("true_photon_count",&true_photon_count,"true_photon_count/I");
   outtree->Branch("photon_times",&photon_times,"photon_times[300]/D");
   outtree->Branch("best_weights",&best_weights,"best_weights[300]/D");
   outtree->Branch("fit_area",&fit_area,"fit_area[300]/D");
   
   // Create variables to store fitter output. 
   double results[20];
   double results_err[20];
   double threshold;
   int n_phot,
       n_fits_near_phot,
       n_other_near_phot,
       fit_area_index,
       pulse_area_index;

   // Create waveform variables.
   std::vector<double> v_waveform;
   std::vector<double> true_times;
   std::vector<double> fit_times;
   std::vector<double> fit_areas;
   std::vector<int> fit_num; 
   WaveformGenerator w;

   // Create random number generator
   TRandom3 rand;
   rand.SetSeed(0);

   // Begin loop.
   for(int ph=10; ph<120; ph++) {
     std::cout << "Pulse size: " << ph << std::endl;
     for(int event=0; event<Nevents; event++) {
        std::cout << "Event: " << event << std::endl;
        // Get true photon count
        w.SetPhotonsInArray( ph );
        true_photon_count = ph;
  
        // Calculate the pulse area
        pulse_area = 0.;
        for(int i=0; i<ph; i++) pulse_area += SampleFromPhdDist();
        pulse_area_index = TMath::Floor( pulse_area / 5. );
  
        // Create waveforms representing photons up to true photon count
        int photons_generated = 0, photons_fit = 0;
        photon_count = 0;
        for(int i=0; i<300; i++) {
           photon_times[i] = -1000.;
           fit_area[i] = -1000.;
           best_weights[i] = -1000.;
        }
  
        while( photons_generated < ph ) {
                  
            w.GenPhotonsInCh();
            //printf("PhotonsInCh = %d\n",w.GetPhotonsInCh());
            if( w.GetPhotonsInCh() > (ph - photons_generated) )  continue;
  
            w.GenPhotonArrivalTimes();
            w.GenRandomPhdArea();
            w.GenerateWaveform();
  
            v_waveform = w.GetWaveVec();
            Fit( v_waveform, w.GetTraceStart()+123.7, 0.01, results, results_err);
            for(int i=0; i<TMath::Floor(results[0]+0.5); i++) {
              int ind = i*2 + 1; 
              fit_area[photon_count + i] = results[ind];
              photon_times[photon_count + i] = results[ind+1]-123.7;
              fit_area_index = TMath::Floor(results[ind]/0.2);
              best_weights[photon_count + i] = area_correction[pulse_area_index][fit_area_index];
  
            }     
            photon_count += TMath::Floor(results[0]+0.5);
            photons_generated += w.GetPhotonsInCh();
        }
        outtree->Fill();
     }
   }

   outfile->Write();
   outfile->Close();

   time(&end_time);
   double seconds = difftime(end_time,start_time);
   std::cout << "Run time: " << seconds/60. << " minutes." << std::endl;


   return 1;
}

/*************************************************************************
    AssignTrueToFitPhotons
**************************************************************************/

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

/*************************************************************************
    SampleFromPhdDist()
**************************************************************************/
double SampleFromPhdDist() {

  TRandom3 r;
  r.SetSeed(0);
  double mu = 0.8327873877969729;
  double spe_sig = 0.35;
  double temp_area = 0.; 

    if(r.Uniform() < 0.2)
      temp_area = r.Gaus( 2 * mu, TMath::Sqrt(2) * spe_sig );
    else
      temp_area = r.Gaus( mu, spe_sig);

     return temp_area;
   
}
