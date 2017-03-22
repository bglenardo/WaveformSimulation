#include "TH1F.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TRandom3.h"
#include "TMath.h"
#include "Area_Corrections_FitCut_0_15_11-18-16.h"

int main() {

  TChain * chain = new TChain("data");
  chain->Add("batch_jobs/bath_11_22_2016/Area*.root");

  double photon_times[300],
         fit_area[300],
         best_weights[300],
         pulse_area;
  int    photon_count,
         true_photon_count;

  chain->SetBranchAddress("photon_count",&photon_count);
  chain->SetBranchAddress("true_photon_count",&true_photon_count);
  chain->SetBranchAddress("fit_area",&fit_area);
  chain->SetBranchAddress("best_weights",&best_weights);
  chain->SetBranchAddress("pulse_area",&pulse_area);
  chain->SetBranchAddress("photon_times",&photon_times); 

  TFile * outfile = new TFile("fit_area_104_sim_histograms_11_22_2016_low_stats.root","RECREATE");
  TH1F * avg_shape[9];
  TH1F * delta_t[9];
  double true_num_photons_in_bin[9];
  for(int i=0; i<9; i++) true_num_photons_in_bin[i] = 0.;
  char histname[100];
  printf("initializing histograms...\n");
  for(int i=0; i<9; i++) {
     sprintf(histname,"h_avg_shape_%d",i+1);
     avg_shape[i] = new TH1F(histname,histname,1400,-50.,130.);
     sprintf(histname,"h_delta_t_%d",i+1);
     delta_t[i] = new TH1F(histname,histname,500,0.,50.);
  }
  TH2F * h_countVsArea_cor = new TH2F("h_countVsArea_cor",
                                      "h_countVsArea_cor",
                                       120,0.,120.,120,0.,120.);
  TH2F * h_countVsArea_fitArea_cor = new TH2F("h_countVsArea_fitArea_cor",
                                              "h_countVsArea_fitArea_cor",
                                              120,0.,120.,120,0.,120.);
  TH2F * h_countVsArea = new TH2F("h_countVsArea",
                                  "h_countVsArea",
                                  120,0.,120.,120,0.,120.);
  TH2F * h_fit_area_vs_avg_shape = new TH2F("h_fit_area_vs_avg_shape",
                                            "h_fit_area_vs_avg_shape",
                                            700,-30.,40.,300,0.,7.5);

  printf("Histograms initialized.\n");

  int pulse_area_index,
      hist_index;
  double corrected_count,
         fit_corrected_count;
         
  TRandom3 r;
  r.SetSeed(0);

  int num_entries = chain->GetEntries();
  printf("Number of events: %d\n",num_entries);
  
  for(int i=0; (double)i<num_entries/10.; i++) {
     corrected_count = 0.;
     fit_corrected_count = 0.;
     if( i % 10000 == 0 ) printf("Processing entry %d\n",i);

     chain->GetEntry(i);
     pulse_area_index = TMath::Floor(pulse_area / 5.);
     hist_index = TMath::Floor(pulse_area/10. - 1.);
     if(pulse_area < 100. && pulse_area > 10.) 
        true_num_photons_in_bin[hist_index] += (double) true_photon_count; 


     for(int j=0; j<photon_count; j++) {

       int fit_area_index = TMath::Floor((fit_area[j])/0.2);
       int fit_area_index_2;


       corrected_count += fit_area[j]*1.04;//area_correction[fit_area_index][pulse_area_index]*1.02;
       fit_corrected_count += fit_area[j]; 
       if(fit_area[j] < 0.25) continue;
       if(pulse_area < 100. && pulse_area > 10.) {
         avg_shape[hist_index]->Fill( r.Uniform() + 
                                             r.Gaus(0.,1.693/TMath::Sqrt((double)true_photon_count)) +
                                             -0.345*TMath::Power(((double)true_photon_count),-0.699) - 0.729 +
//                                             photon_times[j], area_correction[fit_area_index][pulse_area_index]*1.02);
                                             photon_times[j], fit_area[j]*1.04);
 //                                            photon_times[j]);
         h_fit_area_vs_avg_shape->Fill(photon_times[j],fit_area[j]);
         for(int jj=j+1; jj<photon_count; jj++) {
             if(fit_area[jj] < 0.2) continue;
             fit_area_index_2 = TMath::Floor((fit_area[jj])/0.2);
             delta_t[hist_index]->Fill(TMath::Abs( photon_times[j]-photon_times[jj] ),
        //                                     (area_correction[fit_area_index][pulse_area_index] + 
        //                                      area_correction[fit_area_index_2][pulse_area_index]));
                                            (fit_area[jj] + fit_area[j])*1.04 );
                                             //2.);
         }

       }
       
     }     
     h_countVsArea_cor->Fill(pulse_area,corrected_count);
     h_countVsArea_fitArea_cor->Fill(pulse_area,fit_corrected_count);
     h_countVsArea->Fill(pulse_area,(float)photon_count);
  }

  for(int i=0; i<9; i++) {
     avg_shape[i]->Sumw2();
//     avg_shape[i]->Scale(1./true_num_photons_in_bin[i]);
  } 
 
  outfile->Write();
  outfile->Close();  

}
