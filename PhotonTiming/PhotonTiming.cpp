
//
//   PhotonTiming.cpp
//    
//   2015-07-30 - Added an RQ for the earliest detected photon in each
//                channel (Brian L.
//   2015-11-06 - Added an RQ for the total reconstructed area based on
//                the fit photons. (Brian L.)
////////////////////////////////////////////////////////////////////////////





#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <algorithm>

#include <boost/property_tree/ptree.hpp>
#include <boost/foreach.hpp>

#include "../Utilities/LCvtChannel.h"
#include "../Utilities/LCvtEvent.h"
#include "../Utilities/LCvtFile.h"
#include "../Utilities/RQFile_IO.hh"

#include "PhotonFits.h"
#include "TMath.h"

#define _DEBUG_ 0


//#include "Quantities.h"

using namespace std;
using boost::property_tree::ptree;
double small_number = 1e-10;
//______________________________________________________________________________
void FindModule(ptree &full_tree, ptree& module, string module_name) {  
  string dp = "data_processing_settings";
  string tag = "module_name";
  BOOST_FOREACH(ptree::value_type &v, full_tree.get_child(dp)) {
    if (v.second.count(tag)==1) {
      if (v.second.get<string>(tag)==module_name){
        module = v.second;
      }
    }
  }  
}
//______________________________________________________________________________
int main(int argc, char* argv[]) {
 
  if( _DEBUG_ == 1 ) cerr << "Starting the PhotonTiming module" << endl;

 
  //parse inputs_______________________________________________________________
  if( argc !=8  ) {
      cerr << "Please specify:\n-evt_filename, evt_dir, rq_filename, rq_dir, module no., dp_xml_full_path, iq_xml_full_path" << endl;
      return -1;
  }
  
  string evt_dir = string(argv[2]);
  string evt_filename = string(argv[1]);
  string rq_dir = string(argv[4]);
  string rq_filename =string(argv[3]);
  string dp_xml_full_path = string(argv[6]);
  string iq_xml_full_path = string(argv[7]);
  
  // Read input xml files into BOOST ptrees____________________________________
  ptree dp_settings, iqs;  
  
  ifstream dp_settings_ifs(dp_xml_full_path.c_str());
  read_xml(dp_settings_ifs,dp_settings);
  dp_settings_ifs.close();
  
  ifstream iq_settings_ifs(iq_xml_full_path.c_str());
  read_xml(iq_settings_ifs,iqs);
  iq_settings_ifs.close();
  
  // Fill variables from xml __________________________________________________
   
  ptree module_settings;
  FindModule(dp_settings, module_settings, "PhotonTiming");
  cout<<module_settings.get<string>("relative_path")<<endl;  

 
  //Load data into LCvtFile object_____________________________________________
  string cvt_filename(evt_filename);
  cvt_filename.replace(cvt_filename.rfind(".evt"), 4, ".cvt");
  {
    // This code quickly checks if a cvt file exist 
    string name = evt_dir + "/" + cvt_filename;
    ifstream tmp(name.c_str());
    if (!tmp.good()) {
      cerr << "Error: Cannot find cvt file.\n";
      tmp.close();
      return 1;
    }
    tmp.close();
  }  
  //cout << "About to load cvt..." << endl;
  LCvtFile* cvt = new LCvtFile(evt_dir, cvt_filename);
  //cout << "done." << endl;
  

  // Get some module parameters
  int phPerCh;

  try { phPerCh = module_settings.get<int>("parameters.photons_per_ch");}
  catch (exception const &ex) { phPerCh = 5; }
 
 
  //Load data into RQFileIO object_____________________________________________  
  //cout << "About to load RQ file." << endl;
  RQFileIO *rqio = new RQFileIO;
  rqio->ReadFile(rq_dir+"/"+rq_filename);
  //cout << "done." << endl;
  
  // Add the RQs that this module is going to generate__________________________
  size_t pulse_dim=dp_settings.get<size_t>("data_processing_settings.global.max_num_pulses");
  string pulse_dim_str=dp_settings.get<string>("data_processing_settings.global.max_num_pulses");
  string photonsPerPulse_str;
  char photonsPerPulse[10];
  sprintf(photonsPerPulse,",%d",phPerCh*122);
  stringstream ss;
  ss << photonsPerPulse;
  ss >> photonsPerPulse_str;
  cerr << pulse_dim_str + photonsPerPulse_str << endl;
  if (!rqio->events.Get("pulse_template_time"))
    rqio->events.AddRQ("pulse_template_time","float",pulse_dim_str);
  if (!rqio->events.Get("pulse_template_time_err"))
    rqio->events.AddRQ("pulse_template_time_err","float",pulse_dim_str);
  if (!rqio->events.Get("photon_count_by_ch"))
    rqio->events.AddRQ("photon_count_by_ch","int",pulse_dim_str+",122"); 
  if (!rqio->events.Get("photon_count_total"))
    rqio->events.AddRQ("photon_count_total","int",pulse_dim_str); 
  if (!rqio->events.Get("photon_ch_total"))
    rqio->events.AddRQ("photon_ch_total","int",pulse_dim_str+photonsPerPulse_str); 
  if (!rqio->events.Get("photon_times_total"))
    rqio->events.AddRQ("photon_times_total","float",pulse_dim_str+photonsPerPulse_str); 
  if (!rqio->events.Get("photon_times_err_total"))
    rqio->events.AddRQ("photon_times_err_total","float",pulse_dim_str+photonsPerPulse_str);
  if (!rqio->events.Get("photon_fit_area_phe"))
    rqio->events.AddRQ("photon_fit_area_phe","float",pulse_dim_str+photonsPerPulse_str);
//  if (!rqio->events.Get("fit_area_ch_phe"))
//    rqio->events.AddRQ("fit_area_ch_phe","float",pulse_dim_str+",122");
  if (!rqio->events.Get("total_fit_area_phe"))
    rqio->events.AddRQ("total_fit_area_phe","float",pulse_dim_str);
  if (!rqio->events.Get("earliest_photon_time"))
    rqio->events.AddRQ("earliest_photon_time","float",pulse_dim_str);
  if (!rqio->events.Get("earliest_ph_time_by_ch"))
    rqio->events.AddRQ("earliest_ph_time_by_ch","float",pulse_dim_str+",122");
  if (!rqio->events.Get("Gatti_S"))
    rqio->events.AddRQ("Gatti_S","float",pulse_dim_str);
  if (!rqio->events.Get("s1_prompt_to_total_ratio"))
    rqio->events.AddRQ("s1_prompt_to_total_ratio","float",pulse_dim_str);
  if (!rqio->events.Get("photons_allowed_per_ch"))
    rqio->events.AddRQ("photons_allowed_per_ch","int","1");
  if (!rqio->events.Get("photon_times_total_ordered"))
    rqio->events.AddRQ("photon_times_total_ordered","float",pulse_dim_str+photonsPerPulse_str);
  if (!rqio->events.Get("avg_5_25_quantile_time"))
    rqio->events.AddRQ("avg_5_25_quantile_time","float",pulse_dim_str);
  if (!rqio->events.Get("avg_10_30_quantile_time"))
    rqio->events.AddRQ("avg_10_30_quantile_time","float",pulse_dim_str);
  if (!rqio->events.Get("ch_num"))
    rqio->events.AddRQ("ch_num","int",pulse_dim_str+",122"); 

  if( _DEBUG_ ==1 ) printf("\nGrabbing the RQ's I created...\n");
  Rq* pulse_template_time = rqio->events.Get("pulse_template_time");
  Rq* pulse_template_time_err = rqio->events.Get("pulse_template_time_err");
  Rq* photon_count_by_ch = rqio->events.Get("photon_count_by_ch");
  Rq* photon_count_total = rqio->events.Get("photon_count_total");
  Rq* photon_ch_total = rqio->events.Get("photon_ch_total");
  Rq* photon_times_total = rqio->events.Get("photon_times_total");
  Rq* photon_fit_area_phe = rqio->events.Get("photon_fit_area_phe");
  Rq* total_fit_area_phe = rqio->events.Get("total_fit_area_phe");
  Rq* earliest_photon_time = rqio->events.Get("earliest_photon_time");
  Rq* Gatti_S = rqio->events.Get("Gatti_S");
  Rq* s1_prompt_to_total_ratio = rqio->events.Get("s1_prompt_to_total_ratio");
  Rq* pstart = rqio->events.Get("pulse_start_samples");
  Rq* pend = rqio->events.Get("pulse_end_samples");  
  Rq* pulse_class = rqio->events.Get("pulse_classification");
  Rq* pulse_area = rqio->events.Get("pulse_area_phe");
  Rq* earliest_ph_time_by_ch = rqio->events.Get("earliest_ph_time_by_ch");
  Rq* photon_times_err_total = rqio->events.Get("photon_times_err_total");
  Rq* photons_allowed_per_ch = rqio->events.Get("photons_allowed_per_ch");
  Rq* photon_times_total_ordered = rqio->events.Get("photon_times_total_ordered");
  Rq* avg_5_25_quantile_time = rqio->events.Get("avg_5_25_quantile_time");
  Rq* avg_10_30_quantile_time = rqio->events.Get("avg_10_30_quantile_time");
  Rq* ch_num = rqio->events.Get("ch_num");

  if( _DEBUG_ == 1 ) printf("\nGrabbing the timing RQ's...\n");
 
  // grap rqs 
//  Rq* pstart = rqio->events.Get("pulse_start_samples");  
//  Rq* pend = rqio->events.Get("pulse_end_samples");  
  Rq* t0 = rqio->events.Get("hft_t0_samples");
  Rq* t2 = rqio->events.Get("hft_t2_samples");
  
  
  //-----EVENT PROCESSING LOOP BEGINS-----
  double results[20], results_err[20], threshold;
  size_t start, end;
  LCvtChannel *ch;
  vector<float> waveform;
  //cout << "About to loop" << endl;
  for (size_t e=0; e<cvt->cvt_events.size(); e++) {
//  for (size_t e=87; e<90; e++) {
//    size_t events[5] = {26, 43, 45, 57, 58};

//  for( int i=0; i<5; i++){
    rqio->GetEvent(e);
    //cout << "Working on event " << e << " of " << cvt->cvt_events.size() << endl;
//    size_t e = events[i];
//    cerr << "Event " << events[i] << "\t LUXStamp: " << luxstamp->GetInt() << endl;
    photons_allowed_per_ch->SetInt(phPerCh,0); 
 
    //-----PULSE PROCESSING LOOP BEGINS---------
    for (size_t p=0; p<pulse_dim; p++) {

      // Initializing all the important arrays   
      photon_count_total->SetInt(0, p);
      Gatti_S->SetDouble(0., p);
      s1_prompt_to_total_ratio->SetDouble(1000., p);  
      for(size_t phot=0; phot<122*size_t(phPerCh); phot++) {
         photon_times_total->SetDouble(1.e7, p, phot);
         photon_times_total_ordered->SetDouble(1.e7, p, phot);
         photon_ch_total->SetInt(1000,p,phot);
      }
      int photonCounter = 0,
          prevPhotonCount = 0;
      float earliestPhotonTime = 1.e7,
            earliestPhotonInCh = 1.e7;
      float total_fit_area = 0.;
      pulse_template_time->SetDouble(0.,p);

      // Skip this pulse if it's not an S1 or is larger than 200phe
      if ( (pulse_class->GetInt(p) !=1) || // && pulse_class->GetInt(p) != 3 && pulse_class->GetInt(p) != 4) || 
            pulse_area->GetDouble(p) >  200.) continue;

      
      start = t0->GetDouble(p) - pstart->GetInt(p);
      end = t2->GetDouble(p) - pstart->GetInt(p);

//      cerr << "Pulse start samples " << pstart->GetInt(p) << endl;
//      cerr << "Pulse end samples " << pend->GetInt(p) << endl;
//      cerr << "My start samples " << start<< endl;
//      cerr << "My end samples " << end  << endl;

      // Fit the full S1 waveform.
      ch = cvt->cvt_events[e]->channels[136];
      waveform = ch->GetPeak( pstart->GetInt(p)-10, pend->GetInt(p)+10);
//      cerr << "Print summed S1 waveform: " << endl;
//      for(size_t ll=0; ll<waveform.size(); ll++) {
//        cerr << waveform[ll] << endl;
//      }

      FitFullS1(&waveform, start, end, threshold, results, results_err);
      pulse_template_time->SetDouble(results[1] + (float) pstart->GetInt(p) ,p);
      pulse_template_time_err->SetDouble(results_err[1],p);                

      //-----PEAK PROCESSING LOOP BEGINS----
      for (size_t pp=0; pp<122; pp++) {
   
        // grab the waveform I want to fit.
        if( _DEBUG_ == 1) cerr << "Event " << e << endl;
        ch = cvt->cvt_events[e]->channels[pp];
        waveform = ch->GetPeak(pstart->GetInt(p)-10, pend->GetInt(p)+10);
        size_t wavesize = waveform.size();
    
        // We need to plan out how the results thing will go.
        // I'll use the following:
        //     results[0] = number of photons best fit in channel
        //     results[1] = norm ph1
        //     results[2] = t0 ph2
        //     And on from there.

        threshold = 0.025;
        if( wavesize > 3 )   
          results[2] = 1.e7;
//        cerr << "Channel: " << pp << "\t" << "starting point: " << pstart->GetInt(p)-10 << endl;
        Fit(&waveform, 0, wavesize, threshold, results, results_err, pp); 
//        cerr << results[0] << " photons returned by the Fit function" << endl;

    
        photon_count_by_ch->SetInt(TMath::Floor(results[0] + 0.2), p, pp);
        ch_num->SetInt(pp,p,pp);

        prevPhotonCount = photonCounter;

        if( results[0] > 0.9 && results[0] < 1.1){
          photon_times_total->SetDouble((float) results[2] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[2],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[1],p,photonCounter);
          total_fit_area += (float) results[1];
          photonCounter++;

          if(results[2] + (float) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[2] + (float) pstart->GetInt(p);

        } else if( results[0] > 1.9 && results[0] < 2.1 ){
          photon_times_total->SetDouble( (float) results[2] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[2],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[1],p,photonCounter);
          total_fit_area += (float) results[1];
          photonCounter++;

          if(results[2] + (float) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[2] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[4] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[4],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[3],p,photonCounter);
          total_fit_area += (float) results[3];
          photonCounter++;

          if(results[4] + (float) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[4] + (float) pstart->GetInt(p);

        } else if( results[0] > 2.9 && results[0] < 3.1 ){
          photon_times_total->SetDouble( (float) results[2] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[2],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[1],p,photonCounter);
          total_fit_area += (float) results[1];
          photonCounter++;

          if(results[2] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[2] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[4] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[4],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[3],p,photonCounter);
          total_fit_area += (float) results[3];
          photonCounter++;          

          if(results[4] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[4] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[6] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[6],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[5],p,photonCounter);
          total_fit_area += (float) results[5];
          photonCounter++;

          if(results[6] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[6] + (float) pstart->GetInt(p);
       
        } else if( results[0] > 3.9 && results[0] < 4.1 ){
          photon_times_total->SetDouble( (float) results[2] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[2],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[1],p,photonCounter);
          total_fit_area += (float) results[1];
          photonCounter++;

          if(results[2] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[2] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[4] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[4],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[3],p,photonCounter);
          total_fit_area += (float) results[3];
          photonCounter++;          

          if(results[4] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[4] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[6] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[6],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[5],p,photonCounter);
          total_fit_area += (float) results[5];
          photonCounter++;

          if(results[6] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[6] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[8] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[8],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[7],p,photonCounter);
          total_fit_area += (float) results[7];
          photonCounter++;

          if(results[8] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[8] + (float) pstart->GetInt(p);

        } else if( results[0] > 4.9 && results[0] < 5.1 ){
          photon_times_total->SetDouble( (float) results[2] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[2],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[1],p,photonCounter);
          total_fit_area += (float) results[1];
          photonCounter++;

          if(results[2] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[2] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[4] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[4],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[3],p,photonCounter);
          total_fit_area += (float) results[3];
          photonCounter++;          

          if(results[4] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[4] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[6] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[6],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[5],p,photonCounter);
          total_fit_area += (float) results[5];
          photonCounter++;

          if(results[6] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[6] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[8] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[8],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[7],p,photonCounter);
          total_fit_area += (float) results[7];
          photonCounter++;

          if(results[8] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[8] + (float) pstart->GetInt(p);

          photon_times_total->SetDouble( (float) results[10] + (float) pstart->GetInt(p), 
                                        p, photonCounter);
          photon_times_err_total->SetDouble((float) results_err[10],p, photonCounter);
          photon_ch_total->SetInt(pp,p,photonCounter);
          photon_fit_area_phe->SetDouble((float) results[9],p,photonCounter);
          total_fit_area += (float) results[9];
          photonCounter++;

          if(results[10] + (double) pstart->GetInt(p) < earliestPhotonTime) 
               earliestPhotonTime = (float) results[10] + (float) pstart->GetInt(p);
       }

      
       earliestPhotonInCh = 1.e7;  
       for(int kk=1; kk<results[0]; kk++) {
          if (results[2*kk] < earliestPhotonInCh) 
             earliestPhotonInCh = (float) results[2*kk] + (float) pstart->GetInt(p);
       }

       earliest_ph_time_by_ch->SetDouble(earliestPhotonInCh,p,pp); 

      } // end loop over channels
      photon_count_total->SetInt(photonCounter,p);
      earliest_photon_time->SetDouble(earliestPhotonTime,p);
      total_fit_area_phe->SetDouble(total_fit_area,p);

      // Calculate the quantile average times
      float temp_times[1220];
      for( int ph_it=0; ph_it < photonCounter; ph_it++) {
        temp_times[ph_it] = photon_times_total->GetDouble(p,ph_it);
      }
      std::sort(temp_times,temp_times + photonCounter);
      int firstInd_5, firstInd_10, secondInd_25, secondInd_30;

      firstInd_5 = TMath::Floor(photonCounter*0.05);
      firstInd_10 = TMath::Floor(photonCounter*0.1);
      secondInd_25 = TMath::Floor(photonCounter*0.25);
      secondInd_30 = TMath::Floor(photonCounter*0.3);
      
      float sum1 = 0.,sum2 = 0.;
      for(int ph_it=firstInd_5; ph_it<secondInd_25; ph_it++)
        sum1 += temp_times[ph_it];
      for(int ph_it=firstInd_10; ph_it<secondInd_30; ph_it++)
        sum2 += temp_times[ph_it];

      avg_5_25_quantile_time->SetDouble(sum1/double(secondInd_25-firstInd_5),p);
      avg_10_30_quantile_time->SetDouble(sum1/double(secondInd_30-firstInd_10),p);

 

    }// end loop over pulses  
  } // end loop over events
  //-----EVENT PROCESSING LOOP ENDS-----
  
  // Write the output rqs
  rqio->WriteFile(rq_dir+"/"+rq_filename);

  //clean up___________________________________________________________________
  delete cvt;
  delete rqio;

  return 0;
}
