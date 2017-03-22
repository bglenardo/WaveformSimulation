{

//char filename_array[14][50] = {"outfile_20.50_3.50.root","outfile_31.20_4.74.root","outfile_47.00_5.68.root","outfile_63.60_4.32.root",
//"outfile_10.70_2.96.root","outfile_23.20_4.19.root","outfile_39.60_5.43.root","outfile_52.60_4.39.root","outfile_7.43_2.82.root",
//"outfile_15.10_3.59.root","outfile_30.80_4.06.root","outfile_41.70_4.36.root","outfile_57.40_8.34.root","outfile_78.30_7.54.root"};

char tt_filename[7][50] = {
"outfile_7.43_2.82.root",
"outfile_15.10_3.59.root",
"outfile_23.20_4.19.root",
"outfile_31.20_4.74.root",
"outfile_39.60_5.43.root",
"outfile_47.00_5.68.root",
"outfile_57.40_8.34.root"};

char dd_filename[7][50] = {
"outfile_10.70_2.96.root",
"outfile_20.50_3.50.root",
"outfile_30.80_4.06.root",
"outfile_41.70_4.36.root",
"outfile_52.60_4.39.root",
"outfile_63.60_4.32.root",
"outfile_78.30_7.54.root"};

bool gaussian = false;

double means[7], stds[7], exp[7];
TF1 * fit = new TF1("fit","[0]/(2*[1])*exp([2]*[2]/(2*[1]))*exp(-(x-[3])/[1])*(1 + TMath::Erf(( (x-[3]) - [2]*[2]/[1])/sqrt(2)/[2]))",-10.,10.);
TF1 * g_fit = new TF1("g_fit","gaus(0)",-10.,10.);

for(int i=0; i<7; i++) {
  printf("Opening %s...\n",tt_filename[i]);
  if(i > 0) break;
  TFile * file = new TFile(tt_filename[i]);
  TH1F * h_aft_25 = (TH1F *)file->Get("h_aft_25");
  fit->SetParameters( h_aft_25->GetEntries(),
                      0.2,
                      h_aft_25->GetRMS(),
                      h_aft_25->GetMean() );
  g_fit->SetParameters(h_aft_25->GetMaximumBin(),
                       h_aft_25->GetMean(),
                       h_aft_25->GetRMS());
  if( i<5)
    h_aft_25->Fit(fit,"","",-5.,5.);
  else
    h_aft_25->Fit(g_fit,"","",-5.,5.);


  char status[50];
  sprintf(status,"%s",gMinuit.fCstatu.Data());
  printf(status);

//  if( status == "CONVERGED\t" ) 
//    printf("Hooray!\n");
//  else
//    printf("Crap!\n");
//  printf("%s",gMinuit.fCstatu.Data());

  if( i < 5 ){
    means[i] = fit->GetParameter(3);
    stds[i]  = fit->GetParameter(2);
    exp[i]   = fit->GetParameter(1);
  } else {
    means[i] = g_fit->GetParameter(1);
    stds[i] = g_fit->GetParameter(2);
    exp[i] = 0.0;
  }
 

}

for(int i=0; i<7; i++) {
   printf("%s\n",dd_filename[i]);
   printf("%f\t%f\t%f\n",means[i],stds[i],exp[i]);

}

}
