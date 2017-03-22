{

char filename_array[7][50] = {"outfile_31.20_4.74.root","outfile_47.00_5.68.root",
"outfile_23.20_4.19.root","outfile_39.60_5.43.root","outfile_7.43_2.82.root",
"outfile_15.10_3.59.root","outfile_57.40_8.34.root"};

double means[14], stds[14];
TF1 * fit = new TF1("fit","gaus(0)",-10.,10.);


for(int i=0; i<7; i++) {
  printf("Opening %s...\n",filename_array[i]);

  TFile * file = new TFile(filename_array[i]);
  TH1F * h_aft_25 = (TH1F *)file->Get("h_aft_25");
  fit->SetParameters( h_aft_25->GetMaximumBin(),
                      h_aft_25->GetMean(),
                      h_aft_25->GetRMS() );

  h_aft_25->Fit(fit);

  means[i] = fit->GetParameter(1);
  stds[i] = fit->GetParameter(2);
 

}

for(int i=0; i<7; i++) {
   printf("%s\n",filename_array[i]);
   printf("%f\t%f\n",means[i],stds[i]);

}

}
