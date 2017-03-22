{

char filename_array[14][50] = {"outfile_20.50_3.50.root","outfile_31.20_4.74.root","outfile_47.00_5.68.root","outfile_63.60_4.32.root",
"outfile_10.70_2.96.root","outfile_23.20_4.19.root","outfile_39.60_5.43.root","outfile_52.60_4.39.root","outfile_7.43_2.82.root",
"outfile_15.10_3.59.root","outfile_30.80_4.06.root","outfile_41.70_4.36.root","outfile_57.40_8.34.root","outfile_78.30_7.54.root"};

double means[14], stds[14];
TF1 * fit = new TF1("fit","gaus(0)",-10.,10.);
fit->SetNpx(5000);
TCanvas * cans[14];

for(int i=0; i<14; i++) {
  printf("Opening %s...\n",filename_array[i]);

  cans[i] = new TCanvas(filename_array[i],filename_array[i]);

  TFile * file = new TFile(filename_array[i]);
  TH1F * h_aft_25 = (TH1F *)file->Get("h_aft_25");
  fit->SetParameters( h_aft_25->GetMaximumBin(),
                      h_aft_25->GetMean(),
                      h_aft_25->GetRMS() );
  

  h_aft_25->Fit(fit);
  h_aft_25->GetXaxis()->SetRangeUser(-5.,5.);
  means[i] = fit->GetParameter(1);
  stds[i] = fit->GetParameter(2);
  h_aft_25->Draw(); 
 

}

for(int i=0; i<14; i++) {
   printf("%s\n",filename_array[i]);
   printf("%f\t%f\n",means[i],stds[i]);

}

}
