{


double means[14], stds[14];
TF1 * fit = new TF1("fit","gaus(0)",-10.,10.);
fit->SetNpx(5000);
TCanvas * cans[14];
TFile * file = new TFile("energy_bin_sigmas_offsets_it02_transit_time_mod.root");


for(int i=0; i<6; i++) {
  printf("Opening file...\n");

  char cname[50];
  sprintf(cname,"c%d",i);
  cans[i] = new TCanvas(cname,cname);

  char hname[50];
  sprintf(hname,"dd_h_aft_05_%02d",i);

  TH1F * h_aft_05 = (TH1F *)file->Get(hname);
  fit->SetParameters( h_aft_05->GetMaximumBin(),
                      h_aft_05->GetMean(),
                      h_aft_05->GetRMS() );
  

  h_aft_05->Fit(fit);
  h_aft_05->GetXaxis()->SetRangeUser(-5.,5.);
  means[i] = fit->GetParameter(1);
  stds[i] = fit->GetParameter(2);
  h_aft_05->Draw(); 
 

}

for(int i=0; i<6; i++) {
   printf("%f\t%f\n",means[i],stds[i]);

}

}
