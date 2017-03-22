{


double means[14], stds[14];
TF1 * fit = new TF1("fit",
                    "[0]/2.*(TMath::Erf((x-[1])/(sqrt(2)*[2])) + TMath::Erf((1-(x-[1]))/sqrt(2)/[2]))",
                    -10.,10.);
fit->SetNpx(5000);
TCanvas * cans[14];
TFile * file = new TFile("tt_outfile_3.root");


for(int i=0; i<6; i++) {
  printf("Opening file...\n");

  char cname[50];
  sprintf(cname,"c%d",i);
  cans[i] = new TCanvas(cname,cname);

  char hname[50];
  sprintf(hname,"tt_h_aft_05_%02d",i);
  TH1F * h_aft_25 = (TH1F *)file->Get(hname);
  fit->SetParameters( h_aft_25->GetMaximumBin()*6.,
                      -0.9,
                      h_aft_25->GetRMS()/2. );
  

  h_aft_25->Fit(fit);
  h_aft_25->GetXaxis()->SetRangeUser(-5.,5.);
  means[i] = fit->GetParameter(1);
  stds[i] = fit->GetParameter(2);
  h_aft_25->Draw(); 
  //fit->Draw("LSAME");

}

for(int i=0; i<6; i++) {
   printf("%f\t%f\n",means[i],stds[i]);

}

}
