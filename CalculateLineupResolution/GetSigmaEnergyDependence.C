{


  TF1 * fit = new TF1("fit","[0]/(2)*( TMath::Erf((x-[2])/sqrt(2*[1]*[1])) + TMath::Erf((1.-(x-[2]))/sqrt(2*[1]*[1])) )",-10.,10.);
  TGraph * g_sigma = new TGraph();  g_sigma->SetName("g_sigma");
  TGraph * g_offset = new TGraph(); g_offset->SetName("g_offset");
  TH1F * hist = new TH1F();

  TFile * infile = new TFile("outfile.root");
  char histname[100];

  for(int i=0; i<60; i++) {
     hist->Reset();
     sprintf(histname,"h_aft_05_%02d",i);
     hist = (TH1F *) infile->Get(histname);
     fit->SetParameters(hist->GetBinContent(hist->GetMaximumBin()),1.6/TMath::Sqrt( (double) (i) * 5. + 2.5 ),-0.8);
     
     hist->Fit(fit,"QM","QM",-5.,5.);
     g_sigma->SetPoint(i,(double) (i) * 5. + 2.5, TMath::Abs(fit->GetParameter(1)) );
     g_offset->SetPoint(i,(double) (i) * 5. + 2.5, fit->GetParameter(2) );


  }

  g_sigma->SetMarkerStyle(21);
  g_sigma->SetMarkerSize(1.5);
  g_sigma->GetXaxis()->SetTitle("Number of photons in pulse");
  g_sigma->GetYaxis()->SetTitle("Sigma");

  g_sigma->Draw("AP");

  TF1 * res_fit = new TF1("res_fit","[0]/sqrt(x-[1])",0.1,500.);
  res_fit->SetParameters(1.3,0.);
  res_fit->FixParameter(1,0.); 
   
  g_sigma->Fit(res_fit,"","",1.,300.);
  res_fit->SetLineColor(2);
  res_fit->SetNpx(3000);
  res_fit->Draw("SAME");

  TCanvas * c2 = new TCanvas();
  g_offset->SetMarkerStyle(21);
  g_offset->SetMarkerSize(1.5);
  g_offset->GetXaxis()->SetTitle("Number of photons in pulse");
  g_offset->GetYaxis()->SetTitle("offset");
  g_offset->Draw("AP");


}
