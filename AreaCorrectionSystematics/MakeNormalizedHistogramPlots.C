{

  TRandom3 r; r.SetSeed(0);

  
  TFile * in = new TFile("fit_area_sim_histograms.root");
  char canname[50];
  char histname[50];

  TH1F * hists[9];
  TH1F * analytic_hists[9];
  TCanvas * cans[9];

  double time;
  double integral;

  for(int i=0; i<9; i++) {

     sprintf(histname,"h_avg_shape_%d",i+1);
     hists[i] = (TH1F *)in->Get(histname);
//     hists[i]->Sumw2();

     sprintf(canname,"can_%d",i+1);
     cans[i] = new TCanvas(canname,canname);
  
     sprintf(histname,"an_%d",i+1);
     analytic_hists[i] = new TH1F(histname,histname,700,-30.,40.);
     analytic_hists[i]->SetLineColor(4);
     analytic_hists[i]->GetXaxis()->SetRangeUser(-10.,30.);

     
     for(int j=0; j<1000000; j++) {
        if( r.Uniform() < 0.4 ){
           time = r.Exp(0.31) + r.Exp(1.08) + r.Gaus(0.,0.75) + r.Gaus(0.,1.693/TMath::Sqrt((double) i+15.)) + 
                  r.Uniform() -0.345*TMath::Power((double) i + 15.,-0.699) - 0.729; 
           analytic_hists[i]->Fill(time);
        } else {
           time = r.Exp(2.7) + r.Exp(1.08) + r.Gaus(0.,0.75) + r.Gaus(0.,1.693/TMath::Sqrt((double) i + 15.)) + 
                 r.Uniform() -0.345*TMath::Power((double) i + 15.,-0.699) - 0.729;
           analytic_hists[i]->Fill(time);
        }
     }

//     analytic_hists[i]->Sumw2();
     analytic_hists[i]->Scale(1./1.e6);
//     analytic_hists[i]->Divide(hists[i]);
     analytic_hists[i]->Draw();

     integral = hists[i]->Integral(1,700);
     hists[i]->Scale(1./integral);

     hists[i]->Draw("SAME");

         

   }
      
}




















}
