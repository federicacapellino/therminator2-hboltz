

Double_t dist(double p) {
   double m = 0.1395699 ;
   double hbarc = 0.1973269631;
   double x = sqrt(p*p + m*m)/0.155 ;
   double V = 30.*30.*30./pow(hbarc,3) ;
   double nevents = 500. ;
   double dp = 2./40.; // bin width

   return 2. * p * p * dp * nevents * V * exp(-x)/(1. - exp(-x)) /(2.*pow(M_PI,2)); 
}

void plot_pitypes() {
   TH1D *h_pi = (TH1D *) gDirectory->Get("pi_dNdp") ;
   TH1D *h_pi_decayed = (TH1D *) gDirectory->Get("pi_decayed_dNdp") ;
   TH1D *h_pi_decayed_short = (TH1D *) gDirectory->Get("pi_decayed_short_25_dNdp") ;
   TH1D *h_pi_all = (TH1D *) gDirectory->Get("pi_all_dNdp") ;

   //  Note that is .155 
   gStyle -> SetPadLeftMargin(.15) ;
   TCanvas *c1 = new TCanvas("c1","A Working Canvas",600,500) ;

   // gROOT->SetStyle("Plain") ;
   gStyle -> SetOptStat(kFALSE) ;
   gStyle -> SetOptTitle(kFALSE) ;

   //  Note that is .155 
   // gStyle -> SetPadLeftMargin(.17) ;

   //gStyle -> SetLabelOffset(0.006, "X") ;
   //gStyle -> SetTitleOffset(1.2,"X") ;

   // No need to conserve Horizontal space usually
   gStyle -> SetTitleOffset(1.0, "Y") ;
   gStyle -> SetPaperSize(TStyle::kUSLetter) ;

   TH1D *hframe = new TH1D("h1","h1", 200, 0., 1.5) ;
   hframe->SetMaximum(350.*pow(10,3)) ;
   hframe->Draw("AXIS") ;
   hframe->GetXaxis()->SetTitle("p (GeV)") ;
   hframe->GetYaxis()->SetTitle("dN/dp (arb units)") ;
   hframe->GetYaxis()->SetLabelSize(0.045) ;
   hframe->GetXaxis()->SetLabelSize(0.045) ;

   hframe->GetXaxis()->SetTitle("p (GeV)");
   hframe->GetXaxis()->SetLabelFont(42);
   hframe->GetXaxis()->SetTitleSize(0.05);
   hframe->GetXaxis()->SetTitleOffset(0.9);
   hframe->GetXaxis()->SetTitleFont(42);
   hframe->GetYaxis()->SetTitle("dN/dp (arb units)");
   hframe->GetYaxis()->SetLabelFont(42);
   hframe->GetYaxis()->SetTitleSize(0.05);
   hframe->GetYaxis()->SetTitleOffset(1.2);

   Int_t cdarkred = TColor::GetFreeColorIndex() ;
   TColor *mydarkred = new TColor(cdarkred, 200./255., 37./255., 6./255.) ;
   Int_t cdarkblue = TColor::GetFreeColorIndex() ;
   TColor *mydarkblue = new TColor(cdarkblue, 1./255., 25./255., 147./255.) ;
   Int_t cdarkgreen = TColor::GetFreeColorIndex() ;
   TColor *mydarkgreen = new TColor(cdarkgreen, 11./255., 93./255., 24./255.) ;
   Int_t cdarkgrey = TColor::GetFreeColorIndex() ;
   TColor *mydarkgrey = new TColor(cdarkgrey, 140./255., 140./255., 140./255.) ;

   gStyle->SetLineStyleString(5,"[40 10]") ;
   gStyle->SetLineStyleString(6,"[20 20]") ;
   gStyle->SetLineStyleString(7,"[40 15 10 15]") ;
   gStyle->SetLineStyleString(8,"[10 10]") ;
   gStyle->SetLineStyleString(9,"[15 15]") ;
   
   h_pi_all->SetLineColor(cdarkred) ;
   h_pi_all->SetLineWidth(3) ;
   h_pi_all->GetXaxis()->SetTitle("p (GeV)") ;
   h_pi_all->GetYaxis()->SetTitle("dN/dp (arb units)") ;
   h_pi_all->Draw("same Chist") ;

   h_pi_decayed->SetLineColor(cdarkred) ;
   h_pi_decayed->SetLineWidth(3) ;
   h_pi_decayed->SetLineStyle(5) ;
   h_pi_decayed->Draw("Chist same") ;

   h_pi_decayed_short->SetLineColor(cdarkgrey) ;
   h_pi_decayed_short->SetLineWidth(3) ;
   h_pi_decayed_short->SetLineStyle(kSolid) ;
   h_pi_decayed_short->Draw("Chist same") ;

   h_pi->SetLineColor(cdarkred) ;
   h_pi->SetLineWidth(3) ;
   h_pi->SetLineStyle(9) ;
   h_pi->Draw("Chist,same") ;

//TF1 f1("f1", "dist(x)",0., 2.) ;
//   f1.Draw("same") ;

   TLegend *leg = new TLegend(0.45,0.55,0.88,0.88,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.045);

   leg->AddEntry(h_pi_all, "all pions") ;
   leg->AddEntry(h_pi_decayed, "decay pions") ;
   leg->AddEntry(h_pi, "thermal pions") ;
   leg->AddEntry(h_pi_decayed_short, "decay pions (#tau<2.5 fm)") ;
   leg->SetBorderSize(0) ;
   leg->Draw() ;

   gPad->Print("c1.pdf") ;
}

