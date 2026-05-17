// ex_plot_o4_pion_spectra.C -- plot pion momentum spectra from ex_o4_pion_spectra.exe
//
// Usage:
//   root -b -q 'ex_plot_o4_pion_spectra.C("o4_pion_spectra.root")'
// or interactively:
//   root o4_pion_spectra.root ex_plot_o4_pion_spectra.C

// Thermal pion spectrum dN/dp for a Boltzmann gas in volume V at temperature T.
// Used as a cross-check against the HRG primordial pions.
Double_t dist(double p) {
   double m      = 0.1395699;
   double hbarc  = 0.1973269631;
   double T      = 0.155;
   double x      = sqrt(p*p + m*m) / T;
   double V      = 30.*30.*30. / pow(hbarc, 3);
   double nevents = 500.;
   double dp     = 2./40.;
   return 2. * p*p * dp * nevents * V * exp(-x) / (1. - exp(-x)) / (2.*pow(M_PI, 2));
}

void ex_plot_o4_pion_spectra(const char* filename = "o4_pion_spectra.root") {
   TFile* f = TFile::Open(filename);
   if (!f || f->IsZombie()) {
      printf("Cannot open %s\n", filename);
      return;
   }

   TH1D* h_pi            = (TH1D*) f->Get("pi_dNdp");
   TH1D* h_pi_decayed    = (TH1D*) f->Get("pi_decayed_dNdp");
   TH1D* h_pi_decayed_short = (TH1D*) f->Get("pi_decayed_short_25_dNdp");
   TH1D* h_pi_all        = (TH1D*) f->Get("pi_all_dNdp");

   gStyle->SetPadLeftMargin(.15);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);
   gStyle->SetTitleOffset(1.0, "Y");
   gStyle->SetPaperSize(TStyle::kUSLetter);

   gStyle->SetLineStyleString(5, "[40 10]");
   gStyle->SetLineStyleString(9, "[15 15]");

   Int_t cdarkred  = TColor::GetFreeColorIndex();
   new TColor(cdarkred,  200./255.,  37./255.,   6./255.);
   Int_t cdarkgrey = TColor::GetFreeColorIndex();
   new TColor(cdarkgrey, 140./255., 140./255., 140./255.);

   TCanvas* c1 = new TCanvas("c1", "O(4) pion spectra", 600, 500);

   TH1D* hframe = new TH1D("hframe", "", 200, 0., 1.5);
   hframe->SetMaximum(350.e3);
   hframe->GetXaxis()->SetTitle("p (GeV)");
   hframe->GetXaxis()->SetLabelFont(42);
   hframe->GetXaxis()->SetTitleSize(0.05);
   hframe->GetXaxis()->SetTitleOffset(0.9);
   hframe->GetXaxis()->SetLabelSize(0.045);
   hframe->GetYaxis()->SetTitle("dN/dp (arb units)");
   hframe->GetYaxis()->SetLabelFont(42);
   hframe->GetYaxis()->SetTitleSize(0.05);
   hframe->GetYaxis()->SetTitleOffset(1.2);
   hframe->GetYaxis()->SetLabelSize(0.045);
   hframe->Draw("AXIS");

   h_pi_all->SetLineColor(cdarkred);
   h_pi_all->SetLineWidth(3);
   h_pi_all->Draw("Chist same");

   h_pi_decayed->SetLineColor(cdarkred);
   h_pi_decayed->SetLineWidth(3);
   h_pi_decayed->SetLineStyle(5);
   h_pi_decayed->Draw("Chist same");

   h_pi_decayed_short->SetLineColor(cdarkgrey);
   h_pi_decayed_short->SetLineWidth(3);
   h_pi_decayed_short->Draw("Chist same");

   h_pi->SetLineColor(cdarkred);
   h_pi->SetLineWidth(3);
   h_pi->SetLineStyle(9);
   h_pi->Draw("Chist same");

   // TF1 f1("f1", "dist(x)", 0., 2.); f1.Draw("same");

   TLegend* leg = new TLegend(0.45, 0.55, 0.88, 0.88, NULL, "brNDC");
   leg->SetBorderSize(0);
   leg->SetTextSize(0.045);
   leg->AddEntry(h_pi_all,           "all pions");
   leg->AddEntry(h_pi_decayed,       "decay pions");
   leg->AddEntry(h_pi,               "thermal pions");
   leg->AddEntry(h_pi_decayed_short, "decay pions (#tau < 2.5 fm)");
   leg->Draw();

   gPad->Print("o4_pion_spectra.pdf");
}
