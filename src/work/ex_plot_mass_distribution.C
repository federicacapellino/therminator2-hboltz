// ex_plot_mass_distribution.C — overlay correlated pair invariant mass for HRG and blastwave
//
// Usage (from repository root):
//   root -b -q src/work/ex_plot_mass_distribution.C
// or with custom file names:
//   root -b -q 'src/work/ex_plot_mass_distribution.C("mass_hrg.root","mass_blastwave.root")'

void ex_plot_mass_distribution(const char* fileHRG  = "mass_hrg.root",
                               const char* fileBW   = "mass_blastwave.root") {
   TFile* fHRG = TFile::Open(fileHRG);
   if (!fHRG || fHRG->IsZombie()) { printf("Cannot open %s\n", fileHRG); return; }

   TFile* fBW = TFile::Open(fileBW);
   if (!fBW  || fBW->IsZombie())  { printf("Cannot open %s\n", fileBW);  return; }

   TH1D* hHRG = (TH1D*) fHRG->Get("hMass");
   TH1D* hBW  = (TH1D*) fBW->Get("hMass");
   if (!hHRG) { printf("Cannot find hMass in %s\n", fileHRG); return; }
   if (!hBW)  { printf("Cannot find hMass in %s\n", fileBW);  return; }

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

   TCanvas* c1 = new TCanvas("c1", "Correlated pair invariant mass", 600, 500);
   c1->SetLogy();

   hHRG->SetLineColor(kBlue + 1);
   hHRG->SetLineWidth(2);
   hBW->SetLineColor(kRed + 1);
   hBW->SetLineWidth(2);

   hHRG->GetXaxis()->SetTitle("M (GeV)");
   hHRG->GetXaxis()->SetTitleSize(0.05);
   hHRG->GetXaxis()->SetTitleOffset(0.9);
   hHRG->GetXaxis()->SetLabelSize(0.045);
   hHRG->GetYaxis()->SetTitle("(1/N_{1}) dN_{pairs}/dM (GeV^{-1})");
   hHRG->GetYaxis()->SetTitleSize(0.05);
   hHRG->GetYaxis()->SetTitleOffset(1.4);
   hHRG->GetYaxis()->SetLabelSize(0.045);

   hHRG->Draw("hist e");
   hBW->Draw("hist e same");

   TLegend* leg = new TLegend(0.55, 0.70, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hHRG, "HRG",        "l");
   leg->AddEntry(hBW,  "Blast wave", "l");
   leg->Draw();

   gPad->Print("mass_distribution.pdf");
}

// Plot hMass sliced by pT (blastwave) or |p| (HRG) of particle 1.
// Each slice is normalized by the corresponding hPid1Counts bin.
//
// Usage:
//   root -b -q 'ex_plot_mass_distribution.C("mass_hrg.root","mass_blastwave.root")'
//   root -l 'ex_plot_mass_distribution.C'; ex_plot_mass_slices("mass_blastwave.root")

void ex_plot_mass_slices(const char* filename = "mass_blastwave.root") {
   TFile* f = TFile::Open(filename);
   if (!f || f->IsZombie()) { printf("Cannot open %s\n", filename); return; }

   TH2D* h2 = (TH2D*) f->Get("hMass2D");
   TH1D* hN = (TH1D*) f->Get("hPid1Counts");
   if (!h2) { printf("Cannot find hMass2D in %s\n", filename); return; }
   if (!hN) { printf("Cannot find hPid1Counts in %s\n", filename); return; }

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

   Int_t colors[] = {kBlue+1, kCyan+2, kGreen+2, kYellow+2, kOrange+1, kRed+1};

   TCanvas* c = new TCanvas("cSlices", "Mass slices by p_{T,1}", 700, 550);
   c->SetLogy();

   TH1D* hFirst = nullptr;
   TLegend* leg = new TLegend(0.55, 0.52, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetHeader(h2->GetXaxis()->GetTitle(), "C");

   int nBins = h2->GetNbinsX();
   for (int i = 1; i <= nBins; i++) {
      TH1D* hSlice = h2->ProjectionY(Form("hSlice_%d", i), i, i);
      double n = hN->GetBinContent(i);
      if (n > 0) hSlice->Scale(1.0 / n);
      hSlice->SetLineColor(colors[i - 1]);
      hSlice->SetLineWidth(2);
      if (!hFirst) {
         hFirst = hSlice;
         hFirst->GetYaxis()->SetTitle("(1/N_{1}) dN_{pairs}/dM (GeV^{-1})");
         hFirst->GetYaxis()->SetTitleSize(0.05);
         hFirst->GetYaxis()->SetTitleOffset(1.4);
         hFirst->GetYaxis()->SetLabelSize(0.045);
         hFirst->GetXaxis()->SetTitle("M (GeV)");
         hFirst->GetXaxis()->SetTitleSize(0.05);
         hFirst->GetXaxis()->SetTitleOffset(0.9);
         hFirst->GetXaxis()->SetLabelSize(0.045);
         hFirst->Draw("hist e");
      } else {
         hSlice->Draw("hist e same");
      }
      leg->AddEntry(hSlice,
          Form("%.2g #minus %.2g GeV",
               h2->GetXaxis()->GetBinLowEdge(i),
               h2->GetXaxis()->GetBinUpEdge(i)), "l");
   }
   leg->Draw();
   gPad->Print("mass_slices.pdf");
}
