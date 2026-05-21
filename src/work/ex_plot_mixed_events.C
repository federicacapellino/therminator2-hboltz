// ex_plot_mixed_events.C — plot histograms from ex_mixed_events.exe output
//
// Usage:
//   root -b -q 'ex_plot_mixed_events.C("mixed_events.root")'
// or interactively:
//   root mixed_events.root ex_plot_mixed_events.C

void ex_plot_mixed_events(const char* filename = "mixed_events.root") {
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    printf("Cannot open %s\n", filename);
    return;
  }

  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);

  // ── Canvas 1: signal+background, background, and connected signal ──────────
  {
    TH1D* hSAB = (TH1D*) f->Get("hSignalAndBackground");
    TH1D* hBkg = (TH1D*) f->Get("hBackground");
    TH1D* hSig = (TH1D*) f->Get("hSignal");

    TCanvas* c1 = new TCanvas("c1", "Invariant mass spectra", 600, 500);
    c1->SetLogy();

    hSAB->SetLineColor(kBlue + 1);
    hSAB->SetLineWidth(2);
    hSAB->GetXaxis()->SetTitle("M (GeV)");
    hSAB->GetXaxis()->SetTitleSize(0.05);
    hSAB->GetXaxis()->SetTitleOffset(0.9);
    hSAB->GetXaxis()->SetLabelSize(0.045);
    hSAB->GetYaxis()->SetTitle("(1/N_{1}) dN_{pairs}/dM (GeV^{-1})");
    hSAB->GetYaxis()->SetTitleSize(0.05);
    hSAB->GetYaxis()->SetTitleOffset(1.4);
    hSAB->GetYaxis()->SetLabelSize(0.045);
    hSAB->Draw("hist e");

    hBkg->SetLineColor(kRed + 1);
    hBkg->SetLineWidth(2);
    hBkg->SetLineStyle(2);
    hBkg->Draw("hist e same");

    hSig->SetLineColor(kGreen + 2);
    hSig->SetLineWidth(2);
    hSig->Draw("hist e same");

    TLegend* leg = new TLegend(0.45, 0.65, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hSAB, "Signal + background", "l");
    leg->AddEntry(hBkg, "Background",          "l");
    leg->AddEntry(hSig, "Signal (difference)", "l");
    leg->Draw();

    gPad->Print("mixed_events_mass.pdf");
  }

  // ── Canvas 2: rapidity distributions ──────────────────────────────────────
  {
    TH1D* hY1 = (TH1D*) f->Get("hRapidityPid1");
    TH1D* hY2 = (TH1D*) f->Get("hRapidityPid2");

    TCanvas* c2 = new TCanvas("c2", "Rapidity distributions", 600, 500);

    hY1->SetLineColor(kBlue + 1);
    hY1->SetLineWidth(2);
    hY1->GetXaxis()->SetTitle("y");
    hY1->GetXaxis()->SetTitleSize(0.05);
    hY1->GetXaxis()->SetTitleOffset(0.9);
    hY1->GetXaxis()->SetLabelSize(0.045);
    hY1->GetYaxis()->SetTitle("dN/dy");
    hY1->GetYaxis()->SetTitleSize(0.05);
    hY1->GetYaxis()->SetTitleOffset(1.4);
    hY1->GetYaxis()->SetLabelSize(0.045);
    hY1->Draw("hist");

    hY2->SetLineColor(kRed + 1);
    hY2->SetLineWidth(2);
    hY2->SetLineStyle(2);
    hY2->Draw("hist same");

    TLegend* leg = new TLegend(0.45, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hY1, "pid1", "l");
    leg->AddEntry(hY2, "pid2", "l");
    leg->Draw();

    gPad->Print("mixed_events_rapidity.pdf");
  }

  // ── Canvas 3: same-rooteid pairs, close vs far in rapidity ────────────────
  {
    TH1D* hClose = (TH1D*) f->Get("hSameRootClose");
    TH1D* hFar   = (TH1D*) f->Get("hSameRootFar");

    TCanvas* c3 = new TCanvas("c3", "Same-rooteid pairs", 600, 500);
    c3->SetLogy();

    hClose->SetLineColor(kBlue + 1);
    hClose->SetLineWidth(2);
    hClose->GetXaxis()->SetTitle("M (GeV)");
    hClose->GetXaxis()->SetTitleSize(0.05);
    hClose->GetXaxis()->SetTitleOffset(0.9);
    hClose->GetXaxis()->SetLabelSize(0.045);
    hClose->GetYaxis()->SetTitle("dN_{pairs}/dM (GeV^{-1})");
    hClose->GetYaxis()->SetTitleSize(0.05);
    hClose->GetYaxis()->SetTitleOffset(1.4);
    hClose->GetYaxis()->SetLabelSize(0.045);
    hClose->Draw("hist e");

    hFar->SetLineColor(kRed + 1);
    hFar->SetLineWidth(2);
    hFar->SetLineStyle(2);
    hFar->Draw("hist e same");

    TLegend* leg = new TLegend(0.45, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(hClose, "Same rooteid, |#Deltay| < DY", "l");
    leg->AddEntry(hFar,   "Same rooteid, |#Deltay| #geq DY", "l");
    leg->Draw();

    gPad->Print("mixed_events_sameroot.pdf");
  }
}
