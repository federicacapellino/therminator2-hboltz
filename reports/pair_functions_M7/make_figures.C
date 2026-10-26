// Figures for the M7 report (L1 table).
// Run: root -l -b -q make_figures.C
const char* T2 = "/home/alice/therminator2-hboltz/build/t2/";

// rho(M): MC at L1 vs the L1 table (pulls below), L0 table for reference.
void rhoPlot(const char* weak, const char* tag, double xlo, double xhi) {
  TFile f1(Form("%sm7/compare_L1_weak%s.root", T2, weak));
  TFile f0(Form("%sm3_weak%s.root", T2, weak));
  TH1D* mc = (TH1D*)f1.Get(Form("rho_L1_%s", tag))->Clone("mc");
  TH1D* p1 = (TH1D*)f1.Get(Form("rho_pred_%s", tag))->Clone("p1");
  TH1D* p0 = (TH1D*)f0.Get(Form("rho_pred_%s", tag))->Clone("p0");
  TCanvas c("c", "", 1400, 1200);
  TPad a("a", "", 0, 0.32, 1, 1); a.SetBottomMargin(0.02); a.SetLeftMargin(0.14); a.SetLogy(); a.Draw();
  TPad b("b", "", 0, 0, 1, 0.32); b.SetTopMargin(0.02); b.SetBottomMargin(0.32); b.SetLeftMargin(0.14); b.Draw();
  a.cd();
  TH1D* pull = (TH1D*)mc->Clone("pull");
  for (int i = 1; i <= pull->GetNbinsX(); i++) {
    double e = mc->GetBinError(i);
    pull->SetBinContent(i, (e > 0 && p1->GetBinContent(i) > 0) ? (mc->GetBinContent(i) - p1->GetBinContent(i)) / e : -99);
    pull->SetBinError(i, 0);
  }
  for (TH1D* h : {mc, p1, p0}) { h->Scale(1., "width"); h->GetXaxis()->SetRangeUser(xlo, xhi); }
  mc->SetTitle(""); mc->GetYaxis()->SetTitle("d#rho/dM  (GeV^{2})"); mc->GetXaxis()->SetLabelSize(0);
  mc->GetYaxis()->SetTitleSize(0.055); mc->GetYaxis()->SetLabelSize(0.045);
  mc->SetMarkerStyle(20); mc->SetMarkerSize(0.4); mc->SetLineColor(kBlack);
  p1->SetLineColor(kRed + 1); p1->SetLineWidth(3);
  p0->SetLineColor(kBlue + 1); p0->SetLineWidth(2); p0->SetLineStyle(2);
  mc->Draw("E1"); p0->Draw("HIST SAME"); p1->Draw("HIST SAME"); mc->Draw("E1 SAME");
  TLegend leg(0.58, 0.72, 0.88, 0.88); leg.SetBorderSize(0);
  leg.AddEntry(mc, "MC, level L1", "lep"); leg.AddEntry(p1, "L1 table", "l");
  leg.AddEntry(p0, "L0 table", "l"); leg.Draw();
  b.cd();
  pull->SetTitle(""); pull->GetXaxis()->SetRangeUser(xlo, xhi); pull->SetMinimum(-4.5); pull->SetMaximum(4.5);
  pull->SetMarkerStyle(20); pull->SetMarkerSize(0.35);
  pull->GetYaxis()->SetTitle("pull (L1)"); pull->GetXaxis()->SetTitle("M (GeV)");
  pull->GetYaxis()->SetNdivisions(505); pull->GetYaxis()->SetTitleSize(0.11); pull->GetYaxis()->SetTitleOffset(0.55);
  pull->GetYaxis()->SetLabelSize(0.09); pull->GetXaxis()->SetLabelSize(0.10); pull->GetXaxis()->SetTitleSize(0.11);
  pull->Draw("P");
  TLine l(xlo, 0, xhi, 0); l.SetLineColor(kRed + 1); l.SetLineStyle(2); l.Draw();
  c.SaveAs(Form("fig/rho_%s_weak%s.png", tag, weak));
}

// dP/dkappa / Boltzmann at fixed M: MC L1, L1 table, L0 table.
void dpdk(const char* weak) {
  TFile f1(Form("%sm7/compare_L1_weak%s.root", T2, weak));
  TFile f0(Form("%sm3_weak%s.root", T2, weak));
  const char* lab[] = {"pole rho770zer", "cont 0.45", "pole Ka0492zer", "pole f00980zer"};
  const char* tit[] = {"#pi^{+}#pi^{-}, #rho^{0} pole bin", "#pi^{+}#pi^{-}, M = 0.45 GeV",
                       "#pi^{+}#pi^{-}, K^{0} pole bin", "#pi^{+}#pi^{-}, f_{0}(980) pole bin"};
  TCanvas c("c", "", 1600, 1200);
  c.Divide(2, 2);
  std::vector<TObject*> keep;
  for (int i = 0; i < 4; i++) {
    c.cd(i + 1); gPad->SetLeftMargin(0.14);
    TH1D* mc = (TH1D*)f1.Get(Form("dPdkRatio_L1_pippim_%s", lab[i]));
    TH1D* t1 = (TH1D*)f1.Get(Form("dPdkRatio_pred_pippim_%s", lab[i]));
    TH1D* t0 = (TH1D*)f0.Get(Form("dPdkRatio_pred_pippim_%s", lab[i]));
    if (!mc || !t1 || !t0) continue;
    mc = (TH1D*)mc->Clone(Form("mc%d", i)); t1 = (TH1D*)t1->Clone(Form("t1%d", i));
    t0 = (TH1D*)t0->Clone(Form("t0%d", i));
    mc->SetTitle(Form("%s;#kappa (GeV);(dP/d#kappa) / Boltzmann", tit[i]));
    mc->GetXaxis()->SetRangeUser(0, 2.0);
    mc->SetMinimum(i == 1 ? 0.3 : 0.6); mc->SetMaximum(i == 1 ? 2.2 : 1.4);
    mc->SetMarkerStyle(20); mc->SetMarkerSize(0.6); mc->SetLineColor(kBlack);
    t1->SetLineColor(kRed + 1); t1->SetLineWidth(2);
    t0->SetLineColor(kBlue + 1); t0->SetLineWidth(2); t0->SetLineStyle(2);
    mc->Draw("E1"); t0->Draw("HIST SAME"); t1->Draw("HIST SAME"); mc->Draw("E1 SAME");
    TLine* l = new TLine(0, 1, 2.0, 1); l->SetLineStyle(3); l->Draw();
    if (i == 0) {
      TLegend* leg = new TLegend(0.55, 0.7, 0.88, 0.88); leg->SetBorderSize(0);
      leg->AddEntry(mc, "MC, level L1", "lep"); leg->AddEntry(t1, "L1 table", "l");
      leg->AddEntry(t0, "L0 table", "l"); leg->Draw(); keep.push_back(leg);
    }
    keep.insert(keep.end(), {mc, t1, t0, l});
  }
  c.SaveAs(Form("fig/dpdk_weak%s.png", weak));
}

// L1 / L0 table ratio of rho(M), pi+pi- and K+K-, both weak settings.
void tableRatio() {
  TCanvas c("c", "", 1600, 700);
  c.Divide(2, 1);
  const char* tags[] = {"pippim", "KpKm"};
  std::vector<TObject*> keep;
  for (int i = 0; i < 2; i++) {
    c.cd(i + 1); gPad->SetLeftMargin(0.13);
    TLegend* leg = new TLegend(0.6, 0.75, 0.88, 0.88); leg->SetBorderSize(0);
    int col[] = {kBlack, kRed + 1};
    const char* ws[] = {"on", "off"};
    for (int w = 0; w < 2; w++) {
      TFile* f1 = TFile::Open(Form("%sm7/compare_L1_weak%s.root", T2, ws[w]));
      TFile* f0 = TFile::Open(Form("%sm3_weak%s.root", T2, ws[w]));
      TH1D* r = (TH1D*)f1->Get(Form("rho_pred_%s", tags[i]))->Clone(Form("r%d%d", i, w));
      r->SetDirectory(nullptr);
      r->Divide((TH1D*)f0->Get(Form("rho_pred_%s", tags[i])));
      r->SetTitle(Form("%s;M (GeV);#rho^{L1} / #rho^{L0} (tables)", i ? "K^{+}K^{-}" : "#pi^{+}#pi^{-}"));
      r->SetLineColor(col[w]); r->SetLineWidth(2); r->SetMinimum(0.9); r->SetMaximum(2.4);
      r->Draw(w ? "HIST SAME" : "HIST");
      leg->AddEntry(r, Form("weak %s", ws[w]), "l");
      keep.push_back(r);
    }
    leg->Draw();
  }
  c.SaveAs("fig/table_ratio.png");
}

void make_figures() {
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  gSystem->mkdir("fig", true);
  for (const char* w : {"on", "off"}) {
    rhoPlot(w, "pippim", 0.27, 2.0);
    rhoPlot(w, "KpKm", 0.98, 2.0);
    dpdk(w);
  }
  tableRatio();
}
