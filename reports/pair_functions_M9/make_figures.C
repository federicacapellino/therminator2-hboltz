// Figures for the M9 report (spherical-harmonic anisotropy).
// Run: root -l -b -q make_figures.C
const char* T2 = "/home/alice/therminator2-hboltz/build/t2/";

// a20 and a40 versus kappa at the test masses of one pair: L2 (filled), L1 a20 (open).
void moments(const char* weak, const char* tag, std::vector<const char*> labels) {
  TFile* f = TFile::Open(Form("%sm9/harm_weak%s_ana.root", T2, weak));
  const int n = labels.size();
  TCanvas c("c", "", 500 * n, 520);
  c.Divide(n, 1);
  std::vector<TObject*> keep;
  for (int i = 0; i < n; i++) {
    c.cd(i + 1);
    gPad->SetLeftMargin(0.17);
    gPad->SetGridy();
    TH1D* a20 = (TH1D*)f->Get(Form("a20_L2_%s_%s", tag, labels[i]));
    TH1D* a40 = (TH1D*)f->Get(Form("a40_L2_%s_%s", tag, labels[i]));
    TH1D* b20 = (TH1D*)f->Get(Form("a20_L1_%s_%s", tag, labels[i]));
    double lo = 1e9, hi = -1e9;
    for (TH1D* h : {a20, a40, b20})
      for (int k = 1; k <= h->GetNbinsX(); k++) {
        if (h->GetBinError(k) <= 0) continue;
        lo = std::min(lo, h->GetBinContent(k) - h->GetBinError(k));
        hi = std::max(hi, h->GetBinContent(k) + h->GetBinError(k));
      }
    double pad = 0.1 * (hi - lo);
    a20->SetTitle(Form("%s, %s;#kappa (GeV);a_{l0}", tag, labels[i]));
    a20->SetMinimum(lo - pad);
    a20->SetMaximum(hi + 3 * pad);
    a20->GetXaxis()->SetRangeUser(0, 2.0);
    a20->SetMarkerStyle(20); a20->SetMarkerColor(kRed + 1); a20->SetLineColor(kRed + 1);
    a40->SetMarkerStyle(21); a40->SetMarkerColor(kBlue + 1); a40->SetLineColor(kBlue + 1);
    b20->SetMarkerStyle(24); b20->SetMarkerColor(kGray + 2); b20->SetLineColor(kGray + 2);
    a20->Draw("E1"); b20->Draw("E1 SAME"); a40->Draw("E1 SAME"); a20->Draw("E1 SAME");
    TLegend* leg = new TLegend(0.45, 0.74, 0.88, 0.89);
    leg->SetBorderSize(0);
    leg->AddEntry(a20, "a_{20}, L2", "lep");
    leg->AddEntry(a40, "a_{40}, L2", "lep");
    leg->AddEntry(b20, "a_{20}, L1", "lep");
    leg->Draw();
    keep.push_back(leg);
  }
  c.SaveAs(Form("fig/moments_%s_weak%s.png", tag, weak));
}

// Pull distributions of T9.2, T9.3, T9.4 (bins with >= 500 entries) with the unit Gaussian.
void pulls() {
  TCanvas c("c", "", 1500, 900);
  c.Divide(3, 2);
  const char* ws[]    = {"on", "off"};
  const char* names[] = {"pulls_T92", "pulls_T93", "pulls_T94"};
  const char* tits[]  = {"T9.2: m #neq 0, L0+L1+L2", "T9.3: odd l, m = 0, unlike-sign",
                         "T9.4: l = 2, 4, m = 0, L0+L1"};
  for (int w = 0; w < 2; w++) {
    TFile* f = TFile::Open(Form("%sm9/harm_weak%s_ana500.root", T2, ws[w])); // >= 500 entries
    for (int t = 0; t < 3; t++) {
      c.cd(3 * w + t + 1);
      gPad->SetLogy();
      TH1D* h = (TH1D*)f->Get(names[t]);
      h->SetTitle(Form("%s (weak %s);a/#sigma;bins", tits[t], ws[w]));
      h->SetMarkerStyle(20); h->SetMarkerSize(0.6);
      h->Draw("E1");
      TF1* g = new TF1(Form("g%d%d", w, t), "[0]*exp(-0.5*x*x)", -5, 5);
      g->SetParameter(0, h->GetEntries() * h->GetBinWidth(1) / std::sqrt(2 * M_PI));
      g->SetLineColor(kRed + 1);
      g->Draw("SAME");
    }
  }
  c.SaveAs("fig/pulls.png");
}

void make_figures() {
  gStyle->SetOptStat(0);
  gSystem->mkdir("fig");
  for (const char* w : {"on", "off"}) {
    moments(w, "pippim", {"pole rho770zer", "cont 0.30", "cont 0.45", "cont 0.60"});
    moments(w, "pippip", {"cont 0.30", "cont 0.45", "cont 0.60"});
    moments(w, "KpKm", {"pole ph1020zer", "cont 1.10", "cont 1.30"});
  }
  pulls();
}
