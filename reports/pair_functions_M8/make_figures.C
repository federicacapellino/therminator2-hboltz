// Figures for the M8 report (chain kernels, L2 table).
// Run: root -l -b -q make_figures.C
const char* T2 = "/home/alice/therminator2-hboltz/build/t2/";

// rho(M): MC at L2 vs the L2 table (pulls below), L1 table for reference.
void rhoPlot(const char* weak, const char* tag, double xlo, double xhi, bool likeSign) {
  TFile f2(Form("%sm8/compare_L2_weak%s.root", T2, weak));
  TFile f1(Form("%sm7/compare_L1_weak%s.root", T2, weak));
  TH1D* mc = (TH1D*)f2.Get(Form("rho_L2_%s", tag))->Clone("mc");
  TH1D* p2 = (TH1D*)f2.Get(Form("rho_pred_%s", tag))->Clone("p2");
  TH1D* p1 = (TH1D*)f1.Get(Form("rho_pred_%s", tag))->Clone("p1");
  TCanvas c("c", "", 1400, 1200);
  TPad a("a", "", 0, 0.32, 1, 1); a.SetBottomMargin(0.02); a.SetLeftMargin(0.14); a.SetLogy(); a.Draw();
  TPad b("b", "", 0, 0, 1, 0.32); b.SetTopMargin(0.02); b.SetBottomMargin(0.32); b.SetLeftMargin(0.14); b.Draw();
  TH1D* pull = (TH1D*)mc->Clone("pull");
  pull->SetTitle("");
  for (int i = 1; i <= pull->GetNbinsX(); i++) {
    double mu = p2->GetBinContent(i);
    double e  = mc->GetBinError(i) * (likeSign ? std::sqrt(2.) : 1.);
    pull->SetBinContent(i, (e > 0 && mu > 0) ? (mc->GetBinContent(i) - mu) / e : -99);
    pull->SetBinError(i, 0);
  }
  a.cd();
  for (TH1D* h : {mc, p2, p1}) { h->Scale(1., "width"); h->GetXaxis()->SetRangeUser(xlo, xhi); }
  mc->SetTitle(""); mc->GetYaxis()->SetTitle("d#rho/dM  (GeV^{2})"); mc->GetXaxis()->SetLabelSize(0);
  mc->GetYaxis()->SetTitleSize(0.055); mc->GetYaxis()->SetLabelSize(0.045);
  mc->SetMarkerStyle(20); mc->SetMarkerSize(0.4); mc->SetLineColor(kBlack);
  p2->SetLineColor(kRed + 1); p2->SetLineWidth(3);
  p1->SetLineColor(kBlue + 1); p1->SetLineWidth(2); p1->SetLineStyle(2);
  mc->Draw("E1"); p1->Draw("HIST SAME"); p2->Draw("HIST SAME"); mc->Draw("E1 SAME");
  TLegend leg(0.58, 0.72, 0.88, 0.88); leg.SetBorderSize(0);
  leg.AddEntry(mc, "MC, level L2", "lep"); leg.AddEntry(p2, "L2 table", "l");
  leg.AddEntry(p1, "L1 table", "l"); leg.Draw();
  b.cd();
  pull->GetXaxis()->SetRangeUser(xlo, xhi); pull->SetMinimum(-4.5); pull->SetMaximum(4.5);
  pull->SetMarkerStyle(20); pull->SetMarkerSize(0.35);
  pull->GetYaxis()->SetTitle("pull (L2)"); pull->GetXaxis()->SetTitle("M (GeV)");
  pull->GetYaxis()->SetNdivisions(505); pull->GetYaxis()->SetTitleSize(0.11); pull->GetYaxis()->SetTitleOffset(0.55);
  pull->GetYaxis()->SetLabelSize(0.09); pull->GetXaxis()->SetLabelSize(0.10); pull->GetXaxis()->SetTitleSize(0.11);
  pull->Draw("P");
  TLine l(xlo, 0, xhi, 0); l.SetLineColor(kRed + 1); l.SetLineStyle(2); l.Draw();
  c.SaveAs(Form("fig/rho_%s_weak%s.png", tag, weak));
}

// Chain fraction of the L2 table versus M.
void chainFraction() {
  TCanvas c("c", "", 1600, 700);
  c.Divide(2, 1);
  const char* ws[] = {"on", "off"};
  const char* tags[] = {"pippim", "pippip", "KpKm"};
  const int col[] = {kBlack, kRed + 1, kBlue + 1};
  std::vector<TObject*> keep;
  for (int w = 0; w < 2; w++) {
    c.cd(w + 1); gPad->SetLeftMargin(0.13);
    TFile* f = TFile::Open(Form("%sm8/compare_L2_weak%s.root", T2, ws[w]));
    TLegend* leg = new TLegend(0.15, 0.15, 0.45, 0.32); leg->SetBorderSize(0);
    for (int i = 0; i < 3; i++) {
      TH1D* r = (TH1D*)f->Get(Form("rho_pred_chain_%s", tags[i]))->Clone(Form("fr%d%d", w, i));
      r->SetDirectory(nullptr);
      TH1D* tot = (TH1D*)f->Get(Form("rho_pred_%s", tags[i]));
      for (int b = 1; b <= r->GetNbinsX(); b++) {
        double t = tot->GetBinContent(b);
        r->SetBinContent(b, t > 0 ? r->GetBinContent(b) / t : -1);
        r->SetBinError(b, 0);
      }
      r->SetTitle(Form("weak %s;M (GeV);chain / L2 table", ws[w]));
      r->SetLineColor(col[i]); r->SetLineWidth(2); r->SetMinimum(0); r->SetMaximum(1.05);
      r->GetXaxis()->SetRangeUser(0.27, 2.0);
      r->Draw(i ? "HIST SAME" : "HIST");
      leg->AddEntry(r, i == 0 ? "#pi^{+}#pi^{-}" : (i == 1 ? "#pi^{+}#pi^{+}" : "K^{+}K^{-}"), "l");
      keep.push_back(r);
    }
    leg->Draw();
  }
  c.SaveAs("fig/chain_fraction.png");
}

// dP/dkappa / Boltzmann: MC L2, L2 table, L1 table.
void dpdk(const char* weak) {
  TFile f2(Form("%sm8/compare_L2_weak%s.root", T2, weak));
  TFile f1(Form("%sm7/compare_L1_weak%s.root", T2, weak));
  const char* tag[] = {"pippim", "pippim", "pippim", "pippip"};
  const char* lab[] = {"pole rho770zer", "cont 0.30", "cont 0.60", "cont 0.45"};
  const char* tit[] = {"#pi^{+}#pi^{-}, #rho^{0} pole bin", "#pi^{+}#pi^{-}, M = 0.30 GeV",
                       "#pi^{+}#pi^{-}, M = 0.60 GeV", "#pi^{+}#pi^{+}, M = 0.45 GeV"};
  TCanvas c("c", "", 1600, 1200);
  c.Divide(2, 2);
  std::vector<TObject*> keep;
  for (int i = 0; i < 4; i++) {
    c.cd(i + 1); gPad->SetLeftMargin(0.14);
    TH1D* mc = (TH1D*)f2.Get(Form("dPdkRatio_L2_%s_%s", tag[i], lab[i]));
    TH1D* t2 = (TH1D*)f2.Get(Form("dPdkRatio_pred_%s_%s", tag[i], lab[i]));
    TH1D* t1 = (TH1D*)f1.Get(Form("dPdkRatio_pred_%s_%s", tag[i], lab[i]));
    if (!mc || !t2) continue;
    mc = (TH1D*)mc->Clone(Form("mc%d", i)); t2 = (TH1D*)t2->Clone(Form("t2%d", i));
    mc->SetTitle(Form("%s;#kappa (GeV);(dP/d#kappa) / Boltzmann", tit[i]));
    mc->GetXaxis()->SetRangeUser(0, 2.0); mc->SetMinimum(0.3); mc->SetMaximum(i == 3 ? 2.2 : 1.8);
    mc->SetMarkerStyle(20); mc->SetMarkerSize(0.6); mc->SetLineColor(kBlack);
    t2->SetLineColor(kRed + 1); t2->SetLineWidth(2);
    mc->Draw("E1");
    if (t1) {
      t1 = (TH1D*)t1->Clone(Form("t1%d", i));
      t1->SetLineColor(kBlue + 1); t1->SetLineWidth(2); t1->SetLineStyle(2); t1->Draw("HIST SAME");
      keep.push_back(t1);
    }
    t2->Draw("HIST SAME"); mc->Draw("E1 SAME");
    TLine* l = new TLine(0, 1, 2.0, 1); l->SetLineStyle(3); l->Draw();
    if (i == 0) {
      TLegend* leg = new TLegend(0.55, 0.7, 0.88, 0.88); leg->SetBorderSize(0);
      leg->AddEntry(mc, "MC, level L2", "lep"); leg->AddEntry(t2, "L2 table", "l");
      if (t1) leg->AddEntry(t1, "L1 table", "l");
      leg->Draw(); keep.push_back(leg);
    }
    keep.insert(keep.end(), {mc, t2, l});
  }
  c.SaveAs(Form("fig/dpdk_weak%s.png", weak));
}

// MC a2(kappa) at L2 (from the M3 output) with the table's a2 = 0.
void a2Plot() {
  TFile f(Form("%sm3_weakon.root", T2));
  const char* lab[] = {"pole rho770zer", "cont 0.45", "cont 0.60"};
  const int col[] = {kBlack, kRed + 1, kBlue + 1};
  TCanvas c("c", "", 1200, 800);
  gPad->SetLeftMargin(0.13);
  TLegend leg(0.6, 0.72, 0.88, 0.88); leg.SetBorderSize(0);
  std::vector<TObject*> keep;
  bool first = true;
  for (int i = 0; i < 3; i++) {
    TH1D* h = (TH1D*)f.Get(Form("a2_L2_pippim_%s", lab[i]));
    if (!h) continue;
    h = (TH1D*)h->Clone(Form("a2%d", i));
    h->SetTitle("#pi^{+}#pi^{-}, MC level L2, weak on;#kappa (GeV);a_{2} = #LTP_{2}(cos#theta*)#GT");
    h->GetXaxis()->SetRangeUser(0, 2.0); h->SetMinimum(-0.25); h->SetMaximum(0.05);
    h->SetLineColor(col[i]); h->SetMarkerColor(col[i]); h->SetMarkerStyle(20); h->SetMarkerSize(0.6);
    h->Draw(first ? "E1" : "E1 SAME");
    first = false;
    leg.AddEntry(h, lab[i], "lep");
    keep.push_back(h);
  }
  TLine l(0, 0, 2.0, 0); l.SetLineStyle(2); l.Draw();
  leg.Draw();
  c.SaveAs("fig/a2_L2.png");
}

void make_figures() {
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  gSystem->mkdir("fig", true);
  for (const char* w : {"on", "off"}) {
    rhoPlot(w, "pippim", 0.27, 2.0, false);
    rhoPlot(w, "pippip", 0.27, 2.0, true);
    rhoPlot(w, "KpKm", 0.98, 2.0, false);
    dpdk(w);
  }
  chainFraction();
  a2Plot();
}
