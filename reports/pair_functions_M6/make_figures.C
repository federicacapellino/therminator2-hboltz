// Figures for the M6 report (single-particle currents).
// Run: root -l -b -q make_figures.C
const char* M6 = "/home/alice/therminator2-hboltz/build/t2/m6/";

// MC dN/dp vs the currents, eight species, with ratio panels.
void spectra(const char* weak) {
  TFile f(Form("%scurrents_weak%s_mc.root", M6, weak));
  const int pdg[8] = {211, 321, 2212, 113, 313, 223, 3122, 2224};
  const char* lab[8] = {"#pi^{+}", "K^{+}", "p", "#rho^{0}", "K^{*0}", "#omega", "#Lambda",
                        "#Delta^{++}"};
  TCanvas c("c", "", 2000, 1100);
  c.Divide(4, 2, 0.002, 0.002);
  std::vector<TObject*> keep;
  for (int i = 0; i < 8; i++) {
    TH1D* mc = (TH1D*)f.Get(Form("hSpecMC_%d", pdg[i]));
    TH1D* pr = (TH1D*)f.Get(Form("hSpecPred_%d", pdg[i]));
    if (!mc || !pr) continue;
    c.cd(i + 1);
    TPad* p1 = new TPad(Form("a%d", i), "", 0, 0.35, 1, 1);
    TPad* p2 = new TPad(Form("b%d", i), "", 0, 0, 1, 0.35);
    p1->SetBottomMargin(0.02); p1->SetLeftMargin(0.16); p1->SetLogy(); p1->Draw();
    p2->SetTopMargin(0.02); p2->SetBottomMargin(0.3); p2->SetLeftMargin(0.16); p2->Draw();
    p1->cd();
    mc->SetTitle(""); mc->SetMarkerStyle(20); mc->SetMarkerSize(0.4);
    mc->GetXaxis()->SetLabelSize(0); mc->GetYaxis()->SetTitle("dN/dp (counts per 10 MeV)");
    mc->GetYaxis()->SetTitleSize(0.06); mc->GetYaxis()->SetLabelSize(0.05);
    mc->SetMinimum(0.5);
    pr->SetLineColor(kRed + 1); pr->SetLineWidth(2);
    mc->Draw("E1"); pr->Draw("HIST SAME"); mc->Draw("E1 SAME");
    TLatex* t = new TLatex(0.62, 0.8, lab[i]); t->SetNDC(); t->SetTextSize(0.09); t->Draw();
    p2->cd();
    TH1D* r = (TH1D*)mc->Clone(Form("r%d", i));
    for (int b = 1; b <= r->GetNbinsX(); b++) {
      double q = pr->GetBinContent(b);
      r->SetBinContent(b, q > 5 ? mc->GetBinContent(b) / q : -99);
      r->SetBinError(b, q > 5 ? mc->GetBinError(b) / q : 0);
    }
    r->SetMinimum(0.9); r->SetMaximum(1.1);
    r->GetYaxis()->SetTitle("MC / S_{c}"); r->GetYaxis()->SetNdivisions(505);
    r->GetXaxis()->SetTitle("p (GeV)");
    r->GetYaxis()->SetTitleSize(0.1); r->GetYaxis()->SetTitleOffset(0.7);
    r->GetYaxis()->SetLabelSize(0.09); r->GetXaxis()->SetLabelSize(0.1);
    r->GetXaxis()->SetTitleSize(0.11);
    r->Draw("E1");
    TLine* l = new TLine(0, 1, 3, 1); l->SetLineColor(kRed + 1); l->SetLineStyle(2); l->Draw();
    keep.insert(keep.end(), {p1, p2, t, r, l});
  }
  c.SaveAs(Form("fig/spectra_weak%s.png", weak));
}

// Pull distribution of the species densities, both weak settings.
void pulls() {
  TCanvas c("c", "", 1200, 800);
  TFile fon(Form("%scurrents_weakon_mc.root", M6)), foff(Form("%scurrents_weakoff_mc.root", M6));
  TH1D* hon  = (TH1D*)fon.Get("hPullSpecies");
  TH1D* hoff = (TH1D*)foff.Get("hPullSpecies");
  hon->SetLineColor(kBlue + 1); hon->SetLineWidth(2);
  hoff->SetLineColor(kRed + 1); hoff->SetLineWidth(2); hoff->SetLineStyle(2);
  hon->SetTitle(";(#LT N_{c}#GT/V - n_{c}) / #sigma;species");
  hon->SetMaximum(1.3 * std::max(hon->GetMaximum(), hoff->GetMaximum()));
  hon->Draw("HIST"); hoff->Draw("HIST SAME");
  TF1 g("g", "[0]*exp(-0.5*x*x)", -5, 5);
  g.SetParameter(0, hon->GetEntries() * hon->GetBinWidth(1) / std::sqrt(2 * M_PI));
  g.SetLineColor(kGray + 2); g.SetLineStyle(3); g.Draw("SAME");
  TLegend leg(0.62, 0.72, 0.88, 0.88); leg.SetBorderSize(0);
  leg.AddEntry(hon, "weak on", "l"); leg.AddEntry(hoff, "weak off", "l");
  leg.AddEntry(&g, "unit Gaussian", "l"); leg.Draw();
  c.SaveAs("fig/pulls.png");
}

// Currents of a few species: A+B and B/(A+B) (zero for a thermal species) versus k.
void currents() {
  TFile f(Form("%scurrents_weakon_quad.root", M6));
  TTree* t = (TTree*)f.Get("currents");
  char name[64]; double m;
  std::vector<double>*k = nullptr, *A = nullptr, *AB = nullptr;
  t->SetBranchAddress("name", name); t->SetBranchAddress("mass", &m);
  t->SetBranchAddress("k", &k); t->SetBranchAddress("A", &A); t->SetBranchAddress("AB", &AB);
  const char* show[] = {"pi0139plu", "Ka0492plu", "pr0938plu", "Lm1115zer", "rho770zer"};
  const int col[] = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1};
  TCanvas c("c", "", 1600, 700);
  c.Divide(2, 1);
  TMultiGraph* mg1 = new TMultiGraph(); TMultiGraph* mg2 = new TMultiGraph();
  TLegend* leg = new TLegend(0.6, 0.6, 0.88, 0.88); leg->SetBorderSize(0);
  for (Long64_t e = 0; e < t->GetEntries(); e++) {
    t->GetEntry(e);
    for (int s = 0; s < 5; s++) {
      if (std::string(name) != show[s]) continue;
      TGraph* g1 = new TGraph(); TGraph* g2 = new TGraph();
      for (size_t i = 0; i < k->size(); i++) {
        if ((*k)[i] > 2.0) break;
        g1->SetPoint(g1->GetN(), (*k)[i], (*AB)[i]);
        g2->SetPoint(g2->GetN(), (*k)[i], ((*AB)[i] - (*A)[i]) / (*AB)[i]);
      }
      g1->SetLineColor(col[s]); g1->SetLineWidth(2);
      g2->SetLineColor(col[s]); g2->SetLineWidth(2);
      mg1->Add(g1, "L"); mg2->Add(g2, "L");
      leg->AddEntry(g1, name, "l");
    }
  }
  c.cd(1); gPad->SetLogy(); gPad->SetLeftMargin(0.14);
  mg1->SetTitle(";k (GeV);A_{c}+B_{c}  (GeV^{-3})"); mg1->Draw("A"); leg->Draw();
  c.cd(2); gPad->SetLeftMargin(0.14);
  mg2->SetTitle(";k (GeV);B_{c} / (A_{c}+B_{c})"); mg2->Draw("A");
  c.SaveAs("fig/currents.png");
}

void make_figures() {
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  gSystem->mkdir("fig", true);
  spectra("on");
  spectra("off");
  pulls();
  currents();
}
