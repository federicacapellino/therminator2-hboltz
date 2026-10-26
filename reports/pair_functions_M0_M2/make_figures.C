// Figures for the M0-M2 report.
const char* T2 = "/home/alice/therminator2-hboltz/build/t2/";
const char* T1 = "/home/alice/therminator2-hboltz/build/t1/";

void ratioPlot(TH1* mc, TH1* pred, const char* xt, const char* yt, const char* out, bool logy,
               double rlo, double rhi, double xlo = 0, double xhi = 0, bool pulls = false) {
  TCanvas c("c", "", 1400, 1200);
  TPad p1("p1", "", 0, 0.32, 1, 1); p1.SetBottomMargin(0.02); p1.SetLeftMargin(0.14); p1.Draw();
  TPad p2("p2", "", 0, 0, 1, 0.32); p2.SetTopMargin(0.02); p2.SetBottomMargin(0.32); p2.SetLeftMargin(0.14); p2.Draw();
  p1.cd(); if (logy) p1.SetLogy();
  pred->SetLineColor(kRed + 1); pred->SetLineWidth(3);
  mc->SetMarkerStyle(20); mc->SetMarkerSize(0.4); mc->SetLineColor(kBlack);
  if (xhi > xlo) { mc->GetXaxis()->SetRangeUser(xlo, xhi); pred->GetXaxis()->SetRangeUser(xlo, xhi); }
  mc->SetTitle(""); mc->GetYaxis()->SetTitle(yt); mc->GetXaxis()->SetLabelSize(0);
  mc->GetYaxis()->SetTitleSize(0.055); mc->GetYaxis()->SetLabelSize(0.045);
  mc->Draw("E1"); pred->Draw("HIST SAME"); mc->Draw("E1 SAME");
  TLegend leg(0.6, 0.75, 0.88, 0.88); leg.SetBorderSize(0);
  leg.AddEntry(mc, "MC, level L0", "lep"); leg.AddEntry(pred, "table", "l"); leg.Draw();
  p2.cd();
  TH1* r = (TH1*)mc->Clone("r");
  if (pulls) {
    for (int b = 1; b <= r->GetNbinsX(); b++) {
      double e = mc->GetBinError(b);
      r->SetBinContent(b, (e > 0 && pred->GetBinContent(b) > 0) ? (mc->GetBinContent(b) - pred->GetBinContent(b)) / e : -99);
      r->SetBinError(b, 0);
    }
    r->SetMarkerSize(0.35);
  } else {
    for (int b = 1; b <= r->GetNbinsX(); b++) {
      double q = pred->GetBinContent(b);
      r->SetBinContent(b, q > 0 ? mc->GetBinContent(b) / q : -99);
      r->SetBinError(b, q > 0 ? mc->GetBinError(b) / q : 0);
    }
  }
  r->SetMinimum(rlo); r->SetMaximum(rhi);
  r->GetYaxis()->SetTitle(pulls ? "pull" : "MC / table"); r->GetXaxis()->SetTitle(xt);
  r->GetYaxis()->SetNdivisions(505); r->GetYaxis()->SetTitleSize(0.11); r->GetYaxis()->SetTitleOffset(0.55);
  r->GetYaxis()->SetLabelSize(0.09); r->GetXaxis()->SetLabelSize(0.10); r->GetXaxis()->SetTitleSize(0.11);
  r->GetXaxis()->SetLabelSize(0.10);
  r->Draw(pulls ? "P" : "E1");
  double x1 = (xhi > xlo) ? xlo : r->GetXaxis()->GetXmin(), x2 = (xhi > xlo) ? xhi : r->GetXaxis()->GetXmax();
  TLine l(x1, pulls ? 0 : 1, x2, pulls ? 0 : 1); l.SetLineColor(kRed + 1); l.SetLineStyle(2); l.Draw();
  gStyle->SetPaperSize(c.GetWw() / 35., c.GetWh() / 35.); c.SaveAs(out);
}

void figs() {
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  TFile f(Form("%scompare_weakon.root", T2));
  // rho(M), divided by bin width -> d rho / dM per unit volume
  for (const char* tag : {"pippim", "KpKm"}) {
    TH1D* mc = (TH1D*)f.Get(Form("rho_L0_%s", tag))->Clone();
    TH1D* pr = (TH1D*)f.Get(Form("rho_pred_%s", tag))->Clone();
    mc->Scale(1., "width"); pr->Scale(1., "width");
    bool pi = std::string(tag) == "pippim";
    ratioPlot(mc, pr, "M (GeV)", "d#rho/dM  (GeV^{2})", Form("fig/rho_%s.png", tag), true,
              -4.5, 4.5, pi ? 0.28 : 0.98, 2.0, true);
  }
  // dP/dkappa
  const char* lab[] = {"pole rho770zer", "cont 0.45", "pole Ka0492zer"};
  const char* out[] = {"fig/dpdk_rho.png", "fig/dpdk_cont045.png", "fig/dpdk_K0.png"};
  for (int i = 0; i < 3; i++) {
    TH1D* mc = (TH1D*)f.Get(Form("dPdk_L0_pippim_%s", lab[i]));
    TH1D* pr = (TH1D*)f.Get(Form("dPdk_pred_pippim_%s", lab[i]));
    ratioPlot(mc, pr, "#kappa (GeV)", "dP/d#kappa  (GeV^{-1})", out[i], true, 0.98, 1.02);
  }
  // a1, a2 pulls, pi+ pi-, weak on
  TFile fm(Form("%smc_weakon.root", T2));
  TH1D h1("h1", ";pull a_{#it{l}} / #sigma;bins", 40, -5, 5), h2("h2", "", 40, -5, 5);
  for (int l = 1; l <= 2; l++) {
    TProfile2D* p = (TProfile2D*)fm.Get(Form("hP%d_L0_pippim", l));
    for (int ix = 1; ix <= p->GetNbinsX(); ix++) for (int iy = 1; iy <= p->GetNbinsY(); iy++) {
      int b = p->GetBin(ix, iy);
      if (p->GetBinEntries(b) < 50 || p->GetBinError(b) <= 0) continue;
      (l == 1 ? h1 : h2).Fill(p->GetBinContent(b) / p->GetBinError(b));
    }
  }
  TCanvas c("c", "", 1400, 1000); c.SetLeftMargin(0.13);
  h1.SetLineColor(kBlue + 1); h1.SetLineWidth(3); h2.SetLineColor(kRed + 1); h2.SetLineWidth(3);
  h1.SetMaximum(1.25 * std::max(h1.GetMaximum(), h2.GetMaximum()));
  h1.Draw("HIST"); h2.Draw("HIST SAME");
  TF1 g("g", "[0]*exp(-0.5*x*x)", -5, 5); g.SetParameter(0, h1.GetEntries() * 0.25 / sqrt(2 * M_PI));
  g.SetLineColor(kBlack); g.SetLineStyle(2); g.Draw("SAME");
  TLegend leg(0.62, 0.72, 0.88, 0.88); leg.SetBorderSize(0);
  leg.AddEntry(&h1, "a_{1}", "l"); leg.AddEntry(&h2, "a_{2}", "l"); leg.AddEntry(&g, "N(0,1)", "l"); leg.Draw();
  gStyle->SetPaperSize(c.GetWw() / 35., c.GetWh() / 35.); c.SaveAs("fig/pulls_pippim.png");
  // C_K, C_u for omega -> pi+ pi- pi0 at fixed M (table, weak on)
  TFile ft(Form("%stable_weakon.root", T1));
  TH2* ck = (TH2*)ft.Get("channels/pippim/CK_om0782zer_pi0139plu_pi0139min_pi0135zer");
  TH2* cu = (TH2*)ft.Get("channels/pippim/Cu_om0782zer_pi0139plu_pi0139min_pi0135zer");
  TCanvas c2("c2", "", 1400, 1000); c2.SetLeftMargin(0.14);
  int bx = ck->GetXaxis()->FindBin(0.30);
  TH1D* hk = ck->ProjectionY("hk", bx, bx); TH1D* hu = cu->ProjectionY("hu", bx, bx);
  TH1D* hs = (TH1D*)hk->Clone("hs"); hs->Add(hu);
  hk->SetLineColor(kBlue + 1); hu->SetLineColor(kGreen + 2); hs->SetLineColor(kBlack);
  for (TH1D* h : {hk, hu, hs}) { h->SetLineWidth(3); h->GetXaxis()->SetRangeUser(0, 1.0); }
  hs->SetTitle(Form(";#kappa (GeV);C (GeV^{-4})  at M = %.3f GeV", ck->GetXaxis()->GetBinCenter(bx)));
  hs->SetMinimum(1.3 * hk->GetMinimum()); hs->SetMaximum(1.15 * hs->GetMaximum());
  hs->Draw("HIST"); hk->Draw("HIST SAME"); hu->Draw("HIST SAME");
  TLine z(0, 0, 1.0, 0); z.SetLineStyle(2); z.Draw();
  TLegend leg2(0.6, 0.7, 0.88, 0.88); leg2.SetBorderSize(0);
  leg2.AddEntry(hk, "C_{K}", "l"); leg2.AddEntry(hu, "C_{u}", "l"); leg2.AddEntry(hs, "C_{K}+C_{u}", "l"); leg2.Draw();
  gStyle->SetPaperSize(c2.GetWw() / 35., c2.GetWh() / 35.); c2.SaveAs("fig/ckcu_omega.png");
}
