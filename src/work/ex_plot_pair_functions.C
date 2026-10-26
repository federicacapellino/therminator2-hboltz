// ex_plot_pair_functions.C — plots of README_pair_functions.md Sec. 5 (milestone M3) from the
// output of ex_pair_functions_compare.exe: MC at the cumulative ancestry levels L0, L1, L2 versus
// the pair-function table.
//   1. rho(M) per unit volume: L0, L1, L2 and the table, with ratio panels; the table split by
//      parent (rho_levels_*, rho_decomp_*)
//   2. dP/dkappa at fixed M divided by the Boltzmann reference of a particle of mass M (dPdk_*)
//   3. Legendre moments a1(kappa), a2(kappa) at fixed M (a1_*, a2_*) and raw cos theta* (cth_*)
//   4. a1, a2 and C_Q/(C_K + C_u) in (M, kappa) cells at L2 (cells_*)
//
// Usage:
//   root -b -q 'ex_plot_pair_functions.C("compare_weakon.root", "weakon", "plots")'

#include <string>
#include <vector>

namespace pfplot {

const char* kLev[3]   = {"L0", "L1", "L2"};
const int kLevCol[3]  = {kBlack, kBlue + 1, kGreen + 2};
const int kLevMark[3] = {20, 21, 33};
const int kPredCol    = kRed + 1;

std::string PairLabel(const std::string& tag) {
  if (tag == "pippim") return "#pi^{+}#pi^{-}";
  if (tag == "pippip") return "#pi^{+}#pi^{+}";
  if (tag == "pimpim") return "#pi^{-}#pi^{-}";
  if (tag == "KpKm") return "K^{+}K^{-}";
  if (tag == "KpKp") return "K^{+}K^{+}";
  if (tag == "KmKm") return "K^{-}K^{-}";
  return tag;
}

void Style(TH1* h, int l) {
  h->SetLineColor(kLevCol[l]);
  h->SetMarkerColor(kLevCol[l]);
  h->SetMarkerStyle(kLevMark[l]);
  h->SetMarkerSize(0.5);
}

void Save(TCanvas& c, const std::string& out) { c.SaveAs(out.c_str()); }

// 1. rho(M): levels and table with ratio to the table; the table by parent.
void PlotRho(TFile& f, const std::string& tag, const std::string& wl, const std::string& dir,
             double xlo) {
  TH1D* pred = (TH1D*)f.Get(("rho_pred_" + tag).c_str());
  if (!pred) return;
  TH1D* h[3];
  bool any = false;
  for (int l = 0; l < 3; l++) {
    h[l] = (TH1D*)f.Get(Form("rho_%s_%s", kLev[l], tag.c_str()));
    if (h[l]) {
      h[l] = (TH1D*)h[l]->Clone();
      h[l]->Scale(1., "width");
      Style(h[l], l);
      any = any || h[l]->Integral() > 0;
    }
  }
  if (!any) return;
  pred = (TH1D*)pred->Clone();
  pred->Scale(1., "width");
  pred->SetLineColor(kPredCol);
  pred->SetLineWidth(2);
  const bool predEmpty = pred->Integral() <= 0;

  TCanvas c("c", "", 1400, 1200);
  TPad p1("p1", "", 0, 0.32, 1, 1);
  p1.SetBottomMargin(0.02);
  p1.SetLeftMargin(0.14);
  p1.Draw();
  TPad p2("p2", "", 0, 0, 1, 0.32);
  p2.SetTopMargin(0.02);
  p2.SetBottomMargin(0.32);
  p2.SetLeftMargin(0.14);
  p2.Draw();
  p1.cd();
  p1.SetLogy();
  TH1D* frame = h[2];
  frame->SetTitle("");
  frame->GetXaxis()->SetRangeUser(xlo, 2.0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetTitle("d#rho/dM  (GeV^{2})");
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetLabelSize(0.045);
  double mn = 1e30;
  for (int b = 1; b <= frame->GetNbinsX(); b++) {
    for (int l = 0; l < 3; l++)
      if (h[l]->GetBinContent(b) > 0) mn = std::min(mn, h[l]->GetBinContent(b));
  }
  frame->SetMinimum(0.5 * mn);
  frame->SetMaximum(5. * frame->GetMaximum());
  frame->Draw("E1");
  for (int l = 0; l < 2; l++) h[l]->Draw("E1 SAME");
  if (!predEmpty) pred->Draw("HIST SAME");
  TLegend leg(0.62, 0.68, 0.88, 0.88);
  leg.SetBorderSize(0);
  leg.SetHeader(Form("%s, %s", PairLabel(tag).c_str(), wl.c_str()));
  for (int l = 0; l < 3; l++) leg.AddEntry(h[l], Form("MC %s", kLev[l]), "lep");
  if (!predEmpty) leg.AddEntry(pred, "table", "l");
  leg.Draw();

  p2.cd();
  if (!predEmpty) {
    p2.SetLogy();
    TH1D* r[3];
    for (int l = 0; l < 3; l++) {
      r[l] = (TH1D*)h[l]->Clone();
      for (int b = 1; b <= r[l]->GetNbinsX(); b++) {
        double q = pred->GetBinContent(b);
        r[l]->SetBinContent(b, q > 0 ? h[l]->GetBinContent(b) / q : 0);
        r[l]->SetBinError(b, q > 0 ? h[l]->GetBinError(b) / q : 0);
      }
    }
    r[0]->SetTitle("");
    r[0]->SetMinimum(0.5);
    r[0]->SetMaximum(tag == "pippim" || tag == "KpKm" ? 30. : 3000.);
    r[0]->GetYaxis()->SetTitle("MC / table");
    r[0]->GetXaxis()->SetTitle("M (GeV)");
    r[0]->GetXaxis()->SetRangeUser(xlo, 2.0);
    r[0]->GetYaxis()->SetNdivisions(505);
    r[0]->GetYaxis()->SetTitleSize(0.11);
    r[0]->GetYaxis()->SetTitleOffset(0.55);
    r[0]->GetYaxis()->SetLabelSize(0.09);
    r[0]->GetXaxis()->SetLabelSize(0.10);
    r[0]->GetXaxis()->SetTitleSize(0.11);
    r[0]->Draw("P");
    r[1]->Draw("P SAME");
    r[2]->Draw("P SAME");
    TLine one(xlo, 1, 2.0, 1);
    one.SetLineColor(kPredCol);
    one.SetLineStyle(2);
    one.DrawClone();
  } else {
    TLatex t;
    t.SetNDC();
    t.SetTextSize(0.12);
    t.DrawLatex(0.3, 0.5, "table empty (no one-step channel)");
  }
  Save(c, dir + "/rho_levels_" + tag + "_" + wl + ".png");

  // decomposition of the table by parent
  TNamed* list = (TNamed*)f.Get(("decomp_" + tag).c_str());
  if (predEmpty || !list) return;
  TCanvas c2("c2", "", 1400, 1000);
  c2.SetLogy();
  c2.SetLeftMargin(0.13);
  pred->SetTitle(Form(";M (GeV);d#rho/dM  (GeV^{2})"));
  pred->GetXaxis()->SetRangeUser(xlo, 2.0);
  pred->SetMinimum(1e-4 * pred->GetMaximum());
  pred->SetMaximum(3. * pred->GetMaximum());
  pred->SetLineColor(kBlack);
  pred->Draw("HIST");
  TLegend leg2(0.66, 0.50, 0.89, 0.89);
  leg2.SetBorderSize(0);
  leg2.SetHeader(Form("table, %s, %s", PairLabel(tag).c_str(), wl.c_str()));
  leg2.AddEntry(pred, "total", "l");
  const int cols[8] = {kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2,
                       kViolet - 5, kGray + 1};
  TObjArray* names = TString(list->GetTitle()).Tokenize(",");
  for (int i = 0; i < names->GetEntries(); i++) {
    std::string pn = ((TObjString*)names->At(i))->GetString().Data();
    TH1D* hp = (TH1D*)f.Get(("rho_pred_" + tag + "_" + pn).c_str());
    if (!hp) continue;
    hp = (TH1D*)hp->Clone();
    hp->Scale(1., "width");
    hp->SetLineColor(cols[i % 8]);
    hp->SetLineWidth(2);
    hp->Draw("HIST SAME");
    leg2.AddEntry(hp, pn.c_str(), "l");
  }
  leg2.Draw();
  Save(c2, dir + "/rho_decomp_" + tag + "_" + wl + ".png");
}

// Labels of the test masses of a pair present in the file (from the dPdkRatio_ keys).
std::vector<std::string> Labels(TFile& f, const std::string& tag) {
  std::vector<std::string> out;
  const char* all[] = {"pole Ka0492zer", "pole rho770zer", "pole om0782zer", "pole f00980zer",
                       "cont 0.30",      "cont 0.45",      "cont 0.60",      "pole ph1020zer",
                       "cont 1.10",      "cont 1.30"};
  for (const char* lab : all)
    for (int l = 0; l < 3; l++)
      if (f.Get(Form("dPdkRatio_%s_%s_%s", kLev[l], tag.c_str(), lab))) {
        out.push_back(lab);
        break;
      }
  return out;
}

void Grid(TCanvas& c, int n) {
  int nx = n <= 2 ? n : (n <= 4 ? 2 : (n <= 6 ? 3 : 4));
  int ny = (n + nx - 1) / nx;
  c.SetCanvasSize(600 * nx, 500 * ny);
  c.Divide(nx, ny);
}

// 2. dP/dkappa divided by the Boltzmann reference; 3. a1(kappa), a2(kappa).
void PlotKappa(TFile& f, const std::string& tag, const std::string& wl, const std::string& dir) {
  std::vector<std::string> labs = Labels(f, tag);
  if (labs.empty()) return;
  const bool likeSign = (tag == "pippip" || tag == "pimpim" || tag == "KpKp" || tag == "KmKm");
  for (const char* what : {"dPdkRatio", "a1", "a2"}) {
    std::string w = what;
    if (w == "a1" && likeSign) continue;
    TCanvas c("c", "", 1200, 1000);
    Grid(c, labs.size());
    std::vector<TObject*> keep;
    for (size_t i = 0; i < labs.size(); i++) {
      TPad* pad = (TPad*)c.cd(i + 1);
      pad->SetLeftMargin(0.15);
      pad->SetBottomMargin(0.12);
      TLegend* leg = new TLegend(0.55, 0.68, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetFillStyle(0);
      keep.push_back(leg);
      bool first = true;
      double lo = 1e30, hi = -1e30;
      std::vector<TH1D*> hs;
      for (int l = -1; l < 3; l++) { // -1: table
        std::string lev = (l < 0) ? "pred" : kLev[l];
        if (l < 0 && w != "dPdkRatio") continue;
        TH1D* h = (TH1D*)f.Get(Form("%s_%s_%s_%s", what, lev.c_str(), tag.c_str(), labs[i].c_str()));
        if (!h) continue;
        h = (TH1D*)h->Clone();
        h->GetXaxis()->SetRangeUser(0., 2.0);
        if (l < 0) {
          h->SetLineColor(kPredCol);
          h->SetLineWidth(3);
        } else
          Style(h, l);
        for (int b = 1; b <= h->FindBin(1.999); b++) {
          if (h->GetBinError(b) == 0 && l >= 0) continue;
          if (w == "dPdkRatio" && h->GetBinContent(b) <= 0) continue;
          lo = std::min(lo, h->GetBinContent(b) - h->GetBinError(b));
          hi = std::max(hi, h->GetBinContent(b) + h->GetBinError(b));
        }
        hs.push_back(h);
        leg->AddEntry(h, l < 0 ? "table" : Form("MC %s", kLev[l]), l < 0 ? "l" : "lep");
      }
      if (hs.empty()) continue;
      if (w == "dPdkRatio") {
        pad->SetLogy(hi / std::max(lo, 1e-6) > 5.);
        lo = std::max(lo, 1e-3);
      }
      double span = hi - lo;
      for (TH1D* h : hs) {
        h->SetTitle(Form("%s, %s;#kappa (GeV);%s", PairLabel(tag).c_str(), labs[i].c_str(),
                         w == "dPdkRatio" ? "(dP/d#kappa) / (dP_{th}/d#kappa)"
                                          : (w == "a1" ? "a_{1}" : "a_{2}")));
        if (pad->GetLogy()) {
          h->SetMinimum(0.7 * lo);
          h->SetMaximum(1.5 * hi);
        } else {
          h->SetMinimum(std::min(lo - 0.1 * span, w == "dPdkRatio" ? 0.95 : -0.01));
          h->SetMaximum(std::max(hi + 0.35 * span, w == "dPdkRatio" ? 1.05 : 0.01));
        }
        h->Draw(first ? (h->GetLineWidth() == 3 ? "HIST" : "E1") : (h->GetLineWidth() == 3 ? "HIST SAME" : "E1 SAME"));
        first = false;
        keep.push_back(h);
      }
      TLine* ln = new TLine(0., w == "dPdkRatio" ? 1. : 0., 2.0, w == "dPdkRatio" ? 1. : 0.);
      ln->SetLineStyle(2);
      ln->Draw();
      keep.push_back(ln);
      leg->Draw();
    }
    Save(c, dir + "/" + (w == "dPdkRatio" ? std::string("dPdk") : w) + "_" + tag + "_" + wl + ".png");
    for (TObject* o : keep) delete o;
  }

  // raw cos theta* (shape), two kappa bands, at the first label with L2 entries
  TCanvas c("c", "", 1400, 600);
  c.Divide(2, 1);
  std::vector<TObject*> keep;
  const std::string lab = labs.size() > 1 ? labs[1] : labs[0];
  for (int ib = 0; ib < 2; ib++) {
    TPad* pad = (TPad*)c.cd(ib + 1);
    pad->SetLeftMargin(0.15);
    TLegend* leg = new TLegend(0.6, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    keep.push_back(leg);
    bool first = true;
    for (int l = 0; l < 3; l++) {
      TH1D* h = (TH1D*)f.Get(Form("cth_%s_%s_%s_k%d", kLev[l], tag.c_str(), lab.c_str(), ib));
      if (!h || h->Integral() <= 0) continue;
      h = (TH1D*)h->Clone();
      h->Sumw2();
      h->Scale(h->GetNbinsX() / h->Integral());
      Style(h, l);
      h->SetMinimum(0.6);
      h->SetMaximum(1.4);
      h->GetYaxis()->SetTitle("normalized");
      h->Draw(first ? "E1" : "E1 SAME");
      first = false;
      leg->AddEntry(h, Form("MC %s", kLev[l]), "lep");
      keep.push_back(h);
    }
    leg->Draw();
  }
  Save(c, dir + "/cth_" + tag + "_" + wl + ".png");
  for (TObject* o : keep) delete o;
}

// 4. cell maps at L2
void PlotCells(TFile& f, const std::string& tag, const std::string& wl, const std::string& dir) {
  std::vector<std::string> names = {"a2_cells_L2_"};
  if (tag == "pippim" || tag == "KpKm") names = {"a1_cells_L2_", "a2_cells_L2_", "CQratio_cells_L2_"};
  TCanvas c("c", "", 700 * names.size(), 600);
  c.Divide(names.size(), 1);
  gStyle->SetPaintTextFormat(".4f");
  for (size_t i = 0; i < names.size(); i++) {
    TH2D* h = (TH2D*)f.Get((names[i] + tag).c_str());
    if (!h) continue;
    TPad* pad = (TPad*)c.cd(i + 1);
    pad->SetRightMargin(0.15);
    h->SetMarkerSize(0.9);
    h->Draw("COLZ TEXT E");
  }
  Save(c, dir + "/cells_" + tag + "_" + wl + ".png");
}

} // namespace pfplot

void ex_plot_pair_functions(const char* filename = "compare_weakon.root",
                            const char* weakLabel = "weakon", const char* outdir = ".") {
  using namespace pfplot;
  gStyle->SetOptStat(0);
  gErrorIgnoreLevel = kWarning;
  TFile f(filename);
  if (f.IsZombie()) {
    printf("Cannot open %s\n", filename);
    return;
  }
  gSystem->mkdir(outdir, true);
  const std::string dir = outdir, wl = weakLabel;
  struct P {
    const char* tag;
    double xlo;
  };
  for (P p : {P{"pippim", 0.27}, P{"pippip", 0.27}, P{"pimpim", 0.27}, P{"KpKm", 0.98},
              P{"KpKp", 0.98}, P{"KmKm", 0.98}}) {
    PlotRho(f, p.tag, wl, dir, p.xlo);
    PlotKappa(f, p.tag, wl, dir);
    PlotCells(f, p.tag, wl, dir);
  }
}
