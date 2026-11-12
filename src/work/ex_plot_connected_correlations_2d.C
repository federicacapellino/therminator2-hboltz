// ex_plot_connected_correlations_2d.C — plot histograms from ex_connected_correlations_2d.exe
// output: the two RHS terms of Eq. (2) (2-body, uncorrelated product) and their difference
// (connected correlation). Draws the 2D (Delta y vs Delta pT, Delta y vs Delta phi) projections as
// COLZ panels, plus 1D overlays (Delta y, Delta pT, Delta phi, and the pair invariant mass M^2)
// with the three terms superimposed on one canvas each, plus the comoving energies eps1/eps2
// (p.u of each paired particle in the local fluid rest frame) overlaid on their own canvas.
//
// Usage:
//   root -b -q 'ex_plot_connected_correlations_2d.C("corr2d_pi_pi.root", "ss1")'
// or interactively:
//   root corr2d_pi_pi.root ex_plot_connected_correlations_2d.C
//
// A second function, ex_plot_compare(), overlays one specific histogram across an arbitrary list
// of input files/labels -- e.g. to compare HRG vs BlastWave, or different species, the same way
// ex_plot_mass_distribution.C overlays hMass from two files:
//   root -l ex_plot_connected_correlations_2d.C
//   ex_plot_compare({"corr2d_pi_pi_hrg.root","corr2d_pi_pi_bw.root"}, {"HRG","BlastWave"}, "conn", "dy", "ss1")
//   ex_plot_compare({"corr2d_pi_pi.root","corr2d_K_K.root","corr2d_p_p.root"}, {"pi pi","K K","p p"}, "conn", "dy", "ss1")
//   ex_plot_compare({"corr2d_pi_pi_hrg.root","corr2d_pi_pi_bw.root"}, {"HRG","BlastWave"}, "eps1", "", "ss1")

#include <string>
#include <vector>

void ex_plot_connected_correlations_2d(const char* filename = "corr2d_pi_pi.root",
                                        const char* category = "ss1") {
  TFile* f = TFile::Open(filename);
  if (!f || f->IsZombie()) {
    printf("Cannot open %s\n", filename);
    return;
  }

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kTRUE);

  auto drawTriplet = [&](const char* varTag, const char* varLabel, const char* outfile) {
    TH2D* h2body  = (TH2D*) f->Get(Form("h2_2body_%s_%s",  category, varTag));
    TH2D* huncorr = (TH2D*) f->Get(Form("h2_uncorr_%s_%s", category, varTag));
    TH2D* hconn   = (TH2D*) f->Get(Form("h2_conn_%s_%s",   category, varTag));

    if (!h2body || !huncorr || !hconn) {
      printf("Missing histograms for category '%s', variable '%s' in %s\n", category, varTag, filename);
      return;
    }

    TCanvas* c = new TCanvas(Form("c_%s", varTag), Form("Connected correlations (%s)", varLabel), 1500, 500);
    c->Divide(3, 1);

    c->cd(1);
    h2body->SetTitle(Form("2-body: #LT E_{1}dN_{s}/d^{3}p_{1} E_{2}dN_{s'}/d^{3}p_{2} #GT;#Delta y;%s", varLabel));
    h2body->Draw("COLZ");

    c->cd(2);
    huncorr->SetTitle(Form("Uncorrelated: #LT E_{1}dN_{s}/d^{3}p_{1}#GT#LT E_{2}dN_{s'}/d^{3}p_{2} #GT;#Delta y;%s", varLabel));
    huncorr->Draw("COLZ");

    c->cd(3);
    hconn->SetTitle(Form("Connected C_{ss',conn};#Delta y;%s", varLabel));
    hconn->Draw("COLZ");

    c->Print(outfile);
  };

  drawTriplet("pt",  "#Delta p_{T} [GeV/c]", Form("connected_correlations_%s_pt.pdf",  category));
  drawTriplet("phi", "#Delta#phi [rad]",     Form("connected_correlations_%s_phi.pdf", category));

  // --- 1D overlays: 2-body / uncorrelated / connected superimposed on one canvas ---
  auto drawOverlay1D = [&](const char* varTag, const char* varLabel, const char* outfile, bool logy) {
    TH1D* h2body  = (TH1D*) f->Get(Form("h1_2body_%s_%s",  category, varTag));
    TH1D* huncorr = (TH1D*) f->Get(Form("h1_uncorr_%s_%s", category, varTag));
    TH1D* hconn   = (TH1D*) f->Get(Form("h1_conn_%s_%s",   category, varTag));

    if (!h2body || !huncorr || !hconn) {
      printf("Missing 1D histograms for category '%s', variable '%s' in %s\n", category, varTag, filename);
      return;
    }

    TCanvas* c = new TCanvas(Form("c1_%s", varTag), Form("Connected correlations 1D (%s)", varLabel), 600, 500);
    if (logy) c->SetLogy();

    h2body->SetLineColor(kBlue + 1);
    h2body->SetLineWidth(2);
    h2body->GetXaxis()->SetTitle(varLabel);
    h2body->GetXaxis()->SetTitleSize(0.05);
    h2body->GetXaxis()->SetTitleOffset(0.9);
    h2body->GetXaxis()->SetLabelSize(0.045);
    h2body->GetYaxis()->SetTitle("entries");
    h2body->GetYaxis()->SetTitleSize(0.05);
    h2body->GetYaxis()->SetTitleOffset(1.4);
    h2body->GetYaxis()->SetLabelSize(0.045);
    h2body->Draw("hist e");

    huncorr->SetLineColor(kRed + 1);
    huncorr->SetLineWidth(2);
    huncorr->SetLineStyle(2);
    huncorr->Draw("hist e same");

    hconn->SetLineColor(kGreen + 2);
    hconn->SetLineWidth(2);
    hconn->Draw("hist e same");

    TLegend* leg = new TLegend(0.45, 0.7, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(h2body,  "2-body", "l");
    leg->AddEntry(huncorr, "Uncorrelated", "l");
    leg->AddEntry(hconn,   "Connected (difference)", "l");
    leg->Draw();

    c->Print(outfile);
  };

  drawOverlay1D("dy",  "#Delta y",                       Form("connected_correlations_%s_dy.pdf",  category), false);
  drawOverlay1D("pt",  "#Delta p_{T} [GeV/c]",            Form("connected_correlations_%s_1d_pt.pdf",  category), false);
  drawOverlay1D("phi", "#Delta#phi [rad]",                Form("connected_correlations_%s_1d_phi.pdf", category), false);
  drawOverlay1D("M2",  "M^{2} [(GeV/c^{2})^{2}]",         Form("connected_correlations_%s_M2.pdf",  category), true);

  // --- eps1/eps2 overlay: comoving energy of each paired particle, no uncorrelated/connected
  // counterpart (single-particle quantity), so just a two-series overlay.
  {
    TH1D* hEps1 = (TH1D*) f->Get(Form("h1_eps1_%s", category));
    TH1D* hEps2 = (TH1D*) f->Get(Form("h1_eps2_%s", category));

    if (!hEps1 || !hEps2) {
      printf("Missing eps histograms for category '%s' in %s\n", category, filename);
    } else {
      TCanvas* c = new TCanvas("c1_eps", "Comoving energy", 600, 500);
      c->SetLogy();

      hEps1->SetLineColor(kBlue + 1);
      hEps1->SetLineWidth(2);
      hEps1->GetXaxis()->SetTitle("#epsilon = p#upoint u [GeV]");
      hEps1->GetXaxis()->SetTitleSize(0.05);
      hEps1->GetXaxis()->SetTitleOffset(0.9);
      hEps1->GetXaxis()->SetLabelSize(0.045);
      hEps1->GetYaxis()->SetTitle("entries");
      hEps1->GetYaxis()->SetTitleSize(0.05);
      hEps1->GetYaxis()->SetTitleOffset(1.4);
      hEps1->GetYaxis()->SetLabelSize(0.045);
      hEps1->Draw("hist e");

      hEps2->SetLineColor(kRed + 1);
      hEps2->SetLineWidth(2);
      hEps2->SetLineStyle(2);
      hEps2->Draw("hist e same");

      TLegend* leg = new TLegend(0.45, 0.75, 0.88, 0.88);
      leg->SetBorderSize(0);
      leg->SetTextSize(0.04);
      leg->AddEntry(hEps1, "#epsilon_{1} (species s)",  "l");
      leg->AddEntry(hEps2, "#epsilon_{2} (species s')", "l");
      leg->Draw();

      c->Print(Form("connected_correlations_%s_eps.pdf", category));
    }
  }
}

// Overlay one specific 1D histogram across an arbitrary number of input files, one color per
// file/label -- generalizes the two-file HRG vs BlastWave overlay in ex_plot_mass_distribution.C
// to an arbitrary file count, so the same function also covers e.g. comparing several species'
// output files.
//   term:     "2body" | "uncorr" | "conn" | "eps1" | "eps2"
//   var:      "dy" | "pt" | "phi" | "M2" (ignored when term is "eps1"/"eps2")
//   category: "ss1" | "s" | "s1"
void ex_plot_compare(std::vector<std::string> files, std::vector<std::string> labels,
                     const char* term = "conn", const char* var = "dy",
                     const char* category = "ss1") {
  if (files.size() != labels.size() || files.empty()) {
    printf("ex_plot_compare: files and labels must be the same non-zero size (got %zu files, %zu labels)\n",
           files.size(), labels.size());
    return;
  }

  std::string termStr(term);
  bool isEps = (termStr.rfind("eps", 0) == 0); // term starts with "eps"

  std::string histName = isEps ? Form("h1_%s_%s", term, category)
                               : Form("h1_%s_%s_%s", term, category, var);

  std::string varLabel;
  if (isEps)                                varLabel = "#epsilon = p#upoint u [GeV]";
  else if (std::string(var) == "dy")         varLabel = "#Delta y";
  else if (std::string(var) == "pt")         varLabel = "#Delta p_{T} [GeV/c]";
  else if (std::string(var) == "phi")        varLabel = "#Delta#phi [rad]";
  else if (std::string(var) == "M2")         varLabel = "M^{2} [(GeV/c^{2})^{2}]";
  else                                       varLabel = var;

  gStyle->SetPadLeftMargin(0.12);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetOptStat(kFALSE);
  gStyle->SetOptTitle(kFALSE);

  static const Int_t colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kOrange + 1, kCyan + 2, kMagenta + 1};
  static const int nColors = sizeof(colors) / sizeof(colors[0]);

  TCanvas* c = new TCanvas(Form("cCompare_%s", histName.c_str()), Form("Compare: %s", histName.c_str()), 600, 500);

  TLegend* leg = new TLegend(0.5, 0.68, 0.88, 0.88);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);

  // Keep the open files alive for the lifetime of the canvas (the histograms are owned by them).
  static std::vector<TFile*> keepOpen;

  TH1D* hFirst = nullptr;
  for (size_t i = 0; i < files.size(); ++i) {
    TFile* f = TFile::Open(files[i].c_str());
    if (!f || f->IsZombie()) {
      printf("ex_plot_compare: cannot open %s, skipping\n", files[i].c_str());
      continue;
    }
    keepOpen.push_back(f);

    TH1D* h = (TH1D*) f->Get(histName.c_str());
    if (!h) {
      printf("ex_plot_compare: cannot find %s in %s, skipping\n", histName.c_str(), files[i].c_str());
      continue;
    }

    int color = colors[i % nColors];
    h->SetLineColor(color);
    h->SetLineWidth(2);

    if (!hFirst) {
      hFirst = h;
      hFirst->GetXaxis()->SetTitle(varLabel.c_str());
      hFirst->GetXaxis()->SetTitleSize(0.05);
      hFirst->GetXaxis()->SetTitleOffset(0.9);
      hFirst->GetXaxis()->SetLabelSize(0.045);
      hFirst->GetYaxis()->SetTitle("entries");
      hFirst->GetYaxis()->SetTitleSize(0.05);
      hFirst->GetYaxis()->SetTitleOffset(1.4);
      hFirst->GetYaxis()->SetLabelSize(0.045);
      hFirst->Draw("hist e");
    } else {
      h->Draw("hist e same");
    }

    leg->AddEntry(h, labels[i].c_str(), "l");
  }

  if (!hFirst) {
    printf("ex_plot_compare: no histograms found across the given files, nothing to plot\n");
    delete c;
    return;
  }

  leg->Draw();

  std::string outfile = isEps ? Form("compare_%s_%s.pdf", term, category)
                              : Form("compare_%s_%s_%s.pdf", term, var, category);
  c->Print(outfile.c_str());
}
