// ex_plot_compare_meanpt.C — overlay hMeanPt1VsM (<pT,1> vs pair invariant mass M,
// from ex_mass_pt_distribution) from several runs on one canvas.
//
// Labels default to "HRG (<run folder>)" for folders starting with "hrg", otherwise
// "blastwave, v_{T}=<VelT> (<run folder>)" with VelT read from <run folder>/blastwave.ini.
// Empty M bins (no pairs) are not drawn.  An optional title is drawn above the frame.
//
// ex_plot_compare_meanpt_all() makes the full set of mean-pT plots from the default
// output names of ex_mass_pt_distribution in each run folder:
//   - one model comparison per pair type, written to the current folder:
//       meanpt_compare_omega_pipm.pdf  <- <run>/mass_pt_meanpt_<N>.root     (--parent 223)
//       meanpt_compare_full_pipm.pdf   <- <run>/mass_pt_full_pipm_<N>.root  (--parent 0)
//       meanpt_compare_full_pipp.pdf   <- <run>/mass_pt_full_pipp_<N>.root
//       meanpt_compare_full_pimm.pdf   <- <run>/mass_pt_full_pimm_<N>.root
//   - one charge-combination comparison per run, <run>/meanpt_full_charges.pdf,
//     overlaying full +-, ++, -- with omega +- as reference.
// For ++ and -- every pair enters twice (ordered pairs), so <pT,1> is the mean over
// both members; the error bars are slightly underestimated.
//
// Usage (from build/work, where the run folders live):
//   root -b -q ex_plot_compare_meanpt.C
//   root -b -q 'ex_plot_compare_meanpt.C("a.root,b.root", "label A,label B", "out.pdf", "title")'
//   root -b -q -e '.L ex_plot_compare_meanpt.C' -e 'ex_plot_compare_meanpt_all()'
//   root -b -q -e '.L ex_plot_compare_meanpt.C' \
//               -e 'ex_plot_compare_meanpt_all("hrg_run,bw_run", "1000")'

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TSystem.h>

// First path component of a file, e.g. "bw_run2" for "bw_run2/x.root".
static TString meanptRunFolderOf(const TString& path) {
   Ssiz_t slash = path.Index("/");
   return (slash > 0) ? TString(path(0, slash)) : TString(".");
}

// VelT from <runFolder>/blastwave.ini, or "" if not found.
static TString meanptReadVelT(const TString& runFolder) {
   std::ifstream in((runFolder + "/blastwave.ini").Data());
   std::string line;
   while (std::getline(in, line)) {
      TString l(line.c_str());
      l = l.Strip(TString::kBoth);
      if (!l.BeginsWith("VelT")) continue;
      Ssiz_t eq = l.Index("=");
      if (eq < 0) continue;
      TString v = l(eq + 1, l.Length());
      return v.Strip(TString::kBoth);
   }
   return "";
}

void ex_plot_compare_meanpt(
    const char* files  = "hrg_run/mass_pt_meanpt_500.root,"
                         "bw_run/mass_pt_meanpt_500.root,"
                         "bw_run2/mass_pt_meanpt_500.root,"
                         "bw_run3/mass_pt_meanpt_500.root",
    const char* labels = "",
    const char* output = "meanpt_compare.pdf",
    const char* title  = "") {
   TObjArray* fileTok  = TString(files).Tokenize(",");
   TObjArray* labelTok = TString(labels).Tokenize(",");

   const int colors[]  = {kBlack, kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1, kOrange + 7};
   const int markers[] = {20, 21, 22, 23, 33, 34};
   const int nColors   = sizeof(colors) / sizeof(colors[0]);

   TMultiGraph* mg  = new TMultiGraph();
   TLegend*     leg = new TLegend(0.18, 0.66, 0.70, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.035);
   int nDrawn = 0;

   for (int i = 0; i < fileTok->GetEntries(); ++i) {
      TString path = ((TObjString*) fileTok->At(i))->GetString().Strip(TString::kBoth);
      if (path.IsNull()) continue;

      TFile* f = TFile::Open(path);
      if (!f || f->IsZombie()) { printf("skip: cannot open %s\n", path.Data()); continue; }
      TProfile* p = (TProfile*) f->Get("hMeanPt1VsM");
      if (!p) { printf("skip: no hMeanPt1VsM in %s\n", path.Data()); f->Close(); continue; }

      // Copy the non-empty bins into a graph so empty bins are not drawn as 0.
      TGraphErrors* g = new TGraphErrors();
      for (int b = 1; b <= p->GetNbinsX(); ++b) {
         if (p->GetBinEntries(b) <= 0) continue;
         int n = g->GetN();
         g->SetPoint(n, p->GetBinCenter(b), p->GetBinContent(b));
         g->SetPointError(n, 0., p->GetBinError(b));
      }
      double meanAll = p->GetMean(2), meanAllErr = p->GetMeanError(2);
      f->Close();

      TString label;
      if (labelTok->GetEntries() > nDrawn) {
         label = ((TObjString*) labelTok->At(nDrawn))->GetString().Strip(TString::kBoth);
      } else {
         TString run  = meanptRunFolderOf(path);
         TString velT = meanptReadVelT(run);
         label = run.BeginsWith("hrg") || velT.IsNull()
                     ? Form("HRG (%s)", run.Data())
                     : Form("blastwave, v_{T}=%s (%s)", velT.Data(), run.Data());
      }
      label += Form(": #LTp_{T,1}#GT = %.4f #pm %.4f", meanAll, meanAllErr);

      int col = colors[nDrawn % nColors];
      g->SetLineColor(col);
      g->SetMarkerColor(col);
      g->SetMarkerStyle(markers[nDrawn % nColors]);
      g->SetMarkerSize(0.9);
      mg->Add(g, "pe");
      leg->AddEntry(g, label, "pe");
      ++nDrawn;
   }

   if (nDrawn == 0) { printf("error: nothing to plot\n"); return; }

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

   TCanvas* c1 = new TCanvas("c1", "<pT,1> vs M comparison", 700, 550);
   mg->Draw("a");
   mg->GetXaxis()->SetTitle("M (GeV)");
   mg->GetYaxis()->SetTitle("#LTp_{T,1}#GT (GeV)");
   mg->GetXaxis()->SetTitleSize(0.05);
   mg->GetYaxis()->SetTitleSize(0.05);
   mg->GetYaxis()->SetTitleOffset(1.3);
   mg->SetMinimum(0.);
   mg->SetMaximum(1.6 * mg->GetHistogram()->GetMaximum());
   leg->Draw();
   if (title && title[0]) {
      TLatex tl;
      tl.SetNDC();
      tl.SetTextSize(0.045);
      tl.DrawLatex(0.15, 0.93, title);
   }
   c1->Print(output);
   TString png = output;
   if (png.EndsWith(".pdf")) {
      png.Replace(png.Length() - 4, 4, ".png");
      c1->Print(png);
   }
}

void ex_plot_compare_meanpt_all(const char* runs = "hrg_run,bw_run,bw_run2,bw_run3",
                                const char* events = "500") {
   TObjArray* runTok = TString(runs).Tokenize(",");

   struct Kind { const char* tag; const char* file; const char* title; const char* label; };
   const Kind kinds[] = {
      {"omega_pipm", "mass_pt_meanpt",    "#omega #rightarrow #pi^{+}#pi^{-}", "#omega: #pi^{+}#pi^{-}"},
      {"full_pipm",  "mass_pt_full_pipm", "full (same ancestor): #pi^{+}#pi^{-}", "full: #pi^{+}#pi^{-}"},
      {"full_pipp",  "mass_pt_full_pipp", "full (same ancestor): #pi^{+}#pi^{+}", "full: #pi^{+}#pi^{+}"},
      {"full_pimm",  "mass_pt_full_pimm", "full (same ancestor): #pi^{-}#pi^{-}", "full: #pi^{-}#pi^{-}"},
   };

   // Model comparison, one plot per pair type.
   for (const Kind& k : kinds) {
      TString files;
      for (int i = 0; i < runTok->GetEntries(); ++i) {
         TString run = ((TObjString*) runTok->At(i))->GetString().Strip(TString::kBoth);
         if (run.IsNull()) continue;
         TString path = Form("%s/%s_%s.root", run.Data(), k.file, events);
         if (gSystem->AccessPathName(path)) { printf("skip: %s not found\n", path.Data()); continue; }
         if (!files.IsNull()) files += ",";
         files += path;
      }
      if (files.IsNull()) { printf("note: no %s files — skipping\n", k.tag); continue; }
      ex_plot_compare_meanpt(files, "", Form("meanpt_compare_%s.pdf", k.tag), k.title);
   }

   // Charge-combination comparison, one plot per run (full +-, ++, --, then omega +-).
   for (int i = 0; i < runTok->GetEntries(); ++i) {
      TString run = ((TObjString*) runTok->At(i))->GetString().Strip(TString::kBoth);
      if (run.IsNull()) continue;
      TString files, labels;
      for (int j : {1, 2, 3, 0}) {
         TString path = Form("%s/%s_%s.root", run.Data(), kinds[j].file, events);
         if (gSystem->AccessPathName(path)) continue;
         if (!files.IsNull()) { files += ","; labels += ","; }
         files  += path;
         labels += kinds[j].label;
      }
      if (files.IsNull()) { printf("note: no mean-pT files in %s — skipping\n", run.Data()); continue; }
      TString velT  = meanptReadVelT(run);
      TString title = run.BeginsWith("hrg") || velT.IsNull()
                          ? Form("HRG (%s)", run.Data())
                          : Form("blastwave, v_{T}=%s (%s)", velT.Data(), run.Data());
      ex_plot_compare_meanpt(files, labels, run + "/meanpt_full_charges.pdf", title);
   }
}
