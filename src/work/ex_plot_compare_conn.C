// ex_plot_compare_conn.C — overlay hMassConn (or, via the 4th argument, any other
// histogram such as the uncut hMass; see ex_plot_compare_full.C) from several runs
// on one canvas.
//
// Works for both ex_mass_distribution_connected (mass_con_*) and
// ex_phi_distribution_connected (phi_con_*) outputs; the kind is detected from the
// file name prefix and all files must be of the same kind.  Mass is drawn as
// |hMassConn| on a log axis, Delta phi signed on a linear axis.
//
// Labels default to "<model>[, v_{T}=<VelT>] (<run folder>)", with VelT read from
// <run folder>/blastwave.ini when it exists.
//
// Usage (from build/work, where the run folders live):
//   root -b -q ex_plot_compare_conn.C
//   root -b -q 'ex_plot_compare_conn.C("bw_run/phi_con/phi_con_blastwave_3_1.5_500.root,hrg_run/phi_con/phi_con_hrg_3_1.5_500.root")'
//   root -b -q 'ex_plot_compare_conn.C("a.root,b.root", "label A,label B", "out.pdf")'
//   root -b -q 'ex_plot_compare_conn.C("a.root,b.root", "", "out.pdf", "hMass")'

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include <TLine.h>
#include <TNamed.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TSystem.h>

// First path component of a file, e.g. "bw_run2" for "bw_run2/mass_con/x.root".
static TString runFolderOf(const TString& path) {
   Ssiz_t slash = path.Index("/");
   return (slash > 0) ? TString(path(0, slash)) : TString(".");
}

// VelT from <runFolder>/blastwave.ini, or "" if not found.
static TString readVelT(const TString& runFolder) {
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

void ex_plot_compare_conn(
    const char* files  = "bw_run/mass_con/mass_con_blastwave_3_1.5_500.root,"
                         "bw_run2/mass_con/mass_con_blastwave_3_1.5_500.root,"
                         "bw_run3/mass_con/mass_con_blastwave_3_1.5_500.root,"
                         "hrg_run/mass_con/mass_con_hrg_3_1.5_10000.root",
    const char* labels = "",
    const char* output = "",
    const char* hist   = "hMassConn") {
   const bool isFull = TString(hist) == "hMass";       // uncut distribution
   const bool isWin  = TString(hist).EndsWith("Win");  // per pid1 inside |y| < Ymax/2
   const char* n1    = isWin ? "N_{1}^{win}" : "N_{1}";
   TObjArray* fileTok  = TString(files).Tokenize(",");
   TObjArray* labelTok = TString(labels).Tokenize(",");

   const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 1,
                         kOrange + 7, kCyan + 2, kBlack, kViolet - 3};
   const int nColors = sizeof(colors) / sizeof(colors[0]);

   std::vector<TH1D*>   hists;
   std::vector<TString> names;
   int kind = -1;  // 0 = mass, 1 = Delta phi
   double refYmax = NAN, refDeltaY = NAN;

   for (int i = 0; i < fileTok->GetEntries(); ++i) {
      TString path = ((TObjString*) fileTok->At(i))->GetString().Strip(TString::kBoth);
      if (path.IsNull()) continue;

      TString base = gSystem->BaseName(path);
      int thisKind = 0;
      if (base.BeginsWith("phi_con_"))       thisKind = 1;
      else if (!base.BeginsWith("mass_con_"))
         printf("warning: %s has an unknown prefix, treating it as mass\n", path.Data());
      if (kind < 0) kind = thisKind;
      if (thisKind != kind) {
         printf("error: cannot mix mass_con and phi_con files (%s)\n", path.Data());
         return;
      }

      TFile* f = TFile::Open(path);
      if (!f || f->IsZombie()) { printf("skip: cannot open %s\n", path.Data()); continue; }
      TH1D* h = (TH1D*) f->Get(hist);
      if (!h) { printf("skip: no %s in %s\n", hist, path.Data()); f->Close(); continue; }
      h = (TH1D*) h->Clone(Form("hConn_%d", i));
      h->SetDirectory(nullptr);

      TNamed*             pModel  = (TNamed*)             f->Get("model");
      TParameter<double>* pYmax   = (TParameter<double>*) f->Get("Ymax");
      TParameter<double>* pDeltaY = (TParameter<double>*) f->Get("DeltaY");
      if (pYmax && pDeltaY) {
         if (std::isnan(refYmax)) {
            refYmax = pYmax->GetVal();
            refDeltaY = pDeltaY->GetVal();
         } else if (pYmax->GetVal() != refYmax || pDeltaY->GetVal() != refDeltaY) {
            printf("warning: %s has Ymax=%g, DeltaY=%g (first file: %g, %g)\n", path.Data(),
                   pYmax->GetVal(), pDeltaY->GetVal(), refYmax, refDeltaY);
         }
      }

      TString label;
      if (labelTok->GetEntries() > (int) hists.size()) {
         label = ((TObjString*) labelTok->At(hists.size()))->GetString().Strip(TString::kBoth);
      } else {
         TString run   = runFolderOf(path);
         TString model = pModel ? pModel->GetTitle() : "?";
         TString velT  = (model == "blastwave") ? readVelT(run) : TString("");
         label = velT.IsNull() ? Form("%s (%s)", model.Data(), run.Data())
                               : Form("%s, v_{T}=%s (%s)", model.Data(), velT.Data(), run.Data());
      }
      f->Close();

      int col = colors[hists.size() % nColors];
      h->SetLineColor(col);
      h->SetMarkerColor(col);
      h->SetLineWidth(2);
      hists.push_back(h);
      names.push_back(label);
   }

   if (hists.empty()) { printf("error: nothing to plot\n"); return; }
   bool isPhi = (kind == 1);

   // Mass: |hMassConn| on a log axis (the difference can be negative).
   if (!isPhi && !isFull)
      for (TH1D* h : hists)
         for (int b = 0; b <= h->GetNbinsX() + 1; ++b)
            h->SetBinContent(b, std::fabs(h->GetBinContent(b)));

   double ymax = -1e30, ymin = isPhi ? 0. : 1e30;
   for (TH1D* h : hists)
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
         double c = h->GetBinContent(b), e = h->GetBinError(b);
         ymax = std::max(ymax, c + e);
         if (isPhi)            ymin = std::min(ymin, c - e);
         else if (c > 0.)      ymin = std::min(ymin, c);
      }
   if (!isPhi && ymin >= 1e30) ymin = 1e-6;

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

   TCanvas* c1 = new TCanvas("c1", Form("%s comparison", hist), 600, 500);
   if (!isPhi) c1->SetLogy();

   TH1D* h0 = hists[0];
   h0->GetXaxis()->SetTitle(isPhi ? "#Delta#phi (rad)" : "M (GeV)");
   h0->GetXaxis()->SetTitleSize(0.05);
   h0->GetXaxis()->SetTitleOffset(0.9);
   h0->GetXaxis()->SetLabelSize(0.045);
   if (isFull)
      h0->GetYaxis()->SetTitle(isPhi ? "(1/N_{1}) dN_{pairs}/d#Delta#phi"
                                     : "(1/N_{1}) dN_{pairs}/dM");
   else
      h0->GetYaxis()->SetTitle(isPhi ? Form("(1/%s) dN^{conn}_{pairs}/d#Delta#phi", n1)
                                     : Form("|(1/%s) dN^{conn}_{pairs}/dM|", n1));
   h0->GetYaxis()->SetTitleSize(0.05);
   h0->GetYaxis()->SetTitleOffset(1.4);
   h0->GetYaxis()->SetLabelSize(0.045);
   if (isPhi) {
      h0->SetMinimum(ymin < 0. ? 1.1 * ymin : 0.);
      h0->SetMaximum(1.5 * ymax);
   } else {
      h0->SetMinimum(0.5 * ymin);
      h0->SetMaximum(20.0 * ymax);
   }

   h0->Draw("hist e");
   for (size_t i = 1; i < hists.size(); ++i) hists[i]->Draw("hist e same");

   if (isPhi) {
      TLine* zero = new TLine(h0->GetXaxis()->GetXmin(), 0., h0->GetXaxis()->GetXmax(), 0.);
      zero->SetLineStyle(2);
      zero->SetLineColor(kGray + 2);
      zero->Draw();
   }

   TLegend* leg = new TLegend(0.45, 0.88 - 0.05 * hists.size(), 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->SetTextSize(0.033);
   for (size_t i = 0; i < hists.size(); ++i) leg->AddEntry(hists[i], names[i], "l");
   leg->Draw();

   TPaveText* pt = new TPaveText(0.18, 0.76, 0.42, 0.88, "NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.035);
   pt->AddText(isFull ? "hMass (no cuts)" : hist);
   if (!isFull && !std::isnan(refYmax)) {
      pt->AddText(Form("Y_{max} = %g", refYmax));
      pt->AddText(Form("#DeltaY = %g", refDeltaY));
   }
   pt->Draw();

   TString pdfName = output;
   if (pdfName.IsNull())
      pdfName = Form("%s_compare_%s.pdf", isFull ? "full" : (isWin ? "connwin" : "conn"),
                     isPhi ? "phi" : "mass");
   gPad->Print(pdfName);
}
