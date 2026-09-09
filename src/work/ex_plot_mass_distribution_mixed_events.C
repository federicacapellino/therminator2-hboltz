// ex_plot_mass_distribution_mixed_events.C — overlay hMass, hMassME and the
// connected difference |hMass - hMassME| from one ex_mass_distribution_mixed_events
// output file on a single canvas.
//
// hMass is the same-event correlated distribution (pairs per pid1); hMassME is
// the mixed-event, same-rooteid baseline (pairs per pid1).  Both are drawn raw on
// a log-y axis.  hMassConn = hMass - hMassME is the connected correlation; it can
// go negative, so |hMassConn| is plotted for the log axis.
//
// The output PDF is named mass_me_<model>_<events>.pdf from the parameters
// stored in the input file (falls back to mass_distribution_mixed_events.pdf).
//
// Usage (from repository root):
//   root -b -q src/work/ex_plot_mass_distribution_mixed_events.C
//   root -b -q 'src/work/ex_plot_mass_distribution_mixed_events.C("mass_me_hrg_500.root")'

#include <algorithm>
#include <cmath>

#include <TNamed.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TSystem.h>

void ex_plot_mass_distribution_mixed_events(const char* filename = "mass_me_hrg_500.root") {
   TFile* f = TFile::Open(filename);
   if (!f || f->IsZombie()) { printf("Cannot open %s\n", filename); return; }

   TH1D* hMass   = (TH1D*) f->Get("hMass");
   TH1D* hMassME = (TH1D*) f->Get("hMassME");
   if (!hMass)   { printf("Cannot find hMass in %s\n", filename);   return; }
   if (!hMassME) { printf("Cannot find hMassME in %s\n", filename); return; }

   // Generation parameters written by ex_mass_distribution_mixed_events (optional).
   TNamed*          pModel = (TNamed*)          f->Get("model");
   TParameter<int>* pNev   = (TParameter<int>*) f->Get("nEvents");

   // Connected correlation hMassConn = hMass - hMassME.  It can be negative, so
   // plot |hMassConn| on the log axis.
   TH1D* hConnAbs = (TH1D*) hMass->Clone("hMassConnAbs");
   hConnAbs->SetDirectory(nullptr);
   hConnAbs->Add(hMassME, -1.0);
   for (int i = 0; i <= hConnAbs->GetNbinsX() + 1; ++i)
      hConnAbs->SetBinContent(i, std::fabs(hConnAbs->GetBinContent(i)));

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

   hMass->SetLineColor(kBlue + 1);    hMass->SetLineWidth(2);
   hMassME->SetLineColor(kRed + 1);   hMassME->SetLineWidth(2);
   hConnAbs->SetLineColor(kGreen + 2); hConnAbs->SetLineWidth(2);

   // Y range spanning all curves on a log scale (positive bin contents only).
   double ymax = std::max({hMass->GetMaximum(), hMassME->GetMaximum(), hConnAbs->GetMaximum()});
   double ymin = 1e30;
   for (TH1D* h : {hMass, hMassME, hConnAbs})
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
         double c = h->GetBinContent(i);
         if (c > 0. && c < ymin) ymin = c;
      }
   if (ymin >= 1e30) ymin = 1e-6;

   TCanvas* c1 = new TCanvas("c1", "Same-event vs mixed-event invariant mass", 600, 500);
   c1->SetLogy();

   hMass->GetXaxis()->SetTitle("M (GeV)");
   hMass->GetXaxis()->SetTitleSize(0.05);
   hMass->GetXaxis()->SetTitleOffset(0.9);
   hMass->GetXaxis()->SetLabelSize(0.045);
   hMass->GetYaxis()->SetTitle("dN_{pairs}/dM  (as stored)");
   hMass->GetYaxis()->SetTitleSize(0.05);
   hMass->GetYaxis()->SetTitleOffset(1.4);
   hMass->GetYaxis()->SetLabelSize(0.045);
   hMass->SetMinimum(0.5 * ymin);
   hMass->SetMaximum(5.0 * ymax);

   hMass->Draw("hist e");
   hMassME->Draw("hist e same");
   hConnAbs->Draw("hist e same");

   TLegend* leg = new TLegend(0.45, 0.70, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hMass,    "hMass (same event)",           "l");
   leg->AddEntry(hMassME,  "hMassME (mixed, same rooteid)", "l");
   leg->AddEntry(hConnAbs, "|hMass - hMassME|",             "l");
   leg->Draw();

   // Generation parameters box (top-left, opposite the legend).
   TPaveText* pt = new TPaveText(0.18, 0.74, 0.46, 0.88, "NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.035);
   pt->AddText(Form("model: %s", pModel ? pModel->GetTitle() : "?"));
   if (pNev) pt->AddText(Form("events: %d", pNev->GetVal()));
   pt->Draw();

   // Write the PDF next to the input ROOT file (same folder, e.g. mass_me/).
   TString outDir = gSystem->GetDirName(filename);
   TString pdfName = outDir + "/mass_distribution_mixed_events.pdf";
   if (pModel && pNev)
      pdfName = TString::Format("%s/mass_me_%s_%d.pdf", outDir.Data(),
                                pModel->GetTitle(), pNev->GetVal());
   gPad->Print(pdfName);
}
