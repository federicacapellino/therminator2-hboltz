// ex_plot_mass_distribution_connected.C — overlay hMass, hMassRapCut and |hMassConn|
// from one ex_mass_distribution_connected output file.
//
// The output PDF is named mass_con_<model>_<ymax>_<deltay>_<events>.pdf from the
// parameters stored in the input file (falls back to mass_distribution_connected.pdf).
//
// Usage (from repository root):
//   root -b -q src/work/ex_plot_mass_distribution_connected.C
//   root -b -q 'src/work/ex_plot_mass_distribution_connected.C("mass_con_blastwave_3_1.5_500.root")'

#include <algorithm>
#include <cmath>

#include <TNamed.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TSystem.h>

void ex_plot_mass_distribution_connected(const char* filename = "mass_con_hrg_2_1_500.root") {
   TFile* f = TFile::Open(filename);
   if (!f || f->IsZombie()) { printf("Cannot open %s\n", filename); return; }

   TH1D* hMass   = (TH1D*) f->Get("hMass");
   TH1D* hRapCut = (TH1D*) f->Get("hMassRapCut");
   TH1D* hConn   = (TH1D*) f->Get("hMassConn");
   if (!hMass)   { printf("Cannot find hMass in %s\n", filename);       return; }
   if (!hRapCut) { printf("Cannot find hMassRapCut in %s\n", filename); return; }
   if (!hConn)   { printf("Cannot find hMassConn in %s\n", filename);   return; }

   // Generation parameters written by ex_mass_distribution_connected (optional).
   TNamed*             pModel  = (TNamed*)             f->Get("model");
   TParameter<int>*    pNev    = (TParameter<int>*)    f->Get("nEvents");
   TParameter<double>* pYmax   = (TParameter<double>*) f->Get("Ymax");
   TParameter<double>* pDeltaY = (TParameter<double>*) f->Get("DeltaY");

   // hMassConn is a difference (can be negative) -> plot |hMassConn| for the log axis.
   TH1D* hConnAbs = (TH1D*) hConn->Clone("hMassConnAbs");
   hConnAbs->SetDirectory(nullptr);
   for (int i = 0; i <= hConnAbs->GetNbinsX() + 1; ++i)
      hConnAbs->SetBinContent(i, std::fabs(hConnAbs->GetBinContent(i)));

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

   hMass->SetLineColor(kBlue + 1);    hMass->SetLineWidth(2);
   hRapCut->SetLineColor(kRed + 1);   hRapCut->SetLineWidth(2);
   hConnAbs->SetLineColor(kGreen + 2); hConnAbs->SetLineWidth(2);

   // Y range spanning all three on a log scale (positive bin contents only).
   double ymax = std::max({hMass->GetMaximum(), hRapCut->GetMaximum(), hConnAbs->GetMaximum()});
   double ymin = 1e30;
   for (TH1D* h : {hMass, hRapCut, hConnAbs})
      for (int i = 1; i <= h->GetNbinsX(); ++i) {
         double c = h->GetBinContent(i);
         if (c > 0. && c < ymin) ymin = c;
      }
   if (ymin >= 1e30) ymin = 1e-6;

   TCanvas* c1 = new TCanvas("c1", "Connected pair invariant mass", 600, 500);
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
   hRapCut->Draw("hist e same");
   hConnAbs->Draw("hist e same");

   TLegend* leg = new TLegend(0.45, 0.70, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hMass,    "hMass (correlated)",           "l");
   leg->AddEntry(hRapCut,  "hMassRapCut (rand. rapidity)", "l");
   leg->AddEntry(hConnAbs, "|hMassConn|",                  "l");
   leg->Draw();

   // Generation parameters box (top-left, opposite the legend).
   TPaveText* pt = new TPaveText(0.18, 0.68, 0.46, 0.88, "NDC");
   pt->SetBorderSize(0);
   pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->SetTextSize(0.035);
   pt->AddText(Form("model: %s", pModel ? pModel->GetTitle() : "?"));
   if (pNev)    pt->AddText(Form("events: %d", pNev->GetVal()));
   if (pYmax)   pt->AddText(Form("Y_{max} = %g", pYmax->GetVal()));
   if (pDeltaY) pt->AddText(Form("#DeltaY = %g", pDeltaY->GetVal()));
   pt->Draw();

   // Write the PDF next to the input ROOT file (same folder, e.g. mass_con/).
   TString outDir = gSystem->GetDirName(filename);
   TString pdfName = outDir + "/mass_distribution_connected.pdf";
   if (pModel && pYmax && pDeltaY && pNev)
      pdfName = TString::Format("%s/mass_con_%s_%g_%g_%d.pdf", outDir.Data(),
                                pModel->GetTitle(), pYmax->GetVal(), pDeltaY->GetVal(),
                                pNev->GetVal());
   gPad->Print(pdfName);
}
