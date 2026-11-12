// ex_plot_phi_distribution_connected.C — overlay hMass, hMassRapCut and hMassConn
// (Delta phi = phi1 - phi2 on [-pi/2, 3pi/2)) from one ex_phi_distribution_connected
// output file.
//
// The output PDF is named phi_con_<model>_<ymax>_<deltay>_<events>.pdf from the
// parameters stored in the input file (falls back to phi_distribution_connected.pdf).
//
// Usage (from repository root):
//   root -b -q src/work/ex_plot_phi_distribution_connected.C
//   root -b -q 'src/work/ex_plot_phi_distribution_connected.C("phi_con/phi_con_blastwave_3_1.5_500.root")'

#include <algorithm>
#include <cmath>

#include <TLine.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TPaveText.h>
#include <TSystem.h>

void ex_plot_phi_distribution_connected(const char* filename = "phi_con/phi_con_hrg_2_1_500.root") {
   TFile* f = TFile::Open(filename);
   if (!f || f->IsZombie()) { printf("Cannot open %s\n", filename); return; }

   TH1D* hMass   = (TH1D*) f->Get("hMass");
   TH1D* hRapCut = (TH1D*) f->Get("hMassRapCut");
   TH1D* hConn   = (TH1D*) f->Get("hMassConn");
   if (!hMass)   { printf("Cannot find hMass in %s\n", filename);       return; }
   if (!hRapCut) { printf("Cannot find hMassRapCut in %s\n", filename); return; }
   if (!hConn)   { printf("Cannot find hMassConn in %s\n", filename);   return; }

   // Generation parameters written by ex_phi_distribution_connected (optional).
   TNamed*             pModel  = (TNamed*)             f->Get("model");
   TParameter<int>*    pNev    = (TParameter<int>*)    f->Get("nEvents");
   TParameter<double>* pYmax   = (TParameter<double>*) f->Get("Ymax");
   TParameter<double>* pDeltaY = (TParameter<double>*) f->Get("DeltaY");

   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetOptStat(kFALSE);
   gStyle->SetOptTitle(kFALSE);

   hMass->SetLineColor(kBlue + 1);    hMass->SetLineWidth(2);
   hRapCut->SetLineColor(kRed + 1);   hRapCut->SetLineWidth(2);
   hConn->SetLineColor(kGreen + 2);   hConn->SetLineWidth(2);

   // Linear y range: hMassConn is a difference and may go negative, so keep it signed.
   double ymax = std::max({hMass->GetMaximum(), hRapCut->GetMaximum(), hConn->GetMaximum()});
   double ymin = 0.;
   for (int i = 1; i <= hConn->GetNbinsX(); ++i)
      ymin = std::min(ymin, hConn->GetBinContent(i) - hConn->GetBinError(i));

   TCanvas* c1 = new TCanvas("c1", "Connected pair Delta phi", 600, 500);

   TAxis* ax = hMass->GetXaxis();
   ax->SetTitle("#Delta#phi (rad)");
   ax->SetTitleSize(0.05);
   ax->SetTitleOffset(0.9);
   ax->SetLabelSize(0.045);
   hMass->GetYaxis()->SetTitle("(1/N_{1}) dN_{pairs}/d#Delta#phi");
   hMass->GetYaxis()->SetTitleSize(0.05);
   hMass->GetYaxis()->SetTitleOffset(1.4);
   hMass->GetYaxis()->SetLabelSize(0.045);
   hMass->SetMinimum(ymin < 0. ? 1.1 * ymin : 0.);
   hMass->SetMaximum(1.5 * ymax);

   hMass->Draw("hist e");
   hRapCut->Draw("hist e same");
   hConn->Draw("hist e same");

   TLine* zero = new TLine(ax->GetXmin(), 0., ax->GetXmax(), 0.);
   zero->SetLineStyle(2);
   zero->SetLineColor(kGray + 2);
   zero->Draw();

   TLegend* leg = new TLegend(0.45, 0.70, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hMass,   "hMass (correlated)",           "l");
   leg->AddEntry(hRapCut, "hMassRapCut (rand. rapidity)", "l");
   leg->AddEntry(hConn,   "hMassConn",                    "l");
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

   // Write the PDF next to the input ROOT file (same folder, e.g. phi_con/).
   TString outDir = gSystem->GetDirName(filename);
   TString pdfName = outDir + "/phi_distribution_connected.pdf";
   if (pModel && pYmax && pDeltaY && pNev)
      pdfName = TString::Format("%s/phi_con_%s_%g_%g_%d.pdf", outDir.Data(),
                                pModel->GetTitle(), pYmax->GetVal(), pDeltaY->GetVal(),
                                pNev->GetVal());
   gPad->Print(pdfName);
}
