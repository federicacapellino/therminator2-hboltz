// ex_plot_mass_pt_v2.C — mass-pT densities and pair v2 from ex_mass_pt_distribution, for the
// four correlated pion-pair channels of one run folder:
//   omega +-  <- <run>/mass_pt_meanpt_<N>.root     (--parent 223 --pid1 211 --pid2 -211)
//   full  +-  <- <run>/mass_pt_full_pipm_<N>.root  (--parent 0   --pid1 211 --pid2 -211)
//   full  ++  <- <run>/mass_pt_full_pipp_<N>.root  (--parent 0   --pid1 211 --pid2 211)
//   full  --  <- <run>/mass_pt_full_pimm_<N>.root  (--parent 0   --pid1 -211 --pid2 -211)
//
// Writes into the run folder:
//   masspt_density.pdf  hMassPt, (1/N_1) dN_pairs per (pT1, M) bin, log z, one pad per channel
//   v2_pt_m.pdf         hV2PtM, pair v2 = <cos 2 phi_pair> per (pT1, M) bin; bins with fewer
//                       than kMinPairs pairs are left blank, z range symmetric around 0
//   v2_vs_m.pdf         hV2VsM, pair v2 vs M integrated over pT1, all channels overlaid;
//                       M bins with fewer than kMinPairs pairs are not drawn
// Missing files are skipped.
//
// Usage (from build/work):
//   root -b -q 'ex_plot_mass_pt_v2.C("bwaniso_run", "500")'

#include <cmath>
#include <vector>

#include <TCanvas.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TStyle.h>

static const int kMinPairs = 20;

struct V2Channel {
   TString file, label;
   int color, marker;
};

// Copy of a TProfile2D as a TH2D, with bins under kMinPairs entries pushed below the z range so
// COLZ leaves them blank.  Returns the largest |v2| among the kept bins in zmax.
static TH2D* v2MaskedCopy(TProfile2D* p, const TString& name, double& zmax) {
   TH2D* h = p->ProjectionXY(name, "B");
   h->Reset();
   for (int ix = 1; ix <= p->GetNbinsX(); ++ix) {
      for (int iy = 1; iy <= p->GetNbinsY(); ++iy) {
         int bin = p->GetBin(ix, iy);
         if (p->GetBinEntries(bin) < kMinPairs) {
            h->SetBinContent(ix, iy, -1e9);
         } else {
            double v = p->GetBinContent(bin);
            h->SetBinContent(ix, iy, v);
            zmax = std::max(zmax, std::fabs(v));
         }
      }
   }
   h->SetTitle(p->GetTitle());
   h->GetXaxis()->SetTitle(p->GetXaxis()->GetTitle());
   h->GetYaxis()->SetTitle(p->GetYaxis()->GetTitle());
   h->GetZaxis()->SetTitle("v_{2}^{pair}");
   return h;
}

// TProfile -> TGraphErrors without bins under kMinPairs pairs.
static TGraphErrors* v2Graph(TProfile* p) {
   TGraphErrors* g = new TGraphErrors();
   for (int i = 1; i <= p->GetNbinsX(); ++i) {
      if (p->GetBinEntries(i) < kMinPairs) continue;
      int n = g->GetN();
      g->SetPoint(n, p->GetBinCenter(i), p->GetBinContent(i));
      g->SetPointError(n, 0., p->GetBinError(i));
   }
   return g;
}

void ex_plot_mass_pt_v2(const char* runFolder = "bwaniso_run", const char* nEvents = "500") {
   gROOT->SetBatch(kTRUE);
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   TH1::AddDirectory(kFALSE);

   const TString run(runFolder), n(nEvents);
   std::vector<V2Channel> channels = {
      {run + "/mass_pt_meanpt_" + n + ".root", "#omega #rightarrow #pi^{+}#pi^{-}", kBlack, 20},
      {run + "/mass_pt_full_pipm_" + n + ".root", "full HRG #pi^{+}#pi^{-}", kRed + 1, 21},
      {run + "/mass_pt_full_pipp_" + n + ".root", "full HRG #pi^{+}#pi^{+}", kBlue + 1, 22},
      {run + "/mass_pt_full_pimm_" + n + ".root", "full HRG #pi^{-}#pi^{-}", kGreen + 2, 23},
   };

   std::vector<TH2D*> dens(channels.size(), nullptr);
   std::vector<TProfile2D*> v2pm(channels.size(), nullptr);
   std::vector<TProfile*> v2m(channels.size(), nullptr);
   for (size_t i = 0; i < channels.size(); ++i) {
      TFile* f = TFile::Open(channels[i].file);
      if (!f || f->IsZombie()) {
         printf("skip: %s not found\n", channels[i].file.Data());
         continue;
      }
      if (auto* h = dynamic_cast<TH2D*>(f->Get("hMassPt"))) dens[i] = (TH2D*)h->Clone();
      if (auto* h = dynamic_cast<TProfile2D*>(f->Get("hV2PtM"))) v2pm[i] = (TProfile2D*)h->Clone();
      if (auto* h = dynamic_cast<TProfile*>(f->Get("hV2VsM"))) v2m[i] = (TProfile*)h->Clone();
      if (!v2m[i]) printf("warning: %s has no v2 histograms\n", channels[i].file.Data());
      else
         printf("%-28s pair v2 (pT1, M integrated) = %.4f +- %.4f\n", channels[i].file.Data(),
                v2m[i]->GetMean(2), v2m[i]->GetMeanError(2));
      f->Close();
   }

   TLatex lab;
   lab.SetNDC();
   lab.SetTextSize(0.05);

   // 1. mass-pT densities
   TCanvas c1("c1", "", 1200, 1000);
   c1.Divide(2, 2);
   for (size_t i = 0; i < channels.size(); ++i) {
      if (!dens[i]) continue;
      TVirtualPad* pad = c1.cd(i + 1);
      pad->SetLogz();
      pad->SetRightMargin(0.19);
      dens[i]->SetTitle("");
      dens[i]->GetZaxis()->SetTitle("(1/N_{1}) dN_{pairs}");
      dens[i]->GetZaxis()->SetTitleOffset(1.6);
      dens[i]->Draw("COLZ");
      lab.DrawLatex(0.15, 0.92, channels[i].label);
   }
   c1.SaveAs(run + "/masspt_density.pdf");

   // 2. pair v2 in (pT1, M)
   TCanvas c2("c2", "", 1200, 1000);
   c2.Divide(2, 2);
   std::vector<TH2D*> v2h(channels.size(), nullptr);
   for (size_t i = 0; i < channels.size(); ++i) {
      if (!v2pm[i]) continue;
      double zmax = 0.;
      v2h[i] = v2MaskedCopy(v2pm[i], Form("v2masked_%zu", i), zmax);
      if (zmax <= 0.) zmax = 0.1;
      TVirtualPad* pad = c2.cd(i + 1);
      pad->SetRightMargin(0.19);
      v2h[i]->SetTitle("");
      v2h[i]->SetMinimum(-zmax);
      v2h[i]->SetMaximum(zmax);
      v2h[i]->GetZaxis()->SetTitleOffset(1.3);
      v2h[i]->Draw("COLZ");
      lab.DrawLatex(0.15, 0.92, channels[i].label + Form("  (bins #geq %d pairs)", kMinPairs));
   }
   c2.SaveAs(run + "/v2_pt_m.pdf");

   // 3. pair v2 vs M, integrated over pT1
   TCanvas c3("c3", "", 900, 700);
   c3.SetLeftMargin(0.13);
   TMultiGraph* mg = new TMultiGraph();
   TLegend* leg = new TLegend(0.15, 0.68, 0.48, 0.88);
   leg->SetBorderSize(0);
   for (size_t i = 0; i < channels.size(); ++i) {
      if (!v2m[i]) continue;
      TGraphErrors* g = v2Graph(v2m[i]);
      g->SetLineColor(channels[i].color);
      g->SetMarkerColor(channels[i].color);
      g->SetMarkerStyle(channels[i].marker);
      g->SetMarkerSize(0.7);
      mg->Add(g, "PE");
      leg->AddEntry(g, channels[i].label, "pe");
   }
   mg->SetTitle(";M_{inv} (GeV);v_{2}^{pair} = #LTcos 2#phi_{pair}#GT");
   mg->Draw("A");
   mg->GetXaxis()->SetLimits(0., 2.);
   c3.Update();
   TLine zero(0., 0., 2., 0.);
   zero.SetLineStyle(2);
   zero.Draw();
   TLine omega(0.7827, c3.GetUymin(), 0.7827, c3.GetUymax());
   omega.SetLineStyle(3);
   omega.SetLineColor(kGray + 1);
   omega.Draw();
   leg->Draw();
   c3.SaveAs(run + "/v2_vs_m.pdf");
}
