// Diagnostic for the M9 report: T9.2 pull width versus bin entries, compared with the Student-t
// expectation sqrt(N/(N-3)) for pulls with an estimated error (like-sign: N/2 independent entries).
// Run: root -l -b -q 'width_vs_entries.C("<merged harmonics file>")'
// T9.2-type pulls (m != 0, non-mirror pairs) pooled in classes of bin entries, per cumulative level.
void width_vs_entries(const char* file) {
  TFile f(file);
  const char* tags[] = {"pippip","pippim","pimpim","KpKp","KpKm","KmKm","pp","ppbar","pbarpbar"};
  const char* ex[] = {"L0","L1x","L2x"};
  double edges[] = {50,100,200,500,1000,1e4,1e5,1e12};
  const int nc = 7;
  for (int lev = 0; lev < 3; lev++) {
    std::vector<double> s1(nc), s2(nc), n(nc), predW(nc);
    for (int ip = 0; ip < 9; ip++) {
      bool like = (ip==0||ip==2||ip==3||ip==5||ip==6||ip==8);
      for (int l = 1; l <= 4; l++) for (int m = 1; m <= l; m++) for (const char* cs : {"c","s"}) {
        if (like && l % 2) continue;
        TProfile2D* h = nullptr;
        for (int k = 0; k <= lev; k++) {
          auto* g = f.Get<TProfile2D>(Form("hY%d%d%s_%s_%s", l, m, cs, ex[k], tags[ip]));
          if (!h) { h = (TProfile2D*)g->Clone(); h->SetDirectory(nullptr);} else h->Add(g);
        }
        double sc = (like && l % 2 == 0) ? sqrt(2.) : 1.;
        for (int b = 0; b < h->GetNcells(); b++) {
          double N = h->GetBinEntries(b);
          if (N < 50 || h->GetBinError(b) <= 0) continue;
          double p = h->GetBinContent(b) / (h->GetBinError(b) * sc);
          int c = 0; while (N >= edges[c+1]) c++;
          s1[c] += p; s2[c] += p*p; n[c]++;
          double Ne = like ? N/2 : N;   // independent entries
          predW[c] += sqrt(Ne/(Ne-1.) * (Ne-1.)/(Ne-3.));
        }
        delete h;
      }
    }
    for (int c = 0; c < nc; c++) if (n[c] > 0) {
      double mean = s1[c]/n[c], w = sqrt(s2[c]/n[c]-mean*mean);
      printf("L%d  N in [%6.0f,%6.0g)  bins %7.0f  mean %+.4f  width %.4f +- %.4f  t-pred %.4f\n",
             lev, edges[c], edges[c+1], n[c], mean, w, 1/sqrt(2*n[c]), predW[c]/n[c]);
    }
  }
}
