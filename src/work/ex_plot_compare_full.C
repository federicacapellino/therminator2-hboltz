// ex_plot_compare_full.C — overlay the full (uncut) hMass from several runs on one
// canvas.  Same inputs, labels and mass / Delta-phi detection as
// ex_plot_compare_conn.C, which does the actual work.
//
// Usage (from build/work, where the run folders live):
//   root -b -q ex_plot_compare_full.C
//   root -b -q 'ex_plot_compare_full.C("bw_run/phi_con/phi_con_blastwave_3_1.5_500.root,hrg_run/phi_con/phi_con_hrg_3_1.5_500.root")'

#include "ex_plot_compare_conn.C"

void ex_plot_compare_full(
    const char* files  = "bw_run/mass_con/mass_con_blastwave_3_1.5_500.root,"
                         "bw_run2/mass_con/mass_con_blastwave_3_1.5_500.root,"
                         "bw_run3/mass_con/mass_con_blastwave_3_1.5_500.root,"
                         "hrg_run/mass_con/mass_con_hrg_3_1.5_10000.root",
    const char* labels = "",
    const char* output = "") {
   ex_plot_compare_conn(files, labels, output, "hMass");
}
