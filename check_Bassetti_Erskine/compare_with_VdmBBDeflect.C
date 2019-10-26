// Cross checks the consistency of the electric fields calculated by
// 1) E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0() from E_field.hh
// and 2) the relevant part of "VdmBBDeflect" class functions (copy-pasted
// to VdmBBDeflect_Bassetti_Erskine_part.C)
// for N points (eg. 100000, configurable) randomly distributed in the ranges
// specified in "compare_with_VdmBBDeflect_config.txt" together with fixed
// Gaussian bunch widths.
//
// Compile it with "make compare_with_VdmBBDeflect"
// and run as "compare_with_VdmBBDeflect compare_with_VdmBBDeflect_config.txt"

// Author V. Balagura, balagura@cern.ch (Oct 2019)

#include <E_field.hh>
#include <config.hh>
#include <fstream>

using namespace std;
struct Gen {
  double min, max;
  Gen(double fMin, double fMax) : min(fMin), max(fMax) {}
  double operator()() { return min + (max - min) * drand48(); }
};

std::pair<double, double> GetExEy(double sepx, double sepy,
				  double sigmax, double sigmay);

void calc_max(double* maximum, double* rel_maximum, double value1, double value2) {
  double d = abs(value1 - value2);
  double d_rel = d * 2 / abs(value1 + value2);
  *maximum = max(*maximum, d);
  *rel_maximum = max(*rel_maximum, d_rel);
}

int main(int argc, char** argv) {
  string config_file;
  if (argc != 2) {
    cout << "Usage: <program name> <configuration file>\n";
    return 0;
  }
  ifstream conf_file(argv[1]);
  Config c(conf_file);
  if (c.defined("seed")) {
    srand48(c("seed"));
  } else {
    srand48(time(NULL));
  }
  Gen
    x_gen(c["xmin"], c["xmax"]),
    y_gen(c["ymin"], c["ymax"]),
    sigx_gen(c["sigx.min"], c["sigx.max"]),
    sigy_gen(c["sigy.min"], c["sigy.max"]);
  int n = c("n");
  int print_E = c.defined("n.times.print.E") ? c("n.times.print.E") : 0;
  double max_discrepancy = 0, max_rel_discrepancy = 0;
  for (int i=0; i<n; ++i) {
    double x = x_gen();
    double y = y_gen();
    double sigx = sigx_gen();
    double sigy = sigy_gen();
    if ( abs(sigy/sigx - 1) < 1e-8 ) continue; // avoid round beam case
    complex<double> e1 = E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y, sigx, sigy);
    //
    std::pair<double, double> e2 = GetExEy(x, y, sigx, sigy);
    calc_max(&max_discrepancy, &max_rel_discrepancy, real(e1), e2.first);
    calc_max(&max_discrepancy, &max_rel_discrepancy, imag(e1), e2.second);
    if (i < print_E) cout << "x,y = (" << x << ", " << y
			  << "), sigx,sigy = (" << sigx << ", " << sigy
			  << "), field 1 = (" << real(e1) << ", " << imag(e1)
			  << "), field 2 = (" << e2.first << ", " << e2.second << ")" << endl;
  }
  cout << "\n"
       << "Maximal absolute mismatch in field components = " << max_discrepancy << "\n"
       << "Maximal relative mismatch in field components = " << max_rel_discrepancy << endl;
  return 0;
}
