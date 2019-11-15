#ifndef multi_gaussian_hh 
#define multi_gaussian_hh 1

#include <vector>
#include <array>
#include <string>
#include <complex>
#include <iostream>
#include <bilinear_interpolator.hh>

using namespace std;

struct Mutli_XY_Gaussian_bunches {
  // In the following "coor" argument is 0 for X, and 1 for Y;
  // "ip = 0, 1, 2, ..." is the sequantial number of the interaction point,
  // "kicked" denotes the bunch perturbed by the beam-beam interaction,
  // "kicker" is the bunch producing an electromagnetic field.
  //
  // Ie. there is only one "kicked" bunch, while the number of "kickers" is
  // equal to the number of the interation points.
  //
  // Bunch densities along X- or Y-axis are represented as a weighted sum of
  // an arbitrary number of Gaussians with the common center.  Their number
  // for representing X- and Y-densities, and for kicker(s) and kicked bunches
  // are not constrained and can be different.
  //
  // Kicked bunch is perturbed by the beam-beam interaction at every
  // interaction point, their number can be chosen arbitrarily.
  // ----------------------------------------------------------------------
  //
  // Functions to create Mutli_XY_Gaussian_bunches object
  //
  // The number of Gaussians is determined by the number of elements (Gaussian
  // widths) in the vector "sigmas".  The vector of "weights" should be either
  //
  // - of the same length as "sigmas",
  // - one less, then the last weight is complemented to the total weight of
  //  1, or
  // - absent in a single Gaussian case when "sigmas" has only one element. In
  // - this case the weight is set to 1.
  //
  // "n_sig_cut" will be explained below
  //
  void reset_kicked_bunch(int coor, // 0 for x, 1 for y
			  const vector<double>& sigmas,
			  const vector<double>& weights,
			  double n_sig_cut); // see explanations below
  //
  // reset number of interaction points (usually followed by the reset of all
  // kicker bunches)
  void reset_ip_number(int n_ip);
  //
  // "phase_advance" is the phase change of the betatron transverse
  // oscillations relative to the previous kick. In case of only one
  // interaction point, it is the full tune of the accelerator
  // ring. Otherwise, the sum of all "phase_advance"s should be equal to the
  // full tune.
  //
  void reset_kicker_bunch(int ip, // can be 0, 1, 2 , ...
			  int coor, // 0 for x, 1 for y
			  const vector<double>& sigmas,
			  const vector<double>& weights,
			  double phase_advance);
  
  //
  // "positions[ip][coor]" is a vector of beam positions: it contains the
  // coordinates of the kicker bunch centers at the interaction point "ip"
  // along "coor" axis in the frame where the kicked bunch center is at (0,0).
  //
  void reset_kicker_positions(vector<array<vector<double>, 2> > positions);
  //
  // Create bilinear interpolators to speed up density and field
  // calculations. The interpolator for one bunch is reused for all identical
  // bunches at other interaction points.
  //
  void reset_bilinear_interpolators(int n_density_cells,
				    int n_field_cells);
  //
  // Functions required in the simulation. They encapsulate the phase advances
  // and everything what dependends on the bunch densities, namely:
  //
  // Dependencies on the density of the "kicked" bunch
  // 1) probability density functions of X-X' and Y-Y' radii, not necessarily
  //    normalized (as the corresponding weights will be renormalized anyway)
  double not_normalized_r_kicked_density(int coor, double r) const;
  //
  // 2) the maximal radii in X-X' and Y-Y' to be simulated, returned as a pair
  // of doubles (.first corresponds to X, .second - to Y).
  //
  // They are chosen to reproduce the efficiency of +/- N sigmas for a single
  // Gaussian bunch, where N = n_sig_cut (the parameter above in
  // reset_kicked_bunch()). Ie. the area of the multi-Gaussian tails beyond
  // maximal radius == the area of a single Gaussian outside +/-"n_sig_cut"
  // sigmas.
  //
  pair<double, double> r_cut() const;
  //
  // Dependent on kicker bunch(es), ie. on the source of the electrostatic
  // field:
  // 1) kicker bunch density function
  double kicker_density(double x, double y, int ip) const;
  //
  // 2) the following function returns 2*epsilon0 * E-field at (x,y) created
  // by the kicker bunch placed at (0,0) assuming that its total charge is
  // equal to one.
  //
  // if "reset_bilinear_interpolators()" was not called, "field" function
  // below calculates the field directly, otherwise it uses bilinear
  // interpolation which is much faster for elliptical or multi-Gaussian
  // bunches
  complex<double> field(double x, double y, int ip) const;
  //
  // 3) phase advances deltaQ, one per interaction point and per X-/Y-axis,
  //    and the corresponding complex factor exp(2*pi*i * deltaQ)
  double deltaQ(int ip, int coor) const;
  complex<double> exp_2pi_i_deltaQ(int ip, int coor) const;
  //
  // 4) number of interaction points
  int n_ip() const;
  //
  // 5) vector of the vdM scanned x- (for coor=0) or y-coordinates (for
  // coor=1) of the kicker w.r.t. to the kicked bunch at the interaction point
  // "ip"
  const vector<double>& kicker_positions(int ip, int coor) const;
  //
  // 6) Whether bilinear interpolation of the field is used?
  //
  bool is_density_interpolated();
  bool is_field_interpolated();
  //
  // 7) Estimates interpolation accuracy for density and field:
  //
  // for every interaction point, returns the maximal (first in pair) and
  // average (second) absolute mismatches between the exact and interpolated
  // values at "n_random_points" uniformly distributed in the interpolation
  // grid, normalized by the maximal absolute value
  //
  vector<pair<double, double> >
  max_and_average_interpolation_mismatches_relative_to_max_density(int n_random_points);
  vector<pair<double, double> >
  max_and_average_interpolation_mismatches_relative_to_max_field(int n_random_points);
  //
  // Dependent on both kicker and kicked bunch densities:
  // 1) same 2*epsilon0*E-field but averaged over the first bunch centered at
  // (x,y).
  //
  // This function is used only once to calculate the orbit shift.
  complex<double> field_averaged_over_kicked_bunch(double x, double y, int ip) const;
  //
  // 2) analytically calculated overlap integral of kicked and kicker bunches
  // at the interaction point "ip" as a function of the x,y-distance between
  // their centers.
  double overlap_integral(double x, double y, int ip) const;
  //

protected:
  struct MultiG {
    vector<double> sig, w;
  protected:
    static double g(double chi) { return exp(-0.5 * chi * chi); };
  };
  friend ostream& operator<<(ostream& os, const MultiG& mg);
  friend bool operator<(const MultiG&, const MultiG&);
  struct MultiG_SigSq : public MultiG {
    // This and the following reset() functions return a string with an
    // "error" message. It is empty if no errors are detected.
    string reset(const vector<double>& sigmas, const vector<double>& weights);
    vector<double> sig_sq;
  };

  struct Kicked_MultiG : public MultiG_SigSq {
    string reset(const vector<double>& sigmas, const vector<double>& weights,
		 double n_sig_cut);
    double cut;
    double not_normalized_r_density(double r) const;
  protected:
    void find_cut(double n_sig_cut);
  };

  struct Kicker_MultiG : public MultiG_SigSq {
    string reset(const vector<double>& sigmas, const vector<double>& weights,
		 double phase_advance);
    void reset_positions(const vector<double>& positions);
    double density(double x) const;
    double deltaQ;
    complex<double> exp_2pi_i_deltaQ;
    vector<double> position;
    // Kicker bunch density: rho(x,y) = rho(x) * rho(y), where
    //  rho(x) = sum_i wx[i] /sqrt(2 pi)/sigx[i] * exp(-0.5*(x/sigx[i])^2)
    // and same for y.
    // exp_w[i] is the normalization part w[i] / sqrt(2 pi) / sig[i]:
    //
    vector<double> exp_w;
  };

  struct Kicker_MultiG_XY : public array<Kicker_MultiG, 2> {
    Kicker_MultiG_XY() : bi_density(nullptr), bi_field(nullptr) {}
    double          density(double x, double y) const;
    complex<double> field  (double x, double y) const;
    const Bilinear_interpolator<double> *bi_density;
    const Bilinear_interpolator<complex<double> > *bi_field;
  };

  vector<Kicker_MultiG_XY> kicker;
  array<Kicked_MultiG, 2> kicked;
  struct BI_with_ips {
    Bilinear_interpolator<double> density;
    Bilinear_interpolator<complex<double> > field;
    vector<int> ips;
  };
  vector<BI_with_ips> bis;
};

#endif
