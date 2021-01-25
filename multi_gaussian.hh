//
// B*B simulation of beam-beam effects in van der Meer scans at LHC.
// Copyright (C) 2019 Vladislav Balagura (balagura@cern.ch)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// ----------------------------------------------------------------------
//
// Calculation of the bunch densities and the kicker electrostatic fields for
// multi-Gaussian round and elliptical bunch densities.
//
#ifndef multi_gaussian_hh 
#define multi_gaussian_hh 1

#include <vector>
#include <array>
#include <string>
#include <complex>
#include <iostream>
#include "interpolators.hh"

using namespace std;

struct Multi_XY_Gaussian_bunches {
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
  // Functions to create Multi_XY_Gaussian_bunches object
  //
  // The number of Gaussians is determined by the number of elements (Gaussian
  // widths) in the vector "sigmas".  The vector of "weights" should be either
  //
  // - of the same length as "sigmas",
  // - one less, then the last weight is complemented to the total weight of
  //  1, or
  // - absent in a single Gaussian case when "sigmas" has only one element. In
  //   this case the weight is set to 1.
  //
  // "n_sig_cut" will be explained below
  //
  void reset_kicked_bunch(int ip, // Gaussian widths at different ip2 will
			  // be calculated as
			  // sigmas * sqrt(beta_at_ips(ip2)/beta_at_ips(ip))
			  int coor, // 0 for x, 1 for y
			  const vector<double>& sigmas, // at the given "ip"
			  const vector<double>& weights,
			  const vector<double>& beta_at_all_ips_for_given_coordinate,
			  double n_sig_cut); // see explanations below
  //
  // reset number of interaction points (usually followed by the reset of all
  // kicker bunches)
  //
  void reset_ip_number(int n_ip);
  //
  void reset_kicker_bunch(int ip, // can be 0, 1, 2 , ...
			  int coor, // 0 for x, 1 for y
			  const vector<double>& sigmas,
			  const vector<double>& weights,
			  double kicked_sigma_z_projection);
  //
  // Create bilinear interpolators to speed up density and field
  // calculations. The interpolator for one bunch is reused for all identical
  // bunches at other interaction points.
  //
  void reset_interpolators(int n_density_cells,
			   int n_field_cells,
			   vector<vector<double> >* position);
  //
  // Functions required in the simulation. They encapsulate everything that
  // dependends on the bunch densities, namely:
  //
  // Dependencies on the density of the "kicked" bunch
  // 1) probability density functions of X-X' and Y-Y' radii, not necessarily
  //    normalized (as the corresponding weights will be renormalized anyway).
  //    Gaussian widths at "ip" scale according to sqrt(beta_ip / beta_0) *
  //    sigmas_0
  //
  double not_normalized_r_kicked_density(int coor, double r, int ip) const;
  //
  // 2) the maximal radii in X-X' and Y-Y' to be simulated at "ip", returned
  // as a pair of doubles (.first corresponds to X, .second - to Y).
  //
  // They are chosen to reproduce the efficiency of +/- N sigmas for a single
  // Gaussian bunch, where N = n_sig_cut (the parameter above in
  // reset_kicked_bunch()). Ie. the area of the multi-Gaussian tails beyond
  // maximal radius == the area of a single Gaussian outside +/-"n_sig_cut"
  // sigmas.
  //
  pair<double, double> r_cut(int ip) const;
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
  // if "reset_interpolators()" was not called, "field" function below
  // calculates the field directly, otherwise it uses bilinear interpolation
  // which is much faster for elliptical or multi-Gaussian bunches
  complex<double> field(double x, double y, int ip) const;
  //
  // 4) number of interaction points (for completeness, but not used)
  int n_ip() const;
  //
  // 5) Whether bilinear interpolation of the field is used?
  //
  bool is_density_interpolated();
  bool is_field_interpolated();
  //
  // 6) Estimates interpolation accuracy for density and field:
  //
  // for every interaction point "ip" and the coordinate "coor" (0 or 1 for X
  // or Y, respectively), returns the maximal (first in pair) and average
  // (second) absolute mismatches between the exact and interpolated values at
  // "n_random_points" uniformly distributed in the interpolation grid,
  // normalized by the maximal absolute value. In the returned structure the
  // first index is "ip", the second is "coor": [ip][coor]
  //
  vector<array<pair<double, double>, 2> >
  max_and_average_interpolation_mismatches_relative_to_max_density(int n_random_points);
  //
  // The same maximal (first in pair) and average (second) absolute mismatches
  // between the exact and interpolated two-dimensional electric E-field at
  // "n_random_points" uniformly distributed in the 2D interpolation grid,
  // normalized by the maximal absolute |E|-value. The returned structure has
  // only one index "ip": [ip]
  //
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
    static double g_sig2(double x, double sig_sq) { return exp(-0.5 * x * x / sig_sq); };
  };
  friend ostream& operator<<(ostream& os, const MultiG& mg);
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
    Kicker_MultiG() : li_density(nullptr) {}
    string reset(const vector<double>& sigmas, const vector<double>& weights,
		 double sigma_z_projection);
    double density(double x) const;
    // Kicker bunch density: rho(x,y) = rho(x) * rho(y), where
    //  rho(x) = sum_i wx[i] /sqrt(2 pi)/sigx[i] * exp(-0.5*(x/sigx[i])^2)
    // and same for y.
    // exp_w[i] is the normalization part w[i] / sqrt(2 pi) / sig[i]:
    //
    vector<double> exp_w;
    const Linear_interpolator *li_density;
  };
  // to make maps with MultiG and Kicker_MultiG keys, this is needed in
  // reset_interpolators().
  friend bool operator<(const MultiG&, const MultiG&);

  struct Kicker_MultiG_XY : public array<Kicker_MultiG, 2> {
    Kicker_MultiG_XY() : bi_field(nullptr) {}
    complex<double> field  (double x, double y) const;
    const Bilinear_interpolator *bi_field;
  };

  vector<Kicker_MultiG_XY> kicker;
  // kicked bunch internally is stored for ip=0 (its sigmas depend on beta*)
  array<Kicked_MultiG, 2> kicked;
  // interpolators
  struct LI {
    Linear_interpolator density;
    struct Ip_Coor {
      int ip, coor;
    };
    vector<Ip_Coor> ip_coor;
  };
  vector<LI> lis;

  struct BI {
    Bilinear_interpolator field;
    vector<int> ips;
  };
  vector<BI> bis;
  vector<array<double, 2> > beta; // beta[ip][0/1 for x/y]
};

#endif
