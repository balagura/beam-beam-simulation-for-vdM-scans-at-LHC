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
// multi-Gaussian round and elliptical bunch shapes
//
#include <algorithm>
#include <numeric>
#include <map>
#include <limits>
#include "multi_gaussian.hh"
#include "E_field.hh"
//#include <functional> // for placeholders

// set Gaussian sigmas, weights and check their consistency
// (empty input "weights" vector means this parameter is not given)
string Multi_XY_Gaussian_bunches::MultiG_SigSq::reset(const vector<double>& sigmas,
						      const vector<double>& weights) {
  if (sigmas.empty()) return "vector of Gaussian sigmas is empty";
  sig = sigmas;
  sig_sq.resize(sig.size());
  transform(sig.begin(), sig.end(), sig_sq.begin(), [](double e) { return e*e; });
  if (weights.size() != 0) {
    w = weights;
    if (sig.size()-1 == w.size()) {
      // complement last Gaussian weight so that sum of all == 1
      w.push_back(1 - accumulate(w.begin(), w.end(), 0.));
    } else if (sig.size() == w.size()) {
      // normalize
      double wsum = accumulate(w.begin(), w.end(), 0.);
      for_each(w.begin(), w.end(), [wsum](double& x) { x /= wsum; });
    } else {
      return "Gaussian weights and sigmas are inconsistent";
    }
  } else {
    if (sig.size() != 1) {
      return "Gaussian weights are not given";
    } else {
      w = {1};
    }
  }
  return "";
}
// Kicked_MultiG functions
string Multi_XY_Gaussian_bunches::Kicked_MultiG::reset(const vector<double>& sigmas,
						       const vector<double>& weights,
						       double n_sig_cut) {
  string err = MultiG_SigSq::reset(sigmas, weights);
  if (err != "") return err;
  find_cut(n_sig_cut);
  return "";
}
void Multi_XY_Gaussian_bunches::Kicked_MultiG::find_cut(double n_sig_cut) {
  //
  // For a single Gaussian case of the first bunch, X-X' and Y-Y' circles
  // with radii beyond n_sig_cut ("N.sigma.cut" in the configuration file)
  // are not simulated. 
  // The distribution of eg. X-X' radius is:
  // 1/2 /pi / sigx^2 * 2pi r * exp(-0.5 * (r/sig)^2) dr =
  // 1/sig^2 * r * exp(-0.5 * (r/sig)^2) dr
  //
  // By integrating this expression one can find that not simulated part of
  // bunch 1 is exp(-n_sig_cut^2 / 2) in x and the same in y, ie. in total
  // 1 - (1 - exp(-n_sig_cut^2/2))^2, eg. for n_sig_cut = 5 it is 7.4e-6.
  //
  //
  // For multi-Gaussian case, select analogous [0, cut_x] range
  // corresponding to [0, n_sig_cut * sigma] for a single Gaussian, such
  // that it contains 1 - exp(-n_sig_cut^2 / 2) of all radii distributed as
  // a sum of N X-X' 2D Gaussians.
  // 
  // Ie. for N Gaussians we require:
  // sum_{i=1}^N integral_0^cut_x (w[i] / sig[i]^2 * r * exp(-0.5 * (r/sig[i])^2) dr =
  // sum_{i=1}^N w[i] * (1 - exp(-0.5 * (cut_x / sig[i])^2)) = 1 - exp(-n_sig_cut^2 / 2)
  // or, since sum_{i=1}^N w[i] = 1:
  //
  // sum_{i=1}^N w[i] * exp(-0.5 * (cut_x / sig[i])^2) = exp(-n_sig_cut^2 / 2).
  //
  // Solve this equation numerically (and same for y):
  //
  if (sig.size() == 1) {
    cut = n_sig_cut * sig[0];
  } else {
    double rhs = g(n_sig_cut);
    auto f = [this, rhs] (double r) {
	       double p = 0;
	       for (size_t i=0; i<sig.size(); ++i) p += w[i] * g( r/sig[i] );
	       return p - rhs;
	     };
    auto sig_range = minmax_element(sig.begin(), sig.end());
    double
      low  = *sig_range.first  * n_sig_cut,
      middle,
      high = *sig_range.second * n_sig_cut;
    // f(r) is monotonically decreasing, so f(low)>=0, f(high)<=0,
    // find the root f() = 0 by the bisection method:
    while (true) {
      middle = (low + high) / 2;
      double f_middle = f(middle);
      // cout << "l,m,h,f(l,m,h): " << low << " " << middle << " " << high << " "
      //      << f(low) << " " << f_middle << " " << f(high) << endl;
      if (abs(f_middle) < 1e-8) break; // require the integral outside [0, cut_x]
      // to be equal to the desired value (exp(-nsig_cut^2/2)) within +/-(1e-8)
      if (f_middle < 0) high = middle; else low = middle;
    }
    cut = middle;
  }
}
double Multi_XY_Gaussian_bunches::Kicked_MultiG::not_normalized_r_density(double r) const {
  // distribution of radius length for eg. X-X' 2D Gaussian:
  //
  // exp(-r^2/2/sig^2) * r / sig^2 dr
  //
  // the absolute normalization can be arbitrary here
  // and independent on the r bin size since the weights
  // will be renormalized anyway.
  //
  // The resulting weight will be a product of two
  // weights for X-X' and Y-Y'.
  //
  double d = 0;
  for (size_t i=0; i<w.size(); ++i) {
    d += w[i] / sig[i] / sig[i] * r * g(r / sig[i]);
  }
  return d;
}
// Kicker_MultiG functions
string Multi_XY_Gaussian_bunches::Kicker_MultiG::reset(const vector<double>& sigmas,
						       const vector<double>& weights,
						       double sigma_z_projection) {
  string err = Multi_XY_Gaussian_bunches::MultiG_SigSq::reset(sigmas, weights);
  if (err != "") return err;
  exp_w.resize(w.size());
  for (size_t i=0; i<exp_w.size(); ++i) {
    exp_w[i] = w[i] / sqrt(2 * M_PI) / sig[i];
  }
  sig_z_projection_sq = sigma_z_projection * sigma_z_projection;
  return "";
}
double Multi_XY_Gaussian_bunches::Kicker_MultiG::density(double x) const {
  double prob = 0;
  for (size_t i=0; i<sig.size(); ++i) prob += exp_w[i] * g_sig2(x, sig_sq[i] + sig_z_projection_sq);
  return prob;
}
//
void Multi_XY_Gaussian_bunches::reset_kicked_bunch(int ip, // Gaussian widths at different ip2 will
						   // be calculated as
						   // sigmas * sqrt(beta_at_ips(ip2)/beta_at_ips(ip))
						   int coor, // 0 for x, 1 for y
						   const vector<double>& sigmas, // at the given "ip"
						   const vector<double>& weights,
						   const vector<double>& beta_at_all_ips_for_given_coordinate,
						   double n_sig_cut) {
  if (ip<0 || ip>=beta_at_all_ips_for_given_coordinate.size()) {
    cerr << "Given ip number of the kicked bunch (" << ip << ") is out of range of the given vector of beta*\n";
    exit(1);
  }
  beta.resize(beta_at_all_ips_for_given_coordinate.size());
  for (size_t ip=0; ip<beta.size(); ++ip)
    beta[ip][coor] = beta_at_all_ips_for_given_coordinate[ip];
  // rescale sigmas so that "kicked" bunch always corresponds to ip=0
  vector<double> sig(sigmas);
  for (auto& s: sig) s *= sqrt(beta_at_all_ips_for_given_coordinate[0] /
			       beta_at_all_ips_for_given_coordinate[ip]);
  string err = kicked[coor].reset(sigmas, weights, n_sig_cut);
  if (err != "") {
    cerr << "Kicked bunch, " << "xy"[coor] << "-density: " << err << endl;
    exit(1);
  }
}
void Multi_XY_Gaussian_bunches::reset_ip_number(int n_ip) { kicker.resize(n_ip); }
void Multi_XY_Gaussian_bunches::reset_kicker_bunch(int ip, // can be 0, 1, 2 , ...
						   int coor, // 0 for x, 1 for y
						   const vector<double>& sigmas,
						   const vector<double>& weights,
						   double kicked_sigma_z_projection) {
  string err = kicker[ip][coor].reset(sigmas, weights, kicked_sigma_z_projection);
  if (err != "") {
    cerr << "Kicker." << ip+1 << "." << "xy"[coor] << ": " << err << endl;
    exit(1);
  }
}
ostream& operator<<(ostream& os, const Multi_XY_Gaussian_bunches::MultiG& x) {
  os << "sigma (weight) = ";
  for (size_t i=0; i<x.sig.size(); ++i) {
    os << x.sig[i] << " (" << x.w[i] << ")";
    if (i != x.sig.size()-1) os << ", ";
  }
  return os;
}
bool operator<(const Multi_XY_Gaussian_bunches::MultiG& a, const Multi_XY_Gaussian_bunches::MultiG& b) {
  return a.sig < b.sig || (a.sig == b.sig && a.w < b.w);
  // vectors and arrays are compared by default lexicographically
}
void Multi_XY_Gaussian_bunches::reset_interpolators(int n_density_cells,
						    int n_field_cells,
						    vector<vector<double> >* position) {
  if (n_density_cells > 0) {
    // group identical X or Y kicker bunch ditributions into a map
    map<MultiG, vector<LI::Ip_Coor> > m;
    // m[MultiG_X + Y] = vector of kicker (ip, coor)
    for (int ip = 0; ip < int(kicker.size()); ++ip) {
      for (int coor=0; coor<2; ++coor) {
	const MultiG& mg_XY = kicker[ip][coor];
	LI::Ip_Coor i; i.ip = ip; i.coor = coor;
	m[mg_XY].push_back(i);
      }
    }
    lis.clear(); // in case it was already initialized
    // loop over groups of identical X+Y bunch distributions: find the range
    // of the kicker positions inside the group and choose wide enough grid
    // for a linear interpolator
    lis.resize(m.size()); // create one per group
    int i_li = 0;
    for (const auto& p: m) {
      pair<double, double> range = make_pair(numeric_limits<double>::max(), -numeric_limits<double>::max());
      // first/second = min/max, initialize min by maximal double and vice versa
      for (auto ip_coor: p.second) {
	int ip = ip_coor.ip, coor = ip_coor.coor;
	vector<double>& p = position[coor][ip];
	auto r = minmax_element(p.begin(), p.end());
	// Positions are specified in the frame where the kicked bunch is at
	// (0,0), while the density map - with the kicker at (0,0). So, for the
	// density map the positions should be subtracted and min/max - swapped.
	// In addition, add +/-max simulated radius kicked[coor].cut:
	range.first  = min(-(*r.second) - 1.1 * kicked[coor].cut, range.first);
	range.second = max(-(*r.first ) + 1.1 * kicked[coor].cut, range.second);
	// 10% margin (1.1) since beam-beam can increase the radius
      }
      auto& li = lis[i_li].density;
      auto ip_coor0 = p.second[0]; // take first density in the group, others are identical
      int ip0 = ip_coor0.ip, coor0 = ip_coor0.coor;
      li.create(range.first, range.second,
		n_density_cells,
		// note, k in lambda is passed by reference, to avoid
		// capturing the copy of kicker[ip0][coor0]. So, the kicker
		// should not be reallocated in memory afterwards.
		[&k = kicker[ip0][coor0]](double x) { return k.density(x); });
      // store the pointers to the created linear interpolator in the kickers
      for (auto ip_coor: p.second) kicker[ip_coor.ip][ip_coor.coor].li_density = &li;
      lis[i_li].ip_coor = p.second;
      ++i_li;
    }
  }
  if (n_field_cells > 0) {
    // group identical bunches into a map
    map<array<MultiG, 2>, vector<int> > m; // m[(MultiG_X, MultiG_Y)] = vector of kicker ips
    for (int ip = 0; ip < int(kicker.size()); ++ip) {
      const auto& k = kicker[ip];
      array<MultiG, 2> mg_XY;
      for (int coor=0; coor<2; ++coor) mg_XY[coor] = k[coor];
      m[mg_XY].push_back(ip);
    }
    bis.clear(); // in case it was already initialized
    // loop over groups of identical bunches: find the range of the beam
    // movements at IPs where bunches are identical and define bilinear
    // interpolator in a wide enough grid
    bis.resize(m.size()); // create one per group
    int i_bi = 0;
    for (const auto& p: m) {
      array<pair<double, double>, 2> range; // [0/1 for x/y].first/second = min/max
      // initialize min by maximal double and vice versa
      range[1] = range[0] = make_pair(numeric_limits<double>::max(), -numeric_limits<double>::max());
      for (int ip: p.second) { // ip = kicker index
	for (int coor=0; coor<2; ++coor) {
	  vector<double>& p = position[coor][ip];
	  auto r = minmax_element(p.begin(), p.end());
	  // Positions are specified in the frame where the kicked bunch is at
	  // (0,0), while the field map - with the kicker at (0,0). So, for the
	  // field map the positions should be subtracted and min/max - swapped:
	  range[coor].first  = min(-(*r.second), range[coor].first);
	  range[coor].second = max(-(*r.first ), range[coor].second);
	}
      }
      for (int coor=0; coor<2; ++coor) {
	range[coor].first  -= 1.1 * kicked[coor].cut;
	range[coor].second += 1.1 * kicked[coor].cut;
	// 10% margin (1.1) since beam-beam can increase the radius
      }
      auto& bi = bis[i_bi].field;
      bi.create(range[0].first, range[0].second,
		range[1].first, range[1].second,
		n_field_cells,
		[&k = kicker[p.second[0]]](double x, double y) { return k.field(x, y); });
	// Other option using binding and placeholders (for a member function the
	// first "bind" arg is *this, here it is again fixed to kicker[p.second[0]]):
	//			 bind(&Kicker_MultiG_XY::density,
	//			      kicker[p.second[0]],
	//			      placeholders::_1, placeholders::_2));
	// 
      for (int ip: p.second) kicker[ip].bi_field = &bi;
      bis[i_bi].ips = p.second;
      ++i_bi;
    }
  }
}
double Multi_XY_Gaussian_bunches::not_normalized_r_kicked_density(int coor, double r, int ip) const {
  // "kicked" is stored internally with sigmas corresponding to "ip"=0, rescale "r" also to ip=0:
  r *= sqrt(beta[0][coor] / beta[ip][coor]);
  return kicked[coor].not_normalized_r_density(r);
}
pair<double, double> Multi_XY_Gaussian_bunches::r_cut(int ip) const {
  // rescale from ip=0
  return make_pair(sqrt(beta[ip][0] / beta[0][0]) * kicked[0].cut,
		   sqrt(beta[ip][1] / beta[0][1]) * kicked[1].cut);
}
double Multi_XY_Gaussian_bunches::kicker_density(double x, double y, int ip) const {
  const auto& k = kicker[ip];
  if (k[0].li_density == nullptr) {
    return k[0].density(x) * k[1].density(y);
  } else {
    return (*k[0].li_density)(x) * (*k[1].li_density)(y);
  }
}
complex<double> Multi_XY_Gaussian_bunches::field(double x, double y, int ip) const {
  const auto& k = kicker[ip];
  if (k.bi_field == nullptr) {
    return k.field(x, y);
  } else {
    return (*k.bi_field)(x, y);
  }
}
complex<double> Multi_XY_Gaussian_bunches::Kicker_MultiG_XY::field(double x, double y) const {
  const MultiG_SigSq &x2 = (*this)[0], &y2 = (*this)[1];
  complex<double> E = 0;
  for (size_t ix2=0; ix2<x2.sig.size(); ++ix2) {
    for (size_t iy2=0; iy2<y2.sig.size(); ++iy2) {
      E += x2.w[ix2] * y2.w[iy2] *
	E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y, x2.sig[ix2], y2.sig[iy2]);
    }
  }
  return E;
}
int Multi_XY_Gaussian_bunches::n_ip() const { return kicker.size(); }
bool Multi_XY_Gaussian_bunches::is_density_interpolated() { return !lis.empty(); }
// check bis[0] only, others should be the same since the interpolation is
// switched on/off for all IPs at once
bool Multi_XY_Gaussian_bunches::is_field_interpolated  () { return !bis.empty(); }
vector<array<pair<double, double>, 2> >
Multi_XY_Gaussian_bunches::max_and_average_interpolation_mismatches_relative_to_max_density(int n_random_points) {
  vector<array<pair<double, double>, 2> > res(kicker.size());
  for (const auto& x: lis) {
    pair<double, double> p =
      x.density.max_and_average_mismatches_relative_to_max(n_random_points);
    for (auto i: x.ip_coor) res[i.ip][i.coor] = p;
  }
  return res;
}
vector<pair<double, double> >
Multi_XY_Gaussian_bunches::max_and_average_interpolation_mismatches_relative_to_max_field(int n_random_points) {
  vector<pair<double, double> > res(kicker.size());
  for (const auto& x: bis) {
    pair<double, double> p =
      x.field.max_and_average_mismatches_relative_to_max(n_random_points);
    for (int i: x.ips) res[i] = p;
  }
  return res;
}
complex<double> Multi_XY_Gaussian_bunches::field_averaged_over_kicked_bunch(double x, double y, int ip) const {
  complex<double> E = 0;
  const MultiG_SigSq &x1 = kicked[0], &y1 = kicked[1], &x2 = kicker[ip][0], &y2 = kicker[ip][1];
  // scale squared sigmas from IP=0 to IP="ip"
  vector<double>
    sig_sq_x1(x1.sig_sq),
    sig_sq_y1(y1.sig_sq);
  double beta_scale_x = beta[ip][0] / beta[0][0];
  double beta_scale_y = beta[ip][1] / beta[0][1];
  for (auto& s: sig_sq_x1) s *= beta_scale_x;
  for (auto& s: sig_sq_y1) s *= beta_scale_y;
  //
  for (size_t ix1=0; ix1<x1.sig.size(); ++ix1) {
    for (size_t iy1=0; iy1<y1.sig.size(); ++iy1) {
      for (size_t ix2=0; ix2<x2.sig.size(); ++ix2) {
	for (size_t iy2=0; iy2<y2.sig.size(); ++iy2) {
	  double vdm_sigx = sqrt(sig_sq_x1[ix1] + x2.sig_sq[ix2]); // use scaled kicked sig_sq here
	  double vdm_sigy = sqrt(sig_sq_y1[iy1] + y2.sig_sq[iy2]);
	  E += x1.w[ix1] * y1.w[iy1] * x2.w[ix2] * y2.w[iy2] *
	    E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y,
							       vdm_sigx, vdm_sigy);
	}
      }
    }
  }
  return E;
}
double Multi_XY_Gaussian_bunches::overlap_integral(double x, double y, int ip) const {
  double o = 0;
  double x_sq = x*x, y_sq = y*y;
  const MultiG_SigSq &x1 = kicked[0], &y1 = kicked[1], &x2 = kicker[ip][0], &y2 = kicker[ip][1]; 
  // scale squared sigmas from IP=0 to IP="ip"
  vector<double>
    sig_sq_x1(x1.sig_sq),
    sig_sq_y1(y1.sig_sq);
  double beta_scale_x = beta[ip][0] / beta[0][0];
  double beta_scale_y = beta[ip][1] / beta[0][1];
  for (auto& s: sig_sq_x1) s *= beta_scale_x;
  for (auto& s: sig_sq_y1) s *= beta_scale_y;
  //
  for (size_t ix1=0; ix1<x1.sig.size(); ++ix1) {
    for (size_t iy1=0; iy1<y1.sig.size(); ++iy1) {
      for (size_t ix2=0; ix2<x2.sig.size(); ++ix2) {
	for (size_t iy2=0; iy2<y2.sig.size(); ++iy2) {
	  double vdm_sigx_sq = sig_sq_x1[ix1] + x2.sig_sq[ix2]; // use scaled kicked sig_sq here
	  double vdm_sigy_sq = sig_sq_y1[iy1] + y2.sig_sq[iy2];
	  o += x1.w[ix1] * y1.w[iy1] * x2.w[ix2] * y2.w[iy2] *
	    0.5 / M_PI / sqrt(vdm_sigx_sq * vdm_sigy_sq) *
	    exp(-0.5 * (x_sq / vdm_sigx_sq + y_sq / vdm_sigy_sq));
	}
      }
    }
  }
  return o;
}
