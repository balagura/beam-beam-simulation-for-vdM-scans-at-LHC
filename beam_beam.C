// B*B simulation of beam-beam effects in van der Meer scans at LHC.
//
// Author V. Balagura, balagura@cern.ch (Jan 2019)
//
// Bassetti-Erskine beam-beam kick calculation for elliptical bunches is based
// on the VdmBBDeflect class adopted to B*B by Timothy Barklow (SLAC) and then
// rewritten as Bassetti_Erskine_angular_kick(..) by Vladislav Balagura.

#include <vector>
#include <string>
#include <complex>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h> // for mkdir
#include <sys/types.h> // for mkdir
#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <thread>
#include <atomic>
#include <mutex>
#include <gzstream.h>
#include <Faddeeva.hh>
#include <config.hh>
#include <E_field.hh>
#include <bilinear_interpolator.hh>

using namespace std;
using namespace std::complex_literals;

ostream& operator<<(ostream& os, complex<double> c) { return os << real(c) << " " << imag(c); }

// Read configuration and encapsulate everything dependent on the distribution
// of the bunch densities, namely:
//
// Dependent on bunch 1, ie. on the one perturbed by beam-beam:
// 1) probability density functions of X-X' and Y-Y' radii:
//    "not_normalized_rx_density1(rx)" and "not_normalized_ry_density1(ry)",
//    respectively (not necessarily normalized since they will be renormalized
//    anyway),
//
// 2) the maximal radii in X-X' and Y-Y' to be simulated: "cut_x" and "cut_y".
//
// Dependent on bunch 2, ie. on the source of the electrostatic field:
// 1) not normalized bunch 2 density function, "not_normalized_density2(x,y)",
// and
// 2) its normalization constant, "density2_normalization" - for a calculation
// of an overlap integral (they are split in 1)+2) for efficiency, to move the
// overall normalization out of the innermost loop),
//
// 3) 2*epsilon0 * E-field(x,y) for unit charge of the second bunch placed at
// (0,0) "field2(x,y)" - to find the kick.
//
// Dependent on both bunch 1 and bunch 2 densities:
// 1) same 2*epsilon0*E-field but averaged over the first bunch placed at
// (x,y), "field2_averaged_over_bunch_1" - to calculate the orbit shift
//
// 2) analytically calculated overlap integral of two bunches unperturbed by
// beam-beam as a function of their x,y-separation, "overlap_integral(x,y)".
//
struct Bunch_density_dependency {
  double cut_x, cut_y;
  function<double (double)> not_normalized_rx_density1, not_normalized_ry_density1;
  function<double (double, double)> not_normalized_density2;
  double density2_normalization;
  // Bilinear_interpolator will need a function "field2" below with x,y
  // arguments returning a 2*epsilon0*E-field at x,y created by the 2d bunch
  // centered at (0,0). More specifically, Bilinear_interpolator expects
  // std::function<complex<double>(double, double)> or something convertible
  // to it (std::function is a wrapper around various function objects, for
  // details, please, see C++ documentation on std::function). "field2" will
  // be used both in Bilinear_interpolator and for the direct field
  // calculation if the interpolation is switched off
  function<complex<double> (double, double)> field2;
  // "field2_averaged_over_bunch_1" below is used only once to predict the
  // orbit shift
  function<complex<double> (double, double)> field2_averaged_over_bunch_1;
  function<double (double, double)> overlap_integral;
};
// implementation class for multi-Gaussian bunch densities
//
struct Mutli_XY_Gaussian_bunches : public Bunch_density_dependency {
  Mutli_XY_Gaussian_bunches(Config& c) {
    vector<double> sig1x = c.vd("sig1.x"), sig1y = c.vd("sig1.y");
    vector<double> sig2x = c.vd("sig2.x"), sig2y = c.vd("sig2.y");
    vector<double>
      sig1x_sq(sig1x.size()), sig1y_sq(sig1y.size()),
      sig2x_sq(sig2x.size()), sig2y_sq(sig2y.size());
    auto mk_square = [] (const vector<double>& v, vector<double>* v_squared) {
		       transform(v.begin(), v.end(), v_squared->begin(), [](double e) { return e*e; });
		     };
    mk_square(sig1x, &sig1x_sq);    mk_square(sig2x, &sig2x_sq);
    mk_square(sig1y, &sig1y_sq);    mk_square(sig2y, &sig2y_sq);
    vector<double> w1x, w1y, w2x, w2y;
    // get Gaussian weights and check consistency of input:
    auto get_w = [&c] (const string& weights_name, size_t n_sig, vector<double>* w) -> void {
		   if (c.defined(weights_name)) {
		     *w = c.vd(weights_name);
		     if (n_sig-1 == w->size()) {
		       // complement last Gaussian weight so that sum of all == 1
		       w->push_back(1 - accumulate(w->begin(), w->end(), 0.));
		     } else if (n_sig != w->size()) {
		       cerr << "N elements in " << weights_name << " and in corresponding sig are different\n";
		       exit(1);
		     } else {
		       // normalize
		       double wsum = accumulate(w->begin(), w->end(), 0.);
		       transform(w->begin(), w->end(), w->begin(), [wsum](double x) { return x/wsum; });
		     }
		   } else {
		     if (n_sig != 1) {
		       cerr << weights_name << " weights are not given\n";
		       exit(1);
		     } else {
		       (*w) = {1};
		     }
		   }
		 };
    get_w("b1x.gaussian.weight", sig1x.size(), &w1x);
    get_w("b1y.gaussian.weight", sig1y.size(), &w1y);
    get_w("b2x.gaussian.weight", sig2x.size(), &w2x);
    get_w("b2y.gaussian.weight", sig2y.size(), &w2y);
    // For a single Gaussian case of the first bunch, X-X' and Y-Y' circles
    // with radii beyond n_sig_cut ("N.sigma.cut" in the configuration file)
    // are not simulated. 
    // The distribution of eg. X-X' radius is:
    // 1/2 pi/sigx^2 * 2pi r * exp(-0.5 * (r/sig)^2) dr =
    // 1/sig^2 * r * exp(-0.5 * (r/sig)^2) dr
    //
    // By integrating this expression one can find that not simulated part of
    // bunch 1 is exp(-n_sig_cut^2 / 2) in x and the same in y, ie. in total
    // 1 - (1 - exp(-n_sig_cut^2/2))^2, eg. for n_sig_cut = 5 it is 7.4e-6.
    //
    double n_sig_cut = c.defined("N.sigma.cut") ? c("N.sigma.cut") : 5;
    // For both multi-Gaussian case, select [0, cut_x] range corresponding to
    // [0, n_sig_cut * sigma] for a single Gaussian, such that it contains
    // 1 - exp(-n_sig_cut^2 / 2) of all radii distributed as a sum of N X-X'
    // 2D Gaussians.
    // 
    // Ie. for N Gaussians we require:
    // sum_{i=1}^N integral_0^cut_x (w[i] / sig[i]^2 * r * exp(-0.5 * (r/sig[i])^2) dr =
    // sum_{i=1}^N w[i] * (1 - exp(-0.5 * (R / sig[i])^2)) = 1 - exp(-n_sig_cut^2 / 2)
    // or, since sum_{i=1}^N w[i] = 1:
    //
    // sum_{i=1}^N w[i] * exp(-0.5 * (R / sig[i])^2) = exp(-n_sig_cut^2 / 2).
    //
    // Solve this equation numerically (and same for y)
    //
    auto g = [](double chi) { return exp(-0.5 * chi * chi); };
    auto root = [&g] (const vector<double>& sig, const vector<double>& w, double n_sig_cut) {
		if (sig.size() == 1) return n_sig_cut * sig[0];
		//
		double rhs = g(n_sig_cut);
		auto f = [&sig, &w, &g, rhs] (double r) {
			   double p = 0;
			   for (size_t i=0; i<sig.size(); ++i) p += w[i] * g( r/sig[i] );
			   return p - rhs;
			 };
		auto sig_range = minmax_element(sig.begin(), sig.end());
		double
		  low  = *sig_range.first  * n_sig_cut,
		  high = *sig_range.second * n_sig_cut, middle;
		// f(r) is monotonically decreasing, so f(low)>=0, f(high)<=0,
		// find the root f() = 0 by the bisection method:
		cout << sig[0] << " " << sig[1] << " " << w[0] << " " << w[1] << endl;
		while (true) {
		  middle = (low + high) / 2;
		  double f_middle = f(middle);
		  double tmp = 0;
		  for (size_t i = 0; i<sig.size(); ++i) tmp += w[i] / sig[i] / sig[i];
		  cout << "l,m,h,f(m): " << low << " " << middle << " " << high << " "
		       << f(low) << " " << f_middle << " " << f(high)
		       << " " << sqrt(n_sig_cut*n_sig_cut / tmp) << endl;
		  if (abs(f_middle) < 1e-8) break; // require the integral outside [0, cut_x]
		  // to be equal to the desired value (exp(-nsig_cut^2/2)) within +/-(1e-8)
		  if (f_middle < 0) high = middle; else low = middle;
		}
		return middle;
	      };
    cut_x = root(sig1x, w1x, n_sig_cut);
    cut_y = root(sig1y, w1y, n_sig_cut);
    auto r_density1 = [&g] (const vector<double>& w, const vector<double>& sig, double r) {
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
		      };
    not_normalized_rx_density1 = [r_density1, w1x, sig1x] (double rx) { return r_density1(w1x, sig1x, rx); };
    not_normalized_ry_density1 = [r_density1, w1y, sig1y] (double ry) { return r_density1(w1y, sig1y, ry); };
    overlap_integral = [sig1x_sq, sig1y_sq,
			sig2x_sq, sig2y_sq,
			w1x, w1y, w2x, w2y] (double x, double y) {
			 double o = 0;
			 double x_sq = x*x, y_sq = y*y;
			 for (size_t ix1=0; ix1<sig1x_sq.size(); ++ix1) {
			   for (size_t iy1=0; iy1<sig1y_sq.size(); ++iy1) {
			     for (size_t ix2=0; ix2<sig2x_sq.size(); ++ix2) {
			       for (size_t iy2=0; iy2<sig2y_sq.size(); ++iy2) {
				 double vdm_sigx_sq = sig1x_sq[ix1] + sig2x_sq[ix2];
				 double vdm_sigy_sq = sig1y_sq[iy1] + sig2y_sq[iy2];
				 o += w1x[ix1] * w1y[iy1] * w2x[ix2] * w2y[iy2] *
				   0.5 / M_PI / sqrt(vdm_sigx_sq * vdm_sigy_sq) *
				   exp(-0.5 * (x_sq / vdm_sigx_sq + y_sq / vdm_sigy_sq));
			       }
			     }
			   }
			 }
			 return o;
		       };
			
    // In case of double Gaussian shape of bunch 2: rho2(x,y) = rho2(x) * rho2(y), where
    //  rho2(x) = wx /sqrt(2 pi)/sig2x * exp(.1.) + (1 - wx) /sqrt(2 pi)/sig2x_2g * exp(.2.)
    // and wx = c["b2.1g.weight.x"] and same for y.
    // To minimize multiplications, rho2(x,y) is equivalently written as
    //  rho2(x,y) = bunch_2_prob_const *
    //                (exp(.1.) + b2_exp_weight_x_g2 * exp(.2.)) *
    //                (exp(.1.) + b2_exp_weight_y_g2 * exp(.2.))
    // where
    //    bunch_2_prob_const = wx * wy / (2 pi) / sig2x / sig2y,
    //    b2_exp_weight_x_g2 = (1 - wx) / wx * (sig2x / sig2x_2g),
    //    b2_exp_weight_y_g2 = (1 - wy) / wy * (sig2y / sig2y_2g)
    //
    // Probability density of the second bunch = bunch_2_prob_const * not_normalized_density(),
    // multiplication by bunch_2_prob_const is moved out of the mostinner loop
    //
    // TODO: change doc for multi-G
    vector<double> exp_w2x(w2x.size()-1), exp_w2y(w2y.size()-1);
    auto fill_exp_w = [] (const vector<double>& sig, const vector<double>& w, vector<double>* exp_w) {
			for (size_t i=0; i<exp_w->size(); ++i) {
			  (*exp_w)[i] = w[i+1] / sig[i+1] * sig[0] / w[0];
			}
		      };
    fill_exp_w(sig2x, w2x, &exp_w2x);
    fill_exp_w(sig2y, w2y, &exp_w2y);
    //
    if (sig2x.size() == 1 && sig2y.size() == 1) {
      not_normalized_density2 = [sigx = sig2x[0], sigy = sig2y[0]] (double x, double y) {
				 double chi_x = x / sigx;
				 double chi_y = y / sigy;
				 return exp(-0.5 * (chi_x * chi_x + chi_y * chi_y));
			       };
      field2 = [sigx = sig2x[0], sigy = sig2y[0]](double x, double y) {
		return E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y, sigx, sigy);
	      };
      field2_averaged_over_bunch_1 =
	[vdm_sigx = sqrt(sig1x_sq[0] + sig2x_sq[0]),
	 vdm_sigy = sqrt(sig1y_sq[0] + sig2y_sq[0])] (double x, double y) {
	  return E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y,
								    vdm_sigx, vdm_sigy);
	};
    } else {
      not_normalized_density2 = [sig2x, sig2y,
				 exp_w2x, exp_w2y, &g] (double x, double y) {
				 auto prob = [&g] (double x, const vector<double>& sig, const vector<double>& exp_w) {
					       double prob = g(x / sig[0]);
					       for (size_t i=1; i<sig.size(); ++i) prob += exp_w[i-1] * g(x / sig[i]);
					       return prob;
					     };
				 return prob(x, sig2x, exp_w2x) * prob(y, sig2y, exp_w2y);
			       };
      field2 = [sig2x, sig2y,
	       w2x, w2y] (double x, double y) {
		complex<double> E = 0;
		for (size_t ix=0; ix<sig2x.size(); ++ix) {
		  for (size_t iy=0; iy<sig2y.size(); ++iy) {
		    E += w2x[ix] * w2y[iy] *
		      E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y, sig2x[ix], sig2y[iy]);
		  }
		}
		return E;
	      };
      field2_averaged_over_bunch_1 =
	[sig1x_sq, sig1y_sq,
	 sig2x_sq, sig2y_sq,
	 w1x, w1y, w2x, w2y] (double x, double y) {
	  complex<double> E = 0;
	  for (size_t ix1=0; ix1<sig1x_sq.size(); ++ix1) {
	    for (size_t iy1=0; iy1<sig1y_sq.size(); ++iy1) {
	      for (size_t ix2=0; ix2<sig2x_sq.size(); ++ix2) {
		for (size_t iy2=0; iy2<sig2y_sq.size(); ++iy2) {
		  E += w1x[ix1] * w1y[iy1] * w2x[ix2] * w2y[iy2] *
		    E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y,
								       sqrt(sig1x_sq[ix1] + sig2x_sq[ix2]),
								       sqrt(sig1y_sq[iy1] + sig2y_sq[iy2]));
		}
	      }
	    }
	  }
	  return E;
	};
    }
    density2_normalization = 0.5 * w2x[0] * w2y[0] / M_PI / sig2x[0] / sig2y[0];
  }
};

int main(int argc, char** argv) {
  // -------------------- Across 2D X-Y plane with weights --------------------
  // Our integral is 4-dimensional: over X,X',Y,Y'. But, one can integrate over
  // X-X', Y-Y' "circles" by averaging over turns, as one point "rotates" and
  // "swipes" those circles. So, one can integrate only over "circle" radii (rX,
  // rY). This is done below using 2-dimensional rX, rY Gaussian distribution,
  // it is covered by the grid of points inside the circle with R = 5 sigma. The
  // 2D Gaussian weight is assigned to every point, then the integral is
  // calculated as the sum of (2nd bunch densities * weight) over the sample of
  // "particles" starting from this grid with random phases and "turning" many
  // times.
  string config_file;
  if (argc != 2) {
    cout << "Usage: <program name> <configuration file>\n";
    return 0;
  }
  ifstream conf_file(argv[1]);
  Config c(conf_file);
  string output_dir;
  if (c.defined("output.directory")) {
    output_dir = c.s("output.directory");
  } else {
    string inp(argv[1]); // drop .txt from the end of configuration file
    size_t n = inp.size();
    if (inp.substr(n-4) == ".txt") {
      output_dir = inp.substr(0, n-4);
    } else {
      cerr << "The name of the configuration file must have the "
	   << "form *.txt (" << inp << " is given)\n";
      return 1;
    }
  }
  if( mkdir(output_dir.c_str(), 0755) ) {
    if ( errno == EEXIST ) {
      cerr << "Output directory " << output_dir << " already exists. Terminating\n";
    } else {
      cerr << "Can not create output directory " << output_dir << ", errno = " << errno << endl;
    }
    return 1;
  }
  cout << "Output directory: " << output_dir << endl;

  double Qx = c["Qx"];
  double Qy = c["Qy"];
  if (!c.defined("exact.Qxy")) {
    // Qx,y are normally given with 2-3 digits after comma. So, after 100 or
    // 1000 turns the points return to their original positions (as 100* or
    // 1000*Qx,y becomes integer). To avoid this and add extra randomness, add
    // "negligible" irrational numbers exp(-9) = 0.0001234098 and sqrt(2)/1e4
    // = 0.0001414214. This also makes Qx/Qy irrational,and opens otherwise
    // closed Lissajous figures of betatron oscillations in X-Y plane.  In
    // reality Qx,y should be irrational.
    Qx += sqrt(2) / 1e4;
    Qy += exp(-9);
    //
    // If you do not want to add anything to Qx, Qy, please, set
    //   exact.Qxy = TRUE
    // in the configuration file
    //
  }
  complex<double> zQx = exp(2i * M_PI * Qx);
  complex<double> zQy = exp(2i * M_PI * Qy);

  int N_steps = c.vd("x2").size();
  vector<double> x2 = c.vd("x2"), y2 = c.vd("y2");
  if (x2.size() != y2.size()) {
    cerr << "Number of given x2 and y2 coordinates are different\n";
    return 1;
  }
  // if "scale_x_by" parameter is given, multiply x2 by that, same for y2
  if (c.defined("scale_x_by"))
    for_each(x2.begin(), x2.end(), [scale = c["scale_x_by"]] (double& x) { x *= scale; });
  if (c.defined("scale_y_by"))
    for_each(y2.begin(), y2.end(), [scale = c["scale_y_by"]] (double& y) { y *= scale; });
  vector<complex<double> > z2(N_steps);
  transform(x2.begin(), x2.end(), y2.begin(), z2.begin(),
	    [](double x, double y) { return complex<double>(x, y); });
  Mutli_XY_Gaussian_bunches multi_g(c);
  // kick
  double alpha = 1. / 137.035;
  double hbar = 0.197327e-15; // in Gev * m
  double beta0 = 1; // velocity/c of the second bunch in the reference frame of the first
  double kick_const = 2 * c("Z1") * c("Z2") * alpha * hbar * c["N2"] /
    beta0 / c["p"] * c["beta"] * 1e12; // in um^2
  //
  int N_points = c("N.points");
  // Phase 1: no beam-beam
  int N_turns_no_bb = c("N.no.beam.beam.turns");
  // Phase 2: adiabatic switch on of the beam-beam force
  int N_transitional_turns =
    c.defined("N.transitional.turns") ? c("N.transitional.turns") : 1000;
  // Phase 3: stabilization
  // Normally, the luminosity is stable in less than 2000 turns (default)
  // after switching beam-beam ON, and much less if the switch was adiabatic
  int N_stabilization_turns =
    c.defined("N.stabilization.turns") ? c("N.stabilization.turns") : 2000;
  // Phase 4: final calculation of the luminosity
  int N_turns_with_beam_beam =
    c.defined("N.turns.with.beam.beam") ? c("N.turns.with.beam.beam") : 5000;
  const int N_phases = 4;
  int N_turns_in_phase[N_phases] = {N_turns_no_bb,
				    N_transitional_turns,
				    N_stabilization_turns,
				    N_turns_with_beam_beam};
  int N_turns = N_turns_no_bb + N_transitional_turns +
    N_stabilization_turns + N_turns_with_beam_beam;
  enum {PRECISE, PRECISE_MINUS_AVERAGE, AVERAGE} kick_model = PRECISE;
  if (c.defined("kick.model") && c.s("kick.model") != "precise") {
    if (c.s("kick.model") == "precise.minus.average") {
      kick_model = PRECISE_MINUS_AVERAGE;
    } else if (c.s("kick.model") == "average") {
      kick_model = AVERAGE;
    } else {
      cerr << "Wrong kick.model\n";
      return 1;
    }
  }
  //
  bool calculate_field_without_interpolation =
    c.defined("calculate.field.without.interpolation") &&
    c.s("calculate.field.without.interpolation") == "TRUE";
  //
  if (c.defined("seed")) {
    srand48(c("seed"));
  } else {
    srand48(time(NULL));
  }
  //
  enum OUTPUT {RX_RY_WEIGHTS=0, INTEGRALS=1, CENTERS=2,
	       INTEGRALS_PER_CIRCLE=3, CENTERS_PER_CIRCLE=4, POINTS=5, N_OUTPUT=6};
  struct Output{
    ogzstream os;
    bool selected;
    Output() : selected(false) {}
  } output[N_OUTPUT];
  if (c.defined("output")) {
    map<string, OUTPUT> m;
    m["rx_ry_weights"] = RX_RY_WEIGHTS;
    m["integrals"] = INTEGRALS;
    m["centers"] = CENTERS;
    m["integrals_per_circle"] = INTEGRALS_PER_CIRCLE;
    m["centers_per_circle"] = CENTERS_PER_CIRCLE;
    m["points"] = POINTS;
    for (const string& name : c.vs("output")) {
      string name_ = name;
      replace(name_.begin(), name_.end(), '.', '_');
      auto it = m.find(name_);
      if (it != m.end()) {
	OUTPUT o = it->second;
	output[o].os.open((output_dir + "/" + name_ + ".txt.gz").c_str());
	output[o].selected = true;
      } else {
	cerr << "Unrecognized output option given: " << name << ". Terminating" << endl;
	return 1;
      }	     
    }
  }
  ofstream output_summary(output_dir + "/summary.txt");
  // store xZ, yZ once per "select_turns" turns
  int select_turns = c.defined("select.one.turn.out.of") ? c("select.one.turn.out.of") : 0;
  int N_points_in_step = 0; // will be calculated later and will be <= than N_points
  // allocate N_points but only N_points_in_step <= N_points will be used in
  // rx, ry, w:
  vector<double> rx(N_points), ry(N_points), w(N_points);
  {
    // sqrt(N.points) will be used to sample intervals
    //
    // rX = [0...cut_x], rY = [0...cut_y], where cut_x,y were chosen to
    // reproduce the efficiency of "N.sigma" configuration cut for a single
    // Gaussian bunch. Ie. the areas of the rejected multi-Gaussian tails ==
    // the area of a single Gaussian outside +/-N.sigma.
    //
    int N_sqrt = int(sqrt(N_points));
    vector<double> rX_bins(N_sqrt), rY_bins(N_sqrt);
    {
      double
	binx = multi_g.cut_x / N_sqrt,
	biny = multi_g.cut_y / N_sqrt;
      for (int i=0; i<N_sqrt; ++i) {
	rX_bins[i] = (i + 0.5) * binx;
	rY_bins[i] = (i + 0.5) * biny;
      }
      // To reduce the number of simulated points and CPU, select only the
      // points from the ellipse inscribed inside the rectangle [0...cut_x] X
      // [0...cut_y] in X and Y.
      double w_sum = 0;
      for (int ix=0; ix<N_sqrt; ++ix) {
	for (int iy=0; iy<N_sqrt; ++iy) {
	  complex<double> z1(rX_bins[ix]/multi_g.cut_x,
			     rY_bins[iy]/multi_g.cut_y);
	  if (norm(z1) <=1 ) {
	    w[N_points_in_step] =
	      multi_g.not_normalized_rx_density1(rX_bins[ix]) *
	      multi_g.not_normalized_ry_density1(rY_bins[iy]);
	    w_sum += w[N_points_in_step];
	    rx[N_points_in_step] = rX_bins[ix];
	    ry[N_points_in_step] = rY_bins[iy];
	    ++N_points_in_step;
	  }
	}
      }
    }
    // normalize w
    double w_sum = accumulate(w.begin(), w.begin() + N_points_in_step, 0.);
    for (int i = 0; i < N_points_in_step; ++i) w[i] = w[i] / w_sum;
  }
  if (output[RX_RY_WEIGHTS].selected) {
    for (size_t i = 0; i < N_points_in_step; ++i) {
      output[RX_RY_WEIGHTS].os << i << " "
			       << rx[i] << " "
			       << ry[i] << " "
			       << w[i]<< '\n';
    }
  }
  //
  Bilinear_interpolator<complex<double> > E_field;
  int N_random_points_to_probe_field_map =
    c.defined("N.random.points.to.probe.field.map") ? c("N.random.points.to.probe.field.map") : 0;
  if (!calculate_field_without_interpolation) {
    cout << "Precalculating the field map ... " << flush;
    // Find X-Y rectangle of the field map in the frame where the 2d (field "source") bunch is at (0,0).
    //
    // Ranges of the 2d bunch coordinates in the frame where the 1st is at (0,0):
    auto x_range = minmax_element(x2.begin(), x2.end());
    auto y_range = minmax_element(y2.begin(), y2.end());
    double
      xmin = -1.1 * multi_g.cut_x - *x_range.second, // *x/y_range.second/.first == max/min x/y-separation;
      xmax =  1.1 * multi_g.cut_x - *x_range.first,  // *x/y.range is subtracted since the field map is
      ymin = -1.1 * multi_g.cut_y - *y_range.second, // in the frame where the 2d bunch is at (0,0), not the 1st;
      ymax =  1.1 * multi_g.cut_y - *y_range.first;  // 10% margin (1.1) since beam-beam can increase the radius
    //
    size_t n_cells_along_grid_side = c.defined("N.cells.along.grid.side") ? c("N.cells.along.grid.side") : 500;
    E_field.create(xmin ,xmax, ymin, ymax, n_cells_along_grid_side, multi_g.field2);
    cout << "done" << endl;
  }
  //
  // print Config when all type conversions are resolved
  cout << c;
  // main loop
  for (int step = 0; step < N_steps; ++step) {
    struct Summary {
      vector<double> integ;
      vector<complex<double> > avr_z;
      double int0_analytic, int0, int0_to_analytic, int0_rel_err,
	int_to_int0_correction;
      vector<double> integ_per_circle[N_phases]; // initialized by zeros as any vector of doubles
      vector<complex<double> > avr_z_per_circle[N_phases]; // initialized by zeros as any vector of doubles
      complex<double> z, analytic_kick, approx_analytic_z;      
      // Unfortunately, at the moment g++ 8.3 does not support atomic<double> += operations.
      // Therefore, two accumulators, integ and avr_z are expressed via atomic<long int>.
      // This is much faster than using mutex locks in threads.
      vector<atomic<long int> > atomic_integ;
      vector<atomic<long int> > atomic_avr_x, atomic_avr_y;
      // Long integer atomic_*'s are used in threads and then converted back
      // to doubles. To not loose precision, the added doubles are first
      // multiplied by atomic_*_converter below (large number) and converted
      // to long int.  Then, when all threads finish, atomic_*'s are divided
      // by atomic_*_converter.
      double atomic_integral_converter;
      double atomic_avr_xy_converter;
      //
      Summary(int n_turns, int N_points_in_step)
	: atomic_integ(n_turns), atomic_avr_x(n_turns), atomic_avr_y(n_turns),
	  integ(n_turns), avr_z(n_turns)
	{
	  for (int i=0; i<N_phases; ++i) {
	    integ_per_circle[i].resize(N_points_in_step);
	    avr_z_per_circle[i].resize(N_points_in_step);
	  }
	  for (size_t i=0; i<atomic_integ.size(); ++i) {
	    // contrary to other vectors, atomic vectors in current C++ are
	    // not initialized by zeros:
	    // http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2018/p0883r0.pdf
	    atomic_integ[i].store(0);
	    atomic_avr_x[i].store(0);
	    atomic_avr_y[i].store(0);
	  }
	}
    };
    Summary summary(N_turns, N_points_in_step);
    { // Fill summary.atomic_*_converters:
      //
      // they should be large enough to keep precision but small enough to not
      // exceed the long int range
      double max_w = *max_element(w.begin(), w.end());
      // integral contributions in the threads will be limited by max(w[i]) * max_w,
      // where max(w[i]) should be less than multi_g.not_normalized_density2(0,0).
      // If max(w[i]) * max_w < 1: take it to be 1:
      summary.atomic_integral_converter = LONG_MAX /
	max((multi_g.not_normalized_density2(0,0) * max_w), 1.) / 1.0001;
      // 1.0001 makes sure that LONG_MAX is not exceeded even with rounding errors
      //
      // For maximal avr_x,y values take maximal initially simulated radii,
      // max(cut_x,y), times 100 for safety (initial circles are distorted by
      // beam-beam, so maximal distance from 0,0 can grow a little, but not by
      // 100 in average).
      summary.atomic_avr_xy_converter = LONG_MAX / (max(multi_g.cut_x, multi_g.cut_y) * 100);
      //
      // cout << "atomic integral conv " << summary.atomic_integral_converter/1e20 << endl;
      // cout << "atomic avr xy conv " << summary.atomic_avr_xy_converter/1e20 << endl;
    }
    cout << "Processing beam separation (" << z2[step] << ") ... " << flush;
    // average field when the first bunch is at -z2[step] and the second is at
    // (0,0):
    complex<double> kick_average = kick_const *
      multi_g.field2_averaged_over_bunch_1(-real(z2[step]), -imag(z2[step]));
    //
    vector<double> kick_transition_weighted(N_transitional_turns);
    vector<complex<double> > kick_average_transition_weighted(N_transitional_turns);
    
    for (int i_turn = 0; i_turn < N_transitional_turns; ++i_turn) {
      double trans_w = double(i_turn + 1) / N_transitional_turns;
      // linearly increasing from 1/N_transitional_turns for i_turn=0,
      //                       to 1 for i_turn = N_transitional_turns-1
      kick_transition_weighted[i_turn] = kick_const * trans_w;
      kick_average_transition_weighted[i_turn] = kick_average * trans_w;
    }
    complex<double> z2_step = z2[step];
    // prepare N_points_in_step threads to run in parallel
    vector<thread> threads(N_points_in_step);
    // store selected points in
    // v_points[ 2 * (N_turns_stored * i + i_turn) + 0/1 for x/y ].
    //
    // This vector will be filled by all threads in parallel and when all
    // finish, will be written to file. Unfortunately, direct writing to the
    // same file by many threads would block them, this dictates storing large
    // v_points in memory.
    vector<complex<double> > v_points;
    int N_turns_stored = N_turns / select_turns;
    if (output[POINTS].selected) v_points.resize(N_points_in_step * N_turns_stored * 2);
    mutex turn_loop_mutex;
    auto loop = [&rx, &ry, &w, &multi_g, z2_step, &E_field, &summary, &v_points,
		 step, select_turns, zQx, zQy,
		 &kick_transition_weighted, kick_average, kick_const, &kick_average_transition_weighted,
		 N_turns, N_turns_stored, N_turns_in_phase, &output,
		 kick_model, calculate_field_without_interpolation]
      (int i, double rnd1, double rnd2) {
		  // For every (rx,ry) pair obtain X-X', Y-Y' 4-dimensional
		  // coordinates by random rotation around X-X' and Y-Y' circles
		 complex<double>
		   xZ = rx[i] * exp(2i * M_PI * rnd1),
		   yZ = ry[i] * exp(2i * M_PI * rnd2);
		 auto v_points_iter = v_points.begin();
		 if (output[POINTS].selected) v_points_iter += 2 * N_turns_stored * i;
		 //
		 auto apply_kick = [zQx, zQy]
		   (complex<double> kick, complex<double>& xZ, complex<double>& yZ) {
		     // The momentum kick is subtracted below because of the minus sign in
		     // the definition of eg. zX = X - iX'. The kick is added to X',
		     // so that i*kick is subtracted from zX.
				     xZ = (xZ - 1i * real( kick )) * zQx;
				     yZ = (yZ - 1i * imag( kick )) * zQy;
				   };
		 auto kick =
		   [kick_model, z2_step, calculate_field_without_interpolation,
		    &E_field, kick_const, kick_average] (complex<double> z_sep, double k_const,
							 complex<double> k_average) {
		     switch (kick_model) {
		     case PRECISE:
		       if (calculate_field_without_interpolation) {
			 return k_const * E_field.f(real(z_sep), imag(z_sep));
		       } else {
			 return k_const * E_field  (real(z_sep), imag(z_sep));
		       }
		     case AVERAGE: return k_average;
		     default: //  PRECISE_MINUS_AVERAGE: 
		       if (calculate_field_without_interpolation) {
			 return k_const * E_field.f(real(z_sep), imag(z_sep)) - k_average;
		       } else {
			 return k_const * E_field  (real(z_sep), imag(z_sep)) - k_average;
		       }
		     }
		   };		     
		 //
		 auto fill_summary_per_turn =
		   [&w, &multi_g, &summary, i, z2_step,
		    select_turns, &output]
		   (int i_turn, complex<double> xZ, complex<double> yZ,
		    vector<complex<double> >::iterator& v_points_iter, int phase)
		   {
		     // overlap integral is calculated as a sum of (2nd bunch profile
		     // density * weight) over the sample of 1st bunch "particles".
		     //
		     // Calculate rho2 at the particle's position
		     // w.r.t. bunch 2 center:
		     //  new_z1 - z2_step = real( xZ ) + 1i*real( yZ ) - z2_step
		     double rho2 = multi_g.not_normalized_density2(real(xZ) - real(z2_step),
								   real(yZ) - imag(z2_step));
		     summary.integ_per_circle[phase][i] += rho2;
		     if (output[CENTERS_PER_CIRCLE].selected) {
		       summary.avr_z_per_circle[phase][i] += complex<double>(real(xZ), real(yZ));
		     }
		     if (output[INTEGRALS].selected) {
		       summary.atomic_integ[i_turn] +=
			 (long int)(rho2 * w[i] * summary.atomic_integral_converter);
		     }
		     if (output[CENTERS].selected) {
		       summary.atomic_avr_x[i_turn] +=
			 (long int)(real(xZ) * w[i] * summary.atomic_avr_xy_converter);
		       summary.atomic_avr_y[i_turn] +=
			 (long int)(real(yZ) * w[i] * summary.atomic_avr_xy_converter);
		     }
		     //
		     if (output[POINTS].selected && (i_turn + 1) % select_turns == 0) {
		       // last i_turn = N_turn-1 is written, first i_turn=0 - not
		       // (assuming N_turn is a multiple of select_runs)
		       //	cout << "  Turn " << i_turn+1 << endl;
		       *(v_points_iter++) = xZ;
		       *(v_points_iter++) = yZ;
		     }
		   };
		 int i_turn = 0;
		 // first N_turns_no_bb turns without kick to measure the initial, possibly
		 // biased, integral; then calculate how much the kick changes
		 // it. Taking the ratio to the initial integral (instead of the exact
		 // integral) should cancel the bias at least partially and improve
		 // precision
		 for (int phase_turn=0; phase_turn<N_turns_in_phase[0]; ++phase_turn, ++i_turn) {
		   complex<double> z1( real( xZ ), real( yZ ));
		   // coordinate in a frame where the 2nd bunch center is at (0,0):
		   xZ *= zQx; // kick = 0
		   yZ *= zQy;
		   // calculate and store integrals and points here, in the
		   // end, after the propagation through the ring (and after
		   // the beam-beam). The integral and the points (both of
		   // which can be written out), therefore, correspond to each
		   // other. One could do it before, then the last propagation
		   // through the ring for i_turn = N_turn-1 could be
		   // dropped. However, in this case the first integral+points
		   // always corresponded to the case without beam-beam, even
		   // if N_turns_no_bb=0. So, for this reason and for simplicity the
		   // above logic is chosen.
		   fill_summary_per_turn(i_turn, xZ, yZ, v_points_iter, 0);
		 }
		 for (int phase_turn=0; phase_turn < N_turns_in_phase[1]; ++phase_turn, ++i_turn) {
		   complex<double> z_sep = complex<double>( real( xZ ), real( yZ )) - z2_step;
		   complex<double> k = kick(z_sep,
					    kick_transition_weighted[phase_turn],
					    kick_average_transition_weighted[phase_turn]);
		   apply_kick( k, xZ, yZ );
		   fill_summary_per_turn(i_turn, xZ, yZ, v_points_iter, 1);
		 }
		 for (int phase_turn=0; phase_turn < N_turns_in_phase[2]; ++phase_turn, ++i_turn) {
		   complex<double> z_sep = complex<double>( real( xZ ), real( yZ )) - z2_step;
		   complex<double> k = kick(z_sep,
					    kick_const,
					    kick_average);
		   apply_kick( k, xZ, yZ );
		   fill_summary_per_turn(i_turn, xZ, yZ, v_points_iter, 2);
		 }
		 for (int phase_turn=0; phase_turn < N_turns_in_phase[3]; ++phase_turn, ++i_turn) {
		   complex<double> z_sep = complex<double>( real( xZ ), real( yZ )) - z2_step;
		   complex<double> k = kick(z_sep,
					    kick_const,
					    kick_average);
		   apply_kick( k, xZ, yZ );
		   fill_summary_per_turn(i_turn, xZ, yZ, v_points_iter, 3);
		 }
		 // w[i] is common to all contributions to integ/avr_z _per_circle:
		 for (int phase = 0; phase < N_phases; ++phase) {
		   double norm = multi_g.density2_normalization / double(N_turns_in_phase[phase]);
		   summary.integ_per_circle[phase][i] *= norm;
		   summary.avr_z_per_circle[phase][i] /= double(N_turns_in_phase[phase]);
		 }
		}; // end of loop over tracked particle-points
    for (int i = 0; i < N_points_in_step; ++i) {
      // supply random numbers here so that they are reproducible: if
      // drand48() were called in threads, it would depend on the thread's
      // execution order (which is not defined)
      threads[i] = thread(loop, i, drand48(), drand48());
    }
    for (int i = 0; i < N_points_in_step; ++i) threads[i].join();
    // 
    if (output[POINTS].selected) {
      for (int i = 0; i < N_points_in_step; ++i) {
	int ind = 2 * N_turns_stored * i;
	for (int i_turn = 0; i_turn < N_turns_stored; ++i_turn) {
	  output[POINTS].os << step << " "
			    << i_turn << " "
			    << i << " "
			    << v_points[ind++] << " " << v_points[ind++] << '\n';
	}
      }
    }
    // calculate and store summary
    for (int i_turn = 0; i_turn < N_turns; ++i_turn) {
      summary.integ[i_turn] = double(summary.atomic_integ[i_turn]) / summary.atomic_integral_converter;
      summary.avr_z[i_turn] = complex<double>(double(summary.atomic_avr_x[i_turn]) / summary.atomic_avr_xy_converter,
					      double(summary.atomic_avr_y[i_turn]) / summary.atomic_avr_xy_converter);
    }    
    //
    for_each(summary.integ.begin(),
    	     summary.integ.end(), [n = multi_g.density2_normalization](double& x) { x *= n; });
    if (output[INTEGRALS_PER_CIRCLE].selected) {
      for (int phase=0; phase<N_phases; ++phase) {
	for (size_t i = 0; i < summary.integ_per_circle[phase].size(); ++i) {
	  output[INTEGRALS_PER_CIRCLE].os << step << " "
					  << z2_step << " "
					  << phase << " "
					  << i << " "
					  << summary.integ_per_circle[phase][i] << '\n';
	}
      }
    }
    if (output[CENTERS_PER_CIRCLE].selected) {
      for (int phase=0; phase<N_phases; ++phase) {
	for (size_t i = 0; i < summary.integ_per_circle[phase].size(); ++i) {
    	output[CENTERS_PER_CIRCLE].os << step << " "
				      << z2_step << " "
				      << phase << " "
				      << i << " "
				      << summary.avr_z_per_circle[phase][i] << '\n';
	}
      }
    }
    if (output[INTEGRALS].selected) {
      for (int i_turn = 0; i_turn < N_turns; ++i_turn) {
    	output[INTEGRALS].os << step << " "
			     << z2_step << " "
			     << i_turn << " "
			     << summary.integ[i_turn] << '\n';
      }
    }
    if (output[CENTERS].selected) {
      for (int i_turn = 0; i_turn < N_turns; ++i_turn) {
    	output[CENTERS].os << step << " "
			   << z2_step << " "
			   << i_turn << " "
			   << summary.avr_z[i_turn] << "\n";
      }
    }
    summary.int0_analytic = multi_g.overlap_integral(real(z2_step), imag(z2_step));
    //    summary.int0 = accumulate(summary.integ.begin(),
    //			      summary.integ.begin()+N_turns_no_bb, 0.) / double(N_turns_no_bb);
    summary.int0 = 0;
      for (int i=0; i<N_points_in_step; ++i) summary.int0 += summary.integ_per_circle[0][i] * w[i];
    summary.int0_to_analytic = summary.int0 / summary.int0_analytic;
    if (output[INTEGRALS].selected) {
      double sd_int0 = 0;
      for (auto i = summary.integ.begin(); i != summary.integ.begin() + N_turns_no_bb; ++i) {
	double dx = *i - summary.int0;
	sd_int0 += dx * dx;
      }    
      summary.int0_rel_err = sqrt(sd_int0 / (N_turns_no_bb - 1.) / double(N_turns_no_bb)) / summary.int0;
    } else summary.int0_rel_err = -9999;
    // assume that after N_stabilization_turns the distribution stabilizes
    //    summary.int_to_int0_correction =
    //      accumulate(summary.integ.end() - N_turns_with_beam_beam,
    //		 summary.integ.end(), 0.) / N_turns_with_beam_beam / summary.int0;
    summary.int_to_int0_correction = 0;
    for (int i=0; i<N_points_in_step; ++i) summary.int_to_int0_correction += summary.integ_per_circle[3][i] * w[i];
    summary.int_to_int0_correction /= summary.int0;
    //    summary.z =
    //      accumulate(summary.avr_z.end() - N_turns_with_beam_beam,
    //		 summary.avr_z.end(), complex<double>(0, 0)) / double(N_turns_with_beam_beam);
    summary.z = 0;
    for (int i=0; i<N_points_in_step; ++i) summary.z += summary.avr_z_per_circle[3][i] * w[i];
    // the sign of analytic kick below (printed in summary) is not inverted,
    // ie. it is for X', not for z = X - iX'
    //
    // note, it has units of X', ie. of length; to convert to angles, one
    // needs to divide by beta*
    //
    // 0-z2_step is the 1st bunch center when the 2nd is at (0,0),
    // vdm_sigx,y: because the bunch-bunch force is derived from the
    // particle-bunch by extra smearing according to sig1x,y, so in the
    // end sig=sqrt(sig1^2 + sig2^2)
    //
    // kick_const is halved because of the accelerator dynamics: after
    // many turns the dynamics is such that every time the scaled angle =
    // (angle / E_from_unit...()) is kicked from initial -kick_const/2 to
    // final +kick_const/2, and their difference is kick_const.
    //
    summary.analytic_kick = kick_average / 2.;
    // the X shift corresponding to X' kick
    summary.approx_analytic_z =
      complex<double>(real(summary.analytic_kick) / tan(M_PI * Qx),
		      imag(summary.analytic_kick) / tan(M_PI * Qy));
    output_summary << step << " "
		   << z2_step << " "
		   << summary.int_to_int0_correction << " "
		   << summary.int0_analytic << " "
		   << summary.int0 << " "
		   << summary.int0_to_analytic << " "
		   << summary.int0_rel_err << " "
		   << summary.z << " "
		   << summary.analytic_kick << " "
		   << summary.approx_analytic_z << "\n";
    cout << "done" << endl;
  } // end of loop over steps
    //
  if (!calculate_field_without_interpolation && N_random_points_to_probe_field_map > 0) {
    // check the precision of the interpolation by measuring the absolute
    // mismatches with the exact field values at randomly distributed points
    // inside the interpolation grid rectangle
    double max_interpolation_discrepancy = 0, avr_interpolation_discrepancy = 0;
    double max_field = 0;
    {
      auto grid = E_field.grid();
      double
	dx = grid.step_x * grid.N,
	dy = grid.step_y * grid.N;
      for (int i=0; i<N_random_points_to_probe_field_map; ++i) {
	double x = grid.xmin + drand48() * dx;
	double y = grid.ymin + drand48() * dy;
	complex<double>
	  precise = E_field.f(x, y),
	  approx  = E_field(x, y);
	double mismatch = abs(approx - precise);
	max_field       = max(max_field, abs(precise));
	max_interpolation_discrepancy = max(max_interpolation_discrepancy, mismatch);
	avr_interpolation_discrepancy += mismatch;
      }
      avr_interpolation_discrepancy /= N_random_points_to_probe_field_map;
    }
    cout << "Max and average |E field - interpolated| mismatches normalized to max |E| = "
	 << max_interpolation_discrepancy / max_field << ", "
	 << avr_interpolation_discrepancy / max_field << "\n\n";
  }
  cout << "Final results are written to " << output_dir << " directory\n"
       << "The summary appears in \"summary.txt\".\n"
       << "The 4th column \"int_to_int0_correction\" in this file contains the\n"
       << "final correction, i.e. the ratio of the luminosities with and\n"
       << "without beam-beam electromagnetic interaction, for the scan step and\n"
       << "the corresponding beam offset specified in the first three columns.\n\n"
    //
       << "The columns 5-8 refer to no-beam-beam case: : the analytic and numeric\n"
       << "integrals, their ratio and the estimated relative statistical\n"
       << "error of the numeric integral.\n\n"
    //
       << "The (9,10) column pair gives the numerically calculated x-y shift of\n"
       << "of the first bunch, while (11,12) and (13,14) pairs give the\n"
       << "analitically calculated values, the first in the coordinates\n"
       << "beta*(dx/dz, dy/dz) and the second - in (x,y) under the assumption\n"
       << "that its ratio to the first pair is -1/tan(pi*Q).\n";
  {
    ofstream config_out((output_dir + "/config.txt").c_str());
    config_out << c;
  }
  return(0);
}
