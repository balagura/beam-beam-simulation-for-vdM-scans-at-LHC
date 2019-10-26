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
#include <memory> // for shared_ptr
#include <gzstream.h>
#include <Faddeeva.hh>
#include <config.hh>
#include <E_field.hh>
#include <bilinear_interpolator.hh>

using namespace std;
using namespace std::complex_literals;

ostream& operator<<(ostream& os, complex<double> c) { return os << real(c) << " " << imag(c); }

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
  // if "scale_x_by" parameter is given, multiply x2 by that, same for y2
  if (c.defined("scale_x_by"))
    for_each(x2.begin(), x2.end(), [scale = c["scale_x_by"]] (double& x) { x *= scale; });
  if (c.defined("scale_y_by"))
    for_each(y2.begin(), y2.end(), [scale = c["scale_y_by"]] (double& y) { y *= scale; });
  vector<complex<double> > z2(N_steps);
  transform(x2.begin(), x2.end(), y2.begin(), z2.begin(),
	    [](double x, double y) { return complex<double>(x, y); });
  double sig1x = c["sig1.x"];
  double sig2x = c["sig2.x"];
  double sig1y = c["sig1.y"];
  double sig2y = c["sig2.y"];
  double sig1x_sq = sig1x * sig1x, sig1y_sq = sig1y * sig1y;
  double sig2x_sq = sig2x * sig2x, sig2y_sq = sig2y * sig2y;
  double vdm_sigx_sq = sig1x_sq + sig2x_sq;
  double vdm_sigy_sq = sig1y_sq + sig2y_sq;
  double vdm_sigx = sqrt(vdm_sigx_sq);
  double vdm_sigy = sqrt(vdm_sigy_sq);
  // points of the first bunch beyond n_sig_cut ("N.sigma.cut" in
  // the configuration file) from the center are not simulated
  double n_sig_cut = c.defined("N.sigma.cut") ? c("N.sigma.cut") : 5;
  double sig_sq_limit = n_sig_cut * n_sig_cut;
  double overlap_integral_const = 0.5 / M_PI / sig2x / sig2y;
  // kick
  double alpha = 1. / 137.035;
  double hbar = 0.197327e-15; // in Gev * m
  double beta0 = 1; // velocity/c of the second bunch in the reference frame of the first
  double kick_const = 2 * c("Z1") * c("Z2") * alpha * hbar * c["N2"] /
    beta0 / c["p"] * c["beta"] * 1e12; // in um^2
  //
  int N_points = c("N.points");
  int N_points_in_step; // less than N_points
  int N_no_bb = c("N.no.beam.beam.turns");
  // Normally, the luminosity is stable in less than 2000 turns (default)
  // after switching beam-beam ON
  int N_stabilization_turns =
    c.defined("N.stabilization.turns") ? c("N.stabilization.turns") : 2000;
  if (c.vd("x2").size() != c.vd("y2").size()) {
    cerr << "Number of given x2 and y2 coordinates are different\n";
    return 1;
  }
  int N_transitional_turns =
    c.defined("N.transitional.turns") ? c("N.transitional.turns") : 1000;
  // Total number of turns with stabilized beam-beam
  int N_turns_with_beam_beam =
    c.defined("N.turns.with.beam.beam") ? c("N.turns.with.beam.beam") : 5000;
  int N_turns = N_no_bb + N_stabilization_turns +
    N_transitional_turns + N_turns_with_beam_beam;
  //
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
  struct Summary {
    vector<double> integ;
    vector<complex<double> > avr_z;
    double int0_analytic, int0, int0_to_analytic, int0_rel_err,
      int_to_int0_correction;
    complex<double> z, analytic_kick, approx_analytic_z;      
    Summary(int n_turns) : integ(n_turns), avr_z(n_turns) {}
    Summary() {}
  };
  //
  ogzstream output_rx_ry_weights, output_integrals, output_centers, output_points;
  ofstream output_summary(output_dir + "/summary.txt");
  bool output_rx_ry_weights_bool = false, output_integrals_bool = false,
    output_centers_bool = false, output_points_bool = false;
  for (auto i = c.vs("output").begin(); i != c.vs("output").end(); ++i) {
    if        (*i == "rx.ry.weights") {
      output_rx_ry_weights.open((output_dir + "/rx_ry_weights.txt.gz").c_str());
      output_rx_ry_weights_bool = true;
    } else if (*i == "integrals"    ) {
      output_integrals    .open((output_dir + "/integrals.txt.gz").c_str());
      output_integrals_bool = true;
    } else if (*i == "centers"      ) {
      output_centers      .open((output_dir + "/centers.txt.gz").c_str());
      output_centers_bool = true;
    } else if (*i == "points"       ) {
      output_points       .open((output_dir + "/points.txt.gz").c_str());
      output_points_bool = true;
    }
  }
  // store xZ, yZ once per "select_turns" turns
  int select_turns = c.defined("select.one.turn.out.of") ? c("select.one.turn.out.of") : 0;
  //
  // sqrt(N.points) will be used to sample rX and rY intervals
  // [0, n_sig_cut * sig1x/y]
  //
  int N_sqrt = int(sqrt(N_points));
  double bin_area;
  vector<double> rXY_bins(N_sqrt);
  { // bin_area and rXY_bins are in a scaled space with sig1x = sig1y = 1,
    // ie. not multipled by sig1x,y
    double bin = n_sig_cut / N_sqrt;
    bin_area = bin * bin;
    for (int i=0; i<N_sqrt; ++i) rXY_bins[i] = (i + 0.5) * bin;
  }
  // only N_points_in_step <= N_points will be used in rx, ry, w:
  vector<double> rx(N_points), ry(N_points), w(N_points);
  // initialization
  N_points_in_step = 0;
  for (int ix=0; ix<N_sqrt; ++ix) {
    for (int iy=0; iy<N_sqrt; ++iy) {
      complex<double> z1(rXY_bins[ix], rXY_bins[iy]);
      // simulate only points inside ellipse with (n_sig_cut * sig1x,
      // n_sig_cut * sig1y) half-axes, ie. inside R = n_sig_cut circle in
      // the scaled space with sig1x = sig1y = 1
      if (norm(z1) < sig_sq_limit) {
	rx[N_points_in_step] = real(z1);
	ry[N_points_in_step] = imag(z1);
	++N_points_in_step;
      }
    }
  }
  // distribution of radius length for two-dimensional Gaussian:
  //
  // exp(-r^2/2/sig^2) * r / sig^2 dr.
  //
  // So, the weight for two two-dimensional Gaussians for dr = rX, rY grid bin:
  {
    for (int i=0; i<N_points_in_step; ++i)
      w[i] = exp(-(rx[i] * rx[i] + ry[i] * ry[i]) / 2) * rx[i] * ry[i] * bin_area;
    double w_sum = accumulate(w.begin(),
			      w.begin() + N_points_in_step, 0.);
    for (int i = 0; i < N_points_in_step; ++i) {
      w[i] = w[i] / w_sum;
      // rescale X,Y so that round bunch with sig1x = sig1y = 1 is transformed into
      // (sig1x, sig1y) elliptical bunch
      rx[i] *= sig1x;
      ry[i] *= sig1y;
    }
  }
  // Note, the integral outside (n_sig_cut * sig1x, n_sig_cut * sig1y)
  // ellipse is exp(-n_sig_cut^2 / 2) which is 3.7e-6 for n_sig_cut = 5
  if (output_rx_ry_weights_bool) {
    for (size_t i = 0; i < N_points_in_step; ++i) {
      output_rx_ry_weights << i << " "
			   << rx[i] << " "
			   << ry[i] << " "
			   << w[i]<< '\n';
    }
  }
  const char* xyC = "xy";
  //
  Bilinear_interpolator<complex<double> > E_field;
  int N_random_points_to_probe_field_map =
    c.defined("N.random.points.to.probe.field.map") ? c("N.random.points.to.probe.field.map") : 0;
  // to create Bilinear_interpolator, it will need a function with x,y
  // arguments returning a field at x,y created by the 2d bunch centered at
  // (0,0). More specifically, Bilinear_interpolator expects
  // std::function<complex<double>(double, double)> or something convertible
  // to it (std::function is a wrapper around various function objects, for
  // details, please, see C++ documentation on std::function). Below it is
  // constructed from a lambda-function made from
  //
  //   E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y, sig2x, sig2y)
  //
  // by fixing last two arguments. Namely, since the arguments sig2x, sig2y
  // are not in the expected signature, they are transferred as "captured"
  // constants in square brackets [..]:
  //
  auto E_field_fun = [&sig2x, &sig2y](double x, double y) {
		       return E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(x, y,
										 sig2x, sig2y);
		     };
  // E_field_fun(x, y) now returns the field at x,y as required. This
  // lambda-function will be used both in Bilinear_interpolator and for the
  // direct field calculation if the interpolation is switched off
  if (!calculate_field_without_interpolation) {
    cout << "Precalculating the field map ... " << flush;
    // Find X-Y rectangle of the field map in the frame where the 2d (field "source") bunch is at (0,0).
    //
    // Ranges of the 2d bunch coordinates in the frame where the 1st is at (0,0):
    auto x_range = minmax_element(x2.begin(), x2.end());
    auto y_range = minmax_element(y2.begin(), y2.end());
    double
      xmin = -1.1 * n_sig_cut * sig1x - *x_range.second, // *x/y_range.second/.first == max/min x/y-separation;
      xmax =  1.1 * n_sig_cut * sig1x - *x_range.first,  // *x/y.range is subtracted since the field map is
      ymin = -1.1 * n_sig_cut * sig1y - *y_range.second, // in the frame where the 2d bunch is at (0,0), not the 1st;
      ymax =  1.1 * n_sig_cut * sig1y - *y_range.first;  // 10% margin (1.1) since beam-beam can increase the radius
    //
    size_t n_cells_along_grid_side = c.defined("N.cells.along.grid.side") ? c("N.cells.along.grid.side") : 500;
    E_field.create(xmin ,xmax, ymin, ymax, n_cells_along_grid_side, E_field_fun);
    // one can write function<complex<double>(double, double)>(E_field_fun)
    // instead of E_field_fun, but this is not necessary as the conversion is
    // implicit
    //    
    cout << "done" << endl;
  }
  //
  // print Config when all type conversions are resolved
  cout << c;

  // main loop
  for (int step = 0; step < N_steps; ++step) {
    Summary summary(N_turns);
    cout << "Processing beam separation (" << z2[step] << ") ... " << flush;
    // For every (rx,ry) pair obtain X-X', Y-Y' 4-dimensional
    // coordinates by random rotation around X and Y circles
    vector<complex<double> > xZ(N_points_in_step), yZ(N_points_in_step);
    for (int i = 0; i < N_points_in_step; ++i) {
      xZ[i] = rx[i] * exp(2i * M_PI * drand48());
      yZ[i] = ry[i] * exp(2i * M_PI * drand48());
    }
    complex<double> kick_average = kick_const *
      E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(-real(z2[step]), -imag(z2[step]),
							 vdm_sigx, vdm_sigy);
    //
    for (int i_turn = 0; i_turn < N_turns; ++i_turn) {
      double transitional_weight = 0;
      {
	int i_turn_w_bb = i_turn - (N_no_bb - 1); // zero for i_turn = 0 ... N_no_bb-1,
	// first positive kick for i_turn = N_no_bb
	if (i_turn_w_bb <= 0) {
	  transitional_weight = 0;
	} else if (i_turn_w_bb < N_transitional_turns) {
	  transitional_weight = double(i_turn_w_bb) / N_transitional_turns;
	} else {
	  transitional_weight = 1;
	}
      }
      double k_transition_weighted = kick_const * transitional_weight;
      complex<double>
	kick_average_transition_weighted = kick_average * transitional_weight;

      // always store integ, avr_z,
      // summary.avr_z[i_turn] and .integ[i_turn] are initialized by zeros
      // as any vector of doubles
      for (int i = 0; i < N_points_in_step; ++i) {
	complex<double> z1( real( xZ[i] ), real( yZ[i] ));
	summary.avr_z[i_turn] += z1 * w[i];
	// coordinate in a frame where the 2nd bunch center is at (0,0):
	complex<double> z_sep = z1 - z2[step];
	double r2 = norm(z_sep);
	double chi_x = real(z_sep) / sig2x;
	double chi_y = imag(z_sep) / sig2y;
	double e = exp(-0.5 * (chi_x * chi_x + chi_y * chi_y));
	// overlap integral is calculated as a sum of (2nd bunch profile
	// density * weight) over the sample of 1st bunch "particles"
	summary.integ[i_turn] += e * w[i];
	//
	if (i_turn < N_no_bb ||
	    // first N_no_bb turns without kick to measure the initial, possibly
	    // biased, integral; then calculate how much the kick changes
	    // it. Taking the ratio to the initial integral (instead of the exact
	    // integral) should cancel the bias at least partially and improve
	    // precision
	    r2 == 0.) { // kick = 0
	  xZ[i] *= zQx;
	  yZ[i] *= zQy;
	} else if (kick_model == PRECISE) {
	  complex<double> kick;
	  if (calculate_field_without_interpolation) {
	    kick = E_field_fun(real(z_sep), imag(z_sep));
	  } else {
	    kick = E_field(real(z_sep), imag(z_sep));
	  }
	  kick *= k_transition_weighted;
	  // The momentum kick is subtracted below because of the minus sign in
	  // the definition of eg. zX = X - iX'. The kick is added to X',
	  // so that i*kick is subtracted from zX.
	  xZ[i] = (xZ[i] - 1i * real( kick )) * zQx;
	  yZ[i] = (yZ[i] - 1i * imag( kick )) * zQy;
	} else if (kick_model == AVERAGE) {
	  // apply kick_average - the same for all z[i]
	  xZ[i] = (xZ[i] - 1i * real( kick_average_transition_weighted )) * zQx;
	  yZ[i] = (yZ[i] - 1i * imag( kick_average_transition_weighted )) * zQy;
	}
      } // end of loop over tracked particle-points
      //
      summary.integ[i_turn] *= overlap_integral_const;
      if (output_points_bool && (i_turn + 1) % select_turns == 0) {
	//	cout << "  Turn " << i_turn+1 << endl;
	// store all particle coordinates for this turn
	// last i_turn = N_turn-1 is written, first i_turn=0 - not
	// (assuming N_turn is a multiple of select_runs)
	for (size_t i = 0; i < N_points_in_step; ++i) {
	  output_points << step << " "
			<< i_turn << " "
			<< i << " "
			<< xZ[i] << " " << yZ[i] << '\n';
	}
      }
    } // end of loop over turns
    if (output_integrals_bool) {
      for (size_t turn = 0; turn < summary.integ.size(); ++turn) {
	output_integrals << step << " "
			 << z2[step] << " "
			 << turn << " "
			 << summary.integ[turn] << '\n';
      }
    }
    if (output_centers_bool) {
      for (size_t turn = 0; turn < summary.avr_z.size(); ++turn) {
	output_centers << step << " "
		       << z2[step] << " "
		       << turn << " "
		       << summary.avr_z[turn] << "\n";
      }
    }
    // store summary
    {
      double chi_x = real(z2[step]) / vdm_sigx;
      double chi_y = imag(z2[step]) / vdm_sigy;
      summary.int0_analytic = 0.5 / M_PI / vdm_sigx / vdm_sigy *
	exp(-0.5 * (chi_x * chi_x + chi_y * chi_y));
    }
    summary.int0 = accumulate(summary.integ.begin(),
			      summary.integ.begin()+N_no_bb, 0.) / double(N_no_bb);
    summary.int0_to_analytic = summary.int0 / summary.int0_analytic;
    double sd_int0 = 0;
    for (auto i = summary.integ.begin(); i != summary.integ.begin() + N_no_bb; ++i) {
      double dx = *i - summary.int0;
      sd_int0 += dx * dx;
    }    
    summary.int0_rel_err = sqrt(sd_int0 / (N_no_bb - 1.) / double(N_no_bb)) / summary.int0;
    // assume that after N_stabilization_turns the distribution stabilizes
    summary.int_to_int0_correction =
      accumulate(summary.integ.end() - N_turns_with_beam_beam,
		 summary.integ.end(), 0.) / N_turns_with_beam_beam / summary.int0;
    summary.z =
      accumulate(summary.avr_z.end() - N_turns_with_beam_beam,
		 summary.avr_z.end(), complex<double>(0, 0)) / double(N_turns_with_beam_beam);
    // the sign of analytic kick below (printed in summary) is not inverted,
    // ie. it is for X', not for z = X - iX'
    //
    // note, it has units of X', ie. of length; to convert to angles, one
    // needs to divide by beta*
    //
    // 0-z2[step] is the 1st bunch center when the 2nd is at (0,0),
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
		   << z2[step] << " "
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
