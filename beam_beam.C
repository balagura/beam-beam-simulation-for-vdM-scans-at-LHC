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
#include <gzstream.h>
#include <config.hh>
#include <multi_gaussian.hh>
#include <n.hh>
#include <output.hh>

using namespace std;
using namespace std::complex_literals;

ostream& operator<<(ostream& os, complex<double> c) { return os << real(c) << " " << imag(c); }
//
// Same as Mutli_XY_Gaussian_bunches but takes all parameters from Config& c
// and also "n_ip"
//
struct Config_Mutli_XY_Gaussian_bunches : public Mutli_XY_Gaussian_bunches {
  Config_Mutli_XY_Gaussian_bunches(Config& c, int n_ip) {
    Mutli_XY_Gaussian_bunches::reset_ip_number(n_ip);
    //
    double n_sig_cut = c.defined("N.sigma.cut") ? c("N.sigma.cut") : 5;
    // accumulate beam separations:
    vector<array<vector<double>, 2> > positions(n_ip); // [ip][coor][step]
    for (size_t coor=0; coor<2; ++coor) {
      const vector<double>& deltaQs = c.vd(string("kicker.") + "xy"[coor] + ".phase.advance");
      if (deltaQs.size() != n_ip) {
	cerr << "The number of given phase advances (" << deltaQs.size()
	     << ") is not equal to the number of IPs (" << n_ip << ")\n";
	exit(1);
      }
      for (int ip=0; ip<n_ip; ++ip) {
	double deltaQ = deltaQs[ip];
	// Phase advances or full tunes might be given with 2-3 digits after
	// comma. So, after 100 or 1000 turns the points return to their original
	// positions (as 100* or 1000*Qx,y becomes integer). To avoid this and add
	// extra randomness, a "negligible" irrational number randomly distributed
	// in the interval -exp(-8.5)...+exp(-8.5) = +/-0.0001017342 is added by
	// default to deltaQ. This makes it irrational and opens otherwise closed
	// Lissajous figures of betatron oscillations in X-Y plane. In reality the
	// phase advances should be irrational.
	//
	// If this is not desired, an option
	//   exact.phase.advances = TRUE
	// should be set in the configuration. Then all phase advances
	// from the configuration file are used "as is".
	//
	if (!c.defined("exact.phase.advances") ||
	    c.s("exact.phase.advances") != "TRUE") deltaQ += exp(-8.5) * (drand48() - 0.5);
	//
	string s = "kicker." + to_string(ip+1) + string(".") + "xy"[coor];
	vector<double> pos = c.vd(s);
	// if "scale_x_by" parameter is given, multiply x2 by that, same for y2
	if (c.defined(s + ".scale.by")) {
	  for_each(pos.begin(),
		   pos.end(), [scale = c[s + ".scale.by"]] (double& xy) { xy *= scale; });
	}
	positions[ip][coor].resize(pos.size());
	copy(pos.begin(), pos.end(), positions[ip][coor].begin());
	//
	vector<double> weights = c.defined(s + ".weight") ? c.vd(s + ".weight") : vector<double>();
	Mutli_XY_Gaussian_bunches::reset_kicker_bunch(ip, coor,
						      c.vd(s + ".sig"),
						      weights,
						      deltaQ);
      }
      string s = string("kicked.") + "xy"[coor];
      vector<double> weights = c.defined(s + ".weight") ? c.vd(s + ".weight") : vector<double>();
      Mutli_XY_Gaussian_bunches::reset_kicked_bunch(coor,
						    c.vd(s + ".sig"),
						    weights,
						    n_sig_cut);
    }
    Mutli_XY_Gaussian_bunches::reset_kicker_positions(positions);
    //
    for (size_t coor=0; coor<2; ++coor) {
      for (int ip=0; ip<n_ip; ++ip) {
	if (kicker_positions(ip, coor).size() != kicker_positions(0, 0).size()) {
	  cerr << "Number of given " << "xy"[coor] << "-coordinates at IP " << ip
	       << "(" << kicker_positions(ip, coor).size()
	       << ") and x-coordinates at IP 0 ("
	       << kicker_positions(0, 0).size() << ") are different\n";
	  exit(1);
	}
      }
    }
    // Bilinear interpolation
    const vector<long int>& n_cells =
      c.vl("Density.and.Field.bilinear.interpolators.N.cells.along.grid.side");
    if (n_cells.size() != 2) {
      cerr << "Density.and.Field.bilinear.interpolators.N.cells.along.grid.side "
	   << "should contain two numbers (for the density and for the field)\n";
      exit(1);
    }
    reset_interpolators(n_cells[0], n_cells[1]);
  }
protected:
  // mask initialization functions from the base class:
  void reset_kicked_bunch(int coor, const vector<double>& sigmas,
			  const vector<double>& weights, double n_sig_cut);
  void reset_ip_number(size_t n_ip);
  void reset_kicker_bunch(int ip, int coor, const vector<double>& sigmas,
			  const vector<double>& weights, const vector<double>& separations,
			  double phase_advance);
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
    if (n > 4 && inp.substr(n-4) == ".txt") {
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

  N n(c);
  Config_Mutli_XY_Gaussian_bunches bb(c, n.ip());
  if (c.vd("beta").size() != n.ip()) {
    cerr << "Given number of \"beta\" values (" << c.vd("beta").size()
	 << ") is not equal to the number of IPs (" << n.ip() << ")\n";
    return 1;
  }
  if (c.vd("N2").size() != n.ip()) {
    cerr << "Given number of \"N2\" values (" << c.vd("N2").size()
	 << ") is not equal to the number of IPs (" << n.ip() << ")\n";
    return 1;
  }
  vector<double> kick_const(n.ip());
  {
    double alpha = 1. / 137.035;
    double hbar = 0.197327e-15; // in Gev * m
    double beta0 = 1; // velocity/c of the second bunch in the reference frame of the first
    transform(c.vd("beta").begin(), c.vd("beta").end(),
	      c.vd("N2"  ).begin(), kick_const.begin(),
	      [factor = 2 * c("Z1") * c("Z2") * alpha * hbar /
	       beta0 / c["p"] * 1e12] // in um^2
	      (double beta, double N2) { return factor * beta * N2; });
  }
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
  // particle's weights, radii in X-X' and in Y-Y' (some of the particles
  // will not be simulated, so size()'s == N_points_in_step <= N_ponts):
  vector<double> rx, ry, w;
  {
    // sqrt(n.points_max()) will be used to sample intervals
    // rX = [0...cut_x], rY = [0...cut_y], where cut_x,y were chosen to
    // reproduce the efficiency of "N.sigma" configuration cut for a single
    // Gaussian bunch. Ie. the areas of the rejected multi-Gaussian tails ==
    // the area of a single Gaussian outside +/-N.sigma.
    //
    int N_sqrt = int(sqrt(n.points_max()));
    vector<double> rX_bins(N_sqrt), rY_bins(N_sqrt);
    double
      binx = bb.r_cut().first  / N_sqrt,
      biny = bb.r_cut().second / N_sqrt;
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
	complex<double> z1(rX_bins[ix]/bb.r_cut().first,
			   rY_bins[iy]/bb.r_cut().second);
	if (norm(z1) <=1 ) {
	  double rho1 =
	    bb.not_normalized_r_kicked_density(0, rX_bins[ix]) *
	    bb.not_normalized_r_kicked_density(1, rY_bins[iy]);
	  rx.push_back(rX_bins[ix]);
	  ry.push_back(rY_bins[iy]);
	  w.push_back(rho1);
	  w_sum += rho1;
	}
      }
    }
    // normalize w
    for_each(w.begin(), w.end(), [w_sum] (double& x) { x /= w_sum; });
  }
  n.points() = int(w.size());
  if (c.defined("output")) {
    if (c.s("output") == "rx.ry.weights") {
      ogzstream os((output_dir + "/rx_ry_weights.txt.gz").c_str());
      for (size_t i = 0; i < n.points(); ++i) {
	os << i << " "
	   << rx[i] << " "
	   << ry[i] << " "
	   << w[i]<< '\n';
      }
    } else {
      cerr << "Unrecognized output:" << c.s("output")
	   << ". Terminating" << endl;
      return 1;
    }
  }
  //
  // print Config when all type conversions are resolved
  cout << c;
  ofstream output_summary(output_dir + "/summary.txt");
  //
  // main loop
  for (int step = 0; step < bb.kicker_positions(0,0).size(); ++step) {
    vector<complex<double> > z2(bb.n_ip());
    for (size_t ip=0; ip<z2.size(); ++ip) {
      z2[ip] = complex<double>(bb.kicker_positions(ip, 0)[step],
			       bb.kicker_positions(ip, 1)[step]);
    }
    cout << "Processing beam separation";
    if (bb.kicker_positions(0, 0).size() > 1) cout << "s";
    for (size_t ip=0; ip<n.ip(); ++ip) {
      cout << "(" << z2[ip] << ") ";
    }
    cout << "..." << flush;
    // average kick when the first bunch is at -z2[ip] and the second is at
    // (0,0):
    vector<complex<double> > kick_average(n.ip());
    for (int ip=0; ip<n.ip(); ++ip) {
      kick_average[ip] = kick_const[ip] *
	bb.field_averaged_over_kicked_bunch(-real(z2[ip]), -imag(z2[ip]), ip);
    }
    //
    vector<vector<double> > kick_adiabatic(n.ip(), vector<double>(n.turns(N::ADIABATIC)));
    vector<vector<complex<double> > >
      kick_average_adiabatic(n.ip(), vector<complex<double> >(n.turns(N::ADIABATIC)));
    //
    Outputs output(n, c, output_dir,
		   // find max in pair<double,double> returned by bb.r_cut()
		   [](auto p) {return max(p.first, p.second); }(bb.r_cut()));
    //
    for (int i_turn = 0; i_turn < n.turns(N::ADIABATIC); ++i_turn) {
      double trans_w = double(i_turn + 1) / n.turns(N::ADIABATIC);
      // linearly increasing from 1/n.turns(N::ADIABATIC) for i_turn=0,
      //                       to 1 for i_turn = n.turns(N::ADIABATIC)-1
      for (int ip=0; ip<n.ip(); ++ip) {
	kick_adiabatic[ip][i_turn] = kick_const[ip] * trans_w;
	kick_average_adiabatic[ip][i_turn] = kick_average[ip] * trans_w;
      }
    }
    // prepare n.points() threads to run in parallel
    vector<thread> threads(n.points());
    auto loop = [&rx, &ry, &w, &bb, z2, step,
		 &kick_const, &kick_average, &kick_adiabatic, &kick_average_adiabatic,
		 &n, &output,
		 kick_model]
      (int i, double rnd1, double rnd2) {
		  // For every (rx,ry) pair obtain X-X', Y-Y' 4-dimensional
		  // coordinates by random rotation around X-X' and Y-Y' circles
		  complex<double>
		   xZ = rx[i] * exp(2i * M_PI * rnd1),
		   yZ = ry[i] * exp(2i * M_PI * rnd2);
		 //
		 auto turn_and_kick =
		   [&bb, kick_model, &z2](int ip,
					       // either trainsitional or nominal:
					       double k_const,
					       complex<double> k_average,
					       complex<double>& xZ, complex<double>& yZ) {
		     complex<double> kick, z1_minus_z2 = complex<double>(real(xZ), real(yZ)) - z2[ip];
		     if (kick_model == PRECISE) {
		       kick = k_const * bb.field(real(z1_minus_z2), imag(z1_minus_z2), ip);
		     } else if (kick_model == AVERAGE) {
		       kick =  k_average;
		     } else { //  PRECISE_MINUS_AVERAGE: 
		       kick = k_const * bb.field(real(z1_minus_z2), imag(z1_minus_z2), ip) - k_average;
		     }
		     // The momentum kick is subtracted below because of the minus sign in
		     // the definition of eg. zX = X - iX'. The kick is added to X',
		     // so that i*kick is subtracted from zX.
		     xZ = (xZ - 1i * real( kick )) * bb.exp_2pi_i_deltaQ(ip, 0);
		     yZ = (yZ - 1i * imag( kick )) * bb.exp_2pi_i_deltaQ(ip, 1);
		   };
		 //
		 auto fill_output =
		   [&w, &bb, &n, i, &z2,
		    &output]
		   (int ip, int i_turn, complex<double> xZ, complex<double> yZ,
		    int phase)
		   {
		     // overlap integral is calculated as a sum of (2nd bunch profile
		     // density * weight) over the sample of 1st bunch "particles".
		     //
		     complex<double> z1 = complex<double>(real(xZ), real(yZ));
		     complex<double> z1_minus_z2 = z1 - z2[ip];
		     double rho2 = bb.kicker_density(real(z1_minus_z2), imag(z1_minus_z2), ip);
		     // Integrals_Per_Circle keeps integrals for bunches with
		     // densities concentrated at one particle's point
		     // (delta-function density), so weight = 1. They will be
		     // reweighted with w[i] in the calculation of the final
		     // overlap integral appearing in summary
		     output.get<Integrals_Per_Circle>(ip).add(phase, i, rho2); // always filled
		     if (output.is_active<Avr_Z_Per_Circle>(ip)) {
		       output.get<Avr_Z_Per_Circle>(ip).add(phase, i, z1);
		     }
		     if (output.is_active<Integrals>(ip)) {
		       // Integrals a reweighted immediately
		       output.get<Integrals>(ip).add(i_turn, rho2 * w[i]);
		     }
		     if (output.is_active<Avr_Z>(ip)) {
		       output.get<Avr_Z>(ip).add(i_turn, z1 * w[i]);
		     }
		     if (output.is_active<Points>(ip) && (i_turn + 1) % n.select_turns() == 0) {
		       // It is more convenient here to count from one,
		       // eg. consider i_turn+1.  i_turn+1 runs in the range
		       // 1...N_turns. Only multiples of n.select_turns() are
		       // written out. Eg. if n.select_turns() > 1 and N_turns
		       // is the multiple of n.select_turns(): the last i_turn
		       // + 1 = N_turn is written, while the first i_turn + 1
		       // = 1 - not.
		       int j_turn = (i_turn + 1) / n.select_turns() - 1;
		       // in other words: j_turn+1 = (i_turn+1) /
		       // n.select_turns(), where i,j_turn + 1 correspond to
		       // counting from one.  Last j_turn+1 =
		       // N_turn/select_runs.
		       output.get<Points>(ip).add(j_turn, i, xZ, yZ);
		     }
		   };
		 int i_turn = 0;
		 // Phase 0:
		 // first N_turns_in_phase[0] turns without kick to measure the initial, possibly
		 // biased, integral; then calculate how much the kick changes
		 // it. Taking the ratio to the initialE_field   integral (instead of the exact
		 // integral) should cancel the bias at least partially and improve
		 // precision
		 for (int phase_turn=0; phase_turn<n.turns(N::NO_BB); ++phase_turn, ++i_turn) {
		   complex<double> z1( real( xZ ), real( yZ ));
		   for (size_t ip=0; ip<n.ip(); ++ip) {
		     xZ *= bb.exp_2pi_i_deltaQ(ip, 0); // kick = 0
		     yZ *= bb.exp_2pi_i_deltaQ(ip, 1);
		     // calculate and store integrals and points here, in the
		     // end, after the propagation through the ring (and after
		     // the beam-beam). The integral and the points (both of
		     // which can be written out), therefore, correspond to each
		     // other. One could do it before, then the last propagation
		     // through the ring for i_turn = N_turn-1 could be
		     // dropped. However, in this case the first integral+points
		     // always corresponded to the case without beam-beam, even
		     // if N_turns_in_phase[0]=0. So, for this reason and for simplicity the
		     // above logic is chosen.
		     fill_output(ip, i_turn, xZ, yZ, 0);
		   }
		 }
		 // Phase 1: adiabatic switch on
		 for (int phase_turn=0; phase_turn < n.turns(N::ADIABATIC); ++phase_turn, ++i_turn) {
		   for (size_t ip=0; ip<n.ip(); ++ip) {
		     turn_and_kick(ip,
				   kick_adiabatic[ip][phase_turn],
				   kick_average_adiabatic[ip][phase_turn],
				   xZ, yZ);
		     fill_output(ip, i_turn, xZ, yZ, 1);
		   }
		 }
		 // Phase 2: stabilization and phase 3: nominal beam-beam
		 for (int phase = N::STABILIZATION; phase < N::PHASES; ++phase) {
		   for (int phase_turn=0; phase_turn < n.turns(phase); ++phase_turn, ++i_turn) {
		     for (size_t ip=0; ip<n.ip(); ++ip) {
		       turn_and_kick(ip,
				     kick_const[ip],
				     kick_average[ip],
				     xZ, yZ);
		       fill_output(ip, i_turn, xZ, yZ, phase);
		     }
		   }
		 }
		 // Do not multiply *_per_circle variables by w[i] even here,
		 // but assume they correspond to rho1 density 100%
		 // concentrated at the simulated particle (or
		 // circle). Normalize only by n.turns()
		 for (size_t ip=0; ip<n.ip(); ++ip) {
		   for (int phase = 0; phase < N::PHASES; ++phase) {
		     output.get<Integrals_Per_Circle>(ip).normalize_by_n_turns(n.turns(phase), phase, i);
		     if (output.is_active<Avr_Z_Per_Circle>(ip)) {
		       output.get<Avr_Z_Per_Circle>(ip).normalize_by_n_turns(n.turns(phase), phase, i);
		     }
		   }
		 }
		}; // end of loop over tracked particle-points
    // submit threads:
    for (int i = 0; i < n.points(); ++i) {
      // supply random numbers here so that they are reproducible: if
      // drand48() were called in threads, it would depend on undefined
      // thread's execution order
      threads[i] = thread(loop, i, drand48(), drand48());
    }
    // wait for completion of all threads:
    for (int i = 0; i < n.points(); ++i) threads[i].join();
    // 
    output.write(step, z2);
    // write summary
    for (int ip=0; ip<n.ip(); ++ip) {
      auto integ_avr = output.get<Integrals_Per_Circle>(ip).weighted_average(w);
      // int0 can also be calculated as a mean of out<Integrals> (per turn)
      double int0 = integ_avr[N::NO_BB];
      double int_to_int0_correction = integ_avr[N::BB] / int0;

      double int0_analytic = bb.overlap_integral(real(z2[ip]), imag(z2[ip]), ip);
      double int0_to_analytic = int0 / int0_analytic;

      complex<double> avr_z = output.get<Avr_Z_Per_Circle>(ip).weighted_average(w)[N::BB];
      // the sign of "analytic_kick" (printed below) is not inverted, ie. it is
      // for X', not for z = X - iX'
      //
      // note, it has units of X', ie. of length; to convert to angles, one needs
      // to divide by beta*
      //
      // kick_const is halved because of the accelerator dynamics: if field() =
      // const = E, after many turns the ratio angle/E is kicked from initial
      // -kick_const/2 to final +kick_const/2, and their difference is kick_const.
      //
      complex<double> analytic_kick = kick_average[ip] / 2.;
      // the X(Y) shift corresponding to X'(Y') angular kick
      complex<double> approx_analytic_avr_z =
	complex<double>(real(analytic_kick) / tan(M_PI * bb.deltaQ(ip, 0)),
			imag(analytic_kick) / tan(M_PI * bb.deltaQ(ip, 1)));
      double int0_rel_err = -9999;
      if (output.is_active<Integrals>(ip)) {
	double sd_int0 = 0;
	for (int turn=0; turn<n.turns(0); ++turn) {
	  double dx = output.get<Integrals>(ip).integ()[turn] - int0;
	  sd_int0 += dx * dx;
	}    
	int0_rel_err = sqrt(sd_int0 / (n.turns(0) - 1.) / double(n.turns(0))) / int0;
      }
      output_summary << step << " "
		     << ip << " "
		     << z2[ip] << " "
		     << int_to_int0_correction << " "
		     << int0_analytic << " "
		     << int0 << " "
		     << int0_to_analytic << " "
		     << int0_rel_err << " "
		     << avr_z << " "
		     << analytic_kick << " "
		     << approx_analytic_avr_z << "\n";
    }
    cout << " done" << endl;
  } // end of loop over steps
    //
  if ( (bb.is_density_interpolated() || bb.is_field_interpolated()) &&
       c.defined("N.random.points.to.check.bilinear.interpolation")) {
    // check the precision of the interpolation by measuring the absolute
    // mismatches with the exact density/field values at randomly distributed
    // points inside the interpolation grid
    int N_random_points = c("N.random.points.to.check.bilinear.interpolation");
    if (bb.is_density_interpolated()) {
      vector<array<pair<double, double>, 2> > v =
	bb.max_and_average_interpolation_mismatches_relative_to_max_density(N_random_points);
      cout << "Max and average |precise - interpolated| density mismatches normalized to max density at\n";
      for (int ip=0; ip<v.size(); ++ip) {
	cout << "    IP " << ip << " along X: "
	     << v[ip][0].first << ", " << v[ip][0].second << "; along Y: "
	     << v[ip][1].first << ", " << v[ip][1].second << "\n";
      }
    }
    if (bb.is_field_interpolated()) {
      vector<pair<double, double> > v =
	bb.max_and_average_interpolation_mismatches_relative_to_max_field(N_random_points);
      cout << "Max and average |precise - interpolated| E-field mismatches normalized to max |E| at\n";
      for (int ip=0; ip<v.size(); ++ip) {
	cout << "    IP " << ip << ": "
	     << v[ip].first << ", " << v[ip].second << "\n";
      }
    }
  }
  cout << "Final results are written to " << output_dir << " directory\n"
       << "The summary appears in \"summary.txt\".\n"
       << "The 5th column \"int_to_int0_correction\" in this file contains the\n"
       << "final correction, i.e. the ratio of the luminosities with and\n"
       << "without beam-beam electromagnetic interaction, for the scan step and\n"
       << "the corresponding beam offset specified in the first three columns.\n\n"
    //
       << "The columns 6-9 refer to no-beam-beam case: : the analytic and numeric\n"
       << "integrals, their ratio and the estimated relative statistical\n"
       << "error of the numeric integral.\n\n"
    //
       << "The (10,11) column pair gives the numerically calculated x-y shift of\n"
       << "of the first bunch, while (12,13) and (14,15) pairs give the\n"
       << "analitically calculated values, the first in the coordinates\n"
       << "beta*(dx/dz, dy/dz) and the second - in (x,y) under the assumption\n"
       << "that its ratio to the first pair is -1/tan(pi*Q).\n";
  {
    ofstream config_out((output_dir + "/config.txt").c_str());
    config_out << c;
  }
  return(0);
}
