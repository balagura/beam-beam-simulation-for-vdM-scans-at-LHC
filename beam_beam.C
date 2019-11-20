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

//
// Same as Mutli_XY_Gaussian_bunches but takes all parameters from Config& c
// and also "n_ip"
//
struct Config_Mutli_XY_Gaussian_bunches : public Mutli_XY_Gaussian_bunches {
  Config_Mutli_XY_Gaussian_bunches(Config& c, int n_ip) {
    Mutli_XY_Gaussian_bunches::reset_ip_number(n_ip);
    //
    double n_sig_cut = c.defined("N.sigma.cut") ? c("N.sigma.cut") : 5;
    //
    int kicked_ip = c("kicked.ip");
    if (kicked_ip <= 0 || kicked_ip > n_ip) {
      cerr << "\"kicked.ip\" is out of range (its counting starts from one)\n";
      exit(1);
     }
    kicked_ip -= 1; // count from zero
    // accumulate beam separations:
    vector<array<vector<double>, 2> > positions(n_ip); // [ip][coor][step]
    vector<array<double, 2> > betatron_phase_over_2pi_at_next_ip(n_ip);
    for (size_t coor=0; coor<2; ++coor) {
      const vector<double>& deltaQs = c.vd(string("kicker.") + "xy"[coor] + ".next.phase.over.2pi");
      if (deltaQs.size() != n_ip) {
	cerr << "The number of given phases (" << deltaQs.size()
	     << ") is not equal to the number of IPs (" << n_ip << ")\n";
	exit(1);
      }
      for (int ip=0; ip<n_ip; ++ip) {
	// Phases or full tunes might be given with 2-3 digits after
	// comma. So, after 100 or 1000 turns the points return to their original
	// positions (as 100* or 1000*Qx,y becomes integer). To avoid this and add
	// extra randomness, a "negligible" irrational number randomly distributed
	// in the interval -exp(-8.5)...+exp(-8.5) = +/-0.0001017342 is added by
	// default to deltaQ. This makes it irrational and opens otherwise closed
	// Lissajous figures of betatron oscillations in X-Y plane. In reality the
	// phases should be irrational.
	//
	// If this is not desired, an option
	//   exact.phases = TRUE
	// should be set in the configuration. Then all phase
	// from the configuration file are used "as is".
	//
	if (!c.defined("exact.phases") ||
	    c.s("exact.phases") != "TRUE") {
	  betatron_phase_over_2pi_at_next_ip[ip][coor] = deltaQs[ip] + exp(-8.5) * (drand48() - 0.5);
	} else {
	  betatron_phase_over_2pi_at_next_ip[ip][coor] = deltaQs[ip];
	}
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
						      weights);
      }
      string s = string("kicked.") + "xy"[coor];
      if (c.vd(s + ".beta").size() != n_ip) {
	cerr << "Given number of " << "xy"[coor] << "-beta values (" << c.vd(s + ".beta").size()
	     << ") is not equal to the number of IPs (" << n_ip << ")\n";
	exit(1);
      }
      vector<double> weights = c.defined(s + ".weight") ? c.vd(s + ".weight") : vector<double>();
      Mutli_XY_Gaussian_bunches::reset_kicked_bunch(kicked_ip,
						    coor,
						    c.vd(s + ".sig"),
						    weights,
						    c.vd(s + ".beta"),
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
      c.vl("Density.and.Field.interpolators.N.cells.along.grid.side");
    if (n_cells.size() != 2) {
      cerr << "Density.and.Field.interpolators.N.cells.along.grid.side "
	   << "should contain two numbers (for the density and for the field)\n";
      exit(1);
    }
    Mutli_XY_Gaussian_bunches::reset_interpolators(n_cells[0], n_cells[1]);
    Mutli_XY_Gaussian_bunches::reset_phases(betatron_phase_over_2pi_at_next_ip);
  }
protected:
  // mask initialization functions from the base class:
  void reset_kicked_bunch(int coor, const vector<double>& sigmas,
			  const vector<double>& weights, double n_sig_cut);
  void reset_ip_number(size_t n_ip);
  void reset_kicker_bunch(int ip, int coor, const vector<double>& sigmas,
			  const vector<double>& weights, const vector<double>& separations);
  void reset_phases(vector<array<double, 2> > betatron_phase_over_2pi_at_next_ip);
  void reset_interpolators(int n_density_cells,
			   int n_field_cells);
};

int main(int argc, char** argv) {
  string config_file;
  if (argc != 2) {
    cout << "Usage: <program name> <configuration file>\n";
    return 0;
  }
  ifstream conf_file(argv[1]);
  Config c(conf_file);
  //
  // Set a random seed here, in the very beginning as it will be
  // required in Config_Mutli_XY_Gaussian_bunches() constructor for
  // an irrational addition to phases
  if (c.defined("seed")) {
    srand48(c("seed"));
  } else {
    srand48(time(NULL));
  }
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
  N n(c);
  Config_Mutli_XY_Gaussian_bunches bb(c, n.ip());
  if (c.vd("N.kicker.particles").size() != n.ip()) {
    cerr << "Given number of \"N.kicker.particles\" values (" << c.vd("N.kicker.particles").size()
	 << ") is not equal to the number of IPs (" << n.ip() << ")\n";
    return 1;
  }
  vector<array<double, 2> > kick_const(n.ip()); // kick_const[ip][coor]
  {
    double alpha = 1. / 137.035;
    double hbar = 0.197327e-15; // in Gev * m
    double beta0 = 1; // velocity/c of the second bunch in the reference frame of the first
    double factor =
      2 * c("Z.kicked") * c("Z.kicker") * alpha * hbar / beta0 / c["p"] * 1e12; // in um^2
    for (int ip=0; ip<n.ip(); ++ip) {
      double n = c.vd("N.kicker.particles")[ip];
      for (int coor=0; coor<2; ++coor) {
	kick_const[ip][coor] = factor * n * bb.accelerator_beta(ip, coor);
      }
    }
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
  // particle's weights, radii in X-X' and in Y-Y' at ip=0, from where the
  // simulation starts (some of the particles will not be simulated, so
  // size()'s == N_points_in_step <= N_ponts):
  vector<double> rx, ry, w;
  {
    // sqrt(n.points_max()) will be used to sample intervals rX = [0...cut_x],
    // rY = [0...cut_y], where cut_x,y will be chosen to reproduce the
    // efficiency of configuration "N.sigma.cut" for a single Gaussian
    // bunch. Ie. the fraction of eg. X-X' 2D multi-Gaussian inside a circle
    // with the radius cut_x will be equal to the fraction of 2D Gaussian
    // inside +/-N.sigma.cut.
    //
    int N_sqrt = int(sqrt(n.points_max()));
    vector<double> rX_bins(N_sqrt), rY_bins(N_sqrt);
    double
      binx = bb.r_cut(0).first  / N_sqrt, // for ip=0, from where the simulation starts
      biny = bb.r_cut(0).second / N_sqrt;
    for (int i=0; i<N_sqrt; ++i) {
      rX_bins[i] = (i + 0.5) * binx;
      rY_bins[i] = (i + 0.5) * biny;
    }
    // To reduce the number of simulated points and CPU, select only the
    // points from the ellipse inscribed inside the rectangle [0...cut_x] X
    // [0...cut_y] in rX and rY.
    double w_sum = 0;
    for (int ix=0; ix<N_sqrt; ++ix) {
      for (int iy=0; iy<N_sqrt; ++iy) {
	complex<double> z1(rX_bins[ix]/bb.r_cut(0).first,
			   rY_bins[iy]/bb.r_cut(0).second);
	if (norm(z1) <=1 ) {
	  double rho1 = // for ip=0, from where the simulation will start
	    bb.not_normalized_r_kicked_density(0, rX_bins[ix], 0) *
	    bb.not_normalized_r_kicked_density(1, rY_bins[iy], 0);
	  rx.push_back(rX_bins[ix]);
	  ry.push_back(rY_bins[iy]);
	  w.push_back(rho1);
	  w_sum += rho1;
	}
      }
    } // normalize w:
    for_each(w.begin(), w.end(), [w_sum] (double& x) { x /= w_sum; });
  }
  n.points() = int(w.size());
  if (c.defined("output")) {
    if (c.s("output") == "rx.ry.weights") {
      ogzstream os((output_dir + "/rx_ry_weights.txt.gz").c_str());
      os << scientific; // change to scientific format for doubles
      for (int ip = 0; ip < n.ip(); ++ip) {
	double beta_scale_x =sqrt(bb.accelerator_beta(ip, 0) / bb.accelerator_beta(0, 0));
	double beta_scale_y =sqrt(bb.accelerator_beta(ip, 1) / bb.accelerator_beta(0, 1));
	for (size_t i = 0; i < n.points(); ++i) {
	  os << i << " "
	     << ip << " "
	     << rx[i] * beta_scale_x << " "
	     << ry[i] * beta_scale_y << " "
	     << w[i]<< '\n';
	}
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
  ofstream output_kicker_positions(output_dir + "/kicker_positions.txt");
  ofstream output_summary         (output_dir + "/summary.txt");
  output_kicker_positions << scientific;  // change to scientific format for doubles
  output_summary << scientific;
  //
  double max_xy_r = 0; // find max possible no beam-beam radius for atomic converter in output
  for (int ip=0; ip<n.ip(); ++ip) {
    auto p = bb.r_cut(ip);
    max_xy_r = max(max_xy_r, max(p.first, p.second));
  }
  Outputs output(n, c, output_dir, max_xy_r);
  //
  // transportation to next IP:
  //
  // x = emittance * sqrt(beta) * cos(phase(z))
  // dx/dz = -emittance / sqrt(beta) * sin(phase(z)),
  // since d phase/dz = 1/beta and at IP alpha = -0.5 d beta/dz = 0
  // (and same for y).
  //
  // In addition to the simulated phase change (circle's rotation), one needs
  // to change the radius proportionally to sqrt(beta), as emittance is
  // conserved. So, at the next IP the x,y scale = this_x,y_scale *
  // sqrt(next_beta_x,y / this_beta_x,y). For the scaled x' (ie. already
  // scaled by beta, x' = dx/dz * beta which has units of length to form the
  // circle in x-x' plane), the scale is the same.
  //
  // Store sqrt(next_beta_x,y / this_beta_x,y) in a container:
  //
  vector<array<double, 2> > sqrt_beta_ratio_next_to_this(n.ip());
  for (int ip=0; ip<n.ip(); ++ip) {
    int ip_next = (ip==n.ip()-1) ? 0 : (ip+1);
    for (int coor=0; coor<2; ++coor) {
      sqrt_beta_ratio_next_to_this[ip][coor] = sqrt(bb.accelerator_beta(ip_next, coor) /
						    bb.accelerator_beta(ip,      coor));
    }
  }
  // Initial point amplitudes (rX, rY) were prepared for ip=0, so simulation
  // loop should also start from ip=0
  //
  // main loop
  for (int step = 0; step < bb.kicker_positions(0,0).size(); ++step) {
    vector<complex<double> > z2(n.ip());
    for (size_t ip=0; ip<z2.size(); ++ip) {
      z2[ip] = complex<double>(bb.kicker_positions(ip, 0)[step],
			       bb.kicker_positions(ip, 1)[step]);
    }
    cout << "Processing beam separation";
    if (bb.kicker_positions(0, 0).size() > 1) cout << "s";
    for (size_t ip=0; ip<n.ip(); ++ip) {
      cout << z2[ip] << " ";
    }
    cout << "..." << flush;
    // average kick when the kicked bunch is at -z2[ip] and the kicker is at
    // (0,0):
    vector<complex<double> > kick_average(n.ip());
    for (int ip=0; ip<n.ip(); ++ip) {
      complex<double> field = bb.field_averaged_over_kicked_bunch(-real(z2[ip]), -imag(z2[ip]), ip);
      kick_average[ip] = complex<double>(kick_const[ip][0] * real(field),
					 kick_const[ip][1] * imag(field));
    }
    vector<vector<array<double, 2> > >
      kick_adiabatic(n.ip(), vector<array<double, 2> >(n.turns(N::ADIABATIC)));
    vector<vector<complex<double> > >
      kick_average_adiabatic(n.ip(), vector<complex<double> >(n.turns(N::ADIABATIC)));
    output.reset(n);
    //
    for (int i_turn = 0; i_turn < n.turns(N::ADIABATIC); ++i_turn) {
      double trans_w = double(i_turn + 1) / n.turns(N::ADIABATIC);
      // linearly increasing from 1/n.turns(N::ADIABATIC) for i_turn=0,
      //                       to 1 for i_turn = n.turns(N::ADIABATIC)-1
      for (int ip=0; ip<n.ip(); ++ip) {
	for (int coor=0; coor<2; ++coor) {
	  kick_adiabatic[ip][i_turn][coor] = kick_const[ip][coor] * trans_w;
	}
	kick_average_adiabatic[ip][i_turn] = kick_average[ip] * trans_w;
      }
    }
    // Phase 0: no beam-beam 
    // Phase 1: adiabatic switch on
    // Phase 2: stabilization
    // Phase 3: nominal beam-beam
    auto kick = [&bb,
		 &kick_const, &kick_average, &kick_adiabatic, &kick_average_adiabatic,
		 &kick_model](int ip,
			      int phase,
			      int phase_turn,
			      complex<double> z1_minus_z2) -> complex<double>
      {
       if (phase == N::NO_BB) return 0;
       array<double, 2> k_const;
       complex<double> k_average;
       if (phase == N::ADIABATIC) {
	 k_const   = kick_adiabatic        [ip][phase_turn];
	 k_average = kick_average_adiabatic[ip][phase_turn];
       } else {
	 k_const   = kick_const  [ip];
	 k_average = kick_average[ip];
       }
       if (kick_model == PRECISE) {
	 complex<double> field = bb.field(real(z1_minus_z2), imag(z1_minus_z2), ip);
	 return complex<double>(k_const[0] * real(field),
				k_const[1] * imag(field));
       } else if (kick_model == AVERAGE) {
	 return k_average;
       } else { //  PRECISE_MINUS_AVERAGE: 
	 complex<double> field = bb.field(real(z1_minus_z2), imag(z1_minus_z2), ip);
	 return complex<double>(k_const[0] * real(field),
				k_const[1] * imag(field)) - k_average;
       }
      };
    // prepare n.points() threads to run in parallel
    vector<thread> threads(n.points());
    auto loop = [&rx, &ry, &w, &bb, z2, step,
		 &sqrt_beta_ratio_next_to_this, &kick,
		 &n, &output,
		 kick_model]
      (int i, double rnd1, double rnd2) {
		  // For every (rx,ry) pair obtain X-X', Y-Y' 4-dimensional
		  // coordinates by random rotation around X-X' and Y-Y' circles
		  complex<double>
		   xZ = rx[i] * exp(2i * M_PI * rnd1),
		   yZ = ry[i] * exp(2i * M_PI * rnd2);
		  int i_turn = 0;
		 //
		 for (int phase = 0; phase < N::PHASES; ++phase) {
		   for (int phase_turn=0; phase_turn<n.turns(phase); ++phase_turn, ++i_turn) {
		     for (size_t ip=0; ip<n.ip(); ++ip) {
		       { // kick + turn
			 complex<double> z1 = complex<double>(real(xZ), real(yZ));
			 complex<double> z1_minus_z2 = z1 - z2[ip];
			 complex<double> k = kick(ip, phase, phase_turn, z1_minus_z2);
			 // The momentum kick is subtracted below because of
			 // the minus sign in the definition of eg. zX = X -
			 // iX'. The kick is added to X', so that i*kick is
			 // subtracted from zX./
			 xZ = (xZ - 1i * real( k )) * bb.exp_i_next_ip_phase_minus_this(ip, 0);
			 yZ = (yZ - 1i * imag( k )) * bb.exp_i_next_ip_phase_minus_this(ip, 1);
		       }
		       // calculate and store output here, in the end, after
		       // the propagation through the ring (and after the
		       // beam-beam). One could do it before, then the last
		       // propagation through the ring could be
		       // dropped. However, in this case the first output
		       // always corresponded to the case without beam-beam,
		       // even if N_turns for no beam-beam phase would be set
		       // to 0. For this reason and for simplicity the above
		       // logic is chosen.
		       //
		       // Overlap integral is calculated as a sum of (2nd bunch profile
		       // density * weight) over the sample of 1st bunch "particles".
		       //
		       complex<double> z1 = complex<double>(real(xZ), real(yZ));
		       complex<double> z1_minus_z2 = z1 - z2[ip];
		       double rho2 = bb.kicker_density(real(z1_minus_z2), imag(z1_minus_z2), ip);
		       // Integrals_Per_Particle keeps integrals for bunches with
		       // densities concentrated at one particle's point
		       // (delta-function density), so weight = 1. They will be
		       // reweighted with w[i] in the calculation of the final
		       // overlap integral appearing in summary
		       output.get<Integrals_Per_Particle>(ip).add(phase, i, rho2); // always filled
		       if (output.is_active<Avr_XY_Per_Particle>(ip)) {
			 output.get<Avr_XY_Per_Particle>(ip).add(phase, i, z1);
		       }
		       if (output.is_active<Integrals_Per_Turn>(ip)) {
			 // Integrals_Per_Turn are reweighted immediately
			 output.get<Integrals_Per_Turn>(ip).add(i_turn, rho2 * w[i]);
		       }
		       if (output.is_active<Avr_XY_Per_Turn>(ip)) {
			 output.get<Avr_XY_Per_Turn>(ip).add(i_turn, z1 * w[i]);
		       }
		       if (output.is_active<Points>(ip) && (i_turn + 1) % n.select_turns() == 0) {
			 // It is more convenient here to count from one,
			 // eg. consider i_turn+1. Then, it runs in the range
			 // 1...N_turns. Only multiples of n.select_turns()
			 // are written out. Eg. if n.select_turns() > 1 and
			 // N_turns is the multiple of n.select_turns(): the
			 // last i_turn + 1 = N_turn is written, while the
			 // first i_turn + 1 = 1 - not.
			 int j_turn = (i_turn + 1) / n.select_turns() - 1;
			 // in other words: j_turn+1 = (i_turn+1) /
			 // n.select_turns(), where i,j_turn + 1 correspond to
			 // counting from one.  Last j_turn+1 =
			 // N_turn/select_runs.
			 output.get<Points>(ip).add(j_turn, i, xZ, yZ);
		       }
		       // scale to next IP according to 
		       // sqrt(beta(next ip)/beta(this))
		       xZ *= sqrt_beta_ratio_next_to_this[ip][0];
		       yZ *= sqrt_beta_ratio_next_to_this[ip][1];
		     }
		   }
		 }
		 // Do not multiply *_per_particle variables by w[i] even
		 // here, but assume they correspond to rho1 density 100%
		 // concentrated at the simulated particle (or
		 // circle). Normalize only by n.turns()
		 for (size_t ip=0; ip<n.ip(); ++ip) {
		   for (int phase = 0; phase < N::PHASES; ++phase) {
		     int n_turns = n.turns(phase);
		     auto& integ = output.get<Integrals_Per_Particle>(ip);
		     integ.normalize_by_n_turns(n_turns, phase, i);
		     if (output.is_active<Avr_XY_Per_Particle>(ip)) {
		       auto& avr_xy = output.get<Avr_XY_Per_Particle>(ip);
		       avr_xy.normalize_by_n_turns(n_turns, phase, i);
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
    output.write(step);
    // write kicker_positions and summary
    for (int ip=0; ip<n.ip(); ++ip) {
      output_kicker_positions << step << " " << ip << " "
			      << real(z2[ip]) << " " << imag(z2[ip]) << "\n";
      // integrals for summary
      auto integ_avr = output.get<Integrals_Per_Particle>(ip).weighted_average(w);
      double int0_analytic = bb.overlap_integral(real(z2[ip]), imag(z2[ip]), ip);
      double int0_rel_err = nan("");
      if (output.is_active<Integrals_Per_Turn>(ip)) {
	double sd_int0 = 0;
	for (int turn=0; turn<n.turns(N::NO_BB); ++turn) {
	  double dx = output.get<Integrals_Per_Turn>(ip).integ()[turn] - integ_avr[N::NO_BB];
	  sd_int0 += dx * dx;
	}
	int0_rel_err = sqrt(sd_int0 / (n.turns(0) - 1.) / double(n.turns(0))) / int0_analytic;
      }
      // X(Y) orbit shift at ip, analytic formula:
      // sqrt(beta_ip) / 2 / sin(pi*tune) *
      // sum( angular_kick_at_ip2 * sqrt(beta_ip2) * cos(|phase_ip - phase_ip2| - pi*tune) )
      array<double, 2> analytic_avr_xy = {0., 0.};
      for (int coor=0; coor<2; ++coor) {
	double pi_tune = M_PI * bb.accelerator_tune(coor);
	double common_factor = sqrt(bb.accelerator_beta(ip, coor)) / 2 / sin(pi_tune);
	for (size_t ip2=0; ip2<kick_average.size(); ++ip2) {
	  double delta_phase = fabs(bb.betatron_ip_phase(ip, coor) - bb.betatron_ip_phase(ip2, coor));
	  double k = (coor == 0) ? real(kick_average[ip2]) : imag(kick_average[ip2]);
	     analytic_avr_xy[coor] +=
	       common_factor *
	       k *
	       cos(delta_phase - pi_tune) /
	       // divide, since kick_average[ip2] = angular_kick * beta_ip2 -
	       // already multiplied by beta_ip2
	       sqrt(bb.accelerator_beta(ip2, coor));
	}
      }
      output_summary << step << " "
		     << ip << " "
	// overlaps can also be calculated as means of <Integrals_Per_Turn> within a given phase
		     << integ_avr[N::BB] / integ_avr[N::NO_BB] << " "
		     << int0_analytic << " "
		     << integ_avr[N::NO_BB] << " "
		     << integ_avr[N::NO_BB] / int0_analytic << " "
		     << int0_rel_err << " ";
      if (output.is_active<Avr_XY_Per_Particle>(ip)) {
	auto avr_xy = output.get<Avr_XY_Per_Particle>(ip).weighted_average(w);
	output_summary << real(avr_xy[N::NO_BB]) << " " << imag(avr_xy[N::NO_BB]) << " "
		       << real(avr_xy[N::   BB]) << " " << imag(avr_xy[N::   BB]) << " ";
      } else {
	output_summary << nan("") << " " << nan("") << " " << nan("") << " " << nan("") << " ";
      }
      output_summary << analytic_avr_xy[0] << " "
		     << analytic_avr_xy[1] << "\n";
    }
    cout << " done" << endl;
  } // end of loop over steps
    //
  if ( (bb.is_density_interpolated() || bb.is_field_interpolated()) &&
       c.defined("N.random.points.to.check.interpolation")) {
    // check the precision of the interpolation by measuring the absolute
    // mismatches with the exact density/field values at randomly distributed
    // points inside the interpolation grid
    int N_random_points = c("N.random.points.to.check.interpolation");
    if (bb.is_density_interpolated()) {
      vector<array<pair<double, double>, 2> > v =
	bb.max_and_average_interpolation_mismatches_relative_to_max_density(N_random_points);
      cout << "Max and average |precise - interpolated| density mismatches normalized to max density at\n";
      for (int ip=0; ip<v.size(); ++ip) {
	cout << "    IP " << ip << " along X: "
	     << scientific  // change to scientific format for doubles
	     << v[ip][0].first << ", " << v[ip][0].second << "; along Y: "
	     << v[ip][1].first << ", " << v[ip][1].second << "\n"
	     << fixed;      // change back to defaul fixed format
      }
    }
    if (bb.is_field_interpolated()) {
      vector<pair<double, double> > v =
	bb.max_and_average_interpolation_mismatches_relative_to_max_field(N_random_points);
      cout << "Max and average |precise - interpolated| E-field mismatches normalized to max |E| at\n";
      for (int ip=0; ip<v.size(); ++ip) {
	cout << "    IP " << ip << ": "
	     << scientific
	     << v[ip].first << ", " << v[ip].second << "\n"
	     << fixed;
      }
    }
  }
  cout << "\nFinal results are written to " << output_dir << " directory.\n\n"
       << "The summary appears in \"summary.txt\" in the format:\n"
       << "step IP <beam-beam/no beam-beam luminosity correction>\n"
       << "<analytic>, <numeric> no beam-beam overlap integrals, <their ratio> and <its error> \n"
       << "(note, error = \"nan\" if \"integrals.per.turn\" option is not requested in \"output\")\n"
       << "numerically calculated <X>, <Y> center-of-mass shift of the kicked bunch without beam-beam\n"
       << "(should be zero) and <X>, <Y> shift with beam-beam (\"nan\" if \"avr.xy.per.particle\" was\n"
       << "not requested in \"output\"), the same shift calculated for <X> and <Y> analytically.\n\n";
  {
    ofstream config_out((output_dir + "/config.txt").c_str());
    config_out << c;
  }
  return(0);
}
