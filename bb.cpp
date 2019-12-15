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
//  Main simulation C++ function beam_beam()
//
#include <complex>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h> // for mkdir
#include <sys/types.h> // for mkdir
#include <fstream>
#include <algorithm>
#include <numeric>
#include <thread>
#include "gzstream.hh"
#include "multi_gaussian.hh"
#include "output.hh"
#include "bb.hh"

using namespace std;
using namespace std::complex_literals;

ostream& operator<<(ostream& os, const vector<Gaussian>& g) {
  os << "sigma/weight = ";
  for (size_t i=0; i<g.size(); ++i)
    os << "(" << g[i].sig << ", " << g[i].w << ") ";
  return os << "\n";
}
void write_vec_in_brackets(ostream& os, const string& prefix, const vector<double>& v, const string& suffix = "") {
  os << prefix << "(";
  for (size_t i=0; i<v.size(); ++i) {
    os << v[i];
    if (i != v.size()-1) os << ", ";
  }
  os << ")" << suffix;
}
bool consistent_n_ip(const Kicked& kicked, const Kickers& kickers, const Sim& sim) {
  size_t n_ip = kicked.beta[0].size();
  if (kickers.n_particles.size() != n_ip) return false;
  for (int coor=0; coor<2; ++coor) {
    if (kicked.beta[coor].size() != n_ip ||
	kicked.next_phase_over_2pi[coor].size() != n_ip ||
	kickers.gaussian[coor].size() != n_ip ||
	kickers.position[coor].size() != n_ip)
      return false;
  }
  return true;
}
void print(ostream& os, const Kicked& kicked, const Kickers& kickers, const Sim& sim) {
  if (!consistent_n_ip(kicked, kickers, sim)) {
    cerr << "Lengths of kicked X,Y-beta, X,Y-next_phase_over_2pi and "
	 << "kickers n_particles, X,Y-gaussian, X,Y-position\n"
	 << "should coincide and should be equal to the number of IPs\n";
  }
  int n_ip = kicked.beta[0].size();
  int n_step = 0; // set to maximal number of given kicker positions
  for (int coor=0; coor<2; ++coor) {
    for (int ip=0; ip<n_ip; ++ip) {
      n_step = max(n_step, int(kickers.position[coor][ip].size()));
    }
  }
  os << "N IP = " << n_ip << ", N steps = " << n_step << "\n"
     << "Kicked: p = " << kicked.momentum << ", Z = " << kicked.Z << ", kicked IP = " << kicked.ip << "\n";
  write_vec_in_brackets(os, "  beta X = ", kicked.beta[0], "\n");
  write_vec_in_brackets(os, "       Y = ", kicked.beta[1], "\n");
  write_vec_in_brackets(os, "  phases/2pi X = ", kicked.next_phase_over_2pi[0], "\n");
  write_vec_in_brackets(os, "             Y = ", kicked.next_phase_over_2pi[1], "\n");
  if (kicked.exact_phases) cout << "  phases will be used exactly as given\n";
  os << "  at IP " << kicked.ip
     << " Gaussian (sigma, weight) in X "
     << kicked.gaussian[0]
     << "                                   in Y "
     << kicked.gaussian[1];
  os << "Kickers: Z = " << kickers.Z;
  write_vec_in_brackets(os, ", N particles = ", kickers.n_particles, "\n");
  for (int ip=0; ip<n_ip; ++ip) {
    os << "  at IP " << ip
       << " Gaussian (sigma, weight) in X "
       << kickers.gaussian[0][ip]
       << "                                   in Y "
       << kickers.gaussian[1][ip];
  }
  os << "  Positions\n";
  for (int step=0; step<n_step; ++step) {
    os << "    step " << step << ": ";
    for (int ip=0; ip<n_ip; ++ip) {
      os << "(";
      if (step < kickers.position[0][ip].size()) os << kickers.position[0][ip][step];
      os << ", ";
      if (step < kickers.position[1][ip].size()) os << kickers.position[1][ip][step];
      os << ") ";
    }
    os << "\n";
  }
  os << "N points = " << sim.n_points << ", ";
  os << "N turns = (";
  for (int phase=0; phase<PHASES; ++phase) {
    os << sim.n_turns[phase];
    if (phase != PHASES-1) os << ", ";
  }
  os << "), select one turn out of " << sim.select_one_turn_out_of << "\n";
  os << "random seed = " << sim.seed << ", cut " << sim.n_sigma_cut << " sigma\n";
  os << "Kick model \"" << sim.kick_model << "\"\n";
  os << "N grid cells to interpolate density / field = ("
     << sim.density_and_field_interpolators_n_cells_along_grid_side[0]
     << ", " << sim.density_and_field_interpolators_n_cells_along_grid_side[1] << ")\n";
  os << "N points to check interpolation = " << sim.n_random_points_to_check_interpolation << "\n";
  os << "Output directory \"" << sim.output_dir << "\"\n";
  os << "output = " << sim.output << "\n";
}

// -------------------- Simulation --------------------
void beam_beam(const Kicked& kicked, const Kickers& kickers, const Sim& sim,
	       vector<vector<BB_Summary_Per_Step_IP> >* summary, bool quiet) {
  if (quiet) {
    if (!consistent_n_ip(kicked, kickers, sim)) {
      cerr << "Lengths of kicked X,Y-beta, X,Y-next_phase_over_2pi and "
	   << "kickers n_particles, X,Y-gaussian, X,Y-position\n"
	   << "should coincide and should be equal to the number of IPs\n";
      exit(1);
    }
  } else {
    print(cout, kicked, kickers, sim);
  }
  int n_ip = kicked.beta[0].size();
  if (kicked.ip < 0 || kicked.ip >= n_ip) {
    cerr << "kicked IP is out of range\n";
    exit(1);
  }
  int n_step = 0;
  for (int coor=0; coor<2; ++coor) {
    for (int ip=0; ip<n_ip; ++ip) {
      n_step = max(n_step, int(kickers.position[coor][ip].size()));
    }
  }
  vector<vector<double> > position[2];
  for (int coor=0; coor<2; ++coor) {
    position[coor].resize(n_ip);
    for (int ip=0; ip<n_ip; ++ip) {
      const auto& from = kickers.position[coor][ip];
      auto&       to   =         position[coor][ip];
      to.resize(n_step);
      if (from.size() == 1) { // recycle one number to full length
	for (int step=0; step<n_step; ++step) to[step] = from[0];
      } else if (from.size() == n_step) {
	for (int step=0; step<n_step; ++step) to[step] = from[step];
      } else {
	cerr << "Given number of kicker positions (scan steps) should be either maximal or one.\n"
	     << "(The latter case implies a stable position across all scan steps)\n";
	exit(1);
      }
    }
  }  
  // Set the random seed here as it might be needed to add random
  // irrationality to phases
  srand48(sim.seed);
  Multi_XY_Gaussian_bunches bb;
  {
    bb.reset_ip_number(n_ip);
    for (int coor=0; coor<2; ++coor) {
      auto& kg = kicked.gaussian[coor];
      vector<double> sig(kg.size()), w(kg.size());
      for (size_t gauss=0; gauss<kg.size(); ++gauss) {
	sig[gauss] = kg[gauss].sig;
	w  [gauss] = kg[gauss].w;
      }
      bb.reset_kicked_bunch(kicked.ip, coor, sig, w, kicked.beta[coor], sim.n_sigma_cut);
      for (int ip=0; ip<n_ip; ++ip) {
	auto& ksg = kickers.gaussian[coor][ip];
	sig.resize(ksg.size());
	w  .resize(ksg.size());
	for (size_t gauss=0; gauss<ksg.size(); ++gauss) {
	  sig[gauss] = ksg[gauss].sig;
	  w  [gauss] = ksg[gauss].w;
	}
	bb.reset_kicker_bunch(ip, coor, sig, w);
      }
    }
    bb.reset_interpolators(sim.density_and_field_interpolators_n_cells_along_grid_side[0],
			   sim.density_and_field_interpolators_n_cells_along_grid_side[1],
			   position);
  }
  int n_total_turns = accumulate(sim.n_turns, sim.n_turns + PHASES, 0.);
  double tune[2];
  vector<double> ip_phase[2]; // starts from 0 at ip=0
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
  // Store sqrt(next_beta_x,y / this_beta_x,y) * exp(2pi*i*(next phase - this)) in a container:
  //
  vector<complex<double> > transportation_factor[2]; // [coor][ip]
  for (size_t coor=0; coor<2; ++coor) {
    auto phase = kicked.next_phase_over_2pi[coor];
    if (kicked.exact_phases == false) {
      // Phases/2pi or full tunes might be given with 2-3 digits after
      // comma. So, after 100 or 1000 turns the points return to their
      // original positions (as 100* or 1000*Qx,y becomes integer). To avoid
      // this and add extra randomness, a "negligible" irrational number
      // randomly distributed in the interval -exp(-8.5)...+exp(-8.5) =
      // +/-0.0001017342 is added by default to all phases/2pi. This makes
      // them irrational and opens otherwise closed Lissajous figures of
      // betatron oscillations in X-Y plane. In reality the phases should be
      // irrational.
      //
      // If this is not desired, "exact_phases" below should be set to TRUE,
      // then they are used "as is".
      //
      for (auto& ph: phase) ph += exp(-8.5) * (drand48() - 0.5);
    }
    tune[coor] = phase[n_ip-1];
    for (auto& ph: phase) ph *= 2 * M_PI;
    ip_phase[coor].resize(n_ip);
    transportation_factor[coor].resize(n_ip);
    ip_phase[coor][0] = 0;
    transportation_factor[coor][0] = exp(1i * phase[0]);
    for (int ip=1; ip<n_ip; ++ip) {
      double curr = phase[ip-1];
      double next = phase[ip];
      ip_phase[coor][ip] = curr;
      transportation_factor[coor][ip] = exp(1i * (next - curr));
    }
    for (int ip=0; ip<n_ip; ++ip) {
      const auto& beta = kicked.beta[coor];
      int next_ip = (ip==n_ip-1) ? 0 : (ip+1);
      transportation_factor[coor][ip] *= sqrt(beta[next_ip] / beta[ip]);
    }
  }
  vector<array<double, 2> > kick_const(n_ip); // kick_const[ip][coor]
  {
    double alpha = 1. / 137.035;
    double hbar = 0.197327e-15; // in Gev * m
    double beta0 = 1; // velocity/c of the second bunch in the reference frame of the first
    double factor =
      2 * kicked.Z * kickers.Z * alpha * hbar / beta0 / kicked.momentum * 1e12; // in um^2
    for (int ip=0; ip<n_ip; ++ip) {
      double n = kickers.n_particles[ip];
      for (int coor=0; coor<2; ++coor) {
	kick_const[ip][coor] = factor * n * kicked.beta[coor][ip];
      }
    }
  }
  enum {PRECISE, PRECISE_MINUS_AVERAGE, AVERAGE} kick_model = PRECISE;
  if (string(sim.kick_model) != "precise") {
    if (string(sim.kick_model) == "precise.minus.average") {
      kick_model = PRECISE_MINUS_AVERAGE;
    } else if (string(sim.kick_model) == "average") {
      kick_model = AVERAGE;
    } else {
      cerr << "Wrong kick.model\n";
      exit(1);
    }
  }
  // particle's weights, radii in X-X' and in Y-Y' at ip=0, from where the
  // simulation starts (some of the particles will not be simulated, so
  // size()'s == n_points <= sim.n_points):
  vector<double> rx, ry, w;
  {
    // sqrt(sim.n_points) will be used to sample intervals rX = [0...cut_x],
    // rY = [0...cut_y], where cut_x,y will be chosen to reproduce the
    // efficiency of configuration "n_sigma_cut" for a single Gaussian
    // bunch. Ie. the fraction of eg. X-X' 2D multi-Gaussian inside a circle
    // with the radius cut_x will be equal to the fraction of 2D Gaussian
    // inside +/-n_sigma_cut.
    //
    int N_sqrt = int(sqrt(sim.n_points));
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
  int n_points = int(w.size());
  string output_dir = sim.output_dir;
  // if output_dir == "": output_options have no effect (not dir to write)
  ofstream output_summary;
  if (output_dir != "") {
    if( mkdir(sim.output_dir.c_str(), 0755) ) {
      if ( errno == EEXIST ) {
	cerr << "Output directory " << sim.output_dir << " already exists. Terminating\n";
      } else {
	cerr << "Can not create output directory " << sim.output_dir << ", errno = " << errno << endl;
      }
      exit(1);
    }
    // save input
    ofstream input(output_dir + "/input.txt");
    print(input, kicked, kickers, sim);
    //
    ogzstream output_rx_ry_weights((output_dir + "/rx_ry_weights.txt.gz").c_str());
    output_rx_ry_weights << scientific; // change to scientific format for doubles
    for (int ip = 0; ip < n_ip; ++ip) {
      double beta_scale_x =sqrt(kicked.beta[0][ip] / kicked.beta[0][0]);
      double beta_scale_y =sqrt(kicked.beta[1][ip] / kicked.beta[1][0]);
      for (size_t i = 0; i < n_points; ++i) {
	output_rx_ry_weights << i << " "
			     << ip << " "
			     << rx[i] * beta_scale_x << " "
			     << ry[i] * beta_scale_y << " "
			     << w[i]<< '\n';
      }
    }
    ofstream output_kicker_positions(output_dir + "/kicker_positions.txt");
    output_kicker_positions << scientific;  // change to scientific format for doubles
    for (int step=0; step<n_step; ++step) {
      for (int ip=0; ip<n_ip; ++ip) {
	output_kicker_positions << step << " " << ip << " "
				<< position[0][ip][step] << " "
				<< position[1][ip][step] << "\n";
      }
    }
    int n_random = sim.n_random_points_to_check_interpolation;
    if ( (bb.is_density_interpolated() || bb.is_field_interpolated()) &&
	 n_random > 0) {
      // check the precision of the interpolation by measuring the absolute
      // mismatches with the exact density/field values at randomly distributed
      // points inside the interpolation grid
      ofstream output_interp(output_dir + "/interpolation_precision.txt");
      output_interp << scientific;  // change to scientific format for doubles
      if (bb.is_density_interpolated()) {
	vector<array<pair<double, double>, 2> > v =
	  bb.max_and_average_interpolation_mismatches_relative_to_max_density(n_random);
	output_interp << "Max and average |precise - interpolated| density mismatches normalized to max density at\n";
	for (int ip=0; ip<v.size(); ++ip) {
	  output_interp << "    IP " << ip << " along X: "
			<< v[ip][0].first << ", " << v[ip][0].second << "; along Y: "
			<< v[ip][1].first << ", " << v[ip][1].second << "\n";
	}
      }
      if (bb.is_field_interpolated()) {
	vector<pair<double, double> > v =
	  bb.max_and_average_interpolation_mismatches_relative_to_max_field(n_random);
	output_interp << "Max and average |precise - interpolated| E-field mismatches normalized to max |E| at\n";
	for (int ip=0; ip<v.size(); ++ip) {
	  output_interp << "    IP " << ip << ": "
			<< v[ip].first << ", " << v[ip].second << "\n";
	}
      }
    }
    output_summary.open(output_dir + "/summary.txt");
    output_summary << scientific;
  }
  double max_xy_r = 0; // find max possible no beam-beam radius for atomic converter in output
  for (int ip=0; ip<n_ip; ++ip) {
    auto p = bb.r_cut(ip);
    max_xy_r = max(max_xy_r, max(p.first, p.second));
  }
  Outputs output(n_ip, n_points, max_xy_r, sim.output, output_dir);
  summary->resize(n_ip, vector<BB_Summary_Per_Step_IP>(n_step));
  // Initial point amplitudes (rX, rY) were prepared for ip=0, so simulation
  // loop should also start from ip=0
  //
  // main loop
  for (int step = 0; step < n_step; ++step) {
    vector<complex<double> > z2(n_ip);
    for (size_t ip=0; ip<z2.size(); ++ip) {
      z2[ip] = complex<double>(position[0][ip][step],
			       position[1][ip][step]);
    }
    if (!quiet) {
      cout << "Simulating kicker position";
      if (n_ip > 1) cout << "s";
      cout << " ";
      for (size_t ip=0; ip<n_ip; ++ip) cout << z2[ip] << " ";
      cout << "..." << flush;
    }
    // average kick when the kicked bunch is at -z2[ip] and the kicker is at
    // (0,0):
    vector<complex<double> > kick_average(n_ip);
    for (int ip=0; ip<n_ip; ++ip) {
      complex<double> field = bb.field_averaged_over_kicked_bunch(-real(z2[ip]), -imag(z2[ip]), ip);
      kick_average[ip] = complex<double>(kick_const[ip][0] * real(field),
					 kick_const[ip][1] * imag(field));
    }
    vector<vector<array<double, 2> > >
      kick_adiabatic(n_ip, vector<array<double, 2> >(sim.n_turns[ADIABATIC]));
    vector<vector<complex<double> > >
      kick_average_adiabatic(n_ip, vector<complex<double> >(sim.n_turns[ADIABATIC]));
    int selected_turns = (sim.select_one_turn_out_of < 1) ? 1 : sim.select_one_turn_out_of;
    output.reset(n_ip,
    		 n_points,
    		 n_total_turns,
    		 selected_turns);
    for (int i_turn = 0; i_turn < sim.n_turns[ADIABATIC]; ++i_turn) {
      double trans_w = double(i_turn + 1) / sim.n_turns[ADIABATIC];
      // linearly increasing from 1/sim.n_turns[ADIABATIC] for i_turn=0,
      //                       to 1 for i_turn = sim.n_turns[ADIABATIC]-1
      for (int ip=0; ip<n_ip; ++ip) {
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
       if (phase == NO_BB) return 0;
       array<double, 2> k_const;
       complex<double> k_average;
       if (phase == ADIABATIC) {
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
    // prepare n_points threads to run in parallel
    vector<thread> threads(n_points);
    auto loop = [&rx, &ry, &w, &bb, z2, step,
		 &transportation_factor, &kick,
		 &sim, &n_ip, &output, &selected_turns,
		 kick_model]
      (int i, double rnd1, double rnd2) {
		  // For every (rx,ry) pair obtain X-X', Y-Y' 4-dimensional
		  // coordinates by random rotation around X-X' and Y-Y' circles
		  complex<double>
		   xZ = rx[i] * exp(2i * M_PI * rnd1),
		   yZ = ry[i] * exp(2i * M_PI * rnd2);
		  int i_turn = 0;
		 //
		 for (int phase = 0; phase < PHASES; ++phase) {
		   for (int phase_turn=0; phase_turn<sim.n_turns[phase]; ++phase_turn, ++i_turn) {
		     for (size_t ip=0; ip<n_ip; ++ip) {
		       // It is more convenient to fill statistics before the
		       // propagation through the ring and before the
		       // beam-beam: then the statistics and the
		       // transportation are performed from the same ip (in
		       // particular, with the same beta) and one can start
		       // both from ip=0.  Note, however, that with the chosen
		       // logic the first turn in ip=0 is always without
		       // beam-beam (one starts from the prepared randomized xZ,
		       // yZ coordinates)
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
		       if (output.is_active<Points>(ip) && (i_turn + 1) %
			   selected_turns == 0) {
			 // It is more convenient here to count from one,
			 // eg. consider i_turn+1. Then, it runs in the range
			 // 1...N_turns. Only multiples of
			 // selected_turns are written out. Eg. if
			 // selected_turns > 1 and N_turns is the
			 // multiple of selected_turns: the last
			 // i_turn + 1 = N_turn is written, while the first
			 // i_turn + 1 = 1 - not.
			 int j_turn = (i_turn + 1) / selected_turns - 1;
			 // in other words: j_turn+1 = (i_turn+1) /
			 // selected_turns, where i,j_turn + 1
			 // correspond to counting from one.  Last j_turn+1 =
			 // N_turn/select_runs.
			 output.get<Points>(ip).add(j_turn, i, xZ, yZ);
		       }
		       // transport to next ip: kick + turn
		       complex<double> k = kick(ip, phase, phase_turn, z1_minus_z2);
		       // The momentum kick is subtracted below because of
		       // the minus sign in the definition of eg. zX = X -
		       // iX'. The kick is added to X', so that i*kick is
		       // subtracted from zX.
		       //
		       // "transportation_factor" changes scale in next IP
		       // according to sqrt(beta(next ip)/beta(this)) and
		       // turns the phase
		       xZ = (xZ - 1i * real( k )) * transportation_factor[0][ip];
		       yZ = (yZ - 1i * imag( k )) * transportation_factor[1][ip];
		     }
		   }
		 }
		 // Do not multiply *_per_particle variables by w[i] even
		 // here, but assume they correspond to rho1 density 100%
		 // concentrated at the simulated particle (or
		 // circle). Normalize only by n_turns
		 for (size_t ip=0; ip<n_ip; ++ip) {
		   for (int phase = 0; phase < PHASES; ++phase) {
		     int n_turns = sim.n_turns[phase];
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
    for (int i = 0; i < n_points; ++i) {
      // supply random numbers here so that they are reproducible: if
      // drand48() were called in threads, it would depend on undefined
      // thread's execution order
      threads[i] = thread(loop, i, drand48(), drand48());
    }
    // wait for completion of all threads:
    for (int i = 0; i < n_points; ++i) threads[i].join();
    //
    
    output.write(step);
    // write kicker_positions and summary
    for (int ip=0; ip<n_ip; ++ip) {
      // integral overlaps are calculated below from
      // <Integrals_Per_Particles>, but can also be calculated by averaging
      // <Integrals_Per_Turn> within a given phase
      auto integ_avr = output.get<Integrals_Per_Particle>(ip).weighted_average(w);
      // estimate numeric/analytic no-beam-beam ratio from the spread of turn integrals
      double int0_analytic = bb.overlap_integral(real(z2[ip]), imag(z2[ip]), ip);
      double int0_rel_err = nan("");
      if (output.is_active<Integrals_Per_Turn>(ip)) {
	double sd_int0 = 0;
	for (int turn=0; turn<sim.n_turns[NO_BB]; ++turn) {
	  double dx = output.get<Integrals_Per_Turn>(ip).integ()[turn] - integ_avr[NO_BB];
	  sd_int0 += dx * dx;
	}
	int0_rel_err = sqrt(sd_int0 / (sim.n_turns[NO_BB] - 1.) / double(sim.n_turns[NO_BB])) / int0_analytic;
      }
      // X(Y) orbit shift at ip, analytic formula:
      // sqrt(beta_ip) / 2 / sin(pi*tune) *
      // sum( angular_kick_at_ip2 * sqrt(beta_ip2) * cos(|phase_ip - phase_ip2| - pi*tune) )
      array<double, 2> analytic_avr_xy = {0., 0.};
      if (kick_model != PRECISE_MINUS_AVERAGE) { // 0 for precise.minus.average
	for (int coor=0; coor<2; ++coor) {
	  double pi_tune = M_PI * tune[coor];
	  double common_factor = sqrt(kicked.beta[coor][ip]) / 2 / sin(pi_tune);
	  for (size_t ip2=0; ip2<n_ip; ++ip2) {
	    double delta_phase = fabs(ip_phase[coor][ip] - ip_phase[coor][ip2]);
	    double k = (coor == 0) ? real(kick_average[ip2]) : imag(kick_average[ip2]);
	    analytic_avr_xy[coor] +=
	      common_factor *
	      k *
	      cos(delta_phase - pi_tune) /
	      // divide, since kick_average[ip2] = angular_kick * beta_ip2 -
	      // already multiplied by beta_ip2
	      sqrt(kicked.beta[coor][ip2]);
	  }
	}
      } // fill summary
      BB_Summary_Per_Step_IP& s = (*summary)[ip][step];
      s.correction = integ_avr[BB] / integ_avr[NO_BB];
	s.no_bb_analytic_integ = int0_analytic;
	s.no_bb_numeric_over_analytic_integ = integ_avr[NO_BB] / int0_analytic;
	s.no_bb_numeric_over_analytic_integ_err = int0_rel_err;
	s.avr_analytic[0] = analytic_avr_xy[0];
	s.avr_analytic[1] = analytic_avr_xy[1];
      if (output.is_active<Avr_XY_Per_Particle>(ip)) {
	auto avr_xy = output.get<Avr_XY_Per_Particle>(ip).weighted_average(w);
	s.      avr_numeric[0] = real(avr_xy[   BB]);
	s.      avr_numeric[1] = imag(avr_xy[   BB]);
	s.no_bb_avr_numeric[0] = real(avr_xy[NO_BB]);
	s.no_bb_avr_numeric[1] = imag(avr_xy[NO_BB]);
      } else {
	s.      avr_numeric[0] = s.      avr_numeric[1] = 
	s.no_bb_avr_numeric[0] = s.no_bb_avr_numeric[1] = nan("");
      }
      if (output_dir != "") {
	output_summary << step << " "
		       << ip << " "
		       << s.correction << " "
		       << s.no_bb_analytic_integ << " "
		       << s.no_bb_numeric_over_analytic_integ << " "
		       << s.no_bb_numeric_over_analytic_integ_err << " "
		       << s.no_bb_avr_numeric[0] << " "
		       << s.no_bb_avr_numeric[1] << " "
		       << s.avr_analytic[0] << " "
		       << s.avr_analytic[1] << " "
		       << s.avr_numeric[0] << " "
		       << s.avr_numeric[1] << "\n";
      }
    }
    if (!quiet) cout << " done" << endl;
  } // end of loop over steps
  if (!quiet && output_dir != "") {
    cout << "\nResults are written to " << output_dir << " directory.\n\n"
	 << "The summary is in \"summary.txt\" in the format:\n"
	 << "<step> <IP> <beam-beam/no beam-beam luminosity correction>\n"
	 << "no beam-beam: <analytic> overlap, <numeric/analytic ratio> and <its error>\n"
	 << "no beam-beam: <X>, <Y> numeric center-of-mass shift\n"
	 << "   beam-beam: <X>, <Y> analytic, <X>, <Y> numeric shift\n\n"
	 << "No beam-beam numeric/analytic error is roughly estimated from turn-by-turn variations,\n"
	 << "available only if \"integrals_per_turn\" option is set in \"output\".\n"
	 << "Otherwise this error is set to \"nan\". Similarly, numeric <X>, <Y> shifts are calculated\n"
	 << "only if \"avr_xy_per_particle\" option is chosen, and set to \"nan\" otherwise.\n"
	 << "Without beam-beam these shifts should be close to zero.\n";
  }
}
