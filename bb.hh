#ifndef bb_hh
#define bb_hh 1

//
//  B*B simulation of beam-beam effects in van der Meer scans at LHC.
//
//  Author V. Balagura, balagura@cern.ch (Dec 2019)
//

#include <vector>
#include <string>
#include <iostream>
#include "bb_C_CPP.h"

using namespace std;

struct Gaussian { // single Gaussian component with
  double sig, w;  // sigma and weight
};
// ------------------------------ Input ------------------------------
struct Kicked { // "kicked" ie. simulated bunch
  double momentum;  // in GeV
  int Z;            // in units of the proton charge
  int ip;           // interaction point (counting from zero) of "gaussian"
  vector<double> beta[2];                // [coor][ip], beta-star
  vector<double> next_phase_over_2pi[2]; // [coor][ip], absolute phases
  bool exact_phases;                     // add irrationality to phases?
  vector<Gaussian> gaussian[2];          // [coor][gauss]
};
struct Kickers { // "kicker" bunches creting the field
  int Z;
  vector<double> n_particles;            //       [ip]
  vector<vector<Gaussian> > gaussian[2]; // [coor][ip][gauss]
  vector<vector<double> >   position[2]; // [coor][ip][step]
};
struct Sim {
  int n_points;        // N simulated macro-particles
  int n_turns[PHASES]; // for NO_BB, ADIABATIC, STABILIZATION, BB phases
  string kick_model;   // "precise", "average" or "precise.minus.average"
  int n_sigma_cut;     // default: 5
  int density_and_field_interpolators_n_cells_along_grid_side[2]; // for X/Y
  // 1D densities and 2D field; default: 500, 500; if 0 - no interpolation
  int n_random_points_to_check_interpolation; // default: 10000, if 0 - skipped
  int select_one_turn_out_of; // for output option "points"
  long int seed;              // random seed
  string output_dir;          // subdirectory name where output will go
  string output;              // white-space separated list of output options
};
// ------------------------------ Main simulation function ------------------------------
void beam_beam(const Kicked& kicked,
	       const Kickers& kickers, const Sim& sim,
	       vector<vector<BB_Summary_Per_Step_IP> >* summary, // results: summary[ip][step]
	       bool quiet = false);      // print input and progress to cout?

void print(ostream& os, const Kicked& kicked, const Kickers& kickers, const Sim& sim);

#endif
