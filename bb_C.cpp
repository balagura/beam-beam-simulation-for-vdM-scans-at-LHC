#include "bb_C.h"
#include "bb.hh"
#include <vector>

using namespace std;

void convert(const Multi_Gaussian_C& from, vector<Gaussian>* to) {
  to->resize(from.n);
  for (int i=0; i<from.n; ++i) {
    (*to)[i].sig = from.sig[i];
    (*to)[i].w   = from.w  [i];
  }
}
  
void beam_beam(const Kicked_C* kicked, const Kickers_C* kickers, const Sim_C* sim,
	       BB_Summary_Per_Step_IP* summary, bool quiet) {
  int n_ip = sim->n_ip;
  int n_step = sim->n_step;
  Kicked k;
  Kickers ks; 
  Sim s;
  // kicked
  k.momentum = kicked->momentum;
  k.Z = kicked->Z;
  k.ip = kicked->ip;
  for (int coor=0; coor<2; ++coor) {
    k.beta[coor].resize(n_ip);
    k.next_phase_over_2pi[coor].resize(n_ip);
    for (int ip=0; ip<n_ip; ++ip) {
      k.beta[coor][ip] = kicked->beta[coor*n_ip + ip];
      k.next_phase_over_2pi[coor][ip] = kicked->next_phase_over_2pi[coor*n_ip + ip];
    }
    convert(kicked->gaussian[coor], &k.gaussian[coor]);
  }
  k.exact_phases = kicked->exact_phases;
  // kickers
  ks.Z = kickers->Z;
  ks.n_particles.resize(n_ip);
  for (int ip=0; ip<n_ip; ++ip) {
    ks.n_particles[ip] = kickers->n_particles[ip];
  }
  for (int coor=0; coor<2; ++coor) {
    ks.gaussian[coor].resize(n_ip);
    ks.position[coor].resize(n_ip, vector<double>(n_step));
    for (int ip=0; ip<n_ip; ++ip) {
      convert(kickers->gaussian[coor*n_ip + ip], &ks.gaussian[coor][ip]);
      for (int step=0; step<n_step; ++step) {
	ks.position[coor][ip][step] = kickers->position[(coor*n_ip + ip)*n_step + step];
      }
    }
  }
  // sim
  s.n_points = sim->n_points;
  for (int phase=0; phase<PHASES; ++phase) {
    s.n_turns[phase] = sim->n_turns[phase];
  }
  s.kick_model = sim->kick_model;
  s.n_sigma_cut = sim->n_sigma_cut;
  s.density_and_field_interpolators_n_cells_along_grid_side[0] =
    sim->density_and_field_interpolators_n_cells_along_grid_side[0];
  s.density_and_field_interpolators_n_cells_along_grid_side[1] =
    sim->density_and_field_interpolators_n_cells_along_grid_side[1];
  s.n_random_points_to_check_interpolation = sim->n_random_points_to_check_interpolation;  
  s.select_one_turn_out_of = sim->select_one_turn_out_of;
  s.seed = sim->seed;
  s.output_dir = sim->output_dir;
  s.output = sim->output;
  vector<vector<BB_Summary_Per_Step_IP> > sum(n_ip, vector<BB_Summary_Per_Step_IP>(n_step));
  beam_beam(k, ks, s, &sum, quiet); // sum[ip][step]
  for (int ip=0; ip<n_ip; ++ip) {
    for (int step=0; step<n_step; ++step) {
      summary[ip * n_step + step] = sum[ip][step];
    }
  }
}
