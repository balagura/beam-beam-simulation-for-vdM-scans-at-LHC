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
// Reads configuration "Config" parameters from istream and creates Kicked,
// Kickers and Sim structures to run B*B simulation
//
#include "bb_config.hh"

void beam_beam_config(Config& c,
		      Kicked* kicked, Kickers* kickers, Sim* sim) {
  int n_ip = c.vd("kicked.beta.x").size();
  kicked->momentum = c["kicked.momentum"];
  kicked->Z = c["kicked.Z"];
  kicked->ip = c["kicked.ip"];
  kicked->exact_phases = c.defined("kicked.exact_phases") ? (c.s("kicked.exact_phases") == "TRUE") : false;
  kickers->Z = c["kickers.Z"];
  kickers->n_particles = c.vd("kickers.n_particles");
  auto fill_g = [&c](const string& name, vector<Gaussian>& g) {
		  vector<double>
		    sig = c.vd(name + ".sig"),
		    w   = c.vd(name + ".w");
		  if (sig.size() != w.size()) {
		    cerr << "given numbers of Gaussian sigmas and weights differ\n";
		    exit(1);
		  }
		  g.resize(sig.size());
		  for (size_t gauss=0; gauss<sig.size(); ++gauss) {
		    g[gauss].sig = sig[gauss];
		    g[gauss].w   = w  [gauss];
		  }
		};
  auto c_name = [](const string& obj, const string& field, int coor) {
		 return obj + "." + field + "." + string("xy")[coor];
	       };
  for (int coor=0; coor<2; ++coor) {
    kicked->beta[coor] = c.vd(c_name("kicked", "beta", coor));
    kicked->next_phase_over_2pi[coor] = c.vd(c_name("kicked", "next_phase_over_2pi", coor));
    fill_g(c_name("kicked", "gaussian", coor), kicked->gaussian[coor]);
    kickers->gaussian[coor].resize(n_ip);
    kickers->position[coor].resize(n_ip);
    for (int ip=0; ip<n_ip; ++ip) {
      fill_g(c_name("kickers", "gaussian", coor) + "." + to_string(ip),
	     kickers->gaussian[coor][ip]);
      kickers->position[coor][ip] = c.vd(c_name("kickers", "position", coor) + "." + to_string(ip));
    }      
  }
  // provide defaults for all sim.* fields
  sim->n_points = c.defined("sim.n_points") ? c("sim.n_points") : 5000;
  if (c.defined("sim.n_turns")) { 
    if (c.vl("sim.n_turns").size() != PHASES) {
      cerr << "The length of \"n_turns\" array must be " << PHASES << "\n";
      exit(1);
    } else {
      copy(c.vl("sim.n_turns").begin(), c.vl("sim.n_turns").end(), sim->n_turns);
    }
  } else {
    sim->n_turns[0] = 1000; sim->n_turns[1] = 1000; sim->n_turns[2] = 0; sim->n_turns[3] = 5000;
  }
  sim->kick_model = c.defined("sim.kick_model") ? c.s("sim.kick_model") : "precise";
  sim->n_sigma_cut = c.defined("sim.n_sigma_cut") ? c("sim.n_sigma_cut") : 5;
  {
    string name = "sim.density_and_field_interpolators_n_cells_along_grid_side";
    if (c.defined(name)) {
      if (c.vl(name).size() != 2) {
	cerr << "The length of \"" << name << "\" array must be 2\n";
	exit(1);
      } else {
	copy(c.vl(name).begin(), c.vl(name).end(),
	     sim->density_and_field_interpolators_n_cells_along_grid_side);
      }
    } else {
      sim->density_and_field_interpolators_n_cells_along_grid_side[0] = 500;
      sim->density_and_field_interpolators_n_cells_along_grid_side[1] = 500;
    }
  }
  {
    string name = "sim.n_random_points_to_check_interpolation";
    sim->n_random_points_to_check_interpolation = c.defined(name) ? c(name) : 10000;
  }
  {
    string name = "sim.select_one_turn_out_of";
    sim->select_one_turn_out_of = c.defined(name) ? c(name) : 1000;
  }
  sim->seed = c.defined("sim.seed") ? c("sim.seed") : time(NULL);
  sim->output_dir = c.defined("sim.output_dir") ? c.s("sim.output_dir") : "";
  sim->output = "";
  if (c.defined("sim.output")) {
    const vector<string>& vs = c.vs("sim.output");
    for (size_t i=0; i<vs.size(); ++i) sim->output += vs[i] + " ";
  }
}
