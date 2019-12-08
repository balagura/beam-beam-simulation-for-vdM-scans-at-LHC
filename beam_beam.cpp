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
// C++ main program running B*B simulation. Reads the input parameters from
// the "Config" file specified as the first argument. The results are written
// to "output_dir". If this parameter is not given, "output_dir" is obtained
// by dropping .txt from the end of the "Config" file name (an error is
// produced if in this case it does not have *.txt extension).
//
#include <fstream>
#include "bb_config.hh"

int main(int argc, char** argv) {
  string config_file;
  if (argc != 2) {
    cout << "Usage: <program name> <configuration file>\n";
    return 0;
  }
  ifstream conf_file(argv[1]);
  Config c(conf_file);
  Kicked kicked;
  Kickers kickers;
  Sim sim;
  beam_beam_config(c, &kicked, &kickers, &sim);
  if (sim.output_dir == "") { // if not given, obtain it by dropping
    string inp(argv[1]);      // .txt from the configuration file name
    size_t n = inp.size();
    if (n > 4 && inp.substr(n-4) == ".txt") {
      sim.output_dir = inp.substr(0, n-4);
    } else {
      cerr << "If \"output_dir\" is not specified, it is obtained by dropping \".txt\" extension\n"
	   << "from the end of the configuration file name which, therefore, must have the form\n"
	   << "*.txt (" << inp << " is given)\n";
      return 1;
    }
  }
  vector<vector<BB_Summary_Per_Step_IP> > summary;
  beam_beam(kicked, kickers, sim, &summary);
  return 0;
}
