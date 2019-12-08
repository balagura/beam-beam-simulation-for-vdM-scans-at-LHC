//
//  B*B simulation of beam-beam effects in van der Meer scans at LHC.
//
//  Author V. Balagura, balagura@cern.ch (Dec 2019)
//
#include <fstream>
#include "bb_config.hh"

// Runs B*B simulation with the input parameters from the "Config" file
// specified as the first argument. The results are always stored in
// "output_dir". If this parameter is not given, "output_dir" is obtained
// by dropping .txt from the end of the "Config" file name (an error is
// produced if it does not have *.txt extension).
//
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
