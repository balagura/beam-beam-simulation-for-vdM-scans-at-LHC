#ifndef bb_config_hh
#define bb_config_hh 1

//
//  B*B simulation of beam-beam effects in van der Meer scans at LHC.
//
//  Author V. Balagura, balagura@cern.ch (Dec 2019)
//
#include "bb.hh"
#include "config.hh"

// Fills the tripple (kicked, kickers, sim) with the input parameters
// specified in "Config" object
//
void beam_beam_config(Config& c,
		      Kicked* kicked, Kickers* kickers, Sim* sim);

#endif
