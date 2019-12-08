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
#ifndef bb_config_hh
#define bb_config_hh 1

#include "bb.hh"
#include "config.hh"

// Fills the tripple (kicked, kickers, sim) with the input parameters
// specified in "Config" object
//
void beam_beam_config(Config& c,
		      Kicked* kicked, Kickers* kickers, Sim* sim);

#endif
