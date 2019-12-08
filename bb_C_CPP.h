/*
  B*B simulation of beam-beam effects in van der Meer scans at LHC.
  Copyright (C) 2019 Vladislav Balagura (balagura@cern.ch)

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------

  Common declarations for C and C++ interfaces:
    summary structure and definition of simulation phases

*/
#ifndef bb_C_CPP_hh
#define bb_C_CPP_hh 1
/*
   B*B simulation of beam-beam effects in van der Meer scans at LHC.
  
   Author V. Balagura, balagura@cern.ch (Dec 2019)
*/

/* ------------------------------ Summary ------------------------------ */
struct BB_Summary_Per_Step_IP {
  double correction,
    no_bb_analytic_integ,
    no_bb_numeric_over_analytic_integ,
    no_bb_numeric_over_analytic_integ_err,
    avr_analytic[2],
    avr_numeric[2],
    no_bb_avr_numeric[2];
};

enum {NO_BB=0, ADIABATIC=1, STABILIZATION=2, BB=3, PHASES = 4};

#endif
