// Simple class to keep various integer parameters of B*B simulation of
// beam-beam effects in van der Meer scans at LHC.
//
// Author V. Balagura, balagura@cern.ch (Nov 2019)
//

#ifndef n_hh 
#define n_hh 1

#include <config.hh>

using namespace std;

struct N {
  N(Config& c);
  int ip() const { return N_ip; }
  int points() const { return N_points; }
  int& points() { return N_points; }
  int points_max() const { return N_points_max; }
  int turns() const { return N_total_turns; }
  int turns(int phase) const { return N_turns[phase]; }
  int select_turns() const { return N_select_turns; }
  int turns_stored() const { return N_turns_stored; }
  //  int phases() const { return PHASES; }
  enum {NO_BB=0, ADIABATIC=1, STABILIZATION=2, BB=3, PHASES = 4};
protected:
  int N_ip, N_points_max, N_points, N_total_turns, N_select_turns, N_turns_stored;
  array<int, PHASES> N_turns;
};

#endif
