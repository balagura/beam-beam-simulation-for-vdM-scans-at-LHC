#include <n.hh>
#include <numeric> // for accumulate

N::N(Config& c) {
  N_ip = c.defined("N.IP") ? c("N.IP") : 1;
  N_points_max = c("N.points");
  if (c.defined("N.turns")) {
    if (c.vl("N.turns").size() != PHASES) {
      cerr << "There must be " << PHASES << " numbers of accelerator turns ("
	   << c.vl("N.turns").size() << " are given). Terminating\n";
      exit(1);
    }
    for (size_t i=0; i<PHASES; ++i) N_turns[i] = c.vl("N.turns")[i];
  } else {
    N_turns = {1000, 1000, 0, 5000};
  }
  N_total_turns = accumulate(N_turns.begin(), N_turns.end(), 0.);
  // for storing xZ, yZ once per "select_turns" turns
  N_select_turns = c.defined("select.one.turn.out.of") ? c("select.one.turn.out.of") : 1000;
  N_turns_stored = N_total_turns / N_select_turns;
  // N_turns_stored = N_total_turns / N_select_turns is always a number of
  // multiples of N_select_turns in the range 1...N_total_turns, regardless
  // of whether N_total_turns % N_select_turns == 0 or not. Every multiple
  // will be written out, so one will need to store "N_turns_stored" turns.
  //
  // Set also a random seed here, in the very beginning. Note, it might be
  // required early, in Config_Mutli_XY_Gaussian_bunches() constructor for
  // an irrational addition to deltaQ
  if (c.defined("seed")) {
    srand48(c("seed"));
  } else {
    srand48(time(NULL));
  }
}
