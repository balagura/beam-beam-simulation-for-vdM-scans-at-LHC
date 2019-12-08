#ifndef bb_C_CPP_hh
#define bb_C_CPP_hh 1

/*
   B*B simulation of beam-beam effects in van der Meer scans at LHC.
   Common declarations for C and C++ interfaces.
  
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
