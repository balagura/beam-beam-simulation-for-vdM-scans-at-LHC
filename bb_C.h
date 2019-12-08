#ifndef bb_c_h 
#define bb_c_h 1

#include <stdbool.h>
#include "bb_C_CPP.h"

/*
  B*B simulation of beam-beam effects in van der Meer scans at LHC.

  Author V. Balagura, balagura@cern.ch (Nov 2019)
*/

struct Multi_Gaussian_C { /* weighted sum of Gaussians with common center */
  int n;
  double *sig, *w; /* arrays of length "n" with sigmas, in um, and the
		      weights. The latter can be given not normalized */
};
typedef struct Multi_Gaussian_C Multi_Gaussian_C; /* to suit both C and C++ */
/* ------------------------------ Input ------------------------------ */
struct Kicked_C {    /* "kicked" ie. simulated bunch */
  double momentum; /* momentum in GeV */
  int Z; /* charge of the particles in units of the proton charge */
  int ip; /* interaction point (counting from zero) for which multi-Gaussian
	     sigmas are given */  
  double *beta; /* [coor*n_ip + ip], ie. first X, then Y,beta-star function at the
		   simulated interaction points, in meters */
  double *next_phase_over_2pi; /* [coor*ip + ip], betatron absolute phases/2pi */
  Multi_Gaussian_C gaussian[2]; /* X,Y-sigmas and weights at "ip" */
  bool exact_phases; /* default = false, see below */
};

struct Kickers_C { /* "kicker" bunches creting the field */
  int Z;
  double *n_particles; /* number of particles in the kickers */
  Multi_Gaussian_C *gaussian; /* [coor*n_ip + ip], at kicker's ip */
  double *position;  /* [coor*n_step*n_ip + step*n_ip + ip] - wrt. kicked bunch  */
};

struct Sim_C { /* Explanations can be found in the end of this file */
  int n_ip; /* number of interaction points */
  int n_step; /* number of kicker positions */
  int n_points; /* number of simulated macro-particles */
  int n_turns[PHASES]; /* for NO_BB, ADIABATIC, STABILIZATION, BB phases */
  char *kick_model; /* "precise", "average" or "precise.minus.average" */
  int n_sigma_cut;  /* default: 5 */
  int density_and_field_interpolators_n_cells_along_grid_side[2]; /* for X/Y
      1D densities and 2D field; default: 500, 500. If 0 - no interpolation */
  int n_random_points_to_check_interpolation; /* default: 10000, if 0 - skipped */
  int select_one_turn_out_of; /* for output option "points" */
  long int seed;     /* random seed */
  char *output_dir;  /* subdirectory name where output will go */
  char* output;      /* white-space separated list of output options */
};

typedef struct Kicked_C Kicked_C;
typedef struct Kickers_C Kickers_C;
typedef struct Sim_C Sim_C;
typedef struct BB_Summary_Per_Step_IP BB_Summary_Per_Step_IP;

#ifdef __cplusplus
extern "C" {
#endif
  void beam_beam(const Kicked_C* kicked, const Kickers_C* kickers, const Sim_C* sim,
		 BB_Summary_Per_Step_IP* summary, bool quiet);
  /* summary[n_ip * n_step] should be pre-allocated,
     summary[ip * n_step + step]" will  contain BB_Summary_Per_Step_IP for (ip, step),
     "quiet" controls printing to cout */

#ifdef __cplusplus
}
#endif


/*
 -------------------- Detailed explanations --------------------

 Kicked_C::next_phase_over_2pi_x,y
 By definition, at the first IP the phase is zero, so the arrays start from
 the second IP, this is why they are called "next" phase over 2pi. The last
 value in the array is the tune. Note, two beams are rotating in opposite
 directions, and the absolute phase should increase in the direction of the
 beam.

 Kicked_C::gaussian_x,y
 multi-Gaussian sigmas and weights along X or Y at the interaction point "ip";
 at other IPs sigmas are recalculated as sigma(ip2) = sigma(ip) *
 sqrt(beta(ip2) / beta(ip))

 Kicked_C::exact_phases
 Phases/2pi or full tunes might be given with 2-3 digits after the comma. So,
 after 100 or 1000 accelerator turns the points return to their original
 positions (as 100* or 1000*phase/2pi become integer). To avoid this and to
 add an extra randomness, a "negligible" irrational number randomly
 distributed in the interval -exp(-8.5)...+exp(-8.5) = +/-0.0001017342 is
 added by default to all phases/2pi. This makes them irrational and opens
 otherwise closed Lissajous figures of betatron oscillations in X-Y plane. In
 reality the phases/2pi are always irrational.

  If this is not desired, "exact_phases = true" should be set. Then all
  phases/2pi are used "as is". Note, this might slightly deteriorate the
  quality of the simulation if phases are given with only two digits, then it
  is preferable to set "exact_phases = false"

  Kickers_C::gaussian_x,y

  contrary to the kicked bunch, the kicker multi-Gaussian sigmas should be
  given directly at the corresponding IP (where the kick is produced),
  ie. without any extrapolation via beta-functions. In other words, such
  kicker sigmas, in general, can differ from the sigmas measured at kicked
  "ip" 

  Kickers_C::n_positions_x,y
  Kickers_C::position_x,y
  the lengths and the arrays of X, Y kicker bunch center coordinates
  w.r.t. the kicked bunch, in um. Every (X, Y) pair corresponds to one van der
  Meer scan step.  The lengths of all vectors should be either one or
  maximal. If it is one, the corresponding coordinate is assumed to be
  constant during the scan.
 
  ----------------------------------------------------------------------
                            Simulation
 ----------------------------------------------------------------------
 The traced points of the "kicked" bunch are selected in the following way.

 The bunch density factorizes in X and Y. Its projections to X and Y are two
 multi-Gaussian distributions. Because of the circular motion in X-X' and
 Y-Y' planes (where X',Y' denote the angular coordinates scaled by the
 accelerator beta-function, eg. X' = dX/dZ * beta), the projections to X'
 (Y') are identical to X (Y). So, eg. a single Gaussian in X makes a
 two-dimensional Gaussian with equal sigmas in X-X' plane, and a
 multi-Gaussian makes a multi two-dimensional Gaussian with the same weights.

 The four-dimensional bunch density is sampled in a two-dimensional rX-rY
 grid of the radii in the X-X', Y-Y' planes at the first interaction point.

 Each grid side (rX or rY) extends from 0 to a maximal radius rX,Y_max. The
 latter is chosen such that the circle with the radius rX_max (or rY_max)
 contains the the same fraction of a multi two-dimensional Gaussian as that
 of a single two-dimensional Gaussian inside the circle with the radius
 "n_sigma_cut".

 The number of grid lines in rX and rY are the same and equal to
 int(sqrt("n_points")). Finally, only rX-rY points inside the ellipse
 inscribed to the rectangle (rX_max X rY_max) are used for the
 simulation. For every selected rX-rY pair, one point at the X-X', Y-Y'
 circles with these radii and random phases is traced in the accelerator.

 There are several "phases" in the simulation.

 First phase is without beam-beam. It allows to calculate numerically the
 undisturbed overlap integral, compare it with the exact analytic formula and
 estimate the bias of the numerical integration. The final correction is then
 calculated as a ratio of the numerical integrals with (last phase) and
 without (first phase) the beam-beam interaction. The bias is cancelled in
 the ratio at least partially and this potentially allows to improve the
 precision.

 During the second phase the beam-beam interaction is switched on
 "adiabatically": linearly with the turn number from zero to its nominal
 value. If this phase is omitted (zero turns) the switch is abrupt: from
 zero to nominal.

 After beam-beam is fully switched on, there might be a need to wait for the
 stabilization. This is the third phase. With the adiabatic switch on
 (eg. during 1000 turns) the stabilization phase can normally be omitted and
 the corresponding number of turns can be set to zero.

 The final, fourth phase is with beam-beam. Only this phase is used to
 calculate the perturbed overlap integral and the luminosity correction.

 To speed up the simulation, the X- and Y-densities of the kicker bunch at a
 given point are determined using linear (one-dimensional) interpolations,
 while the generated field E - using a bilinear (two-dimensional)
 interpolation. For that, the precise density and the field are precalculated
 at the grid of points and interpolated in between.

 The interpolation grid is chosen to minimally cover the ranges +/- 1.1 *
 rX_max, +/- 1.1 * rY_max around the kicked bunch center but viewed from the
 kicker bunch center for all set kicker separations. Here, rX,Y_max are the
 maximal simulated radii in X-X' and Y-Y' planes.  Therefore, the grid is
 sufficiently wide to cover all simulated without beam-beam particle
 trajectories (circles), while with 10% margins (with the factor 1.1) it
 likely covers also the trajectories deformed by the beam-beam
 interaction. If not, ie. if some point goes beyond the interpolation grid,
 the program falls back to the calculation of the exact field instead of the
 interpolation.

 The linear (bilinear) grid is a raw (a matrix) of N (NxN) cells and has the
 total number of points N+1 ((N+1)*(N+1)). N is defined separately for the
 one-dimensional (X,Y)-density and for the two-dimensional field
 interpolators by the following parameter
 "density_and_field_interpolators_n_cells_along_grid_side". If the
 interpolation is not desired, the corresponding N should be set to zero,
 like "0 500". The first number is for the density, the second - for the
 field, so the density will be calculated using exact formulas while the
 field will be interpolated with N=500.

 The interpolators are created not for each kicker bunch but for each group
 of identical bunches. Ie. if there are identical kickers, the interpolators
 are reused to save memory.

 The accuracy of the interpolation can be estimated if the parameter
 "n_random_points_to_check_interpolation" below is given. This is performed
 by measuring the maximal and average absolute mismatches between the
 interpolated and the exact values at
 "n_random_points_to_check_interpolation" distributed inside the
 interpolation grid. The mismatches are then normalized to the maximal
 absolute exact value and printed out.

 If the interpolation was switched off by setting one or both
 "density_and_field_interpolators_n_cells_along_grid_side" values to zero,
 the corresponding "n_random_points_to_check_interpolation" has no
 effect. The check is not performed either if
 "n_random_points_to_check_interpolation" is set to zero.

 If the kicker density and the field are interpolated, the simulation of
 multi-Gaussian (round or elliptical) bunches takes about the same time
 as the simplest single-Gaussian bunches.

 The results of the simulation will be fully reproducible if the random seed
 is specified. Otherwise, the current time will be used as a seed, and the
 results will be not reproducible.

 For the debugging purposes one can redefine the kick formula by setting
 "kick_model" parameter to one of the following:

  precise
  precise.minus.average
  average

 "precise" is the default. In this case the exact kick formula is used.

 "average" means constant X,Y-independent kick equal to the value averaged
 over the kicked bunch, ie. to the sum of the kicks of all bunch particles
 divided by their number. This average kick depends only on the bunch
 separation, but not on X,Y. For the Gaussian bunches it can be calculated as
 the action of the kicker bunch with (Capital) sigma = sqrt(sigma1^2 +
 sigma2^2) on one particle placed at the kicked bunch center. The constant
 kick only shifts the first bunch center ("orbit shift"), but do not modify
 its shape, so the resulting luminosity change can be computed using analytic
 formula.

 "precise.minus.average" means that the kick is chosen to be the difference
 "precise" - "average".  This model is implemented to demonstrate that the
 "precise" beam-beam luminosity corrections coincide with (and can be
 decoupled to) the sum of the corrections obtained with
 "precise.minus.average" and "average".

 ----------------------------------------------------------------------
                            Output
 ----------------------------------------------------------------------
 "output" controls what should be printed to the "output_dir". It is a
 white-space separated list of options. For every option XXXX.XXXX one
 compressed file will be generated under the name "XXXX_XXXX.txt.gz".

 Options:
   "integrals.per.turn" overlap integrals per turn (averaged over particles).
                        Format: step ip ip.kicker.x ip.kicker.y i.turn integral

   "avr.xy.per.turn" average kicked bunch center per turn (averaged over
                     particles)
                     Format: step ip i.turn average.x average.y


   "integrals.per.particle" overlap integrals per particle, per phase
                            (averaged over turns in the given phase),
                            assuming that the kicked bunch is squashed to
                            this particle, ie. the bunch density is a
                            delta-function at its position.
                            Format: step ip phase i.particle integral

   "avr.xy.per.particle"   average kicked bunch center per particle, per phase
                           as above.
                           Format: step ip phase i.particle avr.x avr.y


   "points"  the traced positions of the kicked bunch particles before
             accelerator "i.turn", note, this file might be very long.
             Format: step ip i.turn i.particle X X' Y Y'

 Storage of the traced particle positions (if the output option "points"
 above is specified) before every accelerator turn might require too much disk
 space. Instead, one can "select_one_turn_out_of" N turns. Eg. if this
 parameter is set to 1000, only the positions before i.turn = 999, 1999, 2999,
 ...  (counting from zero) will be stored. If "points" option is not
 requested in "output", "select_one_turn_out_of" has no effect.

 The integers i.turn and i.particle are counted from zero.

 If "output_dir" is empty ("") or this subdirectory exists, nothing will be
 printed regardless of the options in "output". Otherwise "output_dir" will
 be created with the files "XXXX_XXXX.txt.gz" and in addition the following
 files will be written regardless of "output"

 "rx_ry_weights.txt.gz" for every simulated particle: its initial rX and rY radii
                   (in X,X' and Y,Y' phase spaces) at ip and its "weight"
                   Format: i.particle ip rx ry weight

 "kicker_positions.txt"  coordinates of the kicker bunch centers at all IPs
                         with respect to the kicked bunch center
                         Format: step ip x y
 "summary.txt"   the main results of the simulation including the overlap integrals
                 and the corresponding beam-beam corrections.
   Format: <step> <IP> <beam-beam/no beam-beam luminosity correction>
    no beam-beam: <analytic> overlap integral, <numeric over anaytic> ratio and <its error> 
    (note, error = "nan" if "integrals.per.turn" option is not requested in "output")
    numerically calculated <X>, <Y> center-of-mass shift of the kicked bunch without beam-beam
    (should be zero) and <X>, <Y> shift with beam-beam ("nan" if "avr.xy.per.particle" was
    not requested in "output"), the same shift calculated for <X> and <Y> analytically
*/


#endif
