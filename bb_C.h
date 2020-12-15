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

  Calls B*B simulation from C instead of C++

*/
#ifndef bb_c_h 
#define bb_c_h 1

#include <stdbool.h>
#include "bb_C_CPP.h"

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
  double *beta; /* [coor*n_ip + ip], ie. first X, then Y, beta-star function at the
		   simulated interaction points, in meters */
  double *next_phase_over_2pi; /* [coor*n_ip + ip], betatron absolute phases/2pi */
  bool exact_phases; /* default = false, see below */
  double *sig_z_projection;    /* [coor*n_ip + ip], first X, then Y, sigmaZ transverse projection */
  Multi_Gaussian_C gaussian[2]; /* X,Y-sigmas and weights at "ip" */
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
  char *kick_model; /* "precise", "average" or "precise_minus_average" */
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
 The phase at IP=0 is zero by definition, therefore, the "x" and "y" vectors
 should start from the phase at IP=1, this is why this parameter is called
 "next" phase over 2pi. The order of IPs = 0,1,2, ... is determined by the
 order in which the kicked bunch collides with the kickers in the accelerator,
 so, eg. it is opposite for two beams rotating in the opposite directions. The
 lengths of "x" and "y" vectors should be equal to the total number of
 simulated IPs. The last value in the array is the tune.

 The LHC phases/2pi in vdM scans are given below for reference. They were
 calculated by Guido Sterbini using MAD-X simulation in 2019 and reported at
 https://indico.cern.ch/event/836681/contributions/3507678/attachments/1884149/3105236/2019_07_22.pdf
 indico: https://indico.cern.ch/event/836681/
           Qx1       Qy1       Qx2       Qy2
  IP1   0.000000  0.000000  0.000000  0.000000
  IP2   8.295954  7.669167  8.272802  7.957720
  IP5  31.975690 29.648602 31.984398 29.761319
  IP8  56.064831 51.017069 55.799012 51.715754
  IP1  64.310000 59.320000 64.310000 59.320000
 
 If the simulated kicked bunch is in beam 1 and the kickers are in beam 2:
 next_phase_over_2pi_x = 8.295954 31.975690 56.064831 64.310000
 next_phase_over_2pi_y = 7.669167 29.648602 51.017069 59.320000

 If the simulated kicked bunch is in beam 1 and the kickers are in beam 2,
 one should take into account that beam2 goes in the opposite direction, ie.
 in the order IP1->8->5->2->1. So, one needs to take (full tune - table column)
 in reverse order,
  ie. reversed 64.31 - (0, 8.272802, 31.984398, 55.799012)
  or  reversed 59.32 - (0, 7.957720, 29.761319, 51.715754):
 next_phase_over_2pi_x = 8.510988 32.325602 56.037198 64.310000
 next_phase_over_2pi_y = 7.604246 29.558681 51.362280 59.320000

 Kicked_C::gaussian_x,y
 multi-Gaussian sigmas in um and the corresponding weights of the kicked bunch
 density in "x" and "y" at the specified Kicked_C::ip number. All Gaussians
 should have a common mean at zero. The weights might be given not normalized.
 Kicked bunch sigmas at other IP2 are calculated via beta-function as
 sigma(ip) * sqrt(beta(ip) / beta(IP2)). Note, that contrary to that, the kicker
 sigmas should be specified directly at the IP of the beam-beam interaction,
 without any extrapolation via beta. In other words, such kicker sigmas, in
 general, can differ from the ones at this Kicked_C::ip.
 
 Kicked_C::exact_phases
 If values in "next_phase_over_2pi" are given with 2 (3) digits after the
 comma, after 100 (1000) turns the points return almost to their original
 positions if the beam-beam effect is small (as 100* or 1000 * phases/2pi
 become integer). The simulation, therefore, probes the same part of the
 phase space over and over again, and extra turns almost do not improve the
 convergence. To avoid this and to add extra randomness in the transverse
 particle trajectories, "negligible" irrational numbers randomly distributed
 in the interval -exp(-8.5)...+exp(-8.5) = +/-0.0001017342 are added to all
 phase / 2pi values if "exact_phases" below is FALSE (default). This makes
 them irrational and opens otherwise closed Lissajous figures of betatron
 oscillations in X-Y plane. If this is not desired, "exact_phases" should be
 set to TRUE. Then all phases/2pi are used exactly as they are given.

 Kicked_C::sig_z_projection
 Specifies the longitudinal bunch sigma projections to "x" or "y"
 (perpendicular to the kicker beam), at all simulated interaction points, in
 microns. If the beam crossing angle is zero at the given point, the
 projections vanish. Otherwise eg. "x"-projection is equal to alpha * sigmaZ *
 cos(beta_x), where alpha is the crossing angle between the kicker and the
 vector difference v1 - v2, where v1,v2 are the beam velocities, and beta_x is
 the angle between x and the projection of the crossing plane to the x-y
 plane. Note, these contributions of the longitudinal spread to the transverse
 widths are important only for calculating luminosities and do not affect the
 transverse dynamics of the particles. For example, if the luminosity at some
 point ip is not needed, sigma_z_projection$x[ip], sigma_z_projection$y[ip]
 can be set arbitrarily, eg. to zeros.

 Kickers_C::gaussian_x,y
 Defines multi-Gaussian sigmas in um and the corresponding weights of the
 kicker bunch densities in "x" and "y" at IP = 0,1,2... .  All Gaussians of
 one kicker should have a common mean. The weights might be given not
 normalized. Note, that contrary to the kicked bunch, the kicker sigmas should
 be specified directly at the IP of the beam-beam interaction, without any
 extrapolation via beta. In other words, such kicker sigmas, in general, can
 differ from the ones at the Kicked_C::ip where the kicked sigmas are
 specified.

 Kickers_C::n_positions_x,y
 Kickers_C::position_x,y
 the lengths and the arrays of kicker "x" and "y" positions at the given IP
 in um in the frame where the bunch center of the kicked bunch before
 beam-beam was at zero. Each vector can have either one or "n_step" kicker
 positions for a given IP, where "n_step" is the number of scan points in van
 der Meer scan. If the vector contains only one position, it is "recycled"
 and the kicker is assumed to be immovable during the scan along the given
 coordinate.
 
 -------------------- Description of the simulation --------------------

 ....... Initial distribution of macro-particles:
 The traced macro-particles of the "kicked" bunch are selected in the
 following way.

 It is assumed that the initial bunch densities unperturbed by beam-beam
 factorize in X and Y and their X,Y-projections are multi-Gaussian
 distributions with common centers. Because of the circular motion in X-X'
 and Y-Y' planes (where X',Y' denote the angular coordinates scaled by the
 accelerator beta-function, eg. X' = dX/dZ * beta), the projections to X'
 (Y') are identical to X (Y). So, eg. a single Gaussian in X makes a
 two-dimensional Gaussian with equal sigmas in X-X' plane, and a
 multi-Gaussian makes a multi two-dimensional Gaussian with the same weights.

 The four-dimensional kicked bunch density is sampled internally in a
 two-dimensional rX-rY grid of the radii in the X-X', Y-Y' planes at the
 interaction point IP = 0. Each grid side (rX or rY) extends from 0 to a
 maximal radius rX,Y_max. The latter is chosen such that the circle with the
 radius rX_max (or rY_max) contains the same fraction of a multi
 two-dimensional Gaussian as that of a single two-dimensional Gaussian inside
 the circle with the radius "n_sigma_cut" (5 by default). The grid sides
 [0...rX_max] and [0...rY_max] are divided into int(sqrt("n_points"))
 intervals and the grid lines are drawn through their centers. Finally, only
 rX-rY points inside the ellipse inscribed to the rectangle [0...rX_max] X
 [0...rY_max] are used for the simulation. The number of selected points is,
 therefore, less than "n_points". For every rX-rY pair selected in this way,
 one four-dimensional point at the X-X', Y-Y' circles with these radii and
 random phases is chosen as the initial position of one macro-particle. Its
 weight is chosen such that the ensemble of all macro-particles represents
 the initial bunch density. The macro-particles are then traced individually
 in the accelerator in separate threads. This allows to achieve significant
 speed-up of the simulation in the computers with multi-core processors, the
 gain is almost proportional to the number of cores.
 
 ....... Four phases in the simulation:
 First phase is without beam-beam. It allows to calculate numerically the
 undisturbed overlap integral, compare it with the exact analytic formula and
 estimate the bias of the numerical integration. The final correction is then
 calculated as a ratio of the numerical integrals with (last phase) and
 without (first phase) the beam-beam interaction. The bias is cancelled in
 the ratio at least partially and this potentially allows to improve the
 precision.

 During the second phase the beam-beam interaction is switched on
 "adiabatically": linearly with the turn number from zero to its nominal
 value. If this phase is omitted (zero turns) the switch is abrupt: from zero
 to nominal.

 After beam-beam is fully switched on, there might be a need to wait for the
 stabilization. This is the third phase. With the adiabatic switch on
 (eg. during 1000 turns) the stabilization phase can normally be omitted and
 the corresponding number of turns can be set to zero.

 The final, fourth phase is with beam-beam. Only this phase is used to
 calculate the perturbed overlap integral and the luminosity correction.

 "n_turns" parameter specifies the number of turns in all 4 phases.
 
 ....... Interpolation:
 To speed up the simulation, the X- and Y-densities of the kicker bunch at a
 given point are determined using linear (one-dimensional) interpolations,
 while the generated E-field - using a bilinear (two-dimensional)
 interpolation. For that, the precise density and the field are precalculated
 at the grid of points and interpolated in between.

 The interpolation grid is chosen to minimally cover the ranges +/- 1.1 *
 rX_max, +/- 1.1 * rY_max around the kicked bunch center but viewed from the
 kicker bunch center for all set kicker separations. Here, rX,Y_max are the
 maximal simulated radii in X-X' and Y-Y' planes.  Therefore, the grid is
 sufficiently wide to cover all simulated without beam-beam particle
 trajectories (circles), while with 10\% margins (with the factor 1.1) it
 likely covers also the trajectories deformed by the beam-beam
 interaction. If not, ie. if some point goes beyond the interpolation grid,
 the program falls back to the calculation of the exact field instead of the
 interpolation.

 The linear (bilinear) grid is a raw (a matrix) of N (NxN) cells and has the
 total number of points N+1 ((N+1)*(N+1)). N is defined separately for the
 one-dimensional (X,Y)-density and for the two-dimensional field
 interpolators by the parameter
 "density_and_field_interpolators_n_cells_along_grid_side". If the
 interpolation is not desired, the corresponding N should be set to zero,
 like "0 500". The first number is for the density, the second - for the
 field, so in this case the density will be calculated using exact formulas
 while the field will be interpolated with N=500.

 If the kicker densities and the fields are interpolated, the simulation of
 complex multi-Gaussian elliptical bunches and simple round single-Gaussian
 ones takes about the same time.

 The interpolators are created not for each kicker bunch but for each group
 of identical bunches. Ie. if there are identical kickers, the interpolators
 are reused to save memory.

 The accuracy of the interpolation is estimated if the parameter
 "n_random_points_to_check_interpolation" is larger than zero. This is
 performed by measuring the maximal and average absolute mismatches between
 the interpolated and the exact values at
 "n_random_points_to_check_interpolation" distributed inside the
 interpolation grid. The mismatches are then normalized to the maximal
 absolute exact value and printed to the file
 "output_dir"/interpolation_precision.txt.

 If the interpolation was switched off by setting one or both
 "density_and_field_interpolators_n_cells_along_grid_side" values to zero,
 the corresponding "n_random_points_to_check_interpolation" has no
 effect. The check is not performed either if
 "n_random_points_to_check_interpolation" is zero or negative.

 ....... Kick model:
 For debugging purposes one can redefine the kick formula by setting
 "kick.model" parameter to one of the following:
      precise
      average
      precise_minus_average
 "precise" is the default. In this case the exact kick formula is used.

 "average" means constant X,Y-independent kick equal to the value averaged
 over the kicked bunch. Ie. it is equal to the sum of the kicks of all bunch
 particles divided by their number. For the Gaussian bunches it can be
 calculated as the kick exerted by the kicker bunch with (Capital) sigma =
 sqrt(sigma1^2 + sigma2^2) on one particle placed at the kicked bunch
 center. The constant kick only shifts the kicked bunch as a whole, but does
 not modify its shape, so the resulting luminosity change can be computed
 using analytic formula.

 "precise_minus_average" kick is simply calculated as the difference
 "precise" - "average". This model allows to compare the "precise" beam-beam
 luminosity correction with the sum of the corrections obtained with
 "precise_minus_average" and "average".

 ....... Output:
 "output" parameter controls what should be calculated and printed after the
 simulation. It is a white-space separated list of options. If "output_dir"
 is not empty (""), for every option "XXXX" one compressed file will be
 created under the name "output_dir"/XXXX.txt.gz. The possible options are
 listed below:

 "integrals_per_turn" - overlap integrals averaged over particles, per
   turn.
 Format: step ip ip_kicker_x ip_kicker_y i_turn integral.

 "avr_xy_per_turn" - kicked bunch center averaged over particles, per
   turn.
 Format: step ip i_turn average_x average_y.

 "integrals_per_particle" - overlap integrals averaged over turns, per
   particle, per phase, assuming that the kicked bunch is squashed to this
   particle, ie. the bunch density is a delta-function at its position.
 Format: step ip phase i_particle integral.

 "avr_xy_per_particle" - kicked bunch center averaged over turns, per
   particle, per phase as above.
 Format: step ip phase i_particle average_x average_y.

 "points" -  the traced positions of the kicked bunch particles before
   accelerator "i_turn", note, the generated file might be very long.
 Format: step ip i_turn i_particle X X' Y Y'.

 The integers i_turn and i_particle are counted from zero.

 In addition, the files
 "output_dir"/input.txt,
 "output_dir"/rx_ry_weights.txt.gz,
 "output_dir"/kicker_positions.txt,
 "output_dir"/interpolation_precision.txt and
 "output_dir"/summary.txt are printed out regardless of options in
 "output" with the following content:

 "output_dir"/input.txt - all input parameters.
 "output_dir"/rx_ry_weights.txt.gz - for every simulated particle: its
             initial rX and rY radii (in X,X' and Y,Y' phase spaces) at ip
             and its "weight".
 Format: i_particle ip rx ry weight.

 "output_dir"/kicker_positions.txt - coordinates of the kicker bunch
     centers with respect to the nominal (not perturbed by the
     beam-beam effect) kicked bunch center.
 Format: step ip x y.

 "output_dir"/interpolation_precision.txt - was explained in
 "Interpolation".
 
  "output_dir"/summary.txt - the main results of the simulation
     including the overlap integrals and the corresponding beam-beam
     corrections. Format:
 <step> <IP> <beam-beam/no beam-beam luminosity correction>
 no beam-beam: <analytic> overlap, <numeric/analytic ratio> and <its error>
 no beam-beam: <X>, <Y> numeric center-of-mass shift
    beam-beam: <X>, <Y> analytic, <X>, <Y> numeric shift.
 
 No beam-beam numeric/analytic error is roughly estimated from turn-by-turn
 variations, available only if "integrals_per_turn" option is set in
 "output". Otherwise this error is set to "nan". Similarly,
 numeric <X>, <Y> shifts are calculated only if "avr_xy_per_particle"
 option is chosen, and set to "nan" otherwise. Without beam-beam these
 shifts should be close to zero.

 The same summary data are returned by "beam_beam" function itself: after

   summary = beam_beam(kicked, kickers, sim, quiet)

 summary[ip][step] contains "Summary" structure (namedtuple) defined in
 "BxB.py" with the results. If only the minimal information on the luminosity
 correction is required from the simulation, it is better to set "output" to
 an empty string ("") or to "integrals_per_particle", as any other option in
 "output" requires extra CPU time. The integrals per particle are calculated
 in any case (but the file "output_dir"/integrals_per_particle.txt.gz is
 printed only if the corresponding option is explicitly set in "output").
 
*/


#endif
