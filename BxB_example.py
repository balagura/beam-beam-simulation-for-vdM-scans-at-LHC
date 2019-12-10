#
# B*B simulation of beam-beam effects in van der Meer scans at LHC.
# Copyright (C) 2019 Vladislav Balagura (balagura@cern.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------
#
# Example how to use python interface of B*B
#
from BxB import Kicked, Kickers, Sim, Summary, beam_beam

# -------------------- Kicked parameters --------------------
kicked = Kicked(
# Momentum in GeV.
  momentum = 3500,
# Charge of the particles in the units of the proton charge.
  Z = 1, 
# Interaction point (IP) number (counting from zero) for which the
# Gaussian sigmas are specified. Kicked bunch sigmas at other IP2 are
# calculated via beta-function as sigma(ip) * sqrt(beta(ip) /
# beta(IP2)). Note, that contrary to that, the kicker sigmas should be
# specified directly at the IP of the beam-beam interaction, without any
# extrapolation via beta. In other words, such kicker sigmas, in general,
# can differ from the ones at this kicked "ip".
  ip = 1,
# "x" and "y" vectors with the beta-function values (beta-star) at all
# simulated interaction points, in meters.
  beta = [[1.5, 1.5, 1.5, 6],
          [1.5, 1.5, 1.5, 6]],
# "x" and "y" vectors specifying transverse betatron oscillation phases
# divided by 2pi. The phase at IP=0 is zero by definition, therefore, the
# vectors should start from the phase at IP=1, this is why this
# parameter is called "next" phase over 2pi. The order of IPs = 0,1,2, ... is
# determined by the order in which the kicked bunch collides with the kickers
# in the accelerator, so, eg. it is opposite for two beams rotating in the
# opposite directions. The lengths of "x" and "y" vectors should be equal to
# the total number of simulated IPs.
  next_phase_over_2pi = [[0.31*(i+1) for i in range(4)],
                         [0.32*(i+1) for i in range(4)]],
# Vector of pairs of Gaussian sigmas in um and the corresponding weights of
# the multi-Gaussian kicked bunch density in "x" and "y". All Gaussians should
# have a common mean at zero. The weights might be given not normalized.
  sigma_weight = [[[40, 0.2], [40.0, 0.8]],
                  [[39.99, 0.3], [40, 0.7]]],
# If values in "next_phase_over_2pi" are given with 2 (3) digits after the
# comma, after 100 (1000) turns the points return almost to their original
# positions if the beam-beam effect is small (as 100* or 1000 * phases/2pi
# become integer). The simulation, therefore, probes the same part of the
# phase space over and over again, and extra turns almost do not improve the
# convergence. To avoid this and to add extra randomness in the transverse
# particle trajectories, "negligible" irrational numbers randomly distributed
# in the interval -exp(-8.5)...+exp(-8.5) = +/-0.0001017342 are added to all
# phase / 2pi values if "exact_phases" below is FALSE (default). This makes
# them irrational and opens otherwise closed Lissajous figures of betatron
# oscillations in X-Y plane. If this is not desired, "exact_phases" should be
# set to TRUE. Then all phases/2pi are used exactly as they are given.
  exact_phases = False)

# -------------------- Kickers parameters --------------------
pos_x = [10*i for i in range(21)]
kickers = Kickers(
# Charge of the particles in the units of the proton charge.
  Z = 1,
# Number of particles in all kickers
  n_particles = [8.5e10, 8.5e10, 8.5e10, 8.5e10],
# Vectors of pairs of Gaussian sigmas in um and the corresponding weights of
# the multi-Gaussian kicker bunch densities in "x" and "y" at IP = 0,1,2... .
# All Gaussians of one kicker should have a common mean. The weights might be
# given not normalized. Note, that contrary to the kicked bunch, the kicker
# sigmas should be specified directly at the IP of the beam-beam interaction,
# without any extrapolation via beta. In other words, such kicker sigmas, in
# general, can differ from the ones at the "kicked.ip" where the kicked sigmas
# are specified.
  sigma_weight = [[[[40, 0.2], [40.0, 0.8]],  # IP = 0, x
                   [[40, 0.3], [40., 0.6], [40.0, 0.1]], # IP = 1, x
                   [[40.002, 2], [39.999, 10], [40.001, 10], [39.998, 2]],
                   [[80,1]]],
                  [[[40, 0.2]], # IP = 0, y
                   [[40, 0.2]],
                   [[40, 0.2]],
                   [[80.001,1], [79.999,1]]]], # IP = 3, y
# Kicker "x" and "y" positions at the given IP in um in the frame where the
# bunch center of the kicked bunch before beam-beam was at zero. Each vector
# can have either one or "n_step" kicker positions for a given IP, where
# "n_step" is the number of scan points in van der Meer scan. If the vector
# contains only one position, it is "recycled" and the kicker is assumed to be
# immovable during the scan along the given coordinate.
  position = [[pos_x, pos_x, pos_x, [x*2 for x in pos_x]],
              [[0], [0], [0], [0]]]) # kickers at all IPs are immovable in Y

# -------------------- Simulation parameters --------------------
# All parameters below are set to their defaults
#
sim = Sim(
# Number of "macro-particles" traced in the simulation will be slightly less
# than this number (eg. by 21% for round single Gaussian kicked bunch). See
# "Initial distribution of macro-particles" below.
  n_points = 5000,
# Vector with 4 integers: number of turns in the "no beam-beam", "adiabatic
# switch on", "stabilization" and "beam-beam" simulation phases.
  n_turns = [1000, 1000, 0, 5000],
# One of "precise", "average" or "precise.minus.average". If "precise"
# (default), the kick is calculated according to the exact formula. "average"
# sets the kick to the constant value everywhere in the X,Y-plane. This is
# equivalent to the field of a dipole magnet. Its strength is chosen to
# reproduce the precise overall force on the kicked bunch which depends on the
# distance between the colliding bunches. "precise.minus.average" sets the
# kick to the difference between the "precise" and "average" values.
  kick_model = "precise",
# Controls the limits of 4D phase space volume from which the macro-particles
# are selected. See "Initial distribution of macro-particles" below.
  n_sigma_cut = 5,
# Two integers specifying the number of cells in the interpolation grids. The
# first is for a one-dimensional interpolation of the kicker bunch densities
# in X and Y, the second is for each side of the two-dimensional X-Y
# interpolation grid of the kicker fields. If one of these numbers is set to
# zero, the corresponding interpolation is not performed and the value each
# time is calculated exactly (CPU consuming). See "Interpolation" below.
  density_and_field_interpolators_n_cells_along_grid_side = [500, 500],
# If the following parameter is larger than zero, the quality of the
# interpolations (if any) will be estimated by comparing the mismatches
# between the exact and the interpolated values at
# "n_random_points_to_check_interpolation" points randomly distributed inside
# the interpolation grid. The results will be reported in a file
# "output_dir"/interpolation_precision.txt.
  n_random_points_to_check_interpolation = 10000,
# If "output" option "points" is specified, all macro-particle positions will
# be stored to disk. Storage per accelerator turn might require too much
# space. Instead, one can "select_one_turn_out_of" N turns. Eg. if this
# parameter is set to 1000, only the positions before the turns = 999, 1999,
# 2999, ...  (counting from zero) will be stored. Zero or negative values of
# "select_one_turn_out_of" are substituted by one. If "points" option is not
# requested in "output", "select_one_turn_out_of" has no effect.
  select_one_turn_out_of = 1000,
# Random seed for the simulation. If specified, the results of the simulation
# will be fully reproducible. Otherwise, the current time will be used as a
# seed leading to not reproducible results.
  seed = 123456789,
# The name of the subdirectory (relative to the current path or absolute)
# where all simulated files will be stored. This subdirectory is created by
# the program and should not exist, otherwise the program will terminate
# without overwriting anything.  '
  output_dir = "tmp",
# A character string with white-space separated options controlling what will
# be calculated during the simulation. If "output_dir" is not empty (""), for
# every option "XXXX.XXXX", one compressed file will be created under the name
# "output_dir/XXXX_XXXX.txt.gz". "output" might include any combination of
# "integrals_per_turn", "avr_xy_per_turn", "integrals_per_particle",
# "avr_xy_per_particle", "points" or be empty, see "Output" below.
  output = "integrals.per.turn avr.xy.per.turn integrals.per.particle avr.xy.per.particle points")

# -------------------- Running the simulation --------------------
summary = beam_beam(kicked, kickers, sim,
                    False) # False/True control whether the input and the
                           # execution progress will be printed to cout

# -------------------- Print corrections --------------------
for step in range(len(summary[0])):
  print "step =", step
  for ip in range(len(summary)):
    print " ", summary[ip][step].correction


# -------------------- Description of the simulation --------------------
#
# ....... Initial distribution of macro-particles:
# The traced macro-particles of the "kicked" bunch are selected in the
# following way.
#
# It is assumed that the initial bunch densities unperturbed by beam-beam
# factorize in X and Y and their X,Y-projections are multi-Gaussian
# distributions with common centers. Because of the circular motion in X-X'
# and Y-Y' planes (where X',Y' denote the angular coordinates scaled by the
# accelerator beta-function, eg. X' = dX/dZ * beta), the projections to X'
# (Y') are identical to X (Y). So, eg. a single Gaussian in X makes a
# two-dimensional Gaussian with equal sigmas in X-X' plane, and a
# multi-Gaussian makes a multi two-dimensional Gaussian with the same weights.
#
# The four-dimensional kicked bunch density is sampled internally in a
# two-dimensional rX-rY grid of the radii in the X-X', Y-Y' planes at the
# interaction point IP = 0. Each grid side (rX or rY) extends from 0 to a
# maximal radius rX,Y_max. The latter is chosen such that the circle with the
# radius rX_max (or rY_max) contains the same fraction of a multi
# two-dimensional Gaussian as that of a single two-dimensional Gaussian inside
# the circle with the radius "n_sigma_cut" (5 by default). The grid sides
# [0...rX_max] and [0...rY_max] are divided into int(sqrt("n_points"))
# intervals and the grid lines are drawn through their centers. Finally, only
# rX-rY points inside the ellipse inscribed to the rectangle [0...rX_max] X
# [0...rY_max] are used for the simulation. The number of selected points is,
# therefore, less than "n_points". For every rX-rY pair selected in this way,
# one four-dimensional point at the X-X', Y-Y' circles with these radii and
# random phases is chosen as the initial position of one macro-particle. Its
# weight is chosen such that the ensemble of all macro-particles represents
# the initial bunch density. The macro-particles are then traced individually
# in the accelerator in separate threads. This allows to achieve significant
# speed-up of the simulation in the computers with multi-core processors, the
# gain is almost proportional to the number of cores.
# 
# ....... Four phases in the simulation:
# First phase is without beam-beam. It allows to calculate numerically the
# undisturbed overlap integral, compare it with the exact analytic formula and
# estimate the bias of the numerical integration. The final correction is then
# calculated as a ratio of the numerical integrals with (last phase) and
# without (first phase) the beam-beam interaction. The bias is cancelled in
# the ratio at least partially and this potentially allows to improve the
# precision.
#
# During the second phase the beam-beam interaction is switched on
# "adiabatically": linearly with the turn number from zero to its nominal
# value. If this phase is omitted (zero turns) the switch is abrupt: from zero
# to nominal.
#
# After beam-beam is fully switched on, there might be a need to wait for the
# stabilization. This is the third phase. With the adiabatic switch on
# (eg. during 1000 turns) the stabilization phase can normally be omitted and
# the corresponding number of turns can be set to zero.
#
# The final, fourth phase is with beam-beam. Only this phase is used to
# calculate the perturbed overlap integral and the luminosity correction.
#
# "n_turns" parameter specifies the number of turns in all 4 phases.
# 
# ....... Interpolation:
# To speed up the simulation, the X- and Y-densities of the kicker bunch at a
# given point are determined using linear (one-dimensional) interpolations,
# while the generated E-field - using a bilinear (two-dimensional)
# interpolation. For that, the precise density and the field are precalculated
# at the grid of points and interpolated in between.
#
# The interpolation grid is chosen to minimally cover the ranges +/- 1.1 *
# rX_max, +/- 1.1 * rY_max around the kicked bunch center but viewed from the
# kicker bunch center for all set kicker separations. Here, rX,Y_max are the
# maximal simulated radii in X-X' and Y-Y' planes.  Therefore, the grid is
# sufficiently wide to cover all simulated without beam-beam particle
# trajectories (circles), while with 10\% margins (with the factor 1.1) it
# likely covers also the trajectories deformed by the beam-beam
# interaction. If not, ie. if some point goes beyond the interpolation grid,
# the program falls back to the calculation of the exact field instead of the
# interpolation.
#
# The linear (bilinear) grid is a raw (a matrix) of N (NxN) cells and has the
# total number of points N+1 ((N+1)*(N+1)). N is defined separately for the
# one-dimensional (X,Y)-density and for the two-dimensional field
# interpolators by the parameter
# "density_and_field_interpolators_n_cells_along_grid_side". If the
# interpolation is not desired, the corresponding N should be set to zero,
# like "0 500". The first number is for the density, the second - for the
# field, so in this case the density will be calculated using exact formulas
# while the field will be interpolated with N=500.
#
# If the kicker densities and the fields are interpolated, the simulation of
# complex multi-Gaussian elliptical bunches and simple round single-Gaussian
# ones takes about the same time.
#
# The interpolators are created not for each kicker bunch but for each group
# of identical bunches. Ie. if there are identical kickers, the interpolators
# are reused to save memory.
#
# The accuracy of the interpolation is estimated if the parameter
# "n_random_points_to_check_interpolation" is larger than zero. This is
# performed by measuring the maximal and average absolute mismatches between
# the interpolated and the exact values at
# "n_random_points_to_check_interpolation" distributed inside the
# interpolation grid. The mismatches are then normalized to the maximal
# absolute exact value and printed to the file
# "output_dir"/interpolation_precision.txt.
#
# If the interpolation was switched off by setting one or both
# "density_and_field_interpolators_n_cells_along_grid_side" values to zero,
# the corresponding "n_random_points_to_check_interpolation" has no
# effect. The check is not performed either if
# "n_random_points_to_check_interpolation" is zero or negative.
#
# ....... Kick model:
# For debugging purposes one can redefine the kick formula by setting
# "kick.model" parameter to one of the following:
#      precise
#      average
#      precise.minus.average
# "precise" is the default. In this case the exact kick formula is used.
#
# "average" means constant X,Y-independent kick equal to the value averaged
# over the kicked bunch. Ie. it is equal to the sum of the kicks of all bunch
# particles divided by their number. For the Gaussian bunches it can be
# calculated as the kick exerted by the kicker bunch with (Capital) sigma =
# sqrt(sigma1^2 + sigma2^2) on one particle placed at the kicked bunch
# center. The constant kick only shifts the kicked bunch as a whole, but does
# not modify its shape, so the resulting luminosity change can be computed
# using analytic formula.
#
# "precise.minus.average" kick is simply calculated as the difference
# "precise" - "average". This model allows to compare the "precise" beam-beam
# luminosity correction with the sum of the corrections obtained with
# "precise.minus.average" and "average".
#
# ....... Output:
# "output" parameter controls what should be calculated and printed after the
# simulation. It is a white-space separated list of options. If "output_dir"
# is not empty (""), for every option "XXXX" one compressed file will be
# created under the name "output_dir"/XXXX.txt.gz. The possible options are
# listed below:
#
# "integrals_per_turn" - overlap integrals averaged over particles, per
#   turn.
# Format: step ip ip_kicker_x ip_kicker_y i_turn integral.
#
# "avr_xy_per_turn" - kicked bunch center averaged over particles, per
#   turn.
# Format: step ip i_turn average_x average_y.
#
# "integrals_per_particle" - overlap integrals averaged over turns, per
#   particle, per phase, assuming that the kicked bunch is squashed to this
#   particle, ie. the bunch density is a delta-function at its position.
# Format: step ip phase i_particle integral.
#
# "avr_xy_per_particle" - kicked bunch center averaged over turns, per
#   particle, per phase as above.
# Format: step ip phase i_particle average_x average_y.
#
# "points" -  the traced positions of the kicked bunch particles before
#   accelerator "i_turn", note, the generated file might be very long.
# Format: step ip i_turn i_particle X X' Y Y'.
#
# The integers i_turn and i_particle are counted from zero.
#
# In addition, the files
# "output_dir"/rx_ry_weights.txt.gz,
# "output_dir"/kicker_positions.txt,
# "output_dir"/interpolation_precision.txt and
# "output_dir"/summary.txt are printed out regardless of options in
# "output" with the following content:
#
# "output_dir"/rx_ry_weights.txt.gz - for every simulated particle: its
#             initial rX and rY radii (in X,X' and Y,Y' phase spaces) at ip
#             and its "weight".
# Format: i_particle ip rx ry weight.
#
# "output_dir"/kicker_positions.txt - coordinates of the kicker bunch
#     centers with respect to the nominal (not perturbed by the
#     beam-beam effect) kicked bunch center.
# Format: step ip x y.
#
# "output_dir"/interpolation_precision.txt - was explained in
# "Interpolation".
# 
#  "output_dir"/summary.txt - the main results of the simulation
#     including the overlap integrals and the corresponding beam-beam
#     corrections. Format:
# <step> <IP> <beam-beam/no beam-beam luminosity correction>
# no beam-beam: <analytic> overlap, <numeric/analytic ratio> and <its error>
# no beam-beam: <X>, <Y> numeric center-of-mass shift
#    beam-beam: <X>, <Y> analytic, <X>, <Y> numeric shift.
# 
# No beam-beam numeric/analytic error is roughly estimated from turn-by-turn
# variations, available only if "integrals_per_turn" option is set in
# "output". Otherwise this error is set to "nan". Similarly,
# numeric <X>, <Y> shifts are calculated only if "avr_xy_per_particle"
# option is chosen, and set to "nan" otherwise. Without beam-beam these
# shifts should be close to zero.
#
# The same summary data are returned by "beam_beam" function itself: after
#
#   summary = beam_beam(kicked, kickers, sim, quiet)
#
# summary[ip][step] contains "Summary" structure (namedtuple) defined in
# "BxB.py" with the results. If only the minimal information on the luminosity
# correction is required from the simulation, it is better to set "output" to
# an empty string ("") or to "integrals_per_particle", as any other option in
# "output" requires extra CPU time. The integrals per particle are calculated
# in any case (but the file "output_dir"/integrals_per_particle.txt.gz is
# printed only if the corresponding option is explicitly set in "output").
# 

