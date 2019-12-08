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
# Python interface of B*B
# Usage: define Kicked, Kickers and Sim structures and call
#
#   summary = beam_beam(kicked, kickers, sim, quiet)
#
# Then, summary[ip][step] will contain the results of the simulation.
#
#
from ctypes import *
from collections import namedtuple

# ---------- Input parameters ----------
Kicked = namedtuple("Kicked",
                    ["momentum", "Z", "ip", "beta",
                     "next_phase_over_2pi",
                     "sigma_weight", "exact_phases"])
Kickers = namedtuple("Kickers",
                     ["Z", "n_particles",
                      "sigma_weight", "position"])
Sim = namedtuple("Sim",
                 ["n_points", "n_turns", "kick_model",
                  "n_sigma_cut",
                  "density_and_field_interpolators_n_cells_along_grid_side",
                  "n_random_points_to_check_interpolation",
                  "select_one_turn_out_of", "seed", "output_dir", "output"])
# ---------- Summary format ----------
Summary = namedtuple("Summary",
                     ["correction",
                      "no_bb_analytic_integ",
                      "no_bb_numeric_over_analytic_integ",
                      "no_bb_numeric_over_analytic_integ_err",
                      "avr_analytic",
                      "avr_numeric",
                      "no_bb_avr_numeric"])

# ---------- C structures ----------
class _C_Multi_Gaussian(Structure) :
  _fields_ = [("n", c_int),
              ("sigma", POINTER(c_double)),
              ("weight", POINTER(c_double))]

class _C_Kicked(Structure) :
 _fields_ = [("momentum", c_double),
             ("Z", c_int),
             ("ip", c_int),
             ("beta", POINTER(c_double)),
             ("next_phase_over_2pi", POINTER(c_double)),
             ("gaussian", _C_Multi_Gaussian*2),
             ("exact_phases", c_bool)]

class _C_Kickers(Structure) :
 _fields_ = [("Z", c_int),
             ("n_particles", POINTER(c_double)),
             ("gaussian", POINTER(_C_Multi_Gaussian)),
             ("position", POINTER(c_double))]

class _C_Sim(Structure) :
 _fields_ = [("n_ip", c_int),
             ("n_step", c_int),
             ("n_points", c_int),
             ("n_turns", c_int*4), # for 4 PHASES
             ("kick_model", c_char_p),
             ("n_sigma_cut", c_int),
             ("density_and_field_interpolators_n_cells_along_grid_side", c_int*2),
             ("n_random_points_to_check_interpolation", c_int),
             ("select_one_turn_out_of", c_int),
             ("seed", c_long),
             ("output_dir", c_char_p),
             ("output", c_char_p)]

class _C_Summary(Structure) :
 _fields_ = [("correction", c_double),
             ("no_bb_analytic_integ", c_double),
             ("no_bb_numeric_over_analytic_integ", c_double),
             ("no_bb_numeric_over_analytic_integ_err", c_double),
             ("avr_analytic", c_double*2),
             ("avr_numeric", c_double*2),
             ("no_bb_avr_numeric", c_double*2)]


# -------------------- Helper functions to fill C structures --------------------
def _to_ints(x): return (c_int   *len(x))(*x)

def _to_doubles(x): return (c_double*len(x))(*x)

def _gaussian_c(sig_w) :
 return _C_Multi_Gaussian(n = len(sig_w),
                          sigma  = (c_double*len(sig_w))(*[c_double(sw[0]) for sw in sig_w]),
                          weight = (c_double*len(sig_w))(*[c_double(sw[1]) for sw in sig_w]))

def _gaussians_c(g):  return (_C_Multi_Gaussian*len(g))(*[_gaussian_c(x) for x in g])

# ---------- Helper function to convert C results to Python ----------
def _python_summary(s) :
 return Summary(s.correction,
                s.no_bb_analytic_integ,
                s.no_bb_numeric_over_analytic_integ,
                s.no_bb_numeric_over_analytic_integ_err,
                s.avr_analytic,
                s.avr_numeric,
                s.no_bb_avr_numeric)

# -------------------- C functions --------------------
_bb_c = CDLL("./BxB_python.so")
_bb_c.beam_beam.restype = None
_bb_c.beam_beam.argtypes = [POINTER(_C_Kicked), POINTER(_C_Kickers), POINTER(_C_Sim), POINTER(_C_Summary), c_bool]

# -------------------- Main function --------------------
def beam_beam(kicked, kickers, sim, quiet = False):
  n_ip = len(kicked.beta[0])
  # recycle stable kicker positions specified as one number
  n_step = max([len(pos_ip)
                for pos_coor in kickers.position
                for pos_ip in pos_coor])
  pos = [0] * (2 * n_ip * n_step)
  for coor in range(2):
    for ip in range(n_ip):
      size = len(kickers.position[coor][ip])
      if (size == 1):
        for step in range(n_step):
          pos[coor*n_ip*n_step + ip*n_step + step] = kickers.position[coor][ip][0]
      elif (size == n_step):
        for step in range(n_step):
          pos[coor*n_ip*n_step + ip*n_step + step] = kickers.position[coor][ip][step]
      else:
        raise Exception("The number of kicker positions should be either one or maximal\n")
  kicked_c = _C_Kicked(momentum = c_double(kicked.momentum),
                       Z = c_int(kicked.Z),
                       ip = c_int(kicked.ip),
                       beta = _to_doubles(kicked.beta[0] + kicked.beta[1]),
                       next_phase_over_2pi = _to_doubles(kicked.next_phase_over_2pi[0] +
                                                         kicked.next_phase_over_2pi[1]),
                       gaussian = _gaussians_c(kicked.sigma_weight),
                       exact_phases = c_bool(kicked.exact_phases))
  kickers_c = _C_Kickers(Z = kickers.Z,
                         n_particles = _to_doubles(kickers.n_particles),
                         gaussian = _gaussians_c(kickers.sigma_weight[0] +
                                                 kickers.sigma_weight[1]),
                         position = _to_doubles(pos))
  sim_c = _C_Sim(n_ip = n_ip,
                 n_step = n_step,
                 n_points = sim.n_points,
                 n_turns = _to_ints(sim.n_turns),
                 kick_model = c_char_p(sim.kick_model),
                 n_sigma_cut = sim.n_sigma_cut,
                 density_and_field_interpolators_n_cells_along_grid_side =
                 _to_ints(sim.density_and_field_interpolators_n_cells_along_grid_side),
                 n_random_points_to_check_interpolation = sim.n_random_points_to_check_interpolation,
                 select_one_turn_out_of = sim.select_one_turn_out_of,
                 seed = sim.seed,
                 output_dir = c_char_p(sim.output_dir),
                 output = c_char_p(sim.output))
  summary = (_C_Summary*(n_ip*n_step))()
  _bb_c.beam_beam(byref(kicked_c), byref(kickers_c), byref(sim_c), summary, c_bool(quiet))
  return [[_python_summary(summary[ip * n_step + step]) for step in range(n_step)] for ip in range(n_ip)]

