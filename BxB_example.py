from BxB import Kicked, Kickers, Sim, Summary, beam_beam

# ---------- Usage example ----------
kicked = Kicked(momentum = 3500,
                Z = 1,
                ip = 1,
                beta = [[1.5, 1.5, 1.5, 6],
                        [1.5, 1.5, 1.5, 6]],
                next_phase_over_2pi = [[0.31*(i+1) for i in range(4)],
                                       [0.32*(i+1) for i in range(4)]],
                sigma_weight = [[[40, 0.2], [40.0, 0.8]],
                                [[39.99, 0.3], [40, 0.7]]],
                exact_phases = False)
pos_x = [10*i for i in range(21)]
kickers = Kickers(Z = 1,
                  n_particles = [8.5e10, 8.5e10, 8.5e10, 8.5e10],
                  sigma_weight = [[[[40, 0.2], [40.0, 0.8]],
                                   [[40, 0.3], [40., 0.6], [40.0, 0.1]],
                                   [[40.002, 2], [39.999, 10], [40.001, 10], [39.998, 2]],
                                   [[80,1]]],
                                  [[[40, 0.2]],
                                   [[40, 0.2]],
                                   [[40, 0.2]],
                                   [[80.001,1], [79.999,1]]]],
                  position = [[pos_x, pos_x, pos_x, [x*2 for x in pos_x]],
                              [[0], [0], [0], [0]]])
sim = Sim(n_points = 5000,
          n_turns = [1000, 1000, 0, 5000],
          kick_model = "precise",
          n_sigma_cut = 5,
          density_and_field_interpolators_n_cells_along_grid_side = [500, 500],
          n_random_points_to_check_interpolation = 10000,
          select_one_turn_out_of = 1000,
          seed = 123456789,
          output_dir = "pyth",
          output = "integrals.per.turn avr.xy.per.turn integrals.per.particle avr.xy.per.particle points")

summary = beam_beam(kicked, kickers, sim, False)

for step in range(len(summary[0])):
  print "step =", step
  for ip in range(len(summary)):
    print " ", summary[ip][step].correction

