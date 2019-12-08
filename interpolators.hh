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
// Linear and bilinear interpolators: linear - to speed up the calculation
// of one-dimensional X and Y kicker densities, bilinear - for the
// two-dimensional field.
//
#ifndef interpolators_hh 
#define interpolators_hh 1

#include <utility> // for pair
#include <functional> // for function
#include <vector>
#include <complex>

using namespace std;

struct Linear_interpolator_grid {
  Linear_interpolator_grid() : N(0) {}
  size_t N; // number of grid cells
  double xmin;
  double step_x;
};

struct Linear_interpolator : private Linear_interpolator_grid {
  void create(double x_min, double x_max,
	      size_t n_cells,
	      const function<double (double)>& function_to_be_interpolated);
  // interpolated value:
  double operator()(double x) const;
  // function without interpolation:
  double f(double x) const { return fun(x); }
  // interpolation grid parameters (from the base class)
  const Linear_interpolator_grid& grid() const { return *this; }
  bool empty() const { return N == 0; }
  // next function returns the maximal and average absolute mismatches between
  // the exact and interpolated values at "n_random_points" uniformly
  // distributed in the grid, normalized by the maximal absolute field value
  pair<double, double> max_and_average_mismatches_relative_to_max(int n_random_points) const;
private:
  vector<double> m; // contains (N+1) function values at grid points
  function<double (double)> fun;
};

// Simple class for the bilinear interpolation of the complex field, see
// eg. https://en.wikipedia.org/wiki/Bilinear_interpolation.
//
// Interpolation is performed using precalculated values at "x_min, x_max,
// y_min, y_max" grid with "n_cells_along_grid_side", for the std::function
// "function_to_be_interpolated" (which can be a normal function, a
// lambda-function etc., please, see the C++ documentation on
// std::function). The total number of grid points is (n_cells_along_grid_side
// + 1)^2.
//
// When constructed, the object "bi" of the type Bilinear_interpolator can
// calculate and return the interpolated field via:
//
//  bi(x, y);
//
// If the point (x,y) is outside the interpolation rectangle, the function
// fun(x, y) is called directly instead of the interpolation.

struct Bilinear_interpolator_grid {
  Bilinear_interpolator_grid() : N(0) {}
  size_t N; // number of cells along one grid size
  double xmin, ymin;
  double step_x, step_y,
    step_xy; // = step_x * step_y
};

struct Bilinear_interpolator : private Bilinear_interpolator_grid {
  void create(double x_min, double x_max, double y_min, double y_max,
	      size_t n_cells_along_grid_side,
	      const function<complex<double> (double, double)>& function_to_be_interpolated);
  // interpolated:
  complex<double> operator()(double x, double y) const;
  // function without interpolation:
  complex<double> f(double x, double y) const { return fun(x, y); }
  // interpolation grid parameters (from the base class)
  const Bilinear_interpolator_grid& grid() const { return *this; }
  bool empty() const { return N == 0; }
  // next function returns the maximal and average absolute mismatches (taken
  // with (double complex<double>::abs())) between the exact and interpolated
  // values at "n_random_points" uniformly distributed in the grid, normalized
  // by the maximal absolute field value
  pair<double, double> max_and_average_mismatches_relative_to_max(int n_random_points) const;
private:
  vector<complex<double>> m; // (N+1) x (N+1) field matrix with indexing: m[ix + iy*(N+1)];
  function<complex<double>(double, double)> fun;
};

#endif
