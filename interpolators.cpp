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
#include "interpolators.hh"
#include <cmath> // for floor()
#include <algorithm> // for max
#include <iostream>

void Linear_interpolator::create(double x_min, double x_max,
				 size_t n_cells,
				 const function<double (double)>& function_to_be_interpolated) {
  xmin = x_min;
  step_x = (x_max-x_min) / n_cells;
  fun = function_to_be_interpolated;
  N = n_cells;
  m.resize(N+1);
  vector<double>::iterator m_iter = m.begin();
  for (size_t ix=0; ix<=N; ++ix) {
    double x = xmin + ix * step_x;
    *(m_iter++) = fun(x);
  }
}

double Linear_interpolator::operator()(double x) const {
  double translated_x = x - xmin;
  int ix = int(floor(translated_x / step_x));
  if (0 <= ix && ix < N) {
    double dx_down = translated_x - ix * step_x;
    double dx_up   = step_x - dx_down;
    vector<double>::const_iterator m_ix = m.begin() + ix;
    return
      ((* m_ix     ) * dx_up +
       (*(m_ix + 1)) * dx_down) / step_x;
  } else {
    //    cout << "calling slow 1D function" << endl;
    return fun(x); // fall back to calling a slow function
  }
}

pair<double, double> Linear_interpolator::
max_and_average_mismatches_relative_to_max(int n_random_points) const {
  double max_interpolation_discrepancy = 0, avr_interpolation_discrepancy = 0;
  double max_field = 0;
  {
    double dx = step_x * N;
    for (int i=0; i<n_random_points; ++i) {
      double x = xmin + drand48() * dx;
      double
	precise = fun    (x),
	approx  = (*this)(x);
      double mismatch = fabs(approx - precise);
      max_field       = max(max_field, fabs(precise));
      max_interpolation_discrepancy = max(max_interpolation_discrepancy, mismatch);
      avr_interpolation_discrepancy += mismatch;
    }
    avr_interpolation_discrepancy /= n_random_points;
  }
  return make_pair(max_interpolation_discrepancy / max_field,
		   avr_interpolation_discrepancy / max_field);
}  
// ------------------------------------------------------------
void Bilinear_interpolator::create(double x_min, double x_max, double y_min, double y_max,
				   size_t n_cells_along_grid_side,
				   const function<complex<double>(double, double)>& function_to_be_interpolated) {
  xmin = x_min; ymin = y_min; 
  step_x = (x_max-x_min) / n_cells_along_grid_side;
  step_y = (y_max-y_min) / n_cells_along_grid_side;
  fun = function_to_be_interpolated;
  step_xy = step_x * step_y;
  N = n_cells_along_grid_side;
  m.resize((N+1) * (N+1));
  vector<complex<double>>::iterator m_iter = m.begin();
  for   (size_t iy=0; iy<=N; ++iy) {
    for (size_t ix=0; ix<=N; ++ix) {
      double x = xmin + ix * step_x;
      double y = ymin + iy * step_y;
      *(m_iter++) = fun(x, y);
    }
  }
}

complex<double> Bilinear_interpolator::operator()(double x, double y) const {
  double translated_x = x - xmin;
  double translated_y = y - ymin;
  int ix = int(floor(translated_x / step_x));
  int iy = int(floor(translated_y / step_y));
  if (0 <= ix && ix < N &&
      0 <= iy && iy < N) {
    double dx_down = translated_x - ix * step_x;
    double dy_down = translated_y - iy * step_y;
    double dx_up   = step_x - dx_down;
    double dy_up   = step_y - dy_down;
    vector<complex<double>>::const_iterator m_ix_iy = m.begin() + ix + iy * (N+1);
    return
      ((* m_ix_iy       ) * dx_up   * dy_up +
       (*(m_ix_iy +   1)) * dx_down * dy_up +
       (*(m_ix_iy + N+1)) * dx_up   * dy_down +
       (*(m_ix_iy + N+2)) * dx_down * dy_down) / step_xy;
  } else {
    //    cout << "calling slow 2D function" << endl;
    return fun(x, y); // fall back to calling a slow function
  }
}

pair<double, double> Bilinear_interpolator::
max_and_average_mismatches_relative_to_max(int n_random_points) const {
  double max_interpolation_discrepancy = 0, avr_interpolation_discrepancy = 0;
  double max_field = 0;
  {
    double
      dx = step_x * N,
      dy = step_y * N;
    for (int i=0; i<n_random_points; ++i) {
      double x = xmin + drand48() * dx;
      double y = ymin + drand48() * dy;
      complex<double>
	precise = fun    (x, y),
	approx  = (*this)(x, y);
      double mismatch = abs(approx - precise);
      max_field       = max(max_field, abs(precise));
      max_interpolation_discrepancy = max(max_interpolation_discrepancy, mismatch);
      avr_interpolation_discrepancy += mismatch;
    }
    avr_interpolation_discrepancy /= n_random_points;
  }
  return make_pair(max_interpolation_discrepancy / max_field,
			avr_interpolation_discrepancy / max_field);
}  

