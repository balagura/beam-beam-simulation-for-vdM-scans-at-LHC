#ifndef bilinear_interpolator_hh 
#define bilinear_interpolator_hh 1

// Simple class for the bilinear interpolation, see
// eg. https://en.wikipedia.org/wiki/Bilinear_interpolation.
//
// Interpolation is performed using precalculated values at "x_min, x_max,
// y_min, y_max" grid with "n_cells_along_grid_side", for the std::function
// "field_fun" (which can be a normal function, a lambda-function etc.,
// please, see the C++ documentation on std::function). The total number of
// grid points is (n_cells_along_grid_side + 1)^2.
//
// The type of the two-dimensional field is a template parameter T.
//
// When constructed, the object "bi" of the type Bilinear_interpolator<T> can
// calculate and return the interpolated field of type T via:
//
//  bi(x, y);
//
// If the point (x,y) is outside the interpolation rectangle, the function
// field_fun(x, y) is called directly instead of the interpolation.

// Author Vladislav Balagura, balagura@cern.ch (24.10.2019)

#include <cmath> // for floor()
#include <functional> // for std::function

struct Bilinear_interpolator_grid {
  size_t N; // number of cells along one grid size
  double xmin, ymin;
  double step_x, step_y, step_xy;
};

template <class T>
struct Bilinear_interpolator : private Bilinear_interpolator_grid {
  Bilinear_interpolator() {}
  void create(double x_min, double x_max, double y_min, double y_max,
	      size_t n_cells_along_grid_side,
	      std::function<T(double, double)> field_fun);
  Bilinear_interpolator(double x_min, double x_max, double y_min, double y_max,
			size_t n_cells_along_grid_side,
			std::function<T(double, double)> field_fun)
  { create(x_min, x_max, y_min, y_max, n_cells_along_grid_side, field_fun); }
  T operator()(double x, double y) const;
  ~Bilinear_interpolator();
  // function value without interpolation
  T f(double x, double y) const { return fun(x, y); }
  // interpolation grid parameters (from the base class)
  const Bilinear_interpolator_grid& grid() const { return *this; }
  //
private:
  T* m; // (N+1) x (N+1) field matrix with indexing: m[ix + iy*(N+1)];
  std::function<T(double, double)> fun;
};

template <class T>
void Bilinear_interpolator<T>::create(double x_min, double x_max, double y_min, double y_max,
				      size_t n_cells_along_grid_side,
				      std::function<T(double, double)> field_fun) {
  xmin = x_min; ymin = y_min; 
  step_x = (x_max-x_min) / n_cells_along_grid_side;
  step_y = (y_max-y_min) / n_cells_along_grid_side;
  fun = field_fun;
  step_xy = step_x * step_y;
  N = n_cells_along_grid_side;
  m = new T[(n_cells_along_grid_side+1) * (n_cells_along_grid_side+1)];
  T* m_ptr = m;
  for   (size_t iy=0; iy<=n_cells_along_grid_side; ++iy) {
    for (size_t ix=0; ix<=n_cells_along_grid_side; ++ix) {
      double x = xmin + ix * step_x;
      double y = ymin + iy * step_y;
      *(m_ptr++) = fun(x, y);
    }
  }
}
template <class T>
Bilinear_interpolator<T>::~Bilinear_interpolator() { delete m; }

template <class T>
T Bilinear_interpolator<T>::operator()(double x, double y) const {
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
    const T* m_ix_iy = m + ix + iy * (N+1);
    return
      ((* m_ix_iy       ) * dx_up   * dy_up +
       (*(m_ix_iy +   1)) * dx_down * dy_up +
       (*(m_ix_iy + N+1)) * dx_up   * dy_down +
       (*(m_ix_iy + N+2)) * dx_down * dy_down) / step_xy;
  } else {
    //  std::cout << "calling slow function" << endl;
    return fun(x, y); // fall back to calling a slow function
  }
}

#endif
