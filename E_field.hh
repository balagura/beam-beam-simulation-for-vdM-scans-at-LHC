#ifndef E_field_hh 
#define E_field_hh 1

// Bassetti-Erskine beam-beam kick calculation for elliptical and round
// bunches.
//
// Originally based on VdmBBDeflect class, then adopted to B*B by Timothy
// Barklow (SLAC). Finally, rewritten from scratch by Vladislav Balagura but
// following the logic from VdmBBDeflect to ensure positiveness of the
// imaginary part of Faddeeva::w() argument (to ensure convergence) and using
// formulas from the Bassetti-Erskine publication.

// The program "check_Bassetti_Erskine/compare_with_VdmBBDeflect.C" cross
// checks that the results returned by the function below are consistent with
// "VdmBBDeflect" for N points (eg. 100000) randomly distributed in the ranges
// specified in "check_Bassetti_Erskine/compare_with_VdmBBDeflect_config.txt".
//
// Compile it with "cd check_Bassetti_Erskine; make compare_with_VdmBBDeflect"
// and run as "compare_with_VdmBBDeflect compare_with_VdmBBDeflect_config.txt"
//
// Author V. Balagura, balagura@cern.ch (Nov 2019)
#include <complex>
#include "Faddeeva.hh"

using namespace std;

inline complex<double> E_from_unit_charge_2d_gaussian_times_2pi_epsilon_0(double x, double y,
									  double sigma_x, double sigma_y) {
  // Calculates
  //
  //  2pi*epsilon_0 * (Ex(x,y) + i Ey(x,y))
  //
  // where Ex(x,y), Ey(x,y) are components of the two-dimensional
  // electrostatic field at the point (x,y), created by the source bunch of
  // unit charge with the two-dimensional Gaussian shape with sigma_x, sigma_y
  // widths along x- and y-axis, respectively, and with the center at (0,0)
  //
  // In the limit sigma_x = sigma_y = sigma the length of the returned vector
  // approaches the expression
  //
  //   1/R * (1 - exp(-R^2/2/sigma^2))
  //
  // corresponding to the round bunch, while its direction approaches the
  // direction of the vector from (0,0) to (x,y)
  //
  if (abs(x / sigma_x) < 1e-8 && abs(y / sigma_y) < 1e-8) return complex<double>(0, 0);
  double r =  sigma_y / sigma_x;
  complex<double> res;
  if ( abs(r-1) < 1e-8 ) {
    // round beams: cancel two singularities
    double mean_sigma = (sigma_x + sigma_y) / 2;
    complex<double> z(x, y);
    double r2 = norm(z);
    return (1 - exp(-0.5 * (x*x + y*y) / mean_sigma / mean_sigma)) / r2 * z;
  }
  // elliptical beams:
  //
  // ensure positiveness of sqrt() argument in sqrt(2*(sigma_x^2-sigma_y^2))
  // below such that the complex arguments z1,z2 of Faddeeva function have
  // Im(z1), Im(z2) > 0. Otherwise, Faddeeva function grows exponentially and
  // the difference between the two functions becomes numerically
  // unstable. For simplicity, consider abs(x,y) instead of x,y, and in case
  // x<0 or y<0: invert the sign of the corresponding field components
  //
  bool xy_swapped = sigma_x < sigma_y;
  bool x_positive = x > 0, y_positive = y > 0;
  x = abs(x);
  y = abs(y);
  double sig_x_sq = sigma_x * sigma_x, sig_y_sq = sigma_y * sigma_y;
  //
  double sxy = sqrt( 2 * abs(sig_x_sq - sig_y_sq) );
  //
  double g   = x*x / sig_x_sq + y*y / sig_y_sq;
  // note, in the original B.-E. paper the sign of y-term was inverted
  // (-y*y/sig_y_sq) by mistake
  //
  // define nested C++ function using C++11 lambda-function
  //
  // relerr argument in Faddeeva::w controls the relative error of the
  // result and trades speed / accuracy, set to 0 by default
  auto f = [](complex<double> z1, complex<double> z2, double g) {
	     return Faddeeva::w( z1 ) - exp(-0.5 * g) * Faddeeva::w( z2 );
	     // example of relerr = 1e-5 use in Faddeeva function
	     // return Faddeeva::w( z1, 1e-5 ) - exp(-0.5 * g) * Faddeeva::w( z2, 1e-5 );
	   };
  complex<double> z1(x,   y);   z1 /= sxy;
  complex<double> z2(x*r, y/r); z2 /= sxy;
  //
  // Basetti-Erskine formula, Eq.6 in their original paper CERN-ISR-TH-80-06
  // for Q=1:
  //
  //   2pi*epsilon_0 * (Ex - iEy) = -i sqrt(pi) / sxy * f(z1, z2, g)
  //
  if (!xy_swapped) {
    res = f(z1, z2, g);
  } else {
    // If sigma_x < sigma_y: swap x and y, calculate the result and swap back
    // its real/imaginary parts (no need to swap "g" as it is symmetric)
    //
    // another nested lambda (to avoid pollution of the global namespace)
    auto swap = [](complex<double> z) { return complex<double>(imag(z), real(z)); };
    res = swap( f(swap(z1), swap(z2), g) );
  }
  res *= sqrt(M_PI) / sxy;
  // After that: res = sqrt(pi) / sxy * f(z1, z2, g) = Ey + iEx,
  // invert Ex(y) component if x(y) was negative:
  return complex<double>(x_positive ? imag(res) : -imag(res),
			 y_positive ? real(res) : -real(res));
}

#endif
