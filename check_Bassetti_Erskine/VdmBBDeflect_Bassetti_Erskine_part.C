// Taken from VdmBBDEflect.C,h with cosmetic changes to drop dependence on full class

#include <complex>
#include <Faddeeva.hh>

template <typename T> T Sign(T val) {
    return (T) ((int) (T(0) < val) - (val < T(0)));
}

std::complex<double> wfun(std::complex<double> z)
{
  std::complex<double> w = Faddeeva::w(z);
    
  double     wx = w.real();
  double     wy = w.imag();
    
  return std::complex<double>(wx, wy);    
}
// sigmax,y are added as arguments (change by VB)
std::pair<std::pair<double,double>, std::pair<double,double> > BaseErsk(double sepx, double sepy,
									double sigmax, double sigmay)
{
  //The BassErsk routine compute the x and y electric field vectors
  //produced from one bunch by the opposite one with Np particles
  //located at a distance (x,y)  
    
  // Convert all dimensions to meters and take absolute value of deflections
  //
  // removed conversion - division by 1000 (VB)
  double x = std::fabs(sepx);  // std::fabs(sepx/1000.);
  double y = std::fabs(sepy);  // std::fabs(sepy/1000.);

  // commented out by VB:
  // double sigmax = _Sigx/1000.;
  // double sigmay = _Sigy/1000.;
  
  //I changed the constant factor in front of the fields to have it in
  //units of rp and gamma (K= 2.*rp/gamma) and not eps0 which I have
  //set to 1 instead of eps0 = 8.854187817620e-12 
  
  static const double eps0 = 1;
  // --- end of changes (VB) ---
  
  double x1 = x;
  double x2 = y;

  double sigma1 = sigmax;
  double sigma2 = sigmay;

  if (sigmax < sigmay) {
    x1 = y;
    x2 = x;
    sigma1 = sigmay;
    sigma2 = sigmax;
  }

  double S = std::sqrt(2*(sigma1*sigma1-sigma2*sigma2));
  double factBE = std::sqrt(M_PI)*2./(2*eps0*S);
  std::complex<double> etaBE = std::complex<double>((sigma2/sigma1)*x1 ,(sigma1/sigma2)*x2);
  std::complex<double> zetaBE(x1,x2);
  
  std::complex<double> w1 = wfun(zetaBE/S);
  std::complex<double> w2 = wfun(etaBE/S);
  
  std::complex<double> val = factBE*(w1-std::exp(-x1*x1/(2.*sigma1*sigma1)-x2*x2/(2.*sigma2*sigma2))*w2);
  
  double Ex = 0, Ey = 0;

  if (sigmax > sigmay) {
    Ex = std::fabs(val.imag())*Sign(sepx);
    Ey = std::fabs(val.real())*Sign(sepy);
  }
  else {
    Ey = std::fabs(val.imag())*Sign(sepy);
    Ex = std::fabs(val.real())*Sign(sepx);
  }
   
  std::pair<double,double> pair1(Ex, Ey);
  std::pair<double,double> pair2(std::exp(-x*x/(2.*sigmax*sigmax)), val.imag());
  return std::pair<std::pair<double,double>, std::pair<double, double> >(pair1, pair2);
}

// sigmax,y are added as arguments (change by VB)
std::pair<double, double> GetExEy(double sepx, double sepy,
				  double sigmax, double sigmay)
{
  std::pair<std::pair<double,double>, std::pair<double,double> > bE = BaseErsk(sepx, sepy, sigmax, sigmay);

  return std::pair<double, double>(bE.first.first, bE.first.second);
}

