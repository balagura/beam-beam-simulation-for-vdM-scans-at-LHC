#include <VdmBBDeflect.h>

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>

template <typename T> T Sign(T val) {
    return (T) ((int) (T(0) < val) - (val < T(0)));
}

double VdmBBDeflect::_betax = 0;
double VdmBBDeflect::_betay = 0;
double VdmBBDeflect::_tunex = 0;
double VdmBBDeflect::_tuney = 0;
double VdmBBDeflect::_EMagnetic = 0;
double VdmBBDeflect::_ABeam = 0;
double VdmBBDeflect::_gamma = 0;
double VdmBBDeflect::_rpOverGamma = 0;
// const double VdmBBDeflect::rp = 1.53e-18;
const double VdmBBDeflect::rp = 1.5346955e-18;

VdmBBDeflect::DynBetaLookupMap VdmBBDeflect::_dynBetaLookupDeltaX;
VdmBBDeflect::DynBetaLookupMap VdmBBDeflect::_dynBetaLookupDeltaY;

double VdmBBDeflect::_EmittanceRef;
double VdmBBDeflect::_NBeamRef;
double VdmBBDeflect::_Qx;
double VdmBBDeflect::_Qy;
double VdmBBDeflect::_deltaxBetaXRef;
double VdmBBDeflect::_deltaxBetaYRef;
double VdmBBDeflect::_deltayBetaXRef;
double VdmBBDeflect::_deltayBetaYRef;

bool VdmBBDeflect::_haveBeamParams = false;
bool VdmBBDeflect::_haveDynBetaTable = false;
bool VdmBBDeflect::_debug = true;

void VdmBBDeflect::CalculateScaleFactors()
{
  double sqrtNomEmittX = std::sqrt(_gamma/2/_betax)*_Sigx;
  double sqrtNomEmittY = std::sqrt(_gamma/2/_betay)*_Sigy;
  double sqrtEmittRef = std::sqrt(_EmittanceRef);
  
  double Navg = 0.5*(_NpBeam1 + _NpBeam2);
  
  double cotQRatiox = std::tan(2*M_PI*_Qx)/std::tan(2*M_PI*_tunex); // cot(2pi*thistune)/cot(2pi*reftune)
  double cotQRatioy = std::tan(2*M_PI*_Qy)/std::tan(2*M_PI*_tuney); // cot(2pi*thistune)/cot(2pi*reftune)

  _dynBetaScaleX = cotQRatiox*Navg/_NBeamRef*sqrtEmittRef/sqrtNomEmittX*(2*sqrtEmittRef)/(sqrtNomEmittX + sqrtNomEmittY);
  _dynBetaScaleY = cotQRatioy*Navg/_NBeamRef*sqrtEmittRef/sqrtNomEmittY*(2*sqrtEmittRef)/(sqrtNomEmittX + sqrtNomEmittY);
  
  if (_debug) {
    std::cout << " VdmBBDeflect::CalculateScaleFactors: " << std::endl;
    std::cout << "gamma, beta = " << _gamma << ", " << _betax << std::endl;
    std::cout << "Sigx, Sigy Sigx-Sigy= " << _Sigx << ", " << _Sigy << ", " << _Sigx-_Sigy << ", Nbeam1, Nbeam 2 = " << _NpBeam1 << ", " << _NpBeam2 << std::endl;
    std::cout << "sqrt(emmittance) x, y, ref = " << sqrtNomEmittX << ", " << sqrtNomEmittY << ", " << sqrtEmittRef << std::endl;
    std::cout << "cotQRatiox, cotQRatioy = " << cotQRatiox << ", " << cotQRatioy << std::endl;
    std::cout << "Scale factors x, y = " << _dynBetaScaleX << ", " << _dynBetaScaleY << std::endl;
  }
}

void VdmBBDeflect::SetBeamParameters(double betax, double betay, double tunex, double tuney, double EMagneticGeV, double ABeam)
{
  _betax = betax;
  _betay = betay;

  _tunex = tunex - std::floor(tunex);
  _tuney = tuney - std::floor(tuney);

  _EMagnetic = EMagneticGeV;
  _ABeam = ABeam;

  // Proton mass in GeV
  //
  static const double mProt = 0.938272;

  //Proton classical radius in meters
  //
  //  static const double rp = 1.53e-18;
  static const double rp = 1.5346955e-18;

  //relativistic gamma factor
  //
  _gamma = _EMagnetic/mProt;

  _rpOverGamma = rp/_gamma;

  _haveBeamParams = true;
  
  std::cout << " VdmBBDeflect::SetBeamParameters: " << std::endl;
  std::cout << "gamma, betax, betay = " << _gamma << ", " << _betax << ", " << _betay << std::endl;
}

void VdmBBDeflect::ChangeBetaStar(float newBeta)
{
  _betax = newBeta;
  _betay = newBeta;
}

void VdmBBDeflect::ChangeTune(float tunex, float tuney)
{
  _tunex = tunex;
  _tuney = tuney;
}

std::pair<double, double> VdmBBDeflect::GetDeflectAngles(double sepx_mm, double sepy_mm, unsigned int beam)
{
  double nBeamOpp = 0;

  if (beam == 1) {
    nBeamOpp = _NpBeam2;
  } 
  else if (beam == 2) {
    nBeamOpp = _NpBeam1;
  } 
  else throw;

  double     K = 2*_rpOverGamma;
    
  std::pair<std::pair<double,double>, std::pair<double,double> > E = BaseErsk(sepx_mm, sepy_mm);

  double     Ex = E.first.first;
  double     Ey = E.first.second;
  double     w1x = E.second.first;
  double     w1y = E.second.second;
  
  double     dfleix = K*nBeamOpp*Ex;
  double     dfleiy = K*nBeamOpp*Ey;

  //  std::cout << " VdmBBDeflect::GetDeflectAngles " << " Ex= " << Ex << " Ey= " << Ey << " w1x= " << w1x << " w1y= " << w1y << " K= " << K << " nBeamOpp= " << nBeamOpp << std::endl;

  return std::pair<double, double>(dfleix, dfleiy);
}

std::pair<double, double> VdmBBDeflect::GetDeflectDist(double sepx_mm, double sepy_mm, unsigned int beam)
{
  std::pair<double, double> deflAng = GetDeflectAngles(sepx_mm, sepy_mm, beam);

  double orbx = deflAng.first*_betax*(1./(2.*std::tan(M_PI*_tunex)));
  double orby = deflAng.second*_betay*(1./(2.*std::tan(M_PI*_tuney)));

  // Convert from m to mm
  //
  return std::pair<double, double>(orbx*1000, orby*1000);
}

std::pair<double, double> VdmBBDeflect::GetDeflectSeparationAdjust(double sepx, double sepy)
{
  std::pair<double, double> deflDistBeam1 = GetDeflectDist(sepx, sepy, 1);
  std::pair<double, double> deflDistBeam2 = GetDeflectDist(sepx, sepy, 2);

  return std::pair<double, double>(deflDistBeam1.first + deflDistBeam2.first, 
				   deflDistBeam1.second + deflDistBeam2.second); 
}

std::pair<double, double> VdmBBDeflect::GetDeflectAdjustedSeparation(const std::pair<double, double>& inSep)
{
  double sepx = inSep.first;
  double sepy = inSep.second;

  std::pair<double, double> deflDistBeam1 = GetDeflectDist(sepx, sepy, 1);
  std::pair<double, double> deflDistBeam2 = GetDeflectDist(sepx, sepy, 2);

  return std::pair<double, double>(sepx + deflDistBeam1.first + deflDistBeam2.first, 
				   sepy + deflDistBeam1.second + deflDistBeam2.second); 
}

std::pair<double,double>  VdmBBDeflect::GetDynamicBetaNomRatio(double sep_mm, int directionXY)
{
  if (!HaveDynBetaTable()) throw;

  double lookup;
  if (directionXY == 1) lookup = std::sqrt(2)*std::fabs(sep_mm)/_Sigx;
  else                  lookup = std::sqrt(2)*std::fabs(sep_mm)/_Sigy;

  std::pair<double, double>  deltaBetaRatio = GetDeltaBetaRatiosDeltaXY(lookup, directionXY);

  double deltaBetaRatioX = deltaBetaRatio.first;
  double deltaBetaRatioY = deltaBetaRatio.second;

  // Now scale the ratios per Witold's instructions
  // 
  deltaBetaRatioX *= _dynBetaScaleX;
  deltaBetaRatioY *= _dynBetaScaleY;

  double betaDynBeta0RatioX = 1 + deltaBetaRatioX;
  double betaDynBeta0RatioY = 1 + deltaBetaRatioY;

  return std::pair<double, double>(betaDynBeta0RatioX, betaDynBeta0RatioY);
}

std::pair<double,double>  VdmBBDeflect::GetDynamicBetaRefRatio(double sep_mm,int directionXY)
{
  if (!HaveDynBetaTable()) throw;

  double lookup;
  if (directionXY == 1) lookup = std::sqrt(2)*std::fabs(sep_mm)/_Sigx;
  else                  lookup = std::sqrt(2)*std::fabs(sep_mm)/_Sigy;

  std::pair<double, double>  deltaBetaRatio = GetDeltaBetaRatiosDeltaXY(lookup, directionXY);

  double deltaBetaRatioX = deltaBetaRatio.first;
  double deltaBetaRatioY = deltaBetaRatio.second;

  // Now scale the ratios per Witold's instructions
  // 
  deltaBetaRatioX *= _dynBetaScaleX;
  deltaBetaRatioY *= _dynBetaScaleY;

  double deltaBetaRatioXRefScaled = (directionXY == 1 ? _deltaxBetaXRef : _deltayBetaXRef) * _dynBetaScaleX;
  double deltaBetaRatioYRefScaled = (directionXY == 1 ? _deltaxBetaYRef : _deltayBetaYRef) * _dynBetaScaleY;

  double betaRatioXRef = 1 + deltaBetaRatioXRefScaled; 
  double betaRatioYRef = 1 + deltaBetaRatioYRefScaled;

  double betaDynBeta0RatioX = 1 + deltaBetaRatioX;
  double betaDynBetaRefRatioX = betaDynBeta0RatioX/betaRatioXRef;

  double betaDynBeta0RatioY = 1 + deltaBetaRatioY;
  double betaDynBetaRefRatioY = betaDynBeta0RatioY/betaRatioYRef;

  return std::pair<double, double>(betaDynBetaRefRatioX, betaDynBetaRefRatioY);
}

double VdmBBDeflect::GetDynamicBetaMuVisCorr(double sepx_mm, double sepy_mm)
{
  if (!HaveDynBetaTable()) throw;

  double xlookup = std::sqrt(2)*std::fabs(sepx_mm)/_Sigx;
  double ylookup = std::sqrt(2)*std::fabs(sepy_mm)/_Sigy;

  std::pair<double, double> betaDynBeta0RatiosDeltaX = GetDynamicBetaRefRatio(sepx_mm, 1);
  std::pair<double, double> betaDynBeta0RatiosDeltaY = GetDynamicBetaRefRatio(sepy_mm, 2);

  double betaDynBetaRefRatioXDeltaX = betaDynBeta0RatiosDeltaX.first;
  double betaDynBetaRefRatioYDeltaX = betaDynBeta0RatiosDeltaX.second;

  double betaDynBetaRefRatioXDeltaY = betaDynBeta0RatiosDeltaY.first;
  double betaDynBetaRefRatioYDeltaY = betaDynBeta0RatiosDeltaY.second;

  double corrFactX = std::sqrt(betaDynBetaRefRatioXDeltaX)*std::sqrt(betaDynBetaRefRatioYDeltaX)*std::exp(-xlookup*xlookup/4*(1 - 1./betaDynBetaRefRatioXDeltaX));
  double corrFactY = std::sqrt(betaDynBetaRefRatioXDeltaY)*std::sqrt(betaDynBetaRefRatioYDeltaY)*std::exp(-ylookup*ylookup/4*(1 - 1./betaDynBetaRefRatioYDeltaY));

  return corrFactX*corrFactY;
}

std::pair<double, double> VdmBBDeflect::GetDeltaBetaRatiosDeltaXY(double xylookup, int XY)
{
  DynBetaLookupMap::const_iterator LBIter;
  if (XY == 1) {
    LBIter = _dynBetaLookupDeltaX.lower_bound(xylookup);
    if (LBIter == _dynBetaLookupDeltaX.end()) return std::pair<double, double>(0, 0);
  }
  else if (XY == 2) {
    LBIter = _dynBetaLookupDeltaY.lower_bound(xylookup);
    if (LBIter == _dynBetaLookupDeltaY.end()) return std::pair<double, double>(0, 0);
  }
  else throw;

  if (LBIter->first == xylookup) return LBIter->second;
  else {
    DynBetaLookupMap::const_iterator LBIterBack = LBIter;
    LBIterBack--;
    
    double xy1 = LBIterBack->first;
    double xy2 = LBIter->first;
    double scale = (xylookup - xy1)/(xy2 - xy1);

    double delBetax1 = LBIterBack->second.first;
    double delBetax2 = LBIter->second.first;

    double delBetay1 = LBIterBack->second.second;
    double delBetay2 = LBIter->second.second;

    double delBetaXInterp = delBetax1 + (delBetax2 - delBetax1)*scale; 
    double delBetaYInterp = delBetay1 + (delBetay2 - delBetay1)*scale; 

    return std::pair<double, double>(delBetaXInterp, delBetaYInterp);
  }
}

bool VdmBBDeflect::ReadDynBetaTable(std::string fileName)
{
  std::ifstream infile(fileName.c_str(), std::ifstream::in);
  if (!infile.is_open()) {
    return false;
  }

  std::string line;
  getline(infile, line); // first header line, text only
  getline(infile, line); // 2nd header line

  std::string scan, comment1, comment2;
  int version = 0, lines = 0;

  std::istringstream line2Stream(line); 

  line2Stream >> version >> lines >> scan >> comment1 >> comment2;

  // For now, only handle header version #1
  //
  if (version != 1 || lines != 5) return false;

  getline(infile, line); // 3rd header line, text only
  getline(infile, line); // 4th header line, scan information

  //  std::cout << "Reading fourth line, result = " << line << std::endl;

  int nsteps = 0;

  std::istringstream line4Stream(line); 
  line4Stream >> nsteps >> _NBeamRef >> _EmittanceRef >> _Qx >> _Qy;

  //  Remove the integer part of the tune values
  //
  _Qx -= std::floor(_Qx);
  _Qy -= std::floor(_Qy);

  getline(infile, line); // 5th header line, text only

  //  Now read in the data from the table
  //
  for (size_t ipoint = 0; ipoint < nsteps; ipoint++) {
    getline(infile, line); // 5th header line, text only
    //    if (!infile.good()) return false;

    //    std::cout << "read " << ipoint << "th line, string = " << line << std::endl;

    double sep, dbxRatioXscan, dbyRatioXscan, dbxRatioYscan, dbyRatioYscan;
    std::istringstream dataLineStream(line);
    
    dataLineStream >> sep >> dbxRatioXscan >> dbyRatioXscan >> dbxRatioYscan >> dbyRatioYscan;
    //    if (!dataLineStream.good()) return false;

    _dynBetaLookupDeltaX.insert(DynBetaLookupMap::value_type(sep, std::pair<double, double>(dbxRatioXscan, dbyRatioXscan)));
    _dynBetaLookupDeltaY.insert(DynBetaLookupMap::value_type(sep, std::pair<double, double>(dbxRatioYscan, dbyRatioYscan)));

    std::cout << "sep = " << sep << ", Xscan dbx, dby ratios = " << dbxRatioXscan << ", " << dbyRatioXscan
	      << ", Yscan dbx, dby ratios = " << dbxRatioYscan << ", " << dbyRatioYscan << std::endl;
  }
  
  // Need to scale up the NbeamRef by 10^11
  //
  _NBeamRef *= 1e11;

  _haveDynBetaTable = true;

  //  The reference beta values are taken to be beta_dyn at zero separation.
  //
  std::pair<double, double> deltaxBetaRatiosOrigin = GetDeltaBetaRatiosDeltaXY(0, 1);
 
  _deltaxBetaXRef = deltaxBetaRatiosOrigin.first;
  _deltaxBetaYRef = deltaxBetaRatiosOrigin.second;

  std::pair<double, double> deltayBetaRatiosOrigin = GetDeltaBetaRatiosDeltaXY(0, 2);
 
  _deltayBetaXRef = deltayBetaRatiosOrigin.first;
  _deltayBetaYRef = deltayBetaRatiosOrigin.second;

  return true;
}

std::pair<std::pair<double,double>, std::pair<double,double> > VdmBBDeflect::BaseErsk(double sepx, double sepy)
{
  //The BassErsk routine compute the x and y electric field vectors
  //produced from one bunch by the opposite one with Np particles
  //located at a distance (x,y)  
    
  // Convert all dimensions to meters and take absolute value of deflections
  //
  double x = std::fabs(sepx/1000.);
  double y = std::fabs(sepy/1000.);
  
  double sigmax = _Sigx/1000.;
  double sigmay = _Sigy/1000.;
  
  //I changed the constant factor in front of the fields to have it in
  //units of rp and gamma (K= 2.*rp/gamma) and not eps0 which I have
  //set to 1 instead of eps0 = 8.854187817620e-12 
  
  static const double eps0 = 1;
  
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

std::complex<double> VdmBBDeflect::wfun(std::complex<double> z)
{
  std::complex<double> w = Faddeeva::w(z);
    
  double     wx = w.real();
  double     wy = w.imag();
    
  return std::complex<double>(wx, wy);    
}
std::pair<double, double> VdmBBDeflect::GetExEy(double sepx_mm, double sepy_mm)
{
  std::pair<std::pair<double,double>, std::pair<double,double> > bE = BaseErsk(sepx_mm, sepy_mm);

  return std::pair<double, double>(bE.first.first, bE.first.second);
}

