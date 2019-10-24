#ifndef _VdmBBDeflect_h
#define _VdmBBDeflect_h

#include <map>
#include <iostream>
#include <utility>
#include <Faddeeva.hh>

class VdmBBDeflect
{
public:
  static const double rp;

private:
  static double _betax;
  static double _betay;
  static double _tunex;
  static double _tuney;

  static double _EMagnetic;
  static double _ABeam;
  static double _gamma;
  static double _rpOverGamma;
  static bool _haveBeamParams;
  static bool _debug;

  //  Data for dynamic beta correction
  //
  static bool _haveDynBetaTable;
 
  typedef std::map<double, std::pair<double, double> > DynBetaLookupMap;

  static DynBetaLookupMap _dynBetaLookupDeltaX;
  static DynBetaLookupMap _dynBetaLookupDeltaY;

  static double _EmittanceRef;
  static double _NBeamRef;
  static double _Qx;
  static double _Qy;
  static double _deltaxBetaXRef;
  static double _deltaxBetaYRef;
  static double _deltayBetaXRef;
  static double _deltayBetaYRef;  

  //  Non-static data
  //
  double _Sigx;
  double _Sigy;
  double _NpBeam1;
  double _NpBeam2;
  double _dynBetaScaleX;
  double _dynBetaScaleY;

  std::pair<std::pair<double,double>, std::pair<double,double> > BaseErsk(double sepx, double sepy);

  static std::complex<double> wfun(std::complex<double> z);
  
  void CalculateScaleFactors();

public:

  VdmBBDeflect(double Sigx, double Sigy, double NpBeam1, double NpBeam2) : 
    _Sigx(Sigx), _Sigy(Sigy), _NpBeam1(NpBeam1), _NpBeam2(NpBeam2)
  {
    CalculateScaleFactors();
  }

  static bool ReadDynBetaTable(std::string filename);
  static bool HaveDynBetaTable() {return _haveDynBetaTable;}

  static void SetBeamParameters(double betax, double betay, double tunex, double tuney, double EMagnetic, double ABeam);
  static bool HaveBeamParams() {return _haveBeamParams;}

  // For use in evaluating systematic uncertainties
  //
  static void ChangeBetaStar(float newBeta);
  static void ChangeTune(float tunex, float tuney);

  std::pair<double, double> GetDeflectAngles(double sepx, double sepy, unsigned int beam = 1);
  std::pair<double, double> GetDeflectDist(double sepx, double sepy, unsigned int beam = 1);
  std::pair<double, double> GetDeflectSeparationAdjust(double sepx, double sepy);
  std::pair<double, double> GetDeflectAdjustedSeparation(const std::pair<double, double>& inSep);

  std::pair<double,double> GetDynamicBetaRefRatio(double sep_mm, int direction);
  std::pair<double,double> GetDynamicBetaNomRatio(double sep_mm, int direction);

  double GetDynamicBetaMuVisCorr(double sepx_mm, double sepy_mm);

  std::pair<double,double> GetDynamicBetaRefRatioX(double sepx_mm) {
    return GetDynamicBetaRefRatio(sepx_mm,  1);
  }

  std::pair<double,double> GetDynamicBetaRefRatioY(double sepy_mm) {
    return GetDynamicBetaRefRatio(sepy_mm, 2);
  }

  std::pair<double,double> GetDynamicBetaNomX(double sepx_mm) {
    return GetDynamicBetaNomRatio(sepx_mm,  1);
  }

  std::pair<double,double> GetDynamicBetaNomY(double sepy_mm) {
    return GetDynamicBetaNomRatio(sepy_mm, 2);
  }

  static std::pair<double, double> GetDeltaBetaRatiosDeltaXY(double xlookup, int XY);
  std::pair<double, double> GetExEy(double sepx, double sepy);
};

#endif //#ifndef _VdmBBDeflect_h
