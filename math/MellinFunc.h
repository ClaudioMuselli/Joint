#ifndef __MELLINFUNC_H__
#define __MELLINFUNC_H__

#include "HSum.h"
#include "complex_def.h"
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_log.h>

using namespace std;

class MellinFunc{
public:
  MellinFunc();
  virtual ~MellinFunc();
  std::complex<long double> Li3zplus(std::complex<long double> N); // Li3(z)/(1+z)
  std::complex<long double> S12z2plus(std::complex<long double> N); // S12(z^2)/(1+z)
  std::complex<long double> S12mzplus(std::complex<long double> N); // S12(-z)/(1+z)
  std::complex<long double> S12zplus(std::complex<long double> N); // S12(z)/(1+z)
  std::complex<long double> Li3mzplus(std::complex<long double> N); //Li3(-z)/(1+z)
  std::complex<long double> Li2zLogminus(std::complex<long double> N); //Li2(z)Log(1-z)
  std::complex<long double> Li2zLogplus(std::complex<long double> N); //Li2(z)Log(1+z)
  std::complex<long double> Li2mzLogplus(std::complex<long double> N); //Li2(-z)Log(1+z)
  std::complex<long double> Li2mzLogz(std::complex<long double> N); //Li2(-z)Log(z)
  std::complex<long double> Li2zLogzplusminus(std::complex<long double> N); //Li2(z)Log(z)/((1+z)(1-z))
  std::complex<long double> Li2mzLogzplus(std::complex<long double> N); //Li2(-z)Log(z)/(1+z)
  std::complex<long double> Li2mzLogminus2plus(std::complex<long double> N); //Li2(-z)Log[1-z]^2/(1+z)
  std::complex<long double> Li2mzLogplusplus(std::complex<long double> N); //Li2(-z)Log[1+z]/(1+z)
  std::complex<long double> Logz2Logminusminus(std::complex<long double> N); //Log(z)^2Log(1-z)/(1-z)
  std::complex<long double> LogzLogminus2minus(std::complex<long double> N); //Log(z)Log(1-z)^2/(1-z)
  std::complex<long double> Logz3minusplus(std::complex<long double> N); // Log(z)^3 /((1-z)(1+z))
  std::complex<long double> Logplusplus(std::complex<long double> N); //Log(1+z)/(1+z)
  std::complex<long double> Li2zLogplusplus(std::complex<long double> N); //Li2(z)Log(1+z)/(1+z)
  std::complex<long double> Logz2Logplusplus(std::complex<long double> N); //Log(1+z)Log(z)^2/(1+z)
  std::complex<long double> LogzLogplus(std::complex<long double> N); //Log(z)Log(1+z)
  std::complex<long double> Logz2Logplus(std::complex<long double> N); //Log(z)^2Log(1+z)
  std::complex<long double> LogzLogplus2(std::complex<long double> N); //Log(z)Log(1+z)^2
  std::complex<long double> LogzLogplus2plus(std::complex<long double> N); //Log(z)Log(1+z)^2/(1+z)
  std::complex<long double> Li2z(std::complex<long double> N); //Li2(z)
  std::complex<long double> Li2mz(std::complex<long double> N); //Li2(-z)
  std::complex<long double> S12z(std::complex<long double> N); //S12(z)
  std::complex<long double> S12mz(std::complex<long double> N); //S12(mz)
  std::complex<long double> S12z2(std::complex<long double> N); //S12(z^2)
  std::complex<long double> Li3z(std::complex<long double> N); //Li3(z)
  std::complex<long double> Li3mz(std::complex<long double> N); //Li3(mz)
  std::complex<long double> Li2zLogz(std::complex<long double> N); //Li2(z)Log(z)
  std::complex<long double> Logminus(std::complex<long double> N); //Log(1-z)
  std::complex<long double> Logplus(std::complex<long double> N); //Log(1+z)
  std::complex<long double> Logz(std::complex<long double> N); //Log(z)
  std::complex<long double> Logz2(std::complex<long double> N); //Log(z)^2
  std::complex<long double> Logz3(std::complex<long double> N); //Log(z)^3
  std::complex<long double> LogzLogminus(std::complex<long double> N); //Log(z)Log(1-z)
  std::complex<long double> LogzLogminus2(std::complex<long double> N); //Log(z)Log(1-z)^2
  std::complex<long double> Logz2Logminus(std::complex<long double> N); //Log(z)^2Log(1-z)
  std::complex<long double> Logz2minus(std::complex<long double> N); //log(z)^2/(1-z)
  std::complex<long double> Logminus2(std::complex<long double> N); //Log(1-z)^2
  std::complex<long double> Logminus3(std::complex<long double> N); //Log(1-z)^3
  std::complex<long double> LogzLogminusminus(std::complex<long double> N); //Log(z)Log(1-z)/(1-z)
  std::complex<long double> Logzminus(std::complex<long double> N); //log(z)/(1-z)
  std::complex<long double> Logzminusplus(std::complex<long double> N); //log(z)/((1-z)(1+z))
  std::complex<long double> plus(std::complex<long double> N); //1/(1+z)
  std::complex<long double> Li2minusminus(std::complex<long double> N); //Li2(1-z)/(1-z)
  std::complex<long double> S12zregminus(std::complex<long double> N); //(S12(z)-Zeta(3))/(1-z)
  std::complex<long double> Li2zregminus(std::complex<long double> N); //(Li2(z)-Zeta(2))/(1-z)
  std::complex<long double> Li3zregminus(std::complex<long double> N); //(Li3(z)-Zeta(3))/(1-z)
  std::complex<long double> Li2zregLogminusminus(std::complex<long double> N); //(Li2(z)-Zeta(2))Log(1-z)/(1-z)
  
  std::complex<long double> Logplus3plus(std::complex<long double > N); //Log(1+z)^3/(1+z)
  std::complex<long double> Li3zregminusplus(std::complex<long double> N); //(Li3-Zeta[3])/(1-z)/(1+z)
  std::complex<long double> Li3zoverplusplus(std::complex<long double> N); //(Li3(z/(1+z))/(1+z)
  std::complex<long double> Logplus3(std::complex<long double> N); // Log[1+z]^3
  std::complex<long double> Li2z2(std::complex<long double> N); //Li2(z^2)
  std::complex<long double> Li2z2Logz(std::complex<long double> N); // Li2(z^2) Log(z)
  std::complex<long double> Li3z2(std::complex<long double> N); //Li3(z^2)
  std::complex<long double> Li3overplus(std::complex<long double> N); //Li3(1/(1+z))
  std::complex<long double> Li2minus(std::complex<long double> N); //Li2(1-z)
  
  std::complex<long double> D0(std::complex<long double> N); //[1/(1-z)]_+
  std::complex<long double> D1(std::complex<long double> N); //[Log(1-z)/(1-z)]_+
  std::complex<long double> D2(std::complex<long double> N); //[Log(1-z)^2/(1-z)]_+
  std::complex<long double> D3(std::complex<long double> N); //[Log(1-z)^3/(1-z)]_+
  
  
  
  
private:
  HSum H;
  long double zeta2;
  long double zeta3;
  long double Li4;
  long double log2;
  long double log2q;
  long double log2c;
  long double zeta2q;
  long double EulerGamma;
  
};



#endif