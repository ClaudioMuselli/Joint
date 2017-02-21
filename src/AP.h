#ifndef _AP_H_
#define _AP_H_


#include <gsl/gsl_math.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <complex>
#include "../math/complex_def.h"
#include "../math/HSum.h"

using namespace std;

struct gamma1sums {
  std::complex<long double> S1, S2, S3, S4;
  std::complex<long double> NS,NT,NFO,NFI,NSI,NSE,NE,NN;
  std::complex<long double> NI,NI2,NI3,NMI, NMI2,N1I,N1I2,N1I3,N2I,N1,N2,NM,NMS,N1S,N1T,N2S,N2T;
  std::complex<long double> N3,N4,N5,N6;
  std::complex<long double> S11,S12,S13,S14,S15,S16,S1M,S21,S31,S2M;
  std::complex<long double> SPMOM,SLC,SLV,SSCHLM,SSTR2M,SSTR3M,SSCHLP;
  std::complex<long double> SSTR2P,SSTR3P;
};


class GammaAP{
public:
  GammaAP(const long double Nf,const long double Nc);
  virtual ~GammaAP();
  void ComputeGamma(std::complex<long double> N, int order);
  std::complex<long double> gg0,SS0,Sg0,gS0,WW0,TT0,VV0,plus0,minus0;
  std::complex<long double> gg1,SS1,Sg1,gS1,WW1,TT1,VV1;
  std::complex<long double> gg2,SS2,Sg2,gS2,WW2,TT2,VV2;
  void sums(std::complex<long double> N);
private:
  long double CA,CF;
  long double NF,NC;
  long double EulerGamma=0.577215664901532860606512090082;
  HSum HAP;
  gamma1sums g1s;
  std::complex<long double> PNPA,PNSB,PNSC,PNMA,PPSA,PQGA,PQGB,PGQA,PGQB,PGQC,PGGA,PGGB,PGGC;
  std::complex<long double> P2PLSN,P2MINN,P2VALN,P2QGN,P2GQN,P2GGN,P2PSN;
  
  //LO anomalous dimensions 
  std::complex<long double> gammagg0(std::complex<long double> N);
  std::complex<long double> gammaSg0(std::complex<long double> N);
  std::complex<long double> gammagS0(std::complex<long double> N);
  std::complex<long double> gammansplus0(std::complex<long double> N);
  //NLO important definitions
  void DEF1(std::complex<long double> N);
  //NLO anomalous dimensions
  std::complex<long double> gammagg1(std::complex<long double> N);
  std::complex<long double> gammaSg1(std::complex<long double> N);
  std::complex<long double> gammagS1(std::complex<long double> N);
  std::complex<long double> gammaps1(std::complex<long double> N);
  std::complex<long double> gammansplus1(std::complex<long double> N);
  std::complex<long double> gammansminus1(std::complex<long double> N);
  //NNLO important definitions
  void DEF2(std::complex<long double> N);
  //NNLO anomalous dimensions
  std::complex<long double> gammagg2(std::complex<long double> N);
  std::complex<long double> gammaSg2(std::complex<long double> N);
  std::complex<long double> gammagS2(std::complex<long double> N);
  std::complex<long double> gammaps2(std::complex<long double> N);
  std::complex<long double> gammansplus2(std::complex<long double> N);
  std::complex<long double> gammansminus2(std::complex<long double> N);
  std::complex<long double> gammansval2(std::complex<long double> N);
  
};



#endif
