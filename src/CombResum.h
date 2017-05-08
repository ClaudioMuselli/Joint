#ifndef __COMBRESUM_H__
#define __COMBRESUM_H__

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <complex>
#include <string>
#include "FixedptResum.h"
#include "Joint.h"
#include "../math/complex_def.h"
#include "../math/integration.h"
#include "Luminosity.h"
#include <functional>
#include <vector>
#include "PhisConst.h"
#include "Borel.h"

using namespace std;
using namespace sp_bessel;
using namespace LHAPDF;
using namespace integration;
using namespace ConstResum;




//WARNING::CONTROLLARE CHE CI SIA N-1 GIUSTO!!!!
//FIXME:: WILSON NOT IMPLEMENTED
 
  
  
  class CombResum {
  public:
    CombResum(int ordres, int ordmatch, int channel, string PDFNAME, bool Wilson, 
	      long double mH=125.09, long double Nc=3., long double Nf=5., long double Mur=125.09, long double Muf=125.09);
    virtual ~CombResum();
    long double alpha_s_muR(long double mu);
    
    void SetMatchingFunction(std::complex<long double> (Func) (std::complex<long double>, long double));
    void SetMUF(long double Muf);
    void SetMUR(long double Mur);
    void SetChannel(int channel);
    void SetOrdRes(int ordres);
    void SetOrdMatch(int ordmatch);
    
    long double  ResummedCrossSection(long double CMS, long double xp, int matching, long double *err);
    std::vector<long double> ResummedCrossSection (long double CMS, std::vector<long double> xp, int matching, std::vector<long double> *err);
    std::vector<std::vector<long double>> ResummedCrossSection (std::vector<long double> CMS, std::vector<long double> xp, int matching
      , std::vector<std::vector<long double>> *err);
    
    long double  MatchedCrossSection(long double CMS, long double xp);
    std::vector<long double> MatchedCrossSection (long double CMS, std::vector<long double> xp);
    std::vector<std::vector<long double>> MatchedCrossSection (std::vector<long double> CMS, std::vector<long double> xp);
    
    //Core Functions
    
    std::complex<long double> JointPartRes(std::complex<long double> N, std::complex<long double> b, long double xp);
    std::complex<long double> FixptPartRes(std::complex<long double> N, long double xp);
    std::complex<long double> JointPartMatch(std::complex<long double> N, std::complex<long double> b, long double xp);
    std::complex<long double> FixptPartMatch(std::complex<long double> N, long double xp);
    
  private:
    Luminosity _Lumi;
    FixedptResum *_Fix;
    Joint *_Joint;
    
    std::function<std::complex<long double>(std::complex<long double>, long double)> _MatchFun;
    
    
    long double _CMS;
    int _ordres,_ordmatch;
    
    int _channel;
    bool _Wilson;
    
    long double tauprime;
    
    //Important constant
    
    long double beta_0, beta_1, beta_2;
    const long double alpha_Mz=0.118;
    const long double Mz=91.1876;
    
    long double _Nc,_Nf;
    long double _Ca, _Cf;
    
    
    
    
    
  };



#endif