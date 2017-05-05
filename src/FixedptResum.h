#ifndef __FIXED_PT__
#define __FIXED_PT__

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <vector>
#include <string>
#include <sstream>
#include "AP.h"
#include "../math/MellinFunc.h"
#include "../math/complex_def.h"
#include "PhisConst.h"

using namespace std;
using namespace ConstResum;

//FIXME IMPLEMENTATION OF RUNNING SCALES


class FixedptResum {
public:
    FixedptResum(const int ordres, const int ordmatch, int channel, bool Wilson, long double Nc=3., long double Nf=5.);
    virtual ~FixedptResum();
    
    std::vector<std::complex<long double>> ComputeFixedptResum(std::complex<long double> N, long double xp);
    std::vector<std::complex<long double>> ComputeMatching(std::complex<long double> N, long double xp);
  
    //Higgs Cross section
    //LO cross section
    std::complex<long double> LOgggH(std::complex<long double> NN, long double xp);
    std::complex<long double> LOgqqH(std::complex<long double> NN, long double xp);
    std::complex<long double> LOqqgH(std::complex<long double> NN, long double xp);
    
    //Sudakov
    std::complex<long double> Sudakov_th_gggH(std::complex<long double> N, long double xp);
    std::complex<long double> Sudakov_th_gqqH(std::complex<long double> N, long double xp);
    std::complex<long double> Sudakov_th_qqgH(std::complex<long double> N, long double xp);
    
    //MatchingConstant - hard part
    long double Hth1gggH(long double xp);
    long double Hth1gqqH(long double xp);
    long double Hth1qqgH(long double xp);
    
    void SetChannel(int channel);
    void SetOrdRes(int ordres);
    void SetOrdMatch(int ordmatch);
  
private:
    int _ordres,_ordmatch;
    //Important Constants
    long double _Nc, _Nf;
    long double _Ca, _Cf;
    long double sigma_0;
    //Channels
    int _channel;
    bool _Wilson;
    
    long double EulerGamma=0.577215664901532860606512090082;
    
    
    //Coefficients anomalous dimensions
    long double Ath1g;
    long double Ath2g;
    long double Ath1q;
    long double Ath2q;
    long double Bth1g;
    long double Bth1q;
    
    long double W1;
    long double W2;
    
    //beta function
    long double beta_0,beta_1,beta_2;
    
  
    
     
};

#endif