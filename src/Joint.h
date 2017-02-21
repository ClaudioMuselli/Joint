#ifndef _JOINT_H_
#define _JOINT_H_

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include <vector>
#include <string>
#include <sstream>
#include "AP.h"
#include "../math/MellinFunc.h"

using namespace std;

class Joint{
  public: 
    Joint(const int ordres, const int ordmatch, string channel, bool Wilson);
    virtual ~Joint();
    
    std::vector<std::complex<long double> > ComputeJointRes(std::complex<long double> N, std::complex<long double> b);
    std::vector<std::complex<long double> > ComputeMatching(std::complex<long double> N, std::complex<long double> b);
    
    //Hard part
    std::complex<long double> Hgggq1(std::complex<long double> N);
    std::complex<long double> Hgggq2(std::complex<long double> N);
    std::complex<long double> Hggggreg2(std::complex<long double> N);
    std::complex<long double> Hggqq2(std::complex<long double> N);
    
  private:
    int _ordres,_ordmatch;
    //Scale and alpha_s
    long double _LF=0.1,_LR=0.1;
    long double _as=1./128.;
    //Important Constants
    const long double _Nc=3., _Nf=5.;
    long double _Ca, _Cf;
    //Channels
    string _channel;
    bool _Wilson;
    
    long double EulerGamma=0.577215664901532860606512090082;
    
    
    GammaAP AP;
    std::complex<long double> **ULL=NULL,**UNLL=NULL,**UNNLL=NULL, **V1=NULL;
    
    //Higgs Cross section
    //Coefficients anomalous dimensions
    long double Apt1g;
    long double Apt2g;
    long double Apt3g;
    long double Bpt1g;
    long double Bpt2g;
    long double Dpt2g;
    long double Hpt1g;
    long double Hpt2g;
    
    long double W1;
    long double W2;
    
    //beta function
    long double beta_0,beta_1,beta_2;
    
    //Sudakov
    std::complex<long double> Sudakov_g(std::complex<long double> N, std::complex<long double> b);
    
    //Evolution
    void ComputeEvolution(std::complex<long double> N);
    
};














































#endif