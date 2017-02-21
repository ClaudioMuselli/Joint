#ifndef _JOINT_H_
#define _JOINT_H_

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include "../math/HSum.h"
#include <vector>
#include <string>

using namespace std;

namespace Joint{
  
  
  
  class Joint{
  public: 
    Joint(const int ordres, const int ordmatch, string channel);
    virtual ~Joint();
    
    std::vector<std::complex<long double> > ComputeJointRes(std::complex<long double> N, std::complex<long double> b);
    std::vector<std::complex<long double> > ComputeMatching(std::complex<long double> N, std::complex<long double> b);
    
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
    
    
    
    std::complex<long double> **U=NULL, **V1=NULL, **V2=NULL;
    
    //Higgs Cross section
    //Coefficients anomalous dimensions
    const long double Apt1g;
    const long double Apt2g;
    const long double Apt3g;
    const long double Bpt1g;
    const long double Bpt2g;
    const long double Dpt2g;
    const long double Hpt1g;
    const long double Hpt2g;
    
    //Sudakov
    std::complex<long double> Sudakov_g(std::complex<long double> N, std::complex<long double> b);
    //Hard part
    std::complex<long double> Hgggq1(std::complex<long double> N);
    std::complex<long double> Hgggq2(std::complex<long double> N);
    std::complex<long double> Hgggg2(std::complex<long double> N);
    std::complex<long double> Hggqq2(std::complex<long double> N);
    //Evolution
    void ComputeEvolution(std::complex<long double> N);
    
  };
  

  
  
  


};














































#endif