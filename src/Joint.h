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
#include "PhisConst.h"

using namespace std;
using namespace ConstResum;

//WARNING: N-1 messo nella parte hard ma non nell'evoluzione... non so cosa sia meglio controllare alla fine... 
//In più matching con il NNLO non implementato ancora


class Joint{
  public: 
    Joint(const int ordres, const int ordmatch, int channel, bool Wilson, long double Nc=3., long double Nf=5.);
    virtual ~Joint();
    
    //Return a vector of three components in the case of Higgs which is flavour blind (gg,gq,qq)
    std::vector<std::complex<long double> > ComputeJointRes(std::complex<long double> N, std::complex<long double> lchi);
    std::vector<std::complex<long double> > ComputeMatching(std::complex<long double> N, std::complex<long double> b);
    
    void SetOrdRes(int ordres);
    void SetOrdMatch(int ordmatch);
    int SetChannel(int channel);
    
   
  private:
    int _ordres,_ordmatch;
    //Important Constants
    long double _Nc,_Nf;
    long double _Ca, _Cf;
    //Channels
    int _channel;
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
    
    
    
    //Hard part
    std::complex<long double> Hgggq1(std::complex<long double> NN);
    std::complex<long double> Hgggq2(std::complex<long double> NN);
    std::complex<long double> Hggggreg2(std::complex<long double> NN);
    std::complex<long double> Hggqq2(std::complex<long double> NN);
    
     //Sudakov
    std::complex<long double> Sudakov_g(std::complex<long double> N, std::complex<long double> lchi);
    
    //Evolution
    void ComputeEvolution(std::complex<long double> N, std::complex<long double> lchi);
    
    
    long double Apt1q;
};














































#endif