#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <gsl/gsl_math.h>
#include "CombResum.h"
#include <vector>
#include <complex>
#include <LHAPDF/LHAPDF.h>
#include <fstream>
#include <sstream>

using namespace std;

int main(){
  long double CMS=13000.;
  long double mH=125.09;
  stringstream ss;
  int orderres=1;
  ss << "graph/PlotResumFinal_" << CMS/1000.;
  switch (orderres){
    case(0):{
      ss <<"_LL_tiny.dat";
      break;
    }
    case (1):{
      ss << "_NLL_tiny.dat";
      break;
    }
    case(2):{
      ss << "_NNLL_tiny.dat" ;
      break;
    }
  }
  ofstream OUT((ss.str()).c_str());
  std::vector<long double> Joint,Fix,Matched,errJoint,errFix,errMatched;
  CombResum Com(orderres,0,0, "PDF4LHC15_nnlo_100",false);
  int NUM=100.;
  //pt grid
  long double ptstart=2.;long double ptend=22.;
  long double rate=(ptend-ptstart)/( 20.);
  std::vector<long double> xp;
  for (int i=0;i<20;i++){
    xp.push_back(std::pow((ptstart+((long double) i)*rate)/mH,2));
  }
  ptstart=22.; ptend=250.;
  rate=(ptend-ptstart)/( 20.);
  for (int i=0;i<20;i++){
    xp.push_back(std::pow((ptstart+((long double) i)*rate)/mH,2));
  }
  Joint=Com.ResummedCrossSection(CMS,xp,0,&errJoint);
  Fix=Com.ResummedCrossSection(CMS,xp,1,&errFix);
  Matched=Com.ResummedCrossSection(CMS,xp,2,&errMatched);
  for (int i=0;i<xp.size();i++){
    long double pt=  std::sqrt(xp[i]*mH*mH);
    long double errstandpt=0.;
    long double standpt=Com.ResummedCrossSectionpt(CMS,xp[i],&errstandpt);
    //OUT << CMS << "\t" << pt << "\t" << Joint[i]*2.*pt/mH/mH << "\t" 
    //<<errJoint[i]*2.*pt/mH/mH
    //<< "\t" << Fix[i]*2.*pt/mH/mH << "\t" 
    //<< errFix[i]*2.*pt/mH/mH << endl;
    //OUT << CMS << "\t" << pt << "\t" << Joint[i]*2.*pt/mH/mH << "\t" 
    //<<errJoint[i]*2.*pt/mH/mH
    //<< "\t" <<<< Fix[i]*2.*pt/mH/mH << "\t" 
    //<< errFix[i]*2.*pt/mH/mH << "\t" << Matched[i]*2.*pt/mH/mH
     //<< "\t" << errMatched[i]*2.*pt/mH/mH << endl;
     OUT << CMS << "\t" << pt << "\t" << Joint[i]*2.*pt/mH/mH << "\t" 
    <<errJoint[i]*2.*pt/mH/mH
    << "\t" << Fix[i]*2.*pt/mH/mH << "\t" 
    << errFix[i]*2.*pt/mH/mH << "\t" << Matched[i]*2.*pt/mH/mH
     << "\t" << errMatched[i]*2.*pt/mH/mH 
     << "\t" << standpt*2.*pt/mH/mH  << "\t"<< errstandpt*2.*pt/mH/mH  << endl;
  }
  OUT.close();
  return 0;
}