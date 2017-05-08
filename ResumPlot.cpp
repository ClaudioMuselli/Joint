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
  int orderres=0;
  ss << "graph/PlotResum_" << CMS/1000.;
  switch (orderres){
    case(0):{
      ss <<"_LL.dat";
      break;
    }
    case (1):{
      ss << "_NLL.dat";
      break;
    }
    case(2):{
      ss << "_NNLL.dat" ;
      break;
    }
  }
  ofstream OUT((ss.str()).c_str());
  std::vector<long double> Joint,Fix,Matched,errJoint,errFix,errMatched;
  CombResum Com(orderres,0,0, "PDF4LHC15_nnlo_100",false);
  int NUM=100;
  //pt grid
  long double ptstart=12.;long double ptend=250.;
  long double rate=(ptend-ptstart)/((long double) NUM);
  std::vector<long double> xp;
  for (int i=0;i<NUM;i++){
    xp.push_back(std::pow((ptstart+((long double) i)*rate)/mH,2));
  }
  Joint=Com.ResummedCrossSection(CMS,xp,0,&errJoint);
  Fix=Com.ResummedCrossSection(CMS,xp,1,&errFix);
  Matched=Com.ResummedCrossSection(CMS,xp,2,&errMatched);
  for (int i=0;i<NUM;i++){
    OUT << CMS << "\t" << xp[i] << "\t" << Joint[i] << "\t" <<errJoint[i]
    << "\t" << Fix[i] << "\t" << errFix[i] << "\t" << Matched[i]
     << "\t" << errMatched[i] << endl;
  }
  OUT.close();
  return 0;
}