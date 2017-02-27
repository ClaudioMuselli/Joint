#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include "math/MellinFunc.h"
#include "src/AP.h"
#include <Joint.h>
#include "src/Luminosity.h"
#include <LHAPDF/LHAPDF.h>


using namespace std;
using namespace LHAPDF;






int main(){
  std::complex<long double> NTest(2.,1.);
  std::complex<long double> NTest2(5.,3.);
  cout.precision(10.);
  /*MellinFunc Fun;

  cout <<"D0 " << Fun.D0(NTest)<< endl;
  cout <<"D1 " << Fun.D1(NTest)<< endl;
  cout <<"D2 " << Fun.D2(NTest)<< endl;
  cout << "Li3zplus " << Fun.Li3zplus(NTest) << endl;
  cout << "S12z2plus " << Fun.S12z2plus(NTest) << endl;
  cout << "S12mzplus " << Fun.S12mzplus(NTest) << endl;
  cout << "S12zplus " << Fun.S12zplus(NTest) << endl;
  cout << "Li3mzplus " << Fun.Li3mzplus(NTest) << endl;
  cout << "Li2zLogzplusminus " << Fun.Li2zLogzplusminus(NTest) << endl;
  cout << "Li2mzLogzplus " << Fun.Li2mzLogzplus(NTest) << endl;
  cout << "LogzLogminus2minus " << Fun.LogzLogminus2minus(NTest) << endl;
  cout << "Li2mzLogplusplus" << Fun.Li2mzLogplusplus(NTest) << endl;
  cout << "Li2zLogminus " << Fun.Li2zLogminus(NTest) << endl;
  cout << "Logz3minusplus " << Fun.Logz3minusplus(NTest) << endl;
  cout << "Logplusplus " << Fun.Logplusplus(NTest) << endl;
  cout << "Li2zLogplusplus " << Fun.Li2zLogplusplus(NTest) << endl;
  cout << "Logz2Logplusplus " << Fun.Logz2Logplusplus(NTest) << endl;
  cout << "LogzLogplus " << Fun.LogzLogplus(NTest) << endl;
  cout << "LogzLogplus2plus " << Fun.LogzLogplus2plus(NTest) << endl;
  cout << "Li2z " << Fun.Li2z(NTest) << endl;
  cout << "Li2mz " << Fun.Li2mz(NTest) << endl;
  cout << "Logminus " << Fun.Logminus(NTest) << endl;
  cout << "Logminus2 " << Fun.Logminus2(NTest) << endl;
  cout << "LogzLogminusminus " << Fun.LogzLogminusminus(NTest) << endl;
  cout << "Logminus3 " << Fun.Logminus3(NTest) << endl;
  cout << "Logzminus " << Fun.Logzminus(NTest) << endl;
  cout << "Logzminusplus " << Fun.Logzminusplus(NTest) << endl;
  cout << "plus " << Fun.plus(NTest) << endl;
  cout << "S12zregminus " << Fun.S12zregminus(NTest) << endl;
  cout << "Li2zregminus " << Fun.Li2zregminus(NTest) << endl;
  cout << "Li3zregminus " << Fun.Li3zregminus(NTest) << endl;
  cout << "Li2zregLogminusminus " << Fun.Li2zregLogminusminus(NTest) << endl;
  cout << "S12z " << Fun.S12z(NTest) << endl;
  cout << "Li3z " << Fun.Li3z(NTest) << endl;
  cout << "Li2zLogz " << Fun.Li2zLogz(NTest) << endl;
  cout << "Logz " << Fun.Logz(NTest) << endl;
  cout << "LogzLogminus2 " << Fun.LogzLogminus2(NTest) << endl;
  cout << "Logz3 " << Fun.Logz3(NTest) << endl;
  cout << "Li2minusminus " << Fun.Li2minusminus(NTest) << endl;
  cout << "Logz2Logminus " << Fun.Logz2Logminus(NTest) << endl;
  cout << "Li3mz " << Fun.Li3mz(NTest) << endl;
  cout << "S12mz " << Fun.S12mz(NTest) << endl;
  cout << "S12z2 " << Fun.S12z2(NTest) << endl;
  cout << "Li2mzLogplus " << Fun.Li2mzLogplus(NTest) << endl;
  cout << "Li2mzLogz " << Fun.Li2mzLogz(NTest) << endl;
  cout << "Logz2Logplus " << Fun.Logz2Logplus(NTest) << endl;
  cout << "LogzLogplus2 " << Fun.LogzLogplus2(NTest) << endl;
  cout << "LogzLogminus " << Fun.LogzLogminus(NTest) << endl;
  cout << "Logplus " << Fun.Logplus(NTest) << endl;
  cout << "Li2Logplus " << Fun.Li2zLogplus(NTest) << endl;
  cout << endl;
  cout << "Logplus3plus " << Fun.Logplus3plus(NTest) << endl;
  cout << "Li3zregminusplus " << Fun.Li3zregminusplus(NTest) << endl;
  cout << "Li3zoverplusplus " << Fun.Li3zoverplusplus(NTest) << endl;
  cout << "Logplus3 " << Fun.Logplus3(NTest) << endl;
  cout << "Li2z2 " << Fun.Li2z2(NTest) << endl;
  cout << "Li2z2Logz " << Fun.Li2z2Logz(NTest) << endl;
  cout << "Li3z2 " << Fun.Li3z2(NTest) << endl;
  cout << "Li3overplus " << Fun.Li3overplus(NTest) << endl;
  cout << "Li2minus " << Fun.Li2minus(NTest) << endl;
  cout << endl;
  cout << "Test finito per funzioni speciale: Test per gamma AP" << endl;
  GammaAP AP(5.,3.);
  AP.ComputeGamma(NTest,2);
  cout << "LO anomalous dimension" << endl;
  cout << "gg0= " << AP.gg0 << endl;
  cout << "gS0= " << AP.gS0 << endl;
  cout << "Sg0= " << AP.Sg0 << endl;
  cout << "SS0= " << AP.SS0 << endl;
  cout << "plus0= " << AP.plus0 << endl;
  cout << "minus0= " << AP.minus0 << endl;
  cout << "NLO anomalous dimension" << endl;
  cout << "gg1= " << AP.gg1 << endl;
  cout << "Sg1= " << AP.Sg1 << endl;
  cout << "gS1= " << AP.gS1 << endl;
  cout << "SS1= " << AP.SS1 << endl;
  cout << "WW1= " << AP.WW1 << endl;
  cout << "TT1= " << AP.TT1 << endl;
  cout << "VV1= " << AP.VV1 << endl;
  cout << "qq1= " << (AP.WW1+AP.TT1)/2.+(AP.SS1-AP.WW1+AP.VV1-AP.TT1)/10. << endl;
  cout << "qqb1= " << (AP.WW1-AP.TT1)/2.+(AP.SS1-AP.WW1-AP.VV1+AP.TT1)/10. << endl;
  cout << "qQ1= " << (AP.SS1-AP.WW1+AP.VV1-AP.TT1)/10. << endl;
  cout << "qQb1= " << (AP.SS1-AP.WW1-AP.VV1+AP.TT1)/10. << endl;
  cout << "NNLO anomaous dimension" << endl;  
  cout << "gg2= " << AP.gg2 << endl;
  cout << "Sg2= " << AP.Sg2 << endl;
  cout << "gS2= " << AP.gS2 << endl;
  cout << "SS2= " << AP.SS2 << endl;
  cout << "WW2= " << AP.WW2 << endl;
  cout << "TT2= " << AP.TT2 << endl;
  cout << "VV2= " << AP.VV2 << endl;
  
  cout << endl;
  
  cout << "Prova DILOG" << endl;
  gsl_sf_result redilog,imdilog;
  long double r=std::abs(NTest*NTest/NTest2);
  long double theta=std::arg(NTest*NTest/NTest2);
  gsl_sf_complex_dilog_e(r,theta,&redilog,&imdilog);
  std::complex<long double> risdilog(redilog.val, imdilog.val);
  cout << risdilog << endl;
  cout << endl;*/
  
  Joint JLL(0,2,0,true);
  Joint JNLL(1,2,0,true);
  Joint JNNLL(2,2,0,true);
  //cout << J.Hgggq1(NTest) << " " << J.Hggggreg2(NTest) << " " << J.Hgggq2(NTest)<< endl;
  
  cout<< endl;
  std::complex<long double> b(3.,5.);
  //cout << "Test Sudakov and evolution" << endl;
  //J.ComputeEvolution(NTest,b);
  cout << "TEST TOTAL" << endl;
  std::vector<std::complex<long double> > resLL,resNLL,resNNLL;
  resLL=JLL.ComputeJointRes(NTest,b);
  resNLL=JNLL.ComputeJointRes(NTest,b);
  resNNLL=JNNLL.ComputeJointRes(NTest,b);
  string name[3]={"GG","GQ","QQ"};
  for (int i=0;i<resLL.size();i++){
    cout << name[i] << "= " << resLL[i] << endl; 
  }
  cout << endl;
  for (int i=0;i<resNLL.size();i++){
    cout << name[i] << "= " << resNLL[i] << endl; 
  }
  cout << endl;
  for (int i=0;i<resNNLL.size();i++){
    cout << name[i] << "= " << resNNLL[i] << endl; 
  }
  cout << endl;
  cout << "TEST MATCHING" << endl;
  std::vector<std::complex<long double> > matchLL,matchNLL,matchNNLL;
  matchLL=JLL.ComputeMatching(NTest,b);
  matchNLL=JNLL.ComputeMatching(NTest,b);
  matchNNLL=JNNLL.ComputeMatching(NTest,b);
  for (int i=0;i<resLL.size();i++){
    cout << name[i] << "= " << matchLL[i] << endl; 
  }
  cout << endl;
  for (int i=0;i<resNLL.size();i++){
    cout << name[i] << "= " << matchNLL[i] << endl; 
  }
  cout << endl;
  for (int i=0;i<resNNLL.size();i++){
    cout << name[i] << "= " << matchNNLL[i] << endl; 
  }
  cout << endl;
  
  cout << "Test Luminosity" << endl;
  Luminosity Lum(LHAPDF::mkPDF("NNPDF30_nnlo_as_0118"),125.09);
  Lum.TestLum("TestLum.dat");
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}