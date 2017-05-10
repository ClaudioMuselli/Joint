#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <Joint.h>
#include "math/integration.h"
#include <LHAPDF/LHAPDF.h>
#include <complex_bessel.h>
#include "src/FixedptResum.h"
#include "math/complex_def.h"
#include "PhisConst.h"
#include <CombResum.h>


using namespace std;
using namespace LHAPDF;
using namespace sp_bessel;
using namespace integration;


std::complex<long double> PROVAN(std::complex<long double> N, void * par){
  return (std::log(N*N+25.));
}

std::complex<long double> PROVAN2(std::complex<long double> N, void * par){
  return (1./N);
}

std::complex<long double> PROVAB(std::complex<long double> b, void * par){
  return (-std::log(b));
}
std::complex<long double> PROVAT(std::complex<long double> N, std::complex<long double> b, void * par){
  return (1./N*1./b);
}
std::complex<long double> PROVAT2(std::complex<long double> N, std::complex<long double> b, void * par){
  return (std::pow(std::log(N*N+b*b/4.),1)/(N*N));
}

double test(double x){
  return ((x+2.)*std::log(1.+x)/x);
}

int main(){
  std::complex<long double> NTest(2.,1.);
  std::complex<long double> NTestd(0.5,-20.);
  std::complex<long double> NTest2(3.,1e20);
   std::complex<double> NTest2d(3.,1e20);
   std::complex<long double> Norder(1.4,0.7);
 
  cout.precision(10.);
  
  //cout << Hyp2F1(NTestd,NTest,NTest+1.,NTest+2.5) << endl;
  
 /*FixedptResum Fix(2.,2.,0.,false);
 cout << "LOgggH= " << Fix.LOgggH(NTest,1.2) << endl;
  cout << "LOgqqH= " << Fix.LOgqqH(NTest,1.2) << endl;
  cout << "LOqqgH= " << Fix.LOqqgH(NTest,1.2) << endl;
  long double xppp=0.3;
  std::complex<long double> xprad(std::pow(std::sqrt(1.+xppp)-std::sqrt(xppp),4.),0.);
  std::complex<long double> half(0.5,0.);
  cout << Hyp2F1(half,NTestd,NTestd+0.5,xprad) << endl;
  cout << Hyp2F1(half,NTestd+1.,NTest+1.5,xprad) << endl;
  cout << Hyp2F1(half,NTestd+2.,NTest+2.5,xprad) << endl;
  cout << Hyp2F1(half,NTest+3.,NTestd+3.5,xprad) << endl;*/
  
  
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
  
  /*Joint JLL(0,2,0,true);
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
  
  //cout << "Test Luminosity" << endl;
  //Luminosity Lum(LHAPDF::mkPDF("NNPDF30_nnlo_as_0118"),125.09);
  //Lum.TestLum("TestLum.dat");*/
  
  
  long double x=0.7, xp=std::pow(5./125.09,2);
  long double tau=x*std::pow(std::sqrt(1.+xp)-std::sqrt(xp),2);
  //cout << xp << " " << tau << endl;
  /*
  cout<< "TEST INVERSE MELLIN" << endl;
  //cout << integration::InverseMellin_path(3,PROVAN,0.1,3.,2.,NULL) << endl;
  //cout << integration::InverseMellin_epsilon(3,PROVAN2,tau,4.5,NULL) << endl;
  //cout << integration::InverseMellin_epsilon(3,PROVAN2,tau,1.,NULL) << endl;
  //cout << integration::InverseMellin_epsilon(3,PROVAN2,tau,2.,NULL) << endl;
  ///cout << integration::InverseMellin_epsilon(3,PROVAN2,tau,3.,NULL) << endl;
  
  //cout << "TEST INVERSE BESSEL" << endl;
  //cout << integration::InverseBessel(3,PROVAB,0.3,2.,15.,2.5,NULL) << endl;
  //cout << integration::InverseBessel(3,PROVAB,0.3,4.,15.,2.5,NULL) << endl;
  
  //cout << "TEST INVERSE TOTAL" << endl;
  //double ris=0.,error=0.,chi=0.;
  //integration::InverseMellinBessel(3,PROVAT2,1.,0.1,NULL,&ris,&error,&chi,true);
  //cout << "integral= " << ris << " +- " << error << "\t chisq= " << chi << endl;
  
  cout << "TEST BESSEL" << endl;
  cout << CBesselK(Norder,NTest)<< endl;
  cout << CBesselK(Norder,NTest2)<< endl;
  cout << CBesselK(Norder,NTestd)<< endl;
  
  //cout << "TEST LOG" << endl;
  //cout << log_c_angle(NTest,-M_PIl/2.) << endl;
  //cout << log_c_angle(-3.-0.000001*II,-M_PIl/2.) << endl;
  //cout << log_c_angle(-3.+0.000001*II,-M_PIl/2.) << endl;
  //cout << log_c_angle(-0.000001-3.*II,-M_PIl/2.) << endl;
  //cout << log_c_angle(0.000001-3.*II,-M_PIl/2.) << endl;
  
  
  
   cout << "TEST TOTAL Fixed pt" << endl;
  FixedptResum FixLL(0,3,0,true);
  FixedptResum FixNLL(1,3,0,true);
  FixedptResum FixNNLL(2,3,0,true);
  long double xpp=1.2;
  cout << "Sudakov" << endl;
  cout << FixLL.Sudakov_th_gggH(NTest,xpp) << endl;
  cout << FixNLL.Sudakov_th_gggH(NTest,xpp) << endl;
  cout << FixNNLL.Sudakov_th_gggH(NTest,xpp) << endl;
  cout << endl;
  cout << FixLL.Sudakov_th_gqqH(NTest,xpp) << endl;
  cout << FixNLL.Sudakov_th_gqqH(NTest,xpp) << endl;
  cout << FixNNLL.Sudakov_th_gqqH(NTest,xpp) << endl;
  cout << endl;
  cout << FixLL.Sudakov_th_qqgH(NTest,xpp) << endl;
  cout << FixNLL.Sudakov_th_qqgH(NTest,xpp) << endl;
  cout << FixNNLL.Sudakov_th_qqgH(NTest,xpp) << endl;
  cout << endl;
  cout << FixNNLL.Hth1gggH(xpp) << endl;
  cout << FixNNLL.Hth1gqqH(xpp) << endl;
  cout << FixNNLL.Hth1qqgH(xpp) << endl;
  std::vector<std::complex<long double> > resLL,resNLL,resNNLL;
  resLL=FixLL.ComputeFixedptResum(NTest,xpp);
  resNLL=FixNLL.ComputeFixedptResum(NTest,xpp);
  resNNLL=FixNNLL.ComputeFixedptResum(NTest,xpp);
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
  matchLL=FixLL.ComputeMatching(NTest,xpp);
  matchNLL=FixNLL.ComputeMatching(NTest,xpp);
  matchNNLL=FixNNLL.ComputeMatching(NTest,xpp);
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
  cout << endl;*/
  
  long double CMS=13000.,err=0.;
  CombResum Final(0,0,1,"PDF4LHC15_nnlo_100",true);
  long double ptstart=1.; long double ptend=20.;
  int NUM=20;
  long double ratexp=(ptend-ptstart)/((long double) NUM);
  for (int i=0;i<=NUM;i++){
    long double pt=(ptstart+i*ratexp);
    long double xxpp=(ptstart+i*ratexp)*(ptstart+i*ratexp)/(125.09*125.09);
    long double FinalRes=Final.ResummedCrossSection(CMS,xxpp,0,&err);
    cout << "pt= " << pt << "ris= " << FinalRes*2.*(ptstart+i*ratexp)/(125.09*125.09) 
    << " +- "<< err*2.*(ptstart+i*ratexp)/(125.09*125.09) << endl;
  }
  
  
  
 
  
  
  
  
  
  
  
  
  
}