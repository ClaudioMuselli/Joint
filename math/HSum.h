#ifndef __HSUM_H__
#define __HSUM_H__

#include <iostream>
#include <cmath>
#include <vector>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_zeta.h>
#include <sstream>
#include <string>
#include <functional>
#include "complex_def.h"

using namespace std;

class HSum {	
public:
  HSum(bool verbose=false,bool testinterfun=false, bool testharmsums=false);
  virtual ~HSum();
  std::complex<long double> HS(int i, std::complex<long double> N);
  std::complex<long double> HS(int i, int j, std::complex<long double> N);
  std::complex<long double> HS(int i, int j, int k, std::complex<long double> N);
  std::complex<long double> HS(int i, int j, int k, int m, std::complex<long double> N);
  
  
  
private:
  bool _testinterpolatedfunction;
  bool _testharmonicsums;
  bool _verbose;
  
  void InizializeConst();
  
  //crosscheck function
  
  void TestInterpolatedFunction();
  void TestHarmonicSums(int n);
  std::complex<long double> HS_int(int i, int N);
  std::complex<long double> HS_int(int i, int j, int N);
  std::complex<long double> HS_int(int i, int j, int k, int N);
  std::complex<long double> HS_int(int i, int j, int k, int m, int N);
  
  //Semi-analitics Mellin transform
  std::complex <long double> g1(std::complex<long double > N);
  std::complex <long double> g2(std::complex<long double > N);
  std::complex <long double> g3(std::complex<long double > N);
  std::complex <long double> g4(std::complex<long double > N);
  std::complex <long double> g5(std::complex<long double > N);
  std::complex <long double> g6(std::complex<long double > N);
  std::complex <long double> g7(std::complex<long double > N);
  std::complex <long double> g8(std::complex<long double > N);
  std::complex <long double> g9(std::complex<long double > N);
  std::complex <long double> g10(std::complex<long double > N);
  std::complex <long double> g11(std::complex<long double > N);
  std::complex <long double> g12(std::complex<long double > N);
  std::complex <long double> g13(std::complex<long double > N);
  std::complex <long double> g14(std::complex<long double > N);
  std::complex <long double> g15(std::complex<long double > N);
  std::complex <long double> g16(std::complex<long double > N);
  std::complex <long double> g17(std::complex<long double > N);
  std::complex <long double> g18(std::complex<long double > N);
  std::complex <long double> g19(std::complex<long double > N);
  std::complex <long double> g20(std::complex<long double > N);
  std::complex <long double> g21(std::complex<long double > N);
  std::complex <long double> g22(std::complex<long double > N);
  std::complex <long double> g23(std::complex<long double > N);
  std::complex <long double> g24(std::complex<long double > N);
  std::complex <long double> g25(std::complex<long double > N);
  std::complex <long double> g26(std::complex<long double > N);
  std::complex <long double> g27(std::complex<long double > N);
  std::complex <long double> g28(std::complex<long double > N);
  std::complex <long double> g29(std::complex<long double > N);
  std::complex <long double> g30(std::complex<long double > N);
  std::complex <long double> g31(std::complex<long double > N);
  std::complex <long double> g32(std::complex<long double > N);
  std::complex <long double> g33(std::complex<long double > N);
  std::complex <long double> g34(std::complex<long double > N);
  std::complex <long double> g35(std::complex<long double > N);
  std::complex <long double> g36(std::complex<long double > N);
  std::complex <long double> g37(std::complex<long double > N);
  std::complex <long double> g38(std::complex<long double > N);
  std::complex <long double> g39(std::complex<long double > N);
  
  
  
  
  //Math important Function
  
  
  
  long double Zeta(int i){
  return gsl_sf_zeta_int(i);
  }

 
  
  std::complex<long double> B(int i,std::complex<long double> z){
    return(1./std::pow(2.,(long double)i+1.)*(PolyGamma(i,(1.+z)/2.)-PolyGamma(i,z/2.)));
  }
  
  
  //Important Costants
  long double zeta2;
  long double zeta3;
  long double Li4;
  long double log2;
  long double log2q;
  long double log2c;
  long double zeta2q;
  long double EulerGamma;
  
  
  
  //Vectors of coefficients
  std::vector<long double> a1;
  std::vector<long double> a2;
  std::vector<long double> a3;
  std::vector<long double> b1;
  std::vector<long double> b2;
  std::vector<long double> b3;
  std::vector<long double> c1;
  std::vector<long double> P21;
  std::vector<long double> c2;
  std::vector<long double> P22;
  std::vector<long double> c3;
  std::vector<long double> P23;
  std::vector<long double> P33;
  std::vector<long double> c4;
  std::vector<long double> P24;
  std::vector<long double> P34;
  std::vector<long double> c5;
  std::vector<long double> d5;
  std::vector<long double> q1;
  std::vector<long double> q2;
  std::vector<long double> q3;
  std::vector<long double> q4;
  std::vector<long double> q5;
  std::vector<long double> q6;
  std::vector<long double> q7;
  
  //Table of Harmonic Sums up to weight 4
  //Weight 1
  std::complex<long double> H_1(std::complex<long double> N);
  std::complex<long double> H_m1(std::complex<long double> N);
  //Weight 2
  std::complex<long double> H_2(std::complex<long double> N);
  std::complex<long double> H_m2(std::complex<long double> N);
  std::complex<long double> H_1_1(std::complex<long double> N);
  std::complex<long double> H_1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_1(std::complex<long double> N);
  std::complex<long double> H_m1_m1(std::complex<long double> N);
  //weight 3
  std::complex<long double> H_3(std::complex<long double> N);
  std::complex<long double> H_m3(std::complex<long double> N);
  std::complex<long double> H_m2_m1(std::complex<long double> N);
  std::complex<long double> H_m2_1(std::complex<long double> N);
  std::complex<long double> H_m1_m2(std::complex<long double> N);
  std::complex<long double> H_m1_2(std::complex<long double> N);
  std::complex<long double> H_1_m2(std::complex<long double> N);
  std::complex<long double> H_1_2(std::complex<long double> N);
  std::complex<long double> H_2_1(std::complex<long double> N);
  std::complex<long double> H_2_m1(std::complex<long double> N);
  std::complex<long double> H_1_1_1(std::complex<long double> N);
  std::complex<long double> H_1_1_m1(std::complex<long double> N);
  std::complex<long double> H_1_m1_1(std::complex<long double> N);
  std::complex<long double> H_1_m1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_1_1(std::complex<long double> N);
  std::complex<long double> H_m1_1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_m1_1(std::complex<long double> N);
  std::complex<long double> H_m1_m1_m1(std::complex<long double> N);
  //Weight 4
  std::complex<long double> H_4(std::complex<long double> N);
  std::complex<long double> H_m4(std::complex<long double> N);
  std::complex<long double> H_m3_m1(std::complex<long double> N);
  std::complex<long double> H_m3_1(std::complex<long double> N);
  std::complex<long double> H_m1_m3(std::complex<long double> N);
  std::complex<long double> H_m1_3(std::complex<long double> N);
  std::complex<long double> H_1_m3(std::complex<long double> N);
  std::complex<long double> H_1_3(std::complex<long double> N);
  std::complex<long double> H_3_1(std::complex<long double> N);
  std::complex<long double> H_3_m1(std::complex<long double> N);
  std::complex<long double> H_m2_m2(std::complex<long double> N);
  std::complex<long double> H_m2_2(std::complex<long double> N);
  std::complex<long double> H_2_m2(std::complex<long double> N);
  std::complex<long double> H_2_2(std::complex<long double> N);
  std::complex<long double> H_m2_1_1(std::complex<long double> N);
  std::complex<long double> H_m2_1_m1(std::complex<long double> N);
  std::complex<long double> H_m2_m1_1(std::complex<long double> N);
  std::complex<long double> H_m2_m1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_1_2(std::complex<long double> N);
  std::complex<long double> H_m1_1_m2(std::complex<long double> N);
  std::complex<long double> H_m1_m1_2(std::complex<long double> N);
  std::complex<long double> H_m1_m1_m2(std::complex<long double> N);
  std::complex<long double> H_m1_2_1(std::complex<long double> N);
  std::complex<long double> H_m1_2_m1(std::complex<long double> N);
  std::complex<long double> H_m1_m2_1(std::complex<long double> N);
  std::complex<long double> H_m1_m2_m1(std::complex<long double> N);
  std::complex<long double> H_1_1_2(std::complex<long double> N);
  std::complex<long double> H_1_1_m2(std::complex<long double> N);
  std::complex<long double> H_1_m1_2(std::complex<long double> N);
  std::complex<long double> H_1_m1_m2(std::complex<long double> N);
  std::complex<long double> H_1_2_1(std::complex<long double> N);
  std::complex<long double> H_1_2_m1(std::complex<long double> N);
  std::complex<long double> H_1_m2_1(std::complex<long double> N);
  std::complex<long double> H_1_m2_m1(std::complex<long double> N);
  std::complex<long double> H_2_1_1(std::complex<long double> N);
  std::complex<long double> H_2_1_m1(std::complex<long double> N);
  std::complex<long double> H_2_m1_1(std::complex<long double> N);
  std::complex<long double> H_2_m1_m1(std::complex<long double> N);
  std::complex<long double> H_1_1_1_1(std::complex<long double> N);
  std::complex<long double> H_1_1_1_m1(std::complex<long double> N);
  std::complex<long double> H_1_1_m1_1(std::complex<long double> N);
  std::complex<long double> H_1_1_m1_m1(std::complex<long double> N);
  std::complex<long double> H_1_m1_1_1(std::complex<long double> N);
  std::complex<long double> H_1_m1_1_m1(std::complex<long double> N);
  std::complex<long double> H_1_m1_m1_1(std::complex<long double> N);
  std::complex<long double> H_1_m1_m1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_1_1_1(std::complex<long double> N);
  std::complex<long double> H_m1_1_1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_1_m1_1(std::complex<long double> N);
  std::complex<long double> H_m1_1_m1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_m1_1_1(std::complex<long double> N);
  std::complex<long double> H_m1_m1_1_m1(std::complex<long double> N);
  std::complex<long double> H_m1_m1_m1_1(std::complex<long double> N);
  std::complex<long double> H_m1_m1_m1_m1(std::complex<long double> N);
  
 std::complex<long double> Log(std::complex<long double> N){
  return(std::log(N));
  }

 //Complex Euler Gamma
  std::complex<long double> LogGamma(std::complex<long double> z){
 int g=7;
 double p[9];
 p[0]=0.99999999999980993;
 p[1]=676.5203681218851;
 p[2]=-1259.1392167224028;
 p[3]=771.32342877765313;
 p[4]=-176.61502916214059;
 p[5]=12.507343278686905;
 p[6]=-0.13857109526572012;
 p[7]=9.9843695780195716e-6;
 p[8]=1.5056327351493116e-7;
 std::complex<long double>sum,ris;
 z -=1;
 sum=p[0];
 for(int i=1;i<(g+2);i++)
   sum += p[i]/(z+i);
 ris=0.5*gsl_sf_log(2*M_PI)+(z+0.5)*log(z+g+0.5)-(z+g+0.5)+log(sum);
 return ris;
  }

  std::complex<long double> CGamma(std::complex<long double> z){
  return( std::exp(LogGamma(z)));
  }

  //Gamma Derivatives of order i complex
  std::complex<long double> PolyGamma(int i, std::complex<long double> z){
  if (i == 0) {
    std::complex<long double> SUB = 0. ;
    std::complex<long double> ZZ = z;
    if(std::abs(std::imag(ZZ))<10.) { // if too close to the real axis...
    label1:
      if(std::real(ZZ)<10.) { // ...use recurrence relation to push real(z) large enough
	SUB = SUB - 1./ ZZ;
	ZZ = ZZ + 1.;
	goto label1;
      }
    }
    std::complex<long double> RZ = 1./ ZZ;
    std::complex<long double> DZ = RZ * RZ;
    // SUB + asympt expansion (Abramowitz, Stengun, 6.3.18)
    return (SUB + std::log(ZZ) - 0.5 * RZ - DZ/5040. * ( 420.+ DZ * ( - 42. + DZ * (20. - 21. * DZ) ) ));
  }
  else {
     int K1, K2;
     std::complex<long double> SUB = 0. , SUBM;
     std::complex<long double> ZZ = z;
     if(std::abs(std::imag(ZZ))<10.) { // if too close to the real axis...
      label2:
      SUBM = -1./ZZ;
      for(K1=1; K1<=i; K1++) {
	SUBM = - SUBM * K1 / ZZ;
      }
      if(std::real(ZZ)<10.) { // ...use recurrence relation to push real(z) large enough
	SUB = SUB + SUBM;
	ZZ = ZZ + 1.;
	goto label2;
      }
    }
    // Expansion coefficients for the first derivative
    long double A1 =  1.;
    long double A2 =  1./2.;
    long double A3 =  1./6.;
    long double A4 = -1./30.;
    long double A5 =  1./42.;
    long double A6 = -1./30.;
    long double A7 =  5./66.;
    // Expansion coefficients for the higher derivatives
    if(i>1) {
      for(K2=2;K2<=i;K2++){
	A1 = A1 * (1.*K2-1.);
	A2 = A2 *  1.*K2;
	A3 = A3 * (1.*K2+1.);
	A4 = A4 * (1.*K2+3.);
	A5 = A5 * (1.*K2+5.);
	A6 = A6 * (1.*K2+7.);
	A7 = A7 * (1.*K2+9.);
      }
    }
    std::complex<long double> RZ = 1./ ZZ;
    std::complex<long double> DZ = RZ * RZ;
    // SUB + asympt expansion (Abramowitz, Stengun, 6.4.11)
    return (SUB + std::pow(-1.,i+1) * std::pow(RZ,i) * ( A1 + RZ * (A2 + RZ * (A3 + DZ * (A4 + DZ * (A5 + DZ * (A6 + A7 * DZ ))))) ));
  }
  }
  
  
};


  
#endif
