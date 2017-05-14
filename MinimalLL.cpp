#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include "Luminosity.h"
#include "complex_def.h"
#include "integration.h"
#include "CombResum.h"
#include <vector>
#include <complex>


using namespace std;

const long double CaLL=3.;
const long double CfLL=4./3.;
const long double NfLL=5.;
const long double beta_0LL=(11.*CaLL-2.*NfLL)/(12.*M_PIl);
const long double Apt1gLL=CaLL/M_PIl;
const long double EulerGammaLL=0.577215664901532860606512090082;
const long double b0LL=2.*std::exp(-EulerGammaLL);

long double asLL;
long double sigma0LL;
const long double GfLL = 0.00001166364;// Fermi constant in GeV^-2
const long double GeVtopbLL = 389379304.;// GeV^-2 to pb conversion factor == (hc)^2 


std::complex<long double> S(std::complex<long double> b, void * p){
  std::complex<long double> y=as*beta_0LL*std::log(1.+b*b/(b0LL*b0LL));
  return(std::exp(1./as*Apt1gLL/beta_0LL/beta_0LL*((std::log(1.-y)+y))));
}

struct MinimalParam {
  long double xp;
  long double phi;
  long double N0=3;
  long double slope=1.5;
  long double v=3.;
  Luminosity *Lum;
  long double x;
  
};

int Minimalintpt(int* ndim, double * x, int* ncomp, double* y, void *p){
  MinimalParam par= *(MinimalParam *)p;
  std::complex<long double> b=-std::exp(-II*par.phi)*std::log(x[0]);
  std::complex<long double> N=par.N0+(par.slope+II)*std::log(x[2]);
  std::complex<long double> theta1=-II*par.v*M_PIl+x[1]*M_PIl*(-1.+2.*II*par.v);
  std::complex<long double> theta2=-II*par.v*M_PIl+x[4]*M_PIl*(1.+2.*II*par.v);
  std::complex<long double> h1=1./M_PIl*(1.-2.*II*par.v)*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta1));
  std::complex<long double> h2=1./M_PIl*(1.+2.*II*par.v)*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta2));
  y[0]=std::imag(std::exp(-N*std::log(par.x))*b/4.*h1*S(b,NULL)*par.Lum->Lum_gg_N(N)/x[0]/x[2]*(II)*std::exp(-II*par.phi)*(par.slope+II))/M_PIl;
  b=-std::exp(II*par.phi)*std::log(x[3]);
  y[0]+=std::imag(std::exp(-N*std::log(par.x))*b/4.*h2*S(b,NULL)*par.Lum->Lum_gg_N(N)/x[0]/x[2]*(II)*std::exp(II*par.phi)*(par.slope+II))/M_PIl;
}

  


int main(){
  long double CMS=13000.;
  long double x= std::pow(125.09/CMS,2);
  MinimalParam Min;
  Luminosity *Lumi;
  Lumi=new Luminosity(LHAPDF::mkPDF("PDF4LHC15_nnlo_100"),125.09,5.);
  asLL=Lumi->get_alphaS(125.09);
  Lumi->TestLum("TestLum.dat");
  ofstream OUT("TestLLCompleto.dat");
  ifstream IN("graph/HqtLL_13.dat");
  CombResum Final(0,0,0,"PDF4LHC15_nnlo_100",true);
  std::vector<long double> hqt;
  long double h;
  while(IN >> h){
    IN >> h;
    hqt.push_back(h);
  }
  int NUM=50.;
  long double ptstart=5.; long double ptend=105.;
  long double rate=2.;
  
  sigma0LL=asLL*asLL*GfLL*std::sqrt(2.)/(576.*M_PIl);
  OUT << "pt \t LL_HqT \t LL_pt_MP \t LL_pt_BP \t LL_Joint" << endl;
  for (int i=0; i<NUM;i++){
    cout << i << endl;
    long double pt=ptstart+i*rate;
    long double xp=std::pow(pt/125.09,2);
    std::complex<long double> err;
    long double error;
    std::complex<long double> ris=integration::InverseBessel_path(3,S,xp,0.,10.,1.5,NULL,&err);
    long double LLptminimal=std::real(ris)*Lumi->Lum_gg_x(x)*x*sigma0LL*GeVtopbLL*2.*std::sqrt(xp)/125.09;
    long double LLptBorel=Final.ResummedCrossSectionpt(CMS,xp,&error)*2.*pt/std::pow(125.09,2);
    long double LLJoint=Final.ResummedCrossSection(CMS,xp,0,&error)*2.*pt/std::pow(125.09,2);
    OUT << pt << "\t" << hqt[i] << "\t" << LLptminimal << "\t" << LLptBorel << "\t" << LLJoint << endl;
    //cout << std::real(ris)*Lumi->Lum_gg_x(x)*x*sigma0*GeVtopb*2.*std::sqrt(xp)/125.09 << " +- " << std::real(err)*Lumi->Lum_gg_x(x)*x*sigma0*GeVtopb*2.*std::sqrt(xp)/125.09<< endl;
  }
  IN.close();
  OUT.close();
  
  
}