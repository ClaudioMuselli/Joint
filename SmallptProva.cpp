#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <complex>
#include <vector>
#include <fstream>
#include <sstream>
#include <LHAPDF/LHAPDF.h>
#include "Luminosity.h"
#include "complex_def.h"
#include "integration.h"

using namespace std;

Luminosity *Lumi;

const long double Ca=3.;
const long double Cf=4./3.;
const long double Nf=5.;
const long double beta_0=(11.*Ca-2.*Nf)/(12.*M_PIl);
const long double Apt1g=Ca/M_PIl;
const long double EulerGamma=0.577215664901532860606512090082;

long double as;
long double sigma0;
const long double Gf = 0.00001166364;// Fermi constant in GeV^-2
const long double GeVtopb = 389379304.;// GeV^-2 to pb conversion factor == (hc)^2 

std::complex<long double> SudsmallptLL(std::complex<long double> y, std::complex<long double> N){
  return(std::exp(1./as*Apt1g/beta_0/beta_0*((std::log(1.+y)-y)-std::log(1.+y)*std::log(N*std::exp(EulerGamma))*2.*as*beta_0)));
}

//std::complex<long double> SudsmallptLL(std::complex<long double> y){
//  return(std::exp(1./as*Apt1g/beta_0/beta_0*((std::log(1.+y)-y))));
//}

std::complex<long double> FunM(std::complex<long double> xi, long double xp){
  return(std::pow(2.,-xi)*std::exp(2.*EulerGamma*xi)*std::exp(-(1.-xi)*std::log(xp))*std::exp(LogGamma(1.-xi)-LogGamma(xi)));
}

std::complex<long double> FunJ(std::complex<long double> xi, long double xp, std::complex<long double> N){
  return(2.*std::pow(xp,-0.5-0.5*xi)*std::pow(N,1.+xi)*CBesselK(1.+xi,2.*N*std::sqrt(xp))
  /std::exp(LogGamma(-xi))*std::exp(2.*EulerGamma*xi));
}

struct Borelptparam{
  long double x;
  long double xp;
  long double Ccutoff;
  long double asbar;
};


int Borelprovapt(int* ndim, double * x, int* ncomp, double* y, void *p){
  Borelptparam par=*(Borelptparam *) p;
  long double N0=3.; long double slope=1.0;
  std::complex<long double> N=N0+(slope+II)*std::log(x[0]);
  std::complex<long double> ris;
  long double  omega=par.Ccutoff*x[2];
  long double theta=2.*M_PIl*x[1];
  long double r=0.5+omega/2.; long double xc=omega/2.;
  std::complex<long double> xi=xc+r*std::exp(II*theta);
  //ris=par.Ccutoff/par.asbar*FunM(xi,par.xp)*SudsmallptLL(omega/xi,N)*std::exp(-omega/par.asbar)*(r*std::exp(II*theta))/xi;
  ris=par.Ccutoff/par.asbar*FunJ(xi,par.xp,N)*SudsmallptLL(-omega/xi,N)*Lumi->Lum_gg_N(N)*std::exp(-omega/par.asbar)*(r*std::exp(II*theta))/xi;
  //cout << ris << endl;
  y[0]=std::imag(std::exp(-N*std::log(par.x*std::pow(std::sqrt(1.+par.xp)-std::sqrt(par.xp),2)))*ris/x[0]*(slope+II))/M_PIl;
  return 0;
}





int main(){
  Lumi=new Luminosity(LHAPDF::mkPDF("PDF4LHC15_nnlo_100"),125.09,5.);
  as=Lumi->get_alphaS(125.09);
  sigma0=std::sqrt(2.)*Gf*as*as/(576.*M_PIl);
  cout.precision(10);
  long double ris=0.,err=0.;
  long double x=std::pow(125.09/13000.,2);
  Borelptparam pp;
  pp.asbar=as*beta_0;
  pp.Ccutoff=2.;
  double *res=NULL,*error=NULL,*prob=NULL;
  res=new double[1];
  error=new double[1];
  prob=new double[1];
  int fail;
  double prec=1e-10;
  int NUM=20.;
  std::complex<long double> NTest(2.670000264,-0.3299997364);
  std::complex<long double> lchiTest(0.6776217648,0.3673330522);

  cout << "ris= " << SudsmallptLL(-lchiTest,NTest)<< endl;
  
  long double ptstart=1.; long double ptend=20.;
  long double rate=(ptend-ptstart)/((long double) NUM);
  for (int i=0;i<NUM;i++){
    long double pt=ptstart+i*rate;
    long double xp=std::pow(pt/125.09,2);
    long double tauprime=x*std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2);
    pp.xp=xp;
    pp.x=tauprime;
    integration::CUBA(3,Borelprovapt,3,1,prec,&res,&fail,&error,&prob,&pp);
    ris=(static_cast<long double> (res[0]))*tauprime*sigma0*GeVtopb*2.*pt/(125.09*125.09);
    err=(static_cast<long double> (error[0]))*tauprime*sigma0*GeVtopb*2.*pt/(125.09*125.09);
    cout << "pt= " << pt << " Ris= " << ris << " +- " << err << endl;
  }
  
}