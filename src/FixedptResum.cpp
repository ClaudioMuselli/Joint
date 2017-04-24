#include "FixedptResum.h"

FixedptResum::FixedptResum(const int ordres, const int ordmatch, int channel, bool Wilson){
  _ordres=ordres; _ordmatch=ordmatch; _channel=channel;_Wilson=Wilson;
  
  //Initialise constants
  _Ca=_Nc; 
  _Cf=(_Nc*_Nc-1.)/(2.*_Nc);
  
  const long double zeta2=gsl_sf_zeta_int(2);
  const long double zeta3=gsl_sf_zeta_int(3);
  const long double zeta4=gsl_sf_zeta_int(4);
  
  
  //REMEMBER!!! ALL in unity as
  beta_0=(11.*_Ca-2.*_Nf)/(12.*M_PIl);
  beta_1=((17.*_Ca*_Ca-5.*_Ca*_Nf-3.*_Cf*_Nf)*2./3.)/(16.*M_PIl*M_PIl);
  beta_2=((2857./54.*_Ca*_Ca*_Ca+(_Cf*_Cf-205./18.*_Cf*_Ca-1415./54.*_Ca*_Ca)*_Nf
  +(11./9.*_Cf+79./54.*_Ca)*_Nf*_Nf))/std::pow(4.*M_PIl,3);
  
  Ath1g=_Ca/M_PIl;
  Ath2g=(_Ca/2.*(_Ca*(67./18.-zeta2)-5./9.*_Nf))/std::pow(M_PIl,2);
  
  Ath1q=_Cf/M_PIl;
  Ath2q=(_Cf/2.*(_Ca*(67./18.-zeta2)-5./9.*_Nf))/std::pow(M_PIl,2);
  
  Bth1g=-beta_0;
  Bth1q=-3./4.*_Cf/M_PIl;  
  
  if (Wilson){
    W1=0;
    W2=0;
  }
  else{
    W1=0;
    W2=0;
  }
  if (_ordres>2){
    cout << "Resummation is implemented only up to NNLL! Set resummation order to NNLL" << endl;
    _ordres=2;
  }
  if (_ordmatch>2){
    cout << "Matching is implemented only up to NNLO! Set matching order to NNLO" << endl;
    _ordmatch=2;
  }
  
}


FixedptResum::~FixedptResum(){
}










//Sudakov Exponent

std::complex<long double> FixedptResum::Sudakov_th_gggg(std::complex<long double> N, long double xp){
  
}
std::complex<long double> FixedptResum::Sudakov_th_gggq(std::complex<long double> N, long double xp){
}
std::complex<long double> FixedptResum::Sudakov_th_ggqq(std::complex<long double> N, long double xp){
}












//Higgs cross section
//LO cross section (WITHOUT sigma_0)

std::complex<long double> FixedptResum::LOgggH(std::complex<long double> N, long double xp){
  std::complex<long double> CLOgggH;
  std::complex<long double> half(0.5,0.);
  std::complex<long double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);
  CLOgggH=2.*_as*_Ca/std::sqrt(M_PIl)*1./xp*std::exp(LogGamma(N)-LogGamma(N+0.5))*(Hyp2F1(half,N,N+0.5,xprad)
  -2.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*N/(N+0.5)*Hyp2F1(half,N+1.,N+1.5,xprad)
  +((1.+xp)*(3.+xp))/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),4.))*N*(N+1.)/((N+0.5)*(N+1.5))*Hyp2F1(half,N+2.,N+2.5,xprad)
  -2.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),6.))*N*(N+1.)*(N+2.)/((N+0.5)*(N+1.5)*(N+2.5))*Hyp2F1(half,N+3.,N+3.5,xprad)
  +1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),8))*N*(N+1.)*(N+2.)*(N+3.)/((N+0.5)*(N+1.5)*(N+2.5)*(N+3.5))
  *Hyp2F1(half,N+4.,N+4.5,xprad));
  return CLOgggH;
}
std::complex<long double> FixedptResum::LOgqqH(std::complex<long double> N, long double xp){
  std::complex<long double> CLOgqqH;
  std::complex<long double> half(0.5,0.);
  std::complex<long double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);
  CLOgqqH=_as*_Cf/std::sqrt(M_PIl)*1./xp*std::exp(LogGamma(N)-LogGamma(N+0.5))*(Hyp2F1(half,N,N+0.5,xprad)
  -(4.+3.*xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*N/(N+0.5)*Hyp2F1(half,N+1.,N+1.5,xprad)
  +3.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),4.))*N*(N+1.)/((N+0.5)*(N+1.5))*Hyp2F1(half,N+2.,N+2.5,xprad)
  -1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),6.))*N*(N+1.)*(N+2.)/((N+0.5)*(N+1.5)*(N+2.5))*Hyp2F1(half,N+3.,N+3.5,xprad));
  return CLOgqqH;
}
std::complex<long double> FixedptResum::LOqqgH(std::complex<long double> N, long double xp){
  std::complex<long double> CLOqqgH;
  std::complex<long double> half(0.5,0.);
  std::complex<long double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);
  CLOqqgH=2.*_as*_Cf*_Cf/M_PIl*1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*(Hyp2F1(half,N,N+0.5,xprad)
  -2.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*N/(N+0.5)*Hyp2F1(half,N+1.,N+1.5,xprad)
  +1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),4.))*N*(N+1.)/((N+0.5)*(N+1.5))*Hyp2F1(half,N+2.,N+2.5,xprad));
  return CLOqqgH;
  
}


































































