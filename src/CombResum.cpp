#include "CombResum.h"

std::complex<long double> Tmatch(std::complex<long double> N, long double xp){
  return (std::pow(N,2)*std::pow(xp,4)/(1.+std::pow(N,2)*std::pow(xp,4)));
}


CombResum::CombResum(int ordres, int ordmatch, int channel, string PDFNAME, bool Wilson,
				long double mH, long double Nc, long double Nf, long double Mur, long double Muf)
:_Lumi(LHAPDF::mkPDF(PDFNAME),Muf,Nf),cc(Muf,Mur,0.,mH,0.){
   cc._sigma0Higgs=std::sqrt(2.)*Gf*cc._as*cc._as/(576.*M_PIl);
   cc._as=alpha_s_muR(Mur);
  _ordres=ordres;_ordmatch=ordmatch;_channel=channel;_Wilson=Wilson;_Nc=Nc;_Nf=Nf;
  _Fix=new FixedptResum(ordres, ordmatch,channel, Wilson,cc, Nc, Nf);
  _Joint=new Joint(ordres, ordmatch, channel, Wilson,cc, Nc, Nf);
  
  
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
  
  
  _MatchFun=Tmatch;
  
 
  
  
  
}

CombResum::~CombResum(){
  delete _Fix;
  delete _Joint;
}


void CombResum::SetMatchingFunction(std::complex<long double> (Func)(std::complex<long double>, long double)){
  _MatchFun=Func;
}

void CombResum::SetMUF(long double Muf){
  cc._MUF=Muf;
  _Lumi.Cheb_Lum(cc._MUF);
}

void CombResum::SetMUR(long double Mur){
  cc._MUR=Mur;
  cc._as=alpha_s_muR(cc._MUR);
  cc._sigma0Higgs=std::sqrt(2.)*Gf*cc._as*cc._as/(576.*M_PIl);
}

void CombResum::SetChannel(int channel){
  _channel=channel;
  _Fix->SetChannel(channel);
  _Joint->SetChannel(channel);
}

void CombResum::SetOrdRes(int ordres){
  _ordres=ordres;
  _Fix->SetOrdRes(ordres);
  _Joint->SetOrdRes(ordres);
}

void CombResum::SetOrdMatch(int ordmatch){
  _ordmatch=ordmatch;
  _Fix->SetOrdMatch(ordmatch);
  _Joint->SetOrdMatch(ordmatch);
}


long double CombResum::alpha_s_muR(long double mu)
{
  double X=1+alpha_Mz*beta_0*std::log(std::pow(mu,2)/std::pow(Mz,2));
  double ris=0.0;
  ris= alpha_Mz/X-beta_1/beta_0*std::pow(alpha_Mz/X,2.)*std::log(X)+std::pow(alpha_Mz/X,3)*(beta_2/beta_0*(1.-X)
  +std::pow(beta_1/beta_0,2)*(std::pow(std::log(X),2.)-std::log(X)-1.+X));
  return ris;
}

long double CombResum::ResummedCrossSection(long double CMS, long double xp){
  
}





//Core Functions

std::complex<long double> CombResum::FixptPartRes(std::complex<long double> N, long double xp){
  std::vector<std::complex<long double>> Fixchannel,Lumi;std::complex<long double> zero(0.,0.);
  Fixchannel=_Fix-> ComputeFixedptResum(N,xp,cc);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return (_MatchFun(N,xp)*std::inner_product(Fixchannel.begin(),Fixchannel.end(),Lumi.begin(),zero));
}

std::complex<long double> CombResum::JointPartRes(std::complex<long double> N, std::complex<long double> b, long double xp){
  std::vector<std::complex<long double>>Jointchannel,Lumi;std::complex<long double> zero(0.,0.);
  Jointchannel=_Joint->ComputeJointRes(N,b,cc);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return ((1.-_MatchFun(N,xp))*std::inner_product(Jointchannel.begin(),Jointchannel.end(),Lumi.begin(),zero));
}

std::complex<long double> CombResum::FixptPartMatch(std::complex<long double> N, long double xp){
  std::vector<std::complex<long double>> Fixchannel,Lumi;std::complex<long double>zero(0.,0.);
  Fixchannel=_Fix->ComputeMatching(N,xp,cc);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return (_MatchFun(N,xp)*std::inner_product(Fixchannel.begin(),Fixchannel.end(),Lumi.begin(),zero));
}

std::complex<long double> CombResum::JointPartMatch(std::complex<long double> N, std::complex<long double> b, long double xp){
  std::vector<std::complex<long double>>Jointchannel,Lumi;std::complex<long double> zero(0.,0.);
  Jointchannel=_Joint->ComputeMatching(N,b,cc);
  Lumi=_Lumi.Higgs_Lum_N(N);
  return ((1.-_MatchFun(N,xp))*std::inner_product(Jointchannel.begin(),Jointchannel.end(),Lumi.begin(),zero));
}




