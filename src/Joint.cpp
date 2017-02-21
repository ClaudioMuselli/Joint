#include "Joint.h"

Joint::Joint(const int ordres, const int ordmatch, string channel, bool Wilson): AP(_Nc,_Nf){
  _ordres=ordres; _ordmatch=ordmatch; _channel=channel;_Wilson=Wilson;
  
  //Initialise constants
  _Ca=_Nc; 
  _Cf=(_Nc*_Nc-1.)/(2.*_Nc);
  
  const long double zeta2=gsl_sf_zeta_int(2);
  const long double zeta3=gsl_sf_zeta_int(3);
  const long double zeta4=gsl_sf_zeta_int(4);
  
  
  //REMEMBER!!! ALL in unity as/Pi
  beta_0=(11.*_Ca-2.*_Nf)/(12.);
  beta_1=
  beta_2=
  
  Apt1g=_Ca;
  Apt2g=(_Ca/2.*(_Ca*(67./18.-zeta2)-5./9.*_Nf));
  Apt3g=_Ca/4.*(_Ca*_Ca*(15503./648.-67./9.*zeta2-11.*zeta3+11./2.*zeta4)+_Cf*_Nf*(-55./24.+2.*zeta3)+_Ca*_Nf*(-2051./324.+10./9.*zeta2)-25./81.*_Nf*_Nf);
  
  Bpt1g=-2.*beta_0;
  Bpt2g=-2.*(-8./3.*_Ca*_Nf+(32./3.+12.*zeta3)*_Ca*_Ca-2.*_Cf*_Nf)/(std::pow(4.,2))+beta_0*_Ca*zeta2;
  
  Dpt2g=_Ca*_Ca*(-101./27.+7./2.*zeta3)+14./27.*_Ca*_Nf;
  
  Hpt1g=3.*_Ca*zeta2;
  Hpt2g=_Ca*_Ca*(93./16.+67./12.*zeta2-55./18.*zeta3+65./8.*zeta4)+_Ca*_Nf*(-5./3.-5./6.*zeta2-4./9.*zeta3);
  
  if (Wilson){
    W1=0;
    W2=0;
  }
  else{
    W1=0;
    W2=0;
  } 
}


std::complex<long double> Joint::Sudakov_g(std::complex<long double> N, std::complex<long double> b){
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> bbar= b/2.*std::exp(EulerGamma);
  
  const std::complex<long double> chi= Nbar*Nbar+bbar*bbar;
  const std::complex<long double> LC= std::log(LC);
  const std::complex<long double> lchi= _as*beta_0/M_PIl*LC;
  std::complex<long double> Sud_g1_g(0.,0.),Sud_g2_g(0.,0.),Sud_g3_g(0.,0.);
  if (_ordres>=0) {
    Sud_g1_g=Apt1g/beta_0/beta_0*(lchi+std::log(1.-lchi));
  }
  if (_ordres>=1) {
    Sud_g2_g=Apt1g*beta_1/std::pow(beta_0,3)*((lchi+std::log(1.-lchi))/(1.-lchi)+0.5*std::pow(std::log(1.-lchi),2))
					    -Apt2g/beta_0/beta_0*(std::log(1.-lchi)+lchi/(1.-lchi))+Bpt1g/beta_0*std::log(1.-lchi);
  }
  if (_ordres>=2) {
    gsl_sf_result redilog,imdilog;
    long double r=std::abs(Nbar*Nbar/chi);
    long double theta=std::arg(Nbar*Nbar/chi);
    gsl_sf_complex_dilog_e(r,theta,&redilog,&imdilog);
    std::complex<long double> risdilog(redilog.val,imdilog.val);
    const std::complex<long double> lN=_as*beta_0/M_PIl*std::log(Nbar*Nbar);
    Sud_g3_g=Apt1g*beta_1*beta_1/(2.*std::pow(beta_0,4))*((lchi+std::log(1.-lchi))/(std::pow(1.-lchi,2))*(lchi+(1.-2.*lchi)*std::log(1.-lchi)))
					    +Apt1g*beta_2/std::pow(beta_0,3)*(((2.-3.*lchi)*lchi)/(2.*std::pow(1.-lchi,2))+std::log(1.-lchi))
					    -Apt2g*beta_1/std::pow(beta_0,3)*(((2.-3.*lchi)*lchi)/(2.*std::pow(1.-lchi,2))+((1.-2.*lchi)*std::log(1.-lchi))/(std::pow(1.-lchi,2)))
					    +Bpt1g*beta_1/beta_0*(lchi+std::log(1.-lchi))/(1.-lchi)-Apt3g*lchi*lchi/(2.*beta_0*beta_0*std::pow(1.-lchi,2))
					    -Bpt2g/beta_0*lchi/(1.-lchi)+Apt1g*lN/(1.-lN)*risdilog;
  }
  return std::exp(M_PIl/_as*Sud_g1_g+Sud_g2_g+_as/M_PIl*Sud_g3_g); 
}
























































































