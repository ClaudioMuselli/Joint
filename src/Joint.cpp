#include "Joint.h"

Joint::Joint(const int ordres, const int ordmatch, int channel, bool Wilson): AP(_Nf,_Nc){
  _ordres=ordres; _ordmatch=ordmatch; _channel=channel;_Wilson=Wilson;
  
  //Initialise constants
  _Ca=_Nc; 
  _Cf=(_Nc*_Nc-1.)/(2.*_Nc);
  
  const long double zeta2=gsl_sf_zeta_int(2);
  const long double zeta3=gsl_sf_zeta_int(3);
  const long double zeta4=gsl_sf_zeta_int(4);
  
  
  //REMEMBER!!! ALL in unity as/Pi
  beta_0=(11.*_Ca-2.*_Nf)/(12.*M_PIl);
  beta_1=((17.*_Ca*_Ca-5.*_Ca*_Nf-3.*_Cf*_Nf)*2./3.)/(16.*M_PIl*M_PIl);
  beta_2=((2857./54.*_Ca*_Ca*_Ca+(_Cf*_Cf-205./18.*_Cf*_Ca-1415./54.*_Ca*_Ca)*_Nf
  +(11./9.*_Cf+79./54.*_Ca)*_Nf*_Nf))/std::pow(4.*M_PIl,3);
  
  Apt1g=_Ca/M_PIl;
  Apt2g=(_Ca/2.*(_Ca*(67./18.-zeta2)-5./9.*_Nf))/std::pow(M_PIl,2);
  Apt3g=(_Ca/4.*(_Ca*_Ca*(15503./648.-67./9.*zeta2-11.*zeta3+11./2.*zeta4)
  +_Cf*_Nf*(-55./24.+2.*zeta3)+_Ca*_Nf*(-2051./324.+10./9.*zeta2)-25./81.*_Nf*_Nf))/std::pow(M_PIl,3);
  
  Bpt1g=-2.*beta_0;
  Bpt2g=(-2.*(-8./3.*_Ca*_Nf+(32./3.+12.*zeta3)*_Ca*_Ca-2.*_Cf*_Nf)/(std::pow(4.,2)))/std::pow(M_PIl,2)+beta_0*_Ca*zeta2/M_PIl;
  
  Dpt2g=(_Ca*_Ca*(-101./27.+7./2.*zeta3)+14./27.*_Ca*_Nf)/std::pow(M_PIl,2);
  
  Hpt1g=3.*_Ca*zeta2/M_PIl;
  Hpt2g=(_Ca*_Ca*(93./16.+67./12.*zeta2-55./18.*zeta3+65./8.*zeta4)+_Ca*_Nf*(-5./3.-5./6.*zeta2-4./9.*zeta3))/std::pow(M_PIl,2);
  
  Apt1q=_Cf/M_PIl;
  
  if (Wilson){
    W1=0;
    W2=0;
  }
  else{
    W1=0;
    W2=0;
  }
  
  ULL=new std::complex<long double>*[2];
  ULL[0]=new std::complex<long double>[2];
  ULL[1]=new std::complex<long double>[2];
  
  UNLL=new std::complex<long double>*[2];
  UNLL[0]=new std::complex<long double>[2];
  UNLL[1]=new std::complex<long double>[2];
  
  UNNLL=new std::complex<long double>*[2];
  UNNLL[0]=new std::complex<long double>[2];
  UNNLL[1]=new std::complex<long double>[2];
  
  V1=new std::complex<long double>*[2];
  V1[0]=new std::complex<long double>[2];
  V1[1]=new std::complex<long double>[2];
  
}

Joint::~Joint(){
}

//Sudakov Higgs
std::complex<long double> Joint::Sudakov_g(std::complex<long double> N, std::complex<long double> b){
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> bbar= b/2.*std::exp(EulerGamma);
  
  const std::complex<long double> chi= Nbar*Nbar+bbar*bbar;
  const std::complex<long double> LC= std::log(chi);
  const std::complex<long double> lchi= _as*beta_0*LC;
  std::complex<long double> Sud_g1_g(0.,0.),Sud_g2_g(0.,0.),Sud_g3_g(0.,0.);
  if (_ordres>=0) {
    Sud_g1_g=Apt1g/beta_0/beta_0*(lchi+std::log(1.-lchi));
  }
  if (_ordres>=1) {
    Sud_g2_g=Apt1g*beta_1/std::pow(beta_0,3)*((lchi+std::log(1.-lchi))/(1.-lchi)+0.5*std::pow(std::log(1.-lchi),2))
					    -Apt2g/beta_0/beta_0*(std::log(1.-lchi)+lchi/(1.-lchi))+Bpt1g/beta_0*std::log(1.-lchi)
					    +Apt1g/beta_0*(std::log(1.-lchi)+lchi/(1.-lchi))*_LR;
  }
  if (_ordres>=2) {
    gsl_sf_result redilog,imdilog;
    long double r=std::abs(Nbar*Nbar/chi);
    long double theta=std::arg(Nbar*Nbar/chi);
    gsl_sf_complex_dilog_e(r,theta,&redilog,&imdilog);
    std::complex<long double> risdilog(redilog.val,imdilog.val);
    const std::complex<long double> lN=_as*beta_0*std::log(Nbar*Nbar);
    Sud_g3_g=Apt1g*beta_1*beta_1/(2.*std::pow(beta_0,4))*((lchi+std::log(1.-lchi))/(std::pow(1.-lchi,2))*(lchi+(1.-2.*lchi)*std::log(1.-lchi)))
					    +Apt1g*beta_2/std::pow(beta_0,3)*(((2.-3.*lchi)*lchi)/(2.*std::pow(1.-lchi,2))+std::log(1.-lchi))
					    -Apt2g*beta_1/std::pow(beta_0,3)*(((2.-3.*lchi)*lchi)/(2.*std::pow(1.-lchi,2))+(1.-2.*lchi)*std::log(1.-lchi)/std::pow(1.-lchi,2))
					    +Bpt1g*beta_1/beta_0*(lchi+std::log(1.-lchi))/(1.-lchi)-Apt3g*lchi*lchi/(2.*beta_0*beta_0*std::pow(1.-lchi,2))
					    -Bpt2g/beta_0*lchi/(1.-lchi)+Apt1g*lN/(1.-lN)*risdilog+Bpt1g*lchi/(1.-lchi)*_LR+Apt2g/beta_0*lchi*lchi/std::pow(1.-lchi,2)*_LR
					    +Apt1g*beta_1/std::pow(beta_0,2)*((lchi*(1.-lchi)+(1.-2.*lchi)*std::log(1.-lchi))/std::pow(1.-lchi,2))*_LR-Apt1g/2.*std::pow(lchi/(1.-lchi)*_LR,2);    
  }
  //cout << "Sud= " << 1./_as*Sud_g1_g+Sud_g2_g+_as*Sud_g3_g << endl;
  return std::exp(1./_as*Sud_g1_g+Sud_g2_g+_as*Sud_g3_g); 
}

//Hard part Higgs
std::complex<long double> Joint::Hgggq1(std::complex<long double> N) {
  return (_Cf/2./(N+1.))/M_PIl;
}

std::complex<long double> Joint::Hggggreg2(std::complex<long double> N) {
  MellinFunc Func;
  std::complex<long double> hggggreg2;
  std::complex<long double> Nm1=N-1.;
  std::complex<long double> N1=N+1.;
  std::complex<long double> N2=N+2.;
  std::complex<long double> N3=N+3.;
  std::complex<long double> N4=N+4.;
  
  hggggreg2=1./12.*_Ca*(_Ca-_Nf)*Func.Logminus(N1)+1./216.*(18.*_Cf*_Nf*(Func.Logz3(N)+Func.Logz3(N1))
	   -36.*_Ca*_Ca*(Func.Logz3minusplus(N)+2.*Func.Logz3minusplus(N1)-Func.Logz3minusplus(N2)
	   -2.*Func.Logz3minusplus(N3)+Func.Logz3minusplus(N4)))+_Ca*_Ca*gsl_sf_zeta_int(2)*
	   (Func.Logplusplus(Nm1)+2.*Func.Logplusplus(N)+3.*Func.Logplusplus(N1)+2.*Func.Logplusplus(N2)+Func.Logplusplus(N3))
	   -_Ca*_Ca/3.*(Func.Logplus3plus(Nm1)+2.*Func.Logplus3plus(N)+3.*Func.Logplus3plus(N1)+2.*Func.Logplus3plus(N2)
	   +Func.Logplus3plus(N3))+1./24.*(2.*_Ca*_Nf*(Func.Logz2(N)+Func.Logz2(N1))+3.*_Cf*_Nf*
	   (3.*Func.Logz2(N)+Func.Logz2(N1))+_Ca*_Ca*(25.*Func.Logz2(N)-11.*Func.Logz2(N1)+44.*Func.Logz2(N2)))
	   +_Ca*_Ca/2.*(Func.Logz2Logminusminus(Nm1)-2.*Func.Logz2Logminusminus(N)+3.*Func.Logz2Logminusminus(N1)
	   -2.*Func.Logz2Logminusminus(N2)+Func.Logz2Logminusminus(N3))-_Ca*_Ca/2.*(Func.Logz2Logplusplus(Nm1)
	   +2.*Func.Logz2Logplusplus(N)+3.*Func.Logz2Logplusplus(N1)+2.*Func.Logz2Logplusplus(N2)+Func.Logz2Logplusplus(N3))
	   +_Ca*_Ca/2.*(Func.LogzLogminus2minus(Nm1)-2.*Func.LogzLogminus2minus(N)+3.*Func.LogzLogminus2minus(N1)
	   -2.*Func.LogzLogminus2minus(N2)+Func.LogzLogminus2minus(N3))
	   +_Ca*_Ca*(Func.LogzLogplus2plus(Nm1)+2.*Func.LogzLogplus2plus(N)+3.*Func.LogzLogplus2plus(N1)
	   +2.*Func.LogzLogplus2plus(N2)+Func.LogzLogplus2plus(N3))+1./216.*((78.*_Ca*_Nf+324.*_Cf*_Nf-2319.*_Ca*_Ca)*Func.Logz(N)
	   -216.*_Ca*_Ca*Func.Logz(Nm1)+(60.*_Ca*_Nf+324.*_Cf*_Nf-447.*_Ca*_Ca)*Func.Logz(N1)-1608.*_Ca*_Ca*Func.Logz(N2)
	   +216.*_Ca*_Ca*(Func.Li2mzLogzplus(Nm1)+2.*Func.Li2mzLogzplus(N)+3.*Func.Li2mzLogzplus(N1)+2.*Func.Li2mzLogzplus(N2)
	   + Func.Li2mzLogzplus(N3)))+_Ca*_Ca*(3.*Func.Li2zLogzplusminus(Nm1)-Func.Li2zLogzplusminus(N)+3.*Func.Li2zLogzplusminus(N1)
	   +Func.Li2zLogzplusminus(N2)-3.*Func.Li2zLogzplusminus(N3)+Func.Li2zLogzplusminus(N4))-_Ca*_Ca*(5.*Func.Li3zregminusplus(Nm1)
	   -Func.Li3zregminusplus(N)+5.*Func.Li3zregminusplus(N1)+Func.Li3zregminusplus(N2)-5.*Func.Li3zregminusplus(N3)+Func.Li3zregminusplus(N4))
	   +1./216.*((_Ca*_Ca*(-8272.+4.*N*(-2145.+N*(-1386.+211.*N)))+108.*_Cf*(-16.+7.*N*(1.+N))*_Nf+2.*_Ca*(332.+N*(327.+(141.-74.*N)*N))*_Nf)
	   /((Nm1)*N*N1*N2)+72.*_Ca*_Ca*(11.*Func.Li2minus(Nm1)-12.*Func.Li2minus(N)+12.*Func.Li2minus(N1)-11.*Func.Li2minus(N2))
	   -216.*_Ca*_Ca*(Func.Li3mzplus(Nm1)+2.*Func.Li3mzplus(N)+3.*Func.Li3mzplus(N1)+2.*Func.Li3mzplus(N2)+Func.Li3mzplus(N3))
	   +432.*_Ca*_Ca*(Func.Li3zoverplusplus(Nm1)+2.*Func.Li3zoverplusplus(N)+3.*Func.Li3zoverplusplus(N1)+2.*Func.Li3zoverplusplus(N2)+Func.Li3zoverplusplus(N3))
	   -108.*_Ca*_Ca*gsl_sf_zeta_int(3)*(-2.*Func.plus(Nm1)+17.*Func.plus(N)+22.*Func.plus(N1)+10.*Func.plus(N2)+12.*Func.plus(N3)));
	   return hggggreg2/std::pow(M_PIl,2);
}

std::complex<long double> Joint::Hgggq2(std::complex<long double> N){
  MellinFunc Func;
  std::complex<long double> hgggq2;
  std::complex<long double> Nm1=N-1.;
  std::complex<long double> N1=N+1.;
  std::complex<long double> N2=N+2.;
  std::complex<long double> N3=N+3.;
  std::complex<long double> N4=N+4.;
  hgggq2=_Cf/72.*(_Ca*(-152.*Func.Logminus(Nm1)+152.*Func.Logminus(N)-43.*Func.Logminus(N1))
	+4.*_Nf*(5.*Func.Logminus(Nm1)-5.*Func.Logminus(N)+Func.Logminus(N1))
	+9.*_Cf*(16.*Func.Logminus(Nm1)-16.*Func.Logminus(N)+5.*Func.Logminus(N1)))
	+_Cf/48.*(_Ca*(-22.*Func.Logminus2(Nm1)+22.*Func.Logminus2(N)-5.*Func.Logminus2(N1))
	+3.*_Cf*(6.*Func.Logminus2(Nm1)-6.*Func.Logminus2(N)+Func.Logminus2(N1))
	+2.*_Nf*(2.*Func.Logminus2(Nm1)-2.*Func.Logminus2(N)+Func.Logminus2(N1)))
	-(_Ca-_Cf)*(_Cf*(2.*Func.Logminus3(Nm1)-2.*Func.Logminus3(N)+Func.Logminus3(N1)))/24.
	-1./48.*_Cf*(_Cf*(-2.*Func.Logz3(N)+Func.Logz3(N1))+2.*_Ca*(2.*Func.Logz3(N)+Func.Logz3(N1)))
	-_Ca*_Cf*gsl_sf_zeta_int(2)/4.*(2.*Func.Logplus(Nm1)+2.*Func.Logplus(N)+Func.Logplus(N1))
	+_Ca*_Cf/12.*(2.*Func.Logplus3(Nm1)+2.*Func.Logplus3(N)+Func.Logplus3(N1))
	-1./96.*_Cf*(3.*_Cf*(4.*Func.Logz2(N)+3.*Func.Logz2(N1))-2.*_Ca*(36.*Func.Logz2(N)+9.*Func.Logz2(N1)
	+8.*Func.Logz2(N2)))+_Ca*_Cf/8.*(2.*Func.Logz2Logminus(Nm1)-2.*Func.Logz2Logminus(N)
	+Func.Logz2Logminus(N1))-_Ca*_Cf/8.*(2.*Func.Logz2Logplus(Nm1)+2.*Func.Logz2Logplus(N)
	+Func.Logz2Logplus(N1))+_Ca*_Cf/12.*(22.*Func.Li2z(Nm1)+24.*Func.Li2z(N)-9.*Func.Li2z(N1)+4.*Func.Li2z(N2))
	+1./8.*_Ca*_Cf*Func.Li2z2(N1)+_Cf/144.*(45.*_Cf*(-3.*Func.Logz(N)+Func.Logz(N1))-2.*_Ca*(72.*Func.Logz(Nm1)
	+321.*Func.Logz(N)-6.*Func.Logz(N1)+88.*Func.Logz(N2)))+_Ca*_Cf/12*(-22.*Func.LogzLogminus(Nm1)
	+24.*Func.LogzLogminus(N)-9.*Func.LogzLogminus(N1)+4.*Func.LogzLogminus(N2))
	+_Ca*_Cf/8.*(2.*Func.LogzLogminus2(Nm1)-2.*Func.LogzLogminus2(N)+Func.LogzLogminus2(N1))+1./4.*_Ca*_Cf*Func.LogzLogplus(N1)
	+_Ca*_Cf/2.*(2.*Func.Li2zLogz(Nm1)-2.*Func.Li2zLogz(N)+Func.Li2zLogz(N1))+_Ca*_Cf/8.*(2.*Func.Li2z2Logz(Nm1)
	+2.*Func.Li2z2Logz(N)+Func.Li2z2Logz(N1))-_Ca*_Cf/2.*(2.*Func.Li3z(Nm1)-4.*Func.Li3z(N)+Func.Li3z(N1))
	-3.*_Ca*_Cf/16.*(2.*Func.Li3z2(Nm1)+2.*Func.Li3z2(N)+Func.Li3z2(N1))-_Ca*_Cf/2.*(2.*Func.Li3overplus(Nm1)
	+2.*Func.Li3overplus(N)+Func.Li3overplus(N1))+(_Cf*(-N2*(27.*_Cf*Nm1*(-10.+3.*N)-4.*(56.+N*(43.+13.*N))*_Nf)
	+4.*_Ca*(-2014+72.*6.*gsl_sf_zeta_int(2)+540.*gsl_sf_zeta_int(3)+N*(-2905.+39.*6.*gsl_sf_zeta_int(2)+702.*gsl_sf_zeta_int(3)
	+N*(-1137.+20.*N+63.*6.*gsl_sf_zeta_int(2)+24.*6.*gsl_sf_zeta_int(2)*N+54.*(18.+7.*N)*gsl_sf_zeta_int(3))))))/(432*Nm1*N*N1*N2);
  return hgggq2/std::pow(M_PIl,2);
  
}

std::complex<long double> Joint::Hggqq2(std::complex<long double> N){
  return (_Cf*_Cf*(4.+N*(1.+N)*(8.+(-3.+N)*N))/(4.*N*N*std::pow(1.-N*N,2)))/std::pow(M_PIl,2);
}


void Joint::ComputeEvolution(std::complex<long double> N, std::complex<long double> b){
  AP.ComputeGamma(N,1);
  const std::complex<long double> rad=std::sqrt(std::pow(AP.gg0-AP.SS0,2)+4.*AP.gS0*AP.Sg0);
  std::complex<long double> lnEPLL(0.,0.),lnEMLL(0.,0.),lnEPNLL(0.,0.),lnEMNLL(0.,0.),lnEPNNLL(0.,0.),lnEMNNLL(0.,0.);
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> bbar= b/2.*std::exp(EulerGamma);
  
  const std::complex<long double> chi= Nbar*Nbar+bbar*bbar;
  const std::complex<long double> LC= std::log(chi);
  const std::complex<long double> lchi= _as*beta_0*LC;
  const std::complex<long double> lN=_as*beta_0*std::log(Nbar*Nbar);
  if (_ordres >=0){
    lnEPLL=-Apt1g/2./beta_0/beta_0*lN*std::log(1.-lchi);
    lnEMLL=-Apt1q/2./beta_0/beta_0*lN*std::log(1.-lchi);
    ULL[0][0]=(std::exp(1./_as*lnEMLL));
    ULL[0][1]=0.+II*0.;
    ULL[1][0]=0.+II*0.;
    ULL[1][1]=(std::exp(1./_as*lnEPLL));
  }
  if (_ordres >=1){
    lnEPNLL=AP.plus0/beta_0*(std::log(1.-lchi))-Apt1g/2./beta_0/beta_0*lN*
    (beta_1/beta_0*(std::log(1.-lchi))/(1.-lchi)+beta_0*lchi/(1.-lchi)*_LR+beta_0*_LF);
    lnEMNLL=AP.minus0/beta_0*(std::log(1.-lchi))-Apt1q/2./beta_0/beta_0*lN*
    (beta_1/beta_0*(std::log(1.-lchi))/(1.-lchi)+beta_0*lchi/(1.-lchi)*_LR+beta_0*_LF);
    UNLL[0][0]=(std::exp(lnEMNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEPNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    UNLL[0][1]=(std::exp(lnEPNLL)-std::exp(lnEMNLL))*AP.Sg0/rad;
    UNLL[1][0]=(std::exp(lnEPNLL)-std::exp(lnEMNLL))*AP.gS0/rad;
    UNLL[1][1]=(std::exp(lnEPNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEMNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    V1[0][0]=(-beta_1*(4.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.SS0+beta_0*beta_0*(AP.gS1*AP.Sg0-AP.gS0*AP.Sg1+beta_1*AP.SS0)
	    -std::pow(beta_0,3)*AP.SS1+beta_0*(2.*AP.gg1*AP.gS0*AP.Sg0-(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)*(AP.gg0-AP.SS0)
	    +(2.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.SS1))/(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[0][1]=(-std::pow(beta_0,3)*AP.Sg1-beta_1*AP.Sg0*(4.*AP.Sg0*AP.gS0+std::pow(AP.gg0-AP.SS0,2))+beta_0*AP.Sg0*(2.*(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)
	    +(AP.gg0-AP.SS0)*(AP.gg1-AP.SS1))+beta_0*beta_0*(AP.Sg1*(-AP.gg0+AP.SS0)+AP.Sg0*(beta_1+AP.gg1-AP.SS1)))
	    /(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[1][0]=(-std::pow(beta_0,3)*AP.gS1-beta_1*AP.gS0*(4.*AP.Sg0*AP.gS0+std::pow(AP.gg0-AP.SS0,2))+beta_0*AP.gS0*(2.*(AP.gS1*AP.Sg0+AP.gS0*AP.Sg1)
	    +(AP.SS0-AP.gg0)*(AP.SS1-AP.gg1))+beta_0*beta_0*(AP.gS1*(+AP.gg0-AP.SS0)+AP.gS0*(beta_1-AP.gg1+AP.SS1)))
	    /(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
    V1[1][1]=(-beta_1*(4.*AP.gS0*AP.Sg0+std::pow(AP.gg0-AP.SS0,2))*AP.gg0+beta_0*beta_0*(-AP.gS1*AP.Sg0+AP.gS0*AP.Sg1+beta_1*AP.gg0)
	    -std::pow(beta_0,3)*AP.gg1+beta_0*(2.*AP.SS1*AP.Sg0*AP.gS0-(AP.Sg1*AP.gS0+AP.Sg0*AP.gS1)*(AP.SS0-AP.gg0)
	    +(2.*AP.Sg0*AP.gS0+std::pow(AP.SS0-AP.gg0,2))*AP.gg1))/(beta_0*beta_0*(beta_0*beta_0-4.*AP.gS0*AP.Sg0-std::pow(AP.gg0-AP.SS0,2)));
  }
  if (_ordres >=2){
    lnEPNNLL=AP.plus0/beta_0*(std::log(1.-lchi)+_as*(beta_1/beta_0*(std::log(1.-lchi)/(1.-lchi))+beta_0*lchi/(1.-lchi)*_LR+beta_0*_LF))
    -_as*Apt1g/2./beta_0/beta_0*lN*((2.*(beta_1*beta_1-beta_0*beta_2)*lchi+2.*beta_1*beta_1*std::log(1.-lchi)
    -beta_1*beta_1*std::pow(std::log(1.-lchi),2))/(2.*beta_0*beta_0*std::pow(1.-lchi,2))
    +_LR*(beta_1*(2.-lchi)*lchi-beta_1*std::log(1.-lchi))/std::pow(1.-lchi,2)+beta_0*beta_0*_LR*_LR*(-2.+lchi)*lchi/(2.*std::pow(1.-lchi,2))
    -0.5*beta_0*beta_0*_LR*_LR+0.5*beta_0*beta_0*std::pow(-_LF+_LR,2)+_LF*beta_1);
    lnEMNNLL=AP.minus0/beta_0*(std::log(1.-lchi)+_as*(beta_1/beta_0*(std::log(1.-lchi)/(1.-lchi))+beta_0*lchi/(1.-lchi)*_LR+beta_0*_LF))
    -_as*Apt1q/2./beta_0/beta_0*lN*((2.*(beta_1*beta_1-beta_0*beta_2)*lchi+2.*beta_1*beta_1*std::log(1.-lchi)
    -beta_1*beta_1*std::pow(std::log(1.-lchi),2))/(2.*beta_0*beta_0*std::pow(1.-lchi,2))
    +_LR*(beta_1*(2.-lchi)*lchi-beta_1*std::log(1.-lchi))/std::pow(1.-lchi,2)+beta_0*beta_0*_LR*_LR*(-2.+lchi)*lchi/(2.*std::pow(1.-lchi,2))
    -0.5*beta_0*beta_0*_LR*_LR+0.5*beta_0*beta_0*std::pow(-_LF+_LR,2)+_LF*beta_1);
    UNNLL[0][0]=(std::exp(lnEMNNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEPNNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
    UNNLL[0][1]=(std::exp(lnEPNNLL)-std::exp(lnEMNNLL))*AP.Sg0/rad;
    UNNLL[1][0]=(std::exp(lnEPNNLL)-std::exp(lnEMNNLL))*AP.gS0/rad;
    UNNLL[1][1]=(std::exp(lnEPNNLL)*(AP.gg0-AP.SS0+rad)+std::exp(lnEMNNLL)*(AP.SS0-AP.gg0+rad))/(2.*rad);
  }
  return;
}

std::vector<std::complex<long double> > Joint::ComputeJointRes(std::complex<long double> N, std::complex<long double> b){
  /* channel legends
   * 0 all channels
   * 1 gg only
   * 2 gq only
   * 3 qq only (means all flavours since Higgs resummation is flavour blind)
   */
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> bbar= b/2.*std::exp(EulerGamma);
  const std::complex<long double> chi= Nbar*Nbar+bbar*bbar;
  const std::complex<long double> LC= std::log(chi);
  const std::complex<long double> lchi= _as*beta_0*LC;
  const std::complex<long double> lN=_as*beta_0*std::log(Nbar*Nbar);
  std::complex<long double> RggNLL,RggNNLL,RgqNLL;
  std::vector<std::complex<long double> > ris;
  ComputeEvolution(N,b);
  gsl_sf_result redilog,imdilog;
  long double r=std::abs(Nbar*Nbar/chi);
  long double theta=std::arg(Nbar*Nbar/chi);
  gsl_sf_complex_dilog_e(r,theta,&redilog,&imdilog);
  std::complex<long double> risdilog(redilog.val,imdilog.val);
  //LL resummation
  switch(_ordres){
    case(0):{
      if ((_channel==0)||(_channel==1)){
	ris.push_back(std::pow(ULL[1][1],2)*Sudakov_g(N,b));
      }
      else ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      break;
    }
    case(1):{
      RggNLL=std::exp(1./2./beta_0/beta_0*lN*(Apt2g-beta_1/beta_0*Apt1g)*lchi/(1.-lchi));
      if((_channel==0)||(_channel==1)){
	ris.push_back((std::pow(RggNLL*UNLL[1][1],2)
	+_as*((Hpt1g+Apt1g*risdilog-2.*beta_0*_LR)*std::pow(ULL[1][1],2)))*Sudakov_g(N,b));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==2)){
	ris.push_back((std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]
	+_as*(Hgggq1(N)*ULL[1][1]*ULL[0][0]-std::pow(ULL[1][1],2)*V1[1][0]
	+ULL[1][1]*ULL[0][0]*V1[1][0]))*Sudakov_g(N,b));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==3)){
	ris.push_back((std::pow(RggNLL*UNLL[1][0],2))*Sudakov_g(N,b));
      }
      else ris.push_back((0.,0.));
      break;
    }
    case(2):{
      RggNLL=std::exp(1./2./beta_0/beta_0*lN*(Apt2g-beta_1/beta_0*Apt1g)*lchi/(1.-lchi));
      RggNNLL=std::exp(_as*(V1[1][1]*lchi/(1.-lchi)
      -1./std::pow(beta_0,3)*(beta_1*beta_1/2./beta_0*Apt1g-beta_2*Apt1g/2.-beta_1*Apt2g/2.+beta_0*Apt3g/2.)*lN
      *(2.-lchi)*lchi/(2.*std::pow(1.-lchi,2))+1./2./beta_0*(Apt2g-beta_1/beta_0*Apt1g)*lN*lchi*(2.-lchi)/std::pow(1.-lchi,2)*_LR
      -beta_1*(1./2./beta_0/beta_0*lN*(Apt2g-beta_1/beta_0*Apt1g))/beta_0*std::log(1.-lchi)/std::pow(1.-lchi,2)
      -1./2./beta_0*lN*(Apt2g-beta_1/beta_0*Apt1g)*_LF-1./4.*Dpt2g/beta_0*lN));
      RgqNLL=std::exp(-std::log(1.-lchi));
      if((_channel==0)||(_channel==1)){
	ris.push_back((std::pow(RggNNLL*UNNLL[1][1],2)
	+_as*(2.*Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][1]+(Hpt1g+Apt1g*risdilog)*std::pow(RggNLL*UNLL[1][1],2)
	+2.*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][1]*V1[1][0]-2.*RggNLL*RggNLL*UNLL[1][1]*UNLL[1][0]*V1[0][1]
	-2.*std::pow(RggNLL*UNLL[1][1],2)*beta_0*_LR)+_as*_as*((Hggggreg2(N)+Hpt2g-3.*beta_0*(Hpt1g+Apt1g*risdilog)*_LR
	-2.*beta_1*_LR+3.*beta_0*beta_0*_LR*_LR)*std::pow(ULL[1][1],2)))*Sudakov_g(N,b));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==2)){
	ris.push_back((std::pow(RggNNLL,2)*UNNLL[1][1]*UNNLL[1][0]+_as*(Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][0]
	+std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[1][1]+(Hpt1g+risdilog*Apt1g)*std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]
	  -std::pow(RggNLL*UNLL[1][1],2)*V1[1][0]+RggNLL*RgqNLL*UNLL[1][1]*UNLL[0][0]*V1[1][0]-std::pow(RggNLL*UNLL[1][0],2)*V1[0][1]
	  -std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[0][0]-2.*beta_0*std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*_LR)
	+_as*_as*(Hgggq2(N)*ULL[1][1]*ULL[0][0]
	-3.*Hgggq1(N)*beta_0*ULL[1][1]*ULL[0][0]*_LR-ULL[1][1]*ULL[0][0]*V1[1][0]*beta_0*_LF))*Sudakov_g(N,b));
      }
      else ris.push_back((0.,0.));
      if((_channel==0)||(_channel==3)){
	ris.push_back((std::pow(RggNNLL*UNNLL[1][0],2)+2.*_as*(Hgggq1(N)*RggNLL*RgqNLL*UNLL[1][0]*UNLL[0][0]
	+(Hpt1g+Apt1g*risdilog)*std::pow(RggNLL*UNLL[1][0],2)+std::pow(RggNLL*UNLL[1][0],2)*V1[1][1]
	  -std::pow(RggNLL,2)*UNLL[1][1]*UNLL[1][0]*V1[1][0]+RggNLL*RgqNLL*UNLL[1][0]*UNLL[0][0]*V1[1][0]
	  -std::pow(RggNLL*UNLL[1][0],2)*V1[0][0]-2.*beta_0*_LR*std::pow(RggNLL*UNLL[1][0],2))
	+_as*_as*(Hggqq2(N)*std::pow(ULL[0][0],2)))*Sudakov_g(N,b));
      }
      else ris.push_back((0.,0.));
      break;
    }
  }
  return ris; 
}















































































