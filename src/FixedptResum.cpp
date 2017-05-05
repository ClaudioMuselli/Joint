#include "FixedptResum.h"

FixedptResum::FixedptResum(const int ordres, const int ordmatch, int channel, bool Wilson, long double Nc, long double Nf){
  _ordres=ordres; _ordmatch=ordmatch; _channel=channel;_Wilson=Wilson;_Nc=Nc;_Nf=Nf;
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
  if (_ordmatch>3){
    cout << "Matching is implemented only up to NNNLO! Set matching order to NNNLO" << endl;
    _ordmatch=3;
  }
  
}


FixedptResum::~FixedptResum(){
}


void FixedptResum::SetChannel(int channel){
  _channel=channel;
}

void FixedptResum::SetOrdRes(int ordres){
  _ordres=ordres;
}

void FixedptResum::SetOrdMatch(int ordmatch){
  _ordmatch=ordmatch;
}

















//Sudakov Exponent

std::complex<long double> FixedptResum::Sudakov_th_gggH(std::complex<long double> N, long double xp){
  std::complex<long double> Sud(0.,0.);
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  const long double LR=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUR,2.));
  const long double LF=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUF,2.));
  if (_ordres==0) return std::exp(Sud);
  if (_ordres>=1){
    Sud+=1./ConstResum::as*(Ath1g/beta_0/beta_0*(lN+0.5*(1.-lN)*std::log(1.-lN)+(1.-lN/2.)*std::log(1.-lN/2.)));//Controllare
  }
  if (_ordres>=2){
    const long double radxp=std::log(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2)/xp);
    Sud+=(Ath1g*beta_1/(std::pow(beta_0,3))*(lN+0.5*std::log(1.-lN)+0.25*std::pow(std::log(1.-lN),2)
    +std::log(1.-lN/2.)+0.5*std::pow(std::log(1.-lN/2.),2))+Ath1g*LR/beta_0*(lN+0.5*std::log(1.-lN)+std::log(1.-lN/2.))
    -lN*LF/beta_0*Ath1g-Ath2g/std::pow(beta_0,2)*(lN+0.5*std::log(1.-lN)+std::log(1.-lN/2.))
    +Bth1g/beta_0*std::log(1.-lN/2.)-Ath1g/(2.*beta_0)*std::log(1.-lN)*radxp);
  }
  return std::exp(Sud);
}
std::complex<long double> FixedptResum::Sudakov_th_gqqH(std::complex<long double> N, long double xp){
 std::complex<long double> Sud(0.,0.);
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  const long double LR=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUR,2.));
  const long double LF=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUF,2.));
  if (_ordres==0) return std::exp(Sud);
  if (_ordres>=1){
    Sud+=1./ConstResum::as*(1./2./beta_0/beta_0*((Ath1g+Ath1q)*lN+Ath1g*(1.-lN)*std::log(1.-lN)+2.*Ath1q*(1.-lN/2.)*std::log(1.-lN/2.)));
  }
  if (_ordres>=2){
    const long double radxp=std::log(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2)/xp);
    Sud+=((Ath1g+Ath1q)/2.*beta_1/(std::pow(beta_0,3))*(lN+std::log(1.-lN)+0.5*std::pow(std::log(1.-lN),2))
    -Ath1q*beta_1/(2.*std::pow(beta_0,3))*(std::log(1.-lN)+0.5*std::pow(std::log(1.-lN),2)-2.*std::log(1.-lN/2.)
    -std::pow(std::log(1.-lN/2.),2))+(Ath1g+Ath1q)/2.*LR/beta_0*(lN+std::log(1.-lN))
    +Ath1q/(2.*beta_0)*LR*(2.*std::log(1.-lN/2.)-std::log(1.-lN))
    -lN*LF/beta_0*(Ath1g+Ath1q)/2.-(Ath2g+Ath2q)/(2.*std::pow(beta_0,2))*(lN+std::log(1.-lN))
    +Ath2q/(2.*beta_0*beta_0)*(std::log(1.-lN)-2.*std::log(1.-lN/2.))
    +Bth1q/beta_0*std::log(1.-lN/2.)-Ath1q/(2.*beta_0)*std::log(1.-lN)*radxp);
  }
  return std::exp(Sud);
}
std::complex<long double> FixedptResum::Sudakov_th_qqgH(std::complex<long double> N, long double xp){
  std::complex<long double> Sud(0.,0.);
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> lN=2.*ConstResum::as*beta_0*std::log(Nbar);
  const long double LR=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUR,2.));
  const long double LF=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUF,2.));
  if (_ordres==0) return std::exp(Sud);
  if (_ordres>=1){
    Sud+=1./ConstResum::as*(1./beta_0/beta_0*(Ath1q*lN-0.5*(Ath1g-2.*Ath1q)*(1.-lN)*std::log(1.-lN)+Ath1g*(1.-lN/2.)*std::log(1.-lN/2.)));
  }
  if (_ordres>=2){
    const long double radxp=std::log(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2)/xp);
    Sud+=(Ath1q*beta_1/(std::pow(beta_0,3))*(lN+std::log(1.-lN)+0.5*std::pow(std::log(1.-lN),2))
    -Ath1g*beta_1/(2.*std::pow(beta_0,3))*(std::log(1.-lN)+0.5*std::pow(std::log(1.-lN),2)-2.*std::log(1.-lN/2.)
    -std::pow(std::log(1.-lN/2.),2))+Ath1q*LR/beta_0*(lN+std::log(1.-lN))
    +Ath1g/(2.*beta_0)*LR*(2.*std::log(1.-lN/2.)-std::log(1.-lN))
    -lN*LF/beta_0*Ath1q-Ath2q/std::pow(beta_0,2)*(lN+std::log(1.-lN))
    +Ath2g/(2.*beta_0*beta_0)*(std::log(1.-lN)-2.*std::log(1.-lN/2.))
    +Bth1g/beta_0*std::log(1.-lN/2.)-Ath1g/(2.*beta_0)*std::log(1.-lN)*radxp);
  }
  return std::exp(Sud);
}


//Higgs cross section
//LO cross section (WITHOUT sigma_0)

std::complex<long double> FixedptResum::LOgggH(std::complex<long double> NN, long double xp){
  std::complex<long double> N=NN-1.;
  std::complex<long double> CLOgggH;
  std::complex<long double> half(0.5,0.);
  std::complex<long double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);
  CLOgggH=2.*ConstResum::as*_Ca/std::sqrt(M_PIl)*1./xp*std::exp(LogGamma(N)-LogGamma(N+0.5))*(Hyp2F1(half,N,N+0.5,xprad)
  -2.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*N/(N+0.5)*Hyp2F1(half,N+1.,N+1.5,xprad)
  +((1.+xp)*(3.+xp))/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),4.))*N*(N+1.)/((N+0.5)*(N+1.5))*Hyp2F1(half,N+2.,N+2.5,xprad)
  -2.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),6.))*N*(N+1.)*(N+2.)/((N+0.5)*(N+1.5)*(N+2.5))*Hyp2F1(half,N+3.,N+3.5,xprad)
  +1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),8))*N*(N+1.)*(N+2.)*(N+3.)/((N+0.5)*(N+1.5)*(N+2.5)*(N+3.5))
  *Hyp2F1(half,N+4.,N+4.5,xprad));
  return CLOgggH;
}
std::complex<long double> FixedptResum::LOgqqH(std::complex<long double> NN, long double xp){
  std::complex<long double> N=NN-1.;
  std::complex<long double> CLOgqqH;
  std::complex<long double> half(0.5,0.);
  std::complex<long double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);
  CLOgqqH=ConstResum::as*_Cf/std::sqrt(M_PIl)*1./xp*std::exp(LogGamma(N)-LogGamma(N+0.5))*(Hyp2F1(half,N,N+0.5,xprad)
  -(4.+3.*xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*N/(N+0.5)*Hyp2F1(half,N+1.,N+1.5,xprad)
  +3.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),4.))*N*(N+1.)/((N+0.5)*(N+1.5))*Hyp2F1(half,N+2.,N+2.5,xprad)
  -1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),6.))*N*(N+1.)*(N+2.)/((N+0.5)*(N+1.5)*(N+2.5))*Hyp2F1(half,N+3.,N+3.5,xprad));
  return CLOgqqH;
}
std::complex<long double> FixedptResum::LOqqgH(std::complex<long double> NN, long double xp){
  std::complex<long double> N=NN-1.;
  std::complex<long double> CLOqqgH;
  std::complex<long double> half(0.5,0.);
  std::complex<long double> xprad(std::pow(std::sqrt(1.+xp)-std::sqrt(xp),4.),0.);
  CLOqqgH=2.*ConstResum::as*_Cf*_Cf/std::sqrt(M_PIl)*1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*std::exp(LogGamma(N)-LogGamma(N+0.5))
  *(Hyp2F1(half,N,N+0.5,xprad)
  -2.*(1.+xp)/(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),2.))*N/(N+0.5)*Hyp2F1(half,N+1.,N+1.5,xprad)
  +1./(std::pow(std::sqrt(1.+xp)+std::sqrt(xp),4.))*N*(N+1.)/((N+0.5)*(N+1.5))*Hyp2F1(half,N+2.,N+2.5,xprad));
  return CLOqqgH;
  
}

//Matching Constants

long double FixedptResum::Hth1gggH(long double xp){
  const long double LR=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUR,2.));
  const long double LF=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUF,2.));
  const long double g1gg= 67./36.*_Ca-5./18.*_Nf+_Ca*gsl_sf_zeta(2.)-beta_0*M_PIl*std::log(xp/(1.+xp))
  -1./8.*_Ca*std::pow(std::log(xp/(1.+xp)),2)+2.*_Ca*gsl_sf_dilog(1.-std::sqrt(xp/(1.+xp)))
  +_Ca*std::log(1.-std::sqrt(xp/(1.+xp)))*std::log(xp/(1.+xp))-0.5*_Ca*std::log(1.+std::sqrt(xp/(1.+xp)))
  *std::log(xp/(1.+xp))+0.5*_Ca*std::pow(std::log(1.+std::sqrt(xp/(1.+xp))),2)
  +2.*beta_0*M_PIl*std::pow(std::log(1.+std::sqrt(xp/(1.+xp))),2)+_Ca*gsl_sf_dilog((2.*std::sqrt(xp))/(std::sqrt(1.+xp)+std::sqrt(xp)))
  -((_Ca-_Nf)*(std::sqrt(xp)*std::sqrt(1.+xp)*(1.+xp)-2.*xp-xp*xp))/(6.*(1.+8.*xp+9.*xp*xp))+2.*beta_0*M_PIl*LF-3.*beta_0*M_PIl*LR;
  return g1gg/M_PIl;
}
long double FixedptResum::Hth1gqqH(long double xp){
  const long double LR=std::log(std::pow(ConstResum::Q/ConstResum::MUR,2.));
  const long double LF=std::log(std::pow(ConstResum::Q/ConstResum::MUF,2.));
  const long double g1gq=-7./4.*_Cf+134./36.*_Ca-20./36.*_Nf-8.*_Cf*gsl_sf_zeta(2.)+12.*_Ca*gsl_sf_zeta(2.)
  -4.*beta_0*M_PIl*std::log(xp/(1.+xp))+3./2.*_Cf*std::log(xp/(1.+xp))-0.5*_Ca*std::pow(std::log(xp/(1.+xp)),2)
  +4.*(_Cf+_Ca)*gsl_sf_dilog(1.-std::sqrt(xp/(1.+xp)))
  +(2.*(_Ca-_Cf)*(1.+3.*xp+3.*std::sqrt(xp*(1.+xp))))/(2.*std::sqrt(xp*(1.+xp))+1.+3.*xp)
  +8.*beta_0*M_PIl*std::log(1.+std::sqrt(xp/(1.+xp)))-3.*_Cf*std::log(1.+std::sqrt(xp/(1.+xp)))
  +2.*_Cf*std::log(1.-std::sqrt(xp/(1.+xp)))*std::log(xp/(1.+xp))
  +2.*_Ca*std::log(1.-std::sqrt(xp/(1.+xp)))*std::log(xp/(1.+xp))
  -2.*_Cf*std::log(1.+std::sqrt(xp/(1.+xp)))*std::log(xp/(1.+xp))
  -2.*_Cf*std::pow(std::log(1.+std::sqrt(xp/(1.+xp))),2)+4.*_Cf*gsl_sf_dilog(2.*std::sqrt(xp)/(std::sqrt(1.+xp)+std::sqrt(xp)))
  +beta_0*M_PIl*LF+3./4.*_Cf*LF-3.*beta_0*M_PIl*LR;
  return g1gq/M_PIl;
}
long double FixedptResum::Hth1qqgH(long double xp){
  const long double LR=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUR,2.));
  const long double LF=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUF,2.));
  const long double g1qq=-9./2.*_Cf+79./12.*_Ca-5./6.*_Nf+12.*_Cf*gsl_sf_zeta(2.)-10.*_Ca*gsl_sf_zeta(2.)
  -(_Cf-_Ca)*std::sqrt(1.+xp)/std::sqrt(xp)+4.*_Cf*gsl_sf_dilog(1.-std::sqrt(xp/(1.+xp)))
  -3./4.*_Cf*std::log(xp/(1.+xp))-beta_0*M_PIl*std::log(xp/(1.+xp))+0.25*_Ca*std::pow(std::log(xp/(1.+xp)),2)
  -0.5*_Cf*std::pow(std::log(xp/(1.+xp)),2)+2.*_Cf*std::log(1.-std::sqrt(xp/(1.+xp)))*std::log(xp/(1.+xp))
  +1.5*_Cf*std::log(1.+std::sqrt(xp/(1.+xp)))+2.*beta_0*M_PIl*std::log(1.+std::sqrt(xp/(1.+xp)))
  +_Ca*std::pow(std::log(1.+std::sqrt(xp/(1.+xp))),2)-_Ca*std::log(1.+std::sqrt(xp/(1.+xp)))*std::log(xp/(1.+xp))
  +2.*_Ca*gsl_sf_dilog(2.*std::sqrt(xp)/(std::sqrt(1.+xp)+std::sqrt(xp)))+3./2.*_Cf*LF-3.*beta_0*M_PIl*LR;
  return g1qq/M_PIl;
}


std::vector<std::complex<long double>> FixedptResum::ComputeFixedptResum(std::complex<long double> N, long double xp){
  /* channel legends
   * 0 all channels
   * 1 gg only
   * 2 gq only
   * 3 qq only (means all flavours since Higgs resummation is flavour blind)
   */
  std::vector<std::complex<long double>> ris;
  switch(_ordres){
    case(0):{
      ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      ris.push_back((0.,0.));
      break;
    }
    case(1):{
      if ((_channel==0)||(_channel==1)) 
	ris.push_back(LOgggH(N,xp)*Sudakov_th_gggH(N,xp));
      else ris.push_back((0.,0.));
      if ((_channel==0)||(_channel==2)) 
	ris.push_back(LOgqqH(N,xp)*Sudakov_th_gqqH(N,xp));
      else ris.push_back((0.,0.));
      if ((_channel==0)||(_channel==3)) 
	ris.push_back(LOqqgH(N,xp)*Sudakov_th_qqgH(N,xp));
      else ris.push_back((0.,0.));
      break;
    }
    case(2):{
      if ((_channel==0)||(_channel==1)) 
	ris.push_back(LOgggH(N,xp)*(1.+ConstResum::as*Hth1gggH(xp))*Sudakov_th_gggH(N,xp));
      else ris.push_back((0.,0.));
      if ((_channel==0)||(_channel==2)) 
	ris.push_back(LOgqqH(N,xp)*(1.+ConstResum::as*Hth1gqqH(xp))*Sudakov_th_gqqH(N,xp));
      else ris.push_back((0.,0.));
      if ((_channel==0)||(_channel==3)) 
	ris.push_back(LOqqgH(N,xp)*(1.+ConstResum::as*Hth1qqgH(xp))*Sudakov_th_qqgH(N,xp));
      else ris.push_back((0.,0.));
      break;
    }
  }
  return ris;  
}

std::vector<std::complex<long double>> FixedptResum::ComputeMatching(std::complex<long double> N, long double xp){
  std::vector<std::complex<long double> > match;
  const std::complex<long double> Nbar= N*std::exp(EulerGamma);
  const std::complex<long double> lN=2.*std::log(Nbar);
  const long double LR=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUR,2.));
  const long double LF=std::log(std::pow(ConstResum::Q*(std::sqrt(1.+xp)+std::sqrt(xp))/ConstResum::MUF,2.));
  switch (_ordres){
    case (0):{
      match.push_back((0.,0.));
      match.push_back((0.,0.));
      match.push_back((0.,0.));
      break;
    }
    case(1):{
      switch(_ordmatch){
	case(0):{
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  break;
	}
	case(1):{
	  if ((_channel==0)||(_channel==1)) 
	    match.push_back(LOgggH(N,xp));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2)) 
	    match.push_back(LOgqqH(N,xp));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3)) 
	    match.push_back(LOqqgH(N,xp));
	  else match.push_back((0.,0.));
	  break;
	}
	case(2):{
	  if ((_channel==0)||(_channel==1)) 
	    match.push_back(LOgggH(N,xp)*(1.+3./8.*Ath1g*lN*lN*ConstResum::as));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2)) 
	    match.push_back(LOgqqH(N,xp)*(1.+1./8.*(2.*Ath1g+Ath1q)*lN*lN*ConstResum::as));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3)) 
	    match.push_back(LOqqgH(N,xp)*(1.+1./8.*(4.*Ath1q-Ath1g)*lN*lN*ConstResum::as));
	  else match.push_back((0.,0.));
	  break;
	}
	case(3):{
	  if ((_channel==0)||(_channel==1)) 
	    match.push_back(LOgggH(N,xp)*(1.+3./8.*Ath1g*lN*lN*ConstResum::as+1./384.*Ath1g*lN*lN*lN*(27.*Ath1g*lN+40.*beta_0)*ConstResum::as*ConstResum::as));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2)) 
	    match.push_back(LOgqqH(N,xp)*(1.+1./8.*(2.*Ath1g+Ath1q)*lN*lN*ConstResum::as
	    +1./384.*lN*lN*lN*(3.*std::pow(2.*Ath1g+Ath1q,2)*lN+8.*(4.*Ath1g+Ath1q)*beta_0)*ConstResum::as*ConstResum::as));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3)) 
	    match.push_back(LOqqgH(N,xp)*(1.+1./8.*(4.*Ath1q-Ath1g)*lN*lN*ConstResum::as
	    +1./384.*lN*lN*lN*(3.*std::pow(Ath1g-4.*Ath1q,2)*lN+8.*(-3.*Ath1g+8.*Ath1q)*beta_0)*ConstResum::as*ConstResum::as));
	  else match.push_back((0.,0.));
	  break;
	}
      }
      break;
    }
    case(2):{
      long double radxp=std::pow(std::sqrt(xp)+std::sqrt(1.+xp),2)/xp;
      switch(_ordmatch){
	case(0):{
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  match.push_back((0.,0.));
	  break;
	}
	case(1):{
	  if ((_channel==0)||(_channel==1)) 
	    match.push_back(LOgggH(N,xp));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2)) 
	    match.push_back(LOgqqH(N,xp));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3)) 
	    match.push_back(LOqqgH(N,xp));
	  else match.push_back((0.,0.));
	  break;
	}
	case(2):{
	  if ((_channel==0)||(_channel==1)) 
	    match.push_back(LOgggH(N,xp)*(1.+ConstResum::as*(Hth1gggH(xp)+1./8.*lN*(-4.*Bth1g-8.*Ath1g*LF
	    +3.*Ath1g*lN+4.*Ath1g*std::log(radxp)))));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2)) 
	    match.push_back(LOgqqH(N,xp)*(1.+ConstResum::as*(Hth1gqqH(xp)+1./8.*lN*(-4.*(Bth1q+(Ath1g+Ath1q)*LF)
	      +(2.*Ath1g+Ath1q)*lN+4.*Ath1q*std::log(radxp)))));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3)) 
	    match.push_back(LOqqgH(N,xp)*(1.+ConstResum::as*(Hth1qqgH(xp)+1./8.*lN*(-4.*Bth1g-8.*Ath1q*LF
	      -Ath1g*lN+4.*Ath1q*lN+4.*Ath1g*std::log(radxp)))));
	  else match.push_back((0.,0.));
	  break;
	}
	case(3):{
	  if ((_channel==0)||(_channel==1)) 
	    match.push_back(LOgggH(N,xp)*(1.+ConstResum::as*(Hth1gggH(xp)+1./8.*lN*(-4.*Bth1g-8.*Ath1g*LF
	    +3.*Ath1g*lN+4.*Ath1g*std::log(radxp)))+ConstResum::as*ConstResum::as*(1./8.*(Hth1gggH(xp)*lN*
	    (-4.*Bth1g-8.*Ath1g*LF+3.*Ath1g*lN+4.*Ath1g*std::log(radxp))+1./48.*lN*lN*(3.*
	      std::pow(4.*Bth1g+8.*Ath1g*LF-3.*Ath1g*lN-4.*Ath1g*std::log(radxp),2.)
	      +8.*(18.*Ath2g+(-6.*Bth1g+5.*Ath1g*lN-18.*Ath1g*LR)*beta_0+12.*Ath1g*beta_0*std::log(radxp)))))));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==2)) 
	    match.push_back(LOgqqH(N,xp)*(1.+ConstResum::as*(Hth1gqqH(xp)+1./8.*lN*(-4.*(Bth1q+(Ath1g+Ath1q)*LF)
	      +(2.*Ath1g+Ath1q)*lN+4.*Ath1q*std::log(radxp)))+ConstResum::as*ConstResum::as*
	      (1./8.*(Hth1gqqH(xp)*lN*(-4.*(Bth1q+(Ath1g+Ath1q)*LF)+(2.*Ath1g+Ath1q)*lN+4.*Ath1q*std::log(radxp))
		+1./48.*lN*lN*(3.*std::pow(-4.*(Bth1q+(Ath1g+Ath1q)*LF)+(2.*Ath1g+Ath1q)*lN+4.*Ath1q*std::log(radxp),2.)
		+8.*(12.*Ath2g+6.*Ath2q+(-6.*Bth1q+4.*Ath1g*lN+Ath1q*lN-6.*(2.*Ath1g+Ath1q)*LR)*beta_0
		+12.*Ath1q*beta_0*std::log(radxp)))))));
	  else match.push_back((0.,0.));
	  if ((_channel==0)||(_channel==3)) 
	    match.push_back(LOqqgH(N,xp)*(1.+ConstResum::as*(Hth1qqgH(xp)+1./8.*lN*(-4.*Bth1g-8.*Ath1q*LF
	      -Ath1g*lN+4.*Ath1q*lN+4.*Ath1g*std::log(radxp)))+ConstResum::as*ConstResum::as*
	      (1./8.*(Hth1qqgH(xp)*lN*(-4.*Bth1g-8.*Ath1q*LF-Ath1g*lN+4.*Ath1q*lN+4.*Ath1g*std::log(radxp))
		+1./48.*lN*lN*(3.*std::pow(4.*Bth1g+8.*Ath1q*LF+Ath1g*lN-4.*Ath1q*lN-4.*Ath1g*std::log(radxp),2.)
		+8.*(-6.*Ath2g+24.*Ath2q+(-6.*Bth1g-3.*Ath1g*lN+8.*Ath1q*lN+6.*(Ath1g-4.*Ath1q)*LR)*beta_0
		+12.*Ath1g*beta_0*std::log(radxp)))))));
	  else match.push_back((0.,0.));
	  break;
	}
      }
      break;
    }
  }
  return match;	
	  
	  
}




























































