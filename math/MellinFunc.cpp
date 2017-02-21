#include "MellinFunc.h"

//Table of Important Mellin of logarithms and harmonic polylogarithms
MellinFunc::MellinFunc():H(false,false,false){
  zeta2=gsl_sf_zeta(2.);
  zeta3=gsl_sf_zeta(3.);
  Li4=0.5174790616738993863;
  log2=std::log(2);
  log2q=log2*log2;
  log2c=log2*log2*log2;
  zeta2q=zeta2*zeta2;
  EulerGamma=0.577215664901532860606512090082; 
}

MellinFunc::~MellinFunc(){
}

std::complex<long double> MellinFunc::Li3zplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-pow(-1.,n)*(H.HS(-3,1,N-1.)-zeta2*H.HS(-2,N-1.)+zeta3*H.HS(-1,N-1.)+3./5.*zeta2q-2.*Li4-3./4.*zeta3*log2+0.5*zeta2*log2q-1./12.*log2q*log2q));
}
std::complex<long double> MellinFunc::S12z2plus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(37./20.*zeta2q+2.*(H.HS(-2,-1,-1,N-1.)+H.HS(-2,1,1,N-1.)+H.HS(2,-1,1,N-1.)+H.HS(2,1,-1,N-1.))
  +2.*(H.HS(-2,-1,N-1.)-H.HS(-2,1,N-1.)-H.HS(2,-1,N-1.)+H.HS(2,1,N-1))*log2-4.*Li4-H.HS(-1,N-1.)*zeta3-H.HS(-2,N-1.)*(zeta2-2.*log2q)
  +H.HS(2,N-1.)*(zeta2-2.*log2q)+zeta2*log2q-1./6.*log2q*log2q-9./2.*zeta3*log2));
}
std::complex<long double> MellinFunc::S12mzplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (pow(-1.,n)*(H.HS(2,1,-1,N-1.)+(H.HS(2,1,N-1.)-H.HS(2,-1,N-1.))*log2-0.5*(H.HS(2,N-1.)-H.HS(-2,N-1.))*log2q-1./8.*zeta3*H.HS(-1,N-1.)-3.*Li4+6./5.*zeta2q
  -11./4.*zeta3*log2+3./4.*zeta2*log2q-1./8.*log2q*log2q));
}

std::complex<long double> MellinFunc::S12zplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(H.HS(-2,1,1,N-1.)-zeta3*H.HS(-1,N-1.)+Li4-1./8.*zeta2q-1./8.*zeta3*log2-1./4.*zeta2*log2q+1./24.*log2q*log2q));
}
std::complex<long double> MellinFunc::Li3mzplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-pow(-1.,n)*(H.HS(3,-1,N-1.)+(H.HS(3,N-1.)-H.HS(-3,N-1.))*log2+0.5*zeta2*H.HS(-2,N-1.)-3./4.*zeta3*H.HS(-1,N-1.)+1./8.*zeta2q-3./4.*zeta3*log2));
}
std::complex<long double> MellinFunc::Li2zLogminus(std::complex<long double> N){
  return(1./N*(-2.*zeta3-zeta2*H.HS(1,N)+1./N*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N))+H.HS(2,1,N)));
}

std::complex<long double> MellinFunc::Li2zLogplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./(2.*N*N)*(-zeta2+H.HS(-1,N)*H.HS(-1,N)+H.HS(2,N)+2.*H.HS(-1,N)*log2-2.*H.HS(1,N)*log2+log2q+2.*N*zeta2*log2
  +pow(-1.,n)*(2.*N*H.HS(-2,1,N)-2.*N*zeta2*(H.HS(-1,N)+log2)+(zeta2+2.*H.HS(-1,1,N)-log2q)+5./4.*N*zeta3)));
}
std::complex<long double> MellinFunc::Li2mzLogplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./(N*N)*(log2*(-0.5*zeta2*N+log2)+pow(-1.,n)*(2.*H.HS(1,-1,N)+N*H.HS(2,-1,N)+H.HS(-1,N)*(0.5*N*zeta2-2.*log2)
  +0.5*N*zeta2*log2-log2*(N*H.HS(-2,N)-2.*H.HS(1,N)-N*H.HS(2,N)+log2)-1./4.*N*zeta3)));
}
std::complex<long double> MellinFunc::Li2mzLogz(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./(N*N)*0.5*zeta2-2./(N*N*N)*log2+pow(-1.,n)*(1./(N*N)*(0.5*zeta2+H.HS(-2,N))+1./(N*N*N)*(2.*H.HS(-1,N)+2.*log2)));
}
std::complex<long double> MellinFunc::Li2zLogzplusminus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./960.*(160.*6.*zeta2*(pow(-1.,n)*H.HS(-2,N-1.)+H.HS(2,N-1))+pow(-1.,n)*(36.*zeta2q-960.*H.HS(-3,1,N-1.)-480.*H.HS(-2,2,N-1.))
  -4.*(36.*zeta2q+60.*(H.HS(2,N-1.)*H.HS(2,N-1.)+H.HS(4,N-1.))+240.*H.HS(3,1,N-1.))));
}
std::complex<long double> MellinFunc::Li2mzLogzplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-pow(-1.,n)*(2.*H.HS(3,-1,N-1.)+H.HS(2,-2,N-1.)-2.*H.HS(-3,N-1.)*log2+2.*H.HS(3,N-1.)*log2+0.5*zeta2*H.HS(2,N-1.)
  +0.5*zeta2*H.HS(-2,N-1.)-4.*Li4+13./8.*zeta2q-7./2.*zeta3*log2+zeta2*log2q-1./6.*log2q*log2q)
  );
}
std::complex<long double> MellinFunc::LogzLogminus2minus(std::complex<long double> N){
  return(H.HS(2,N-1.)*H.HS(2,N-1.)-zeta2*H.HS(2,N-1.)+2.*H.HS(4,N-1.)+H.HS(1,N-1.)*H.HS(1,N-1.)*(H.HS(2,N-1.)-zeta2)+2.*H.HS(1,N-1.)*(H.HS(3,N-1.)-zeta3)-4./5.*zeta2q);
}
std::complex<long double> MellinFunc::Li2mzLogplusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-pow(-1.,n)*(H.HS(1,2,-1,N-1.)+2.*H.HS(2,1,-1,N-1.)+(H.HS(2,1,N-1.)-H.HS(1,-2,N-1.)-2.*H.HS(2,-1,N-1.)+H.HS(1,N-1.)*H.HS(2,N-1.)+H.HS(3,N-1.)-0.5*zeta2*H.HS(-1,N-1.))*log2
    +0.5*zeta2*H.HS(1,-1,N-1.)-(H.HS(2,N-1.)-H.HS(-2,N-1.))*log2q-(1./4.*zeta3-0.5*zeta2*log2)*H.HS(1,N-1.)-3.*Li4+6./5.*zeta2q-21./8.*zeta3*log2+0.5*zeta2*log2q-1./8.*log2q*log2q)
  );
}
std::complex<long double> MellinFunc::Logz2Logminusminus(std::complex<long double> N){
  return(-2.*zeta3*H.HS(1,N-1.)+2.*H.HS(1,N-1.)*H.HS(3,N-1.)-2.*zeta2*H.HS(2,N-1.)+H.HS(2,N-1.)*H.HS(2,N-1.)+3.*H.HS(4,N-1.)-zeta2q/5.);
}
std::complex<long double> MellinFunc::Logz3minusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(21./20.*zeta2q+3.*H.HS(-4,N-1.))-6./5.*zeta2q+3.*H.HS(4,N-1.));
}
std::complex<long double> MellinFunc::Logplusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(H.HS(1,-1,N-1.)-0.5*log2q+(-H.HS(-1,N-1.)+H.HS(1,N-1.))*log2));
}
std::complex<long double> MellinFunc::Li2zLogplusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(3./40.*zeta2q-H.HS(-2,-1,-1,N-1.)-H.HS(1,-2,1,N-1.)-H.HS(2,-1,1,N-1.)+0.5*H.HS(-2,N-1.)*(zeta2-log2q)
    -0.5*H.HS(2,N-1.)*(zeta2-log2q)-(H.HS(-2,-1,N-1.)-H.HS(-2,1,N-1.))*log2-1./4.*zeta2*(-4.*H.HS(1,-1,N-1.)+2.*log2q
    +4.*(H.HS(-1,N-1.)-H.HS(1,N-1.))*log2)-5./8.*H.HS(1,N-1.)*zeta3)
  );
}
std::complex<long double> MellinFunc::Logz2Logplusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(2.*H.HS(3,-1,N-1.)+2.*H.HS(1,-3,N-1.)+2.*H.HS(2,-2,N-1.)+zeta2*H.HS(2,N-1.)+3./2.*zeta3*H.HS(1,N-1.)+2.*H.HS(3,N-1)*log2
	-2.*H.HS(-3,N-1.)*log2-4.*Li4+3./2.*zeta2q-7./2.*zeta3*log2+zeta2*log2q-1./6.*log2q*log2q));
}
std::complex<long double> MellinFunc::LogzLogplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./12./(N*N)*(-12.*log2+pow(-1.,n)*(6.*N*zeta2+12.*N*H.HS(-2,N)+12.*H.HS(-1,N)+12.*log2)));
}
std::complex<long double> MellinFunc::Logz2Logplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./(N*N*N)*2.*log2-pow(-1.,n)*(1./(N*N)*(zeta2+2.*H.HS(-2,N))+1./(N*N*N)*(2.*H.HS(-1,N)+2.*log2)+1./N*(2.*H.HS(-3,N)+3./2.*zeta3)));
}
std::complex<long double> MellinFunc::LogzLogplus2(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-1./(12.*N*N)*(12.*log2q+pow(-1.,n)*(2.*H.HS(1,N)*(6.*N*zeta2+12.*log2)+3.*(8.*N*H.HS(1,-2,N)+8.*H.HS(1,-1,N)+8.*N*H.HS(2,-1,N)
  -4.*log2*(2.*N*H.HS(-2,N)+2.*H.HS(-1,N)-2.*N*H.HS(2,N)+log2)-N*zeta3))));
}
std::complex<long double> MellinFunc::LogzLogplus2plus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(2.*pow(-1.,n)*(H.HS(1,1,-2,N-1.)+H.HS(1,2,-1,N-1.)+H.HS(2,1,-1,N-1.)-H.HS(1,-2,N-1.)*log2
    -H.HS(2,-1,N-1)*log2+H.HS(1,N-1.)*H.HS(2,N-1.)*log2+0.25*zeta2*H.HS(1,N-1.)*H.HS(1,N-1.)+H.HS(3,N-1.)*log2
    +0.25*zeta2*H.HS(2,N-1.)-0.5*H.HS(2,N-1.)*log2q+0.5*H.HS(-2,N-1.)*log2q-1./8.*zeta3*H.HS(1,N-1.)-Li4+
    2./5.*zeta2q-7./8.*zeta3*log2+0.25*zeta2*log2q-1./24.*log2q*log2q)
  );
}
std::complex<long double> MellinFunc::Li2z(std::complex<long double> N){
  return(1./N*zeta2-1./(N*N)*H.HS(1,N));
}
std::complex<long double> MellinFunc::Li2mz(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-0.5*zeta2/N-1./(N*N)*(pow(-1.,n)*(H.HS(-1,N)+log2)-log2));
}
std::complex<long double> MellinFunc::S12z(std::complex<long double> N){
  return(-1./(2.*N*N)*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N)-2.*N*zeta3));
}
std::complex<long double> MellinFunc::S12mz(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./(8.*N*N)*(-8.*pow(-1.,n)*H.HS(1,-1,N)+8.*pow(-1.,n)*H.HS(-1,N)*log2
  -8.*pow(-1.,n)*H.HS(1,N)*log2-4.*log2q+4.*pow(-1.,n)*log2q+N*zeta3));
}
std::complex<long double> MellinFunc::S12z2(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-1./(6.*N*N)*((-1.+pow(-1.,n))*6.*zeta2+12.*pow(-1.,n)*H.HS(-2,N)+12.*H.HS(2,N)+6.*(H.HS(-1,N)*H.HS(-1,N)+H.HS(1,N)*H.HS(1,N)+2.*log2q
  +2.*pow(-1.,n)*(H.HS(1,N)-log2)*(H.HS(-1,N)+log2)+(H.HS(-1,N)-H.HS(1,N))*2.*log2)-6.*N*zeta3));
}
std::complex<long double> MellinFunc::Li3z(std::complex<long double> N){
  return(1./(N*N*N)*H.HS(1,N)-1./(N*N)*zeta2+1./N*zeta3);
}
std::complex<long double> MellinFunc::Li3mz(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*1./(N*N*N)*(log2+H.HS(-1,N))-log2/(N*N*N)+0.5/(N*N)*zeta2-3./4./N*zeta3);
}
std::complex<long double> MellinFunc::Li2zLogz(std::complex<long double> N){
  return(2./(N*N*N)*H.HS(1,N)+1./(N*N)*H.HS(2,N)-2./(N*N)*zeta2);
}
std::complex<long double> MellinFunc::Logminus(std::complex<long double> N){
  return(-H.HS(1,N)/N);
}
std::complex<long double> MellinFunc::Logplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(1./N*(log2-pow(-1.,n)*(H.HS(-1,N)+log2)));
}
std::complex<long double> MellinFunc::Logz(std::complex<long double> N){
  return(-1./(N*N));
}
std::complex<long double> MellinFunc::Logz2(std::complex<long double> N){
  return(2./(N*N*N));
}
std::complex<long double> MellinFunc::Logz3(std::complex<long double> N){
  return(-6./(N*N*N*N));
}
std::complex<long double> MellinFunc::LogzLogminus(std::complex<long double> N){
  return(-1./N*zeta2+1./N*H.HS(2,N)+1./(N*N)*H.HS(1,N));
}
std::complex<long double> MellinFunc::LogzLogminus2(std::complex<long double> N){
  return(1./(N*N)*(-H.HS(1,N)*H.HS(1,N)+2.*N*H.HS(1,N)*(zeta2-H.HS(2,N))-H.HS(2,N)-2.*N*H.HS(3,N)+2.*N*zeta3));
}
std::complex<long double> MellinFunc::Logz2Logminus(std::complex<long double> N){
  return(-2./(N*N*N)*H.HS(1,N)+1./(N*N)*(2.*zeta2-2.*H.HS(2,N))-2./N*(H.HS(3,N)-zeta3));
}
std::complex<long double> MellinFunc::Logz2minus(std::complex<long double> N){
  return(2.*(-H.HS(3,N-1.)+zeta3));
}
std::complex<long double> MellinFunc::Logminus2(std::complex<long double> N){
  return(1/N*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N)));
}
std::complex<long double> MellinFunc::Logminus3(std::complex<long double> N){
  return(-1./N*(H.HS(1,N)*H.HS(1,N)*H.HS(1,N)+3.*H.HS(1,N)*H.HS(2,N)+2.*H.HS(3,N)));
}
std::complex<long double> MellinFunc::LogzLogminusminus(std::complex<long double> N){
  return(zeta3+zeta2*H.HS(1,N-1.)-H.HS(1,N-1.)*H.HS(2,N-1.)-H.HS(3,N-1.));
}
std::complex<long double> MellinFunc::Logzminus(std::complex<long double> N){
  return(-zeta2+H.HS(2,N-1.));
}
std::complex<long double> MellinFunc::Logzminusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(0.25*zeta2+0.5*H.HS(-2,N-1.))+0.5*H.HS(2,N-1.)-0.5*zeta2);
}
std::complex<long double> MellinFunc::plus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(-pow(-1.,n)*(H.HS(-1,N-1.)+log2));
}
std::complex<long double> MellinFunc::Li2minusminus(std::complex<long double> N){
  return(-zeta2*H.HS(1,N-1.)+H.HS(1,2,N-1.)+zeta3);
}
std::complex<long double> MellinFunc::S12zregminus(std::complex<long double> N){
  return(-6./5.*zeta2q+H.HS(2,1,1,N-1.));
}
std::complex<long double> MellinFunc::Li2zregminus(std::complex<long double> N){
  return(H.HS(2,1,N-1.)-2.*zeta3);
}
std::complex<long double> MellinFunc::Li3zregminus(std::complex<long double> N){
  return(-0.5*zeta2q+zeta2*H.HS(2,N-1.)-H.HS(3,1,N-1.));
}
std::complex<long double> MellinFunc::Li2zregLogminusminus(std::complex<long double> N){
  return(6./5.*zeta2q-H.HS(1,2,1,N-1.)-2.*H.HS(2,1,1,N-1.)+2.*H.HS(1,N-1.)*zeta3);
}

std::complex<long double> MellinFunc::Logplus3plus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (1./4.*pow(-1.,n)*(24.*H.HS(1,1,1,-1,N-1.)-log2*(-4.*(pow(H.HS(1,N-1.),3)+3.*H.HS(1,N-1.)*H.HS(2,N-1.)+2.*H.HS(3,N-1.))
  +24.*H.HS(1,1,-1,N-1.)+log2*(6.*(pow(H.HS(1,N-1.),2)+H.HS(2,N-1.))-12.*H.HS(1,-1,N-1.)-4.*H.HS(1,N-1.)*log2
  +log2q+4.*H.HS(-1,N-1.)*log2))));
}
std::complex<long double> MellinFunc::Li3zregminusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (1./720.*(-5.*(pow(M_PIl,4)-12.*pow(M_PIl,2)*H.HS(2,N-1.)+72.*H.HS(3,1,N-1.))-6.*pow(-1.,n)*(pow(M_PIl,4)
  -10.*pow(M_PIl,2)*H.HS(-2,N-1.)+60.*H.HS(-3,1,N-1.)+5.*M_PIl*M_PIl*log2q-5.*(log2q*log2q+24.*Li4+21.*log2*zeta3))));
}
std::complex<long double> MellinFunc::Li3zoverplusplus(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (1./288.*pow(-1.,n)*(pow(M_PIl,4)+288.*H.HS(3,-1,N-1.)-288.*H.HS(1,2,-1,N-1.)-288.*H.HS(2,1,-1,N-1.)+288.*H.HS(1,1,1,-1,N-1.)
  +24.*H.HS(-2,N-1.)*(M_PIl*M_PIl-6.*log2q)-24.*H.HS(1,-1,N-1.)*(M_PIl*M_PIl-6.*log2q)-12.*log2*(24.*H.HS(-3,N-1.)
  -4.*pow(H.HS(1,N-1.),3)-8.*(H.HS(3,N-1.)+3.*(H.HS(1,-2,N-1.)+H.HS(2,-1,N-1.)-H.HS(1,1,-1,N-1.)))
  -(M_PIl*M_PIl+6.*H.HS(2,N-1.))*log2+log2c-2.*H.HS(-1,N-1.)*(M_PIl*M_PIl-2.*log2q)
  +2.*H.HS(1,N-1.)*(M_PIl*M_PIl+6.*H.HS(2,N-1.)-2.*log2q)+H.HS(1,N-1.)*H.HS(1,N-1.)*6.*log2)
  -36.*(7.*H.HS(-1,N-1.)-2.*H.HS(1,N-1.)+7.*log2)*zeta3));
}
std::complex<long double> MellinFunc::Logplus3(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (log2c-pow(-1.,n)*(6.*H.HS(1,1,-1,N)+log2*(3.*(H.HS(1,N)*H.HS(1,N)+H.HS(2,N))
  -6.*H.HS(1,-1,N)+log2q+(H.HS(-1,N)-H.HS(1,N))*3.*log2)))/N;
}
std::complex<long double> MellinFunc::Li2z2(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (N*M_PIl*M_PIl-12.*H.HS(1,N)-12.*pow(-1.,n)*(H.HS(-1,N)+log2)+12.*log2)/(6.*N*N);
}
std::complex<long double> MellinFunc::Li2z2Logz(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*(N*M_PIl*M_PIl+12.*N*H.HS(-2,N)+24.*H.HS(-1,N)+24.*log2)-3.*(N*M_PIl*M_PIl-8.*H.HS(1,N)
  -4.*N*H.HS(2,N)+8.*log2))/(6.*N*N*N);
}
std::complex<long double> MellinFunc::Li3z2(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (-N*M_PIl*M_PIl+12.*H.HS(1,N)-12.*log2+12.*pow(-1.,n)*(H.HS(-1,N)+log2)+3.*N*N*zeta3)/(3.*N*N*N);
}
std::complex< long double > MellinFunc::Li3overplus(std::complex< long double > N) {
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return (1./(24.*N)*(4.*log2c-12.*zeta2*log2+21.*zeta3+pow(-1.,n)*(-4.*log2c+2.*H.HS(-1,N)*(6.*zeta2-6.*log2q)
  +2.*H.HS(1,N)*(6.*zeta2+6.*log2q)+12.*zeta2*log2+3.*(8.*H.HS(1,-2,N)-8.*H.HS(1,1,-1,N)
  +8.*(1./2.*(-H.HS(1,N)*H.HS(1,N)-H.HS(2,N))+H.HS(1,-1,N))*log2+zeta3))));
}
std::complex<long double> MellinFunc::Li2minus(std::complex<long double> N){
  return( zeta2/N-H.HS(2,N)/N);
}


std::complex<long double> MellinFunc::D0(std::complex<long double> N){
  return(-H.HS(1,N-1.));
}
std::complex<long double> MellinFunc::D1(std::complex<long double> N){
  return(H.HS(1,1,N-1.));
}
std::complex<long double> MellinFunc::D2(std::complex<long double> N){
  return(-2.*H.HS(1,1,1,N-1.));
}
std::complex<long double> MellinFunc::D3(std::complex<long double> N){
  return(6.*H.HS(1,1,1,1,N-1.));
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  