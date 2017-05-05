#include "Borel.h"

struct BorelParam{
  std::function<std::complex<long double> (std::complex<long double>, std::complex<long double>, void*)> FF;
  long double x;
  long double xp;
  long double N0;
  long double slope;
  long double Ccutoff;
  long double asbar;
  void *param;
};

int BorelintC(int* ndim, double * x, int* ncomp, double* y, void *p){
  BorelParam par=*(BorelParam *)p;
  long double  theta=2.*M_PIl*x[0];
  long double  omega=par.Ccutoff*x[1];
  long double  r=3.*omega;
  
  std::complex<long double> N=par.N0+(par.slope+II)*std::log(x[2]);
  std::complex<long double> xi=r*std::exp(II*theta);
  y[0]=1./M_PIl*std::imag(par.Ccutoff*(par.slope+II)/x[2]*std::exp(-N*std::log(par.x*std::pow(std::sqrt(1.+par.xp)-std::sqrt(par.xp),2)))
       *std::exp(2.*std::log(N)*xi)*2.*std::pow(par.xp,-0.5-0.5*xi)*std::pow(N,1.-xi)*CBesselK(1.+xi,2.*std::sqrt(par.xp)*N)
       /std::exp(LogGamma(-xi))*(1./par.asbar*std::exp(-omega/par.asbar)*par.FF(N,omega/xi,par.param)));
  r=3.*par.Ccutoff;
  xi= r*std::exp(II*theta);
  y[0]+=1./M_PIl*std::imag(par.Ccutoff*(par.slope+II)/x[2]*std::exp(-N*std::log(par.x*std::pow(std::sqrt(1.+par.xp)-std::sqrt(par.xp),2)))
       *std::exp(2.*std::log(N)*xi)*2.*std::pow(par.xp,-0.5-0.5*xi)*std::pow(N,1.-xi)*CBesselK(1.+xi,2.*std::sqrt(par.xp)*N)
       /std::exp(LogGamma(-xi))*(std::exp(-par.Ccutoff/par.asbar)*par.FF(N,par.Ccutoff/xi,par.param)));
  return 0;
}

int BorelintCp(int* ndim, double * x, int* ncomp, double* y, void *p){
  BorelParam par= *(BorelParam *)p;
  long double N0=3.;
  long double slope=1.0;
  std::complex<long double> N=N0+(slope+II)*std::log(x[0]);
  std::complex<long double> ris;
  long double  omega=par.Ccutoff*x[2];
  long double theta=2.*M_PIl*x[1];
  long double r=2.9; long double xc=-0.3;
  std::complex<long double> xi=xc+r*std::exp(II*theta);
  ris=par.Ccutoff/par.asbar*std::exp(-omega/par.asbar)/xi*r*std::exp(II*theta)*par.FF(N,1./xi,par.param)*
  2.*std::pow(par.xp,-0.5-0.5*xi*omega)*std::pow(N,1.+xi*omega)*CBesselK(1.+xi*omega,2.*N*std::sqrt(par.xp))
  /std::exp(LogGamma(-xi*omega));
  //cout << ris << endl;
  y[0]=std::imag(std::exp(-N*std::log(par.x*std::pow(std::sqrt(1.+par.xp)-std::sqrt(par.xp),2)))*ris/x[0]*(slope+II))/M_PIl;
  return 0;
  
}




int Borel::BorelJointC(int method,long double x, long double xp, long double asbar,
		      std::complex<long double> (Func)(std::complex<long double>, std::complex<long double>, void *), void *params,
		      long double* ris, long double* err, long double Ccutoff, long double N0,long double slope){
  BorelParam pp;
  pp.x=x;
  pp.xp=xp;
  pp.asbar=asbar;
  pp.FF=Func;
  pp.Ccutoff=Ccutoff;
  pp.N0=N0;
  pp.slope=slope;
  pp.param=params;
  double *res=NULL,*error=NULL,*prob=NULL;
  res=new double[1];
  error=new double[1];
  prob=new double[1];
  int fail;
  double prec=1e-12;
  integration::CUBA(method,BorelintC,3,1,prec,&res,&fail,&error,&prob,&pp);
  *ris=res[0];
  *err=error[0];  
  return 0;
}

int Borel::BorelJointCp(int method,long double x, long double xp, long double asbar,
		      std::complex<long double> (Func)(std::complex<long double>, std::complex<long double>, void *), void *params,
		      long double* ris, long double* err, long double Ccutoff, long double N0,long double slope){
  BorelParam pp;
  pp.x=x;
  pp.xp=xp;
  pp.asbar=asbar;
  pp.FF=Func;
  pp.Ccutoff=Ccutoff;
  pp.N0=N0;
  pp.slope=slope;
  pp.param=params;
  double *res=NULL,*error=NULL,*prob=NULL;
  res=new double[1];
  error=new double[1];
  prob=new double[1];
  int fail;
  double prec=1e-10;
  integration::CUBA(method,BorelintCp,3,1,prec,&res,&fail,&error,&prob,&pp);
  *ris=static_cast<long double> (res[0]);
  *err=static_cast<long double> (error[0]);  
  return 0;
}
