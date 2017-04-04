#include "integration.h"
 
 
double integration::GSL (double (func)(double, void *), double min, double max, double prec, double *err, void *par){
  gsl_integration_workspace *w=gsl_integration_workspace_alloc(10000000);
  gsl_function F;
  F.function=func;
  F.params=par;
  double ris=0.0,error=0.0;
  gsl_integration_qags(&F,min,max,0,prec,10000000, w,&ris,&error);
  gsl_integration_workspace_free(w);
  if(err) *err=error;
  return ris;  
}


int integration::CUBA( int method, int (Func)(int*, double *, int *, double *, void *), int ndim, int ncomp, double prec, double **res,int *fail, double **error, double **prob, void * par){
  static const double epsabs=1e-15;
  static const int verbose=0;
  static const int last=4;
  static const int seed=0;
  static const int mineval=0;
  static const int maxeval=100000;
  
  static const int nstart=1000;
  static const int nincrease=500;
  static const int nbatch=1000;
  static const int gridno=0;
  const char *statefile=NULL;
  
  static const int nnew=1000;
  static const int nmin=2;
  static const double flatness=25.;
  
  static const int key1=13;
  static const int key2=13;
  static const int key3=1;
  static const int maxpass=5;
  static const double border=0.1;
  static const double maxchisq=10.;
  static const double mindeviation=0.25;
  static const int ngiven=0;
  static const int ldxgiven=ndim;
  static const int nextra=0;
  static const int key=13;
  static const int nvec=1;
  static const int spin=-1;
  
  int neval, nregions;
  double *ris=NULL, *err=NULL, *proba=NULL;
  ris=new double[ncomp];
  err=new double[ncomp];
  proba=new double[ncomp];
  integrand_t func = (integrand_t) Func;

  
  switch (method){
    case 0:
      std::cout << "Call to VEGAS" << std::endl;
      Vegas(ndim,ncomp,func,par,nvec,prec, epsabs,verbose,seed,mineval,maxeval,nstart,
	    nincrease,nbatch,gridno,statefile,NULL,&neval, fail,ris,err,proba);
      break;
    case 1: 
      std::cout << "Call to SUAVE" << std::endl;
      Suave(ndim, ncomp, func, par,nvec, prec, epsabs, verbose, seed, mineval, maxeval, nnew, 5,
	    flatness,statefile,NULL, &nregions, &neval, fail, ris, err, proba);
      break;
    case 2:
      std::cout << "Call to DIVONNE" << std::endl;
      Divonne(ndim, ncomp, func, par,nvec, prec, epsabs, verbose, seed, mineval, maxeval, key1, 
	      key2, key3, maxpass, border, maxchisq, mindeviation, ngiven, ldxgiven, NULL, 
	      nextra, NULL,statefile,NULL, &nregions, &neval, fail, ris, err, proba);
      break;
    case 3: 
      std::cout << "Call to Cuhre" << std::endl;
      Cuhre(ndim, ncomp, func, par,nvec, prec, epsabs, verbose | last, mineval, maxeval, key, statefile, NULL,
	    &nregions, &neval, fail, ris, err, proba);
      break;
  }
  *error=err;
  *prob=proba;
  *res=ris;
  return 0;
}

//Use For INTEGRATION
struct InverMellin{
  std::function<std::complex<long double> (std::complex<long double>, void*)> PP;
  long double N0;
  long double Nslope;
  long double tau;
  void *param;
};

struct InverBessel{
  std::function<std::complex<long double> (std::complex<long double>, void*)> BB;
  long double bc;
  long double bslope;
  long double xp;
  long double v;
  void *param;
};

struct InverTotal{
  std::function<std::complex<long double> (std::complex<long double>,  std::complex<long double>, void* )> TT;
  long double bc;
  long double N0;
  long double Nslope;
  long double bslope;
  long double xp;
  long double v;
  long double tau;
  void *param;
};

//MELLIN INVERSE
int IM1(int* ndim, double * x, int* ncomp, double* y, void *p){
  InverMellin par= *(InverMellin *)p;
  std::complex<long double> N=par.N0+par.Nslope*(1.+II)*std::log(x[0]);
  
  y[0]=std::imag(par.Nslope*(1.+II)*std::exp(-N*std::log(par.tau))/M_PI/x[0]*par.PP(N,&par.param));
  return 0;
}
int IM2(int* ndim, double * x, int* ncomp, double* y, void *p){
  InverMellin par= *(InverMellin *)p;
  std::complex<long double> N=par.N0+par.Nslope*(1.-II)*std::log(x[0]);
  
  y[0]=std::imag(-par.Nslope*(1.-II)*std::exp(-N*std::log(par.tau))/M_PI/x[0]*par.PP(N,&par.param));
  return 0;
}

double integration::InverseMellin(int method, std::complex<long double> (Func)(std::complex<long double> , void* ), 
				  long double tau, long double N0, long double slope, bool sign, void *pp) {
  double *res=NULL, prec=1e-8;
  int fail; double *err=NULL, *prob=NULL;
  res=new double[1];
  err=new double[1];
  prob=new double[1];
  InverMellin IMC;
  IMC.PP=Func;
  IMC.N0=N0;
  IMC.Nslope=slope;
  IMC.tau=tau;
  IMC.param=pp;
  //sign==true sign=-1, sign==false sign=+1
  if (sign)
    CUBA(method,IM2,2,1,prec,&res,&fail,&err,&prob,&IMC);
  else 
    CUBA(method,IM1,2,1,prec,&res,&fail,&err,&prob,&IMC);
  return res[0];
}

int IB3(int* ndim, double *x, int* ncomp, double* y,void* p){
  InverBessel par=*(InverBessel *)p;
  std::complex<long double> b=par.bc*x[0];
  y[0]=std::real(par.bc*b/2.*(std::complex<long double>)(sp_bessel::besselJ(0.,(std::complex<double>)(b*std::sqrt(par.xp))))*par.BB(b,par.param));
  y[1]=std::imag(par.bc*b/2.*(std::complex<long double>)(sp_bessel::besselJ(0.,(std::complex<double>)(b*std::sqrt(par.xp))))*par.BB(b,par.param));
}

int IB1(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverBessel par=*(InverBessel *) p;
  std::complex<long double> b=par.bc-par.bslope*(1.+II)*std::log(x[0]);
  std::complex<long double> theta=-II*par.v*M_PIl+x[1]*M_PIl*(-1.+2.*II*par.v);
  
  y[0]=std::real(-b/4.*(-1.+2.*II*par.v)*par.bslope*(1.+II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  y[1]=std::imag(-b/4.*(-1.+2.*II*par.v)*par.bslope*(1.+II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  return 0;
  
}

int IB2(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverBessel par=*(InverBessel *) p;
  std::complex<long double> b=par.bc-par.bslope*(1.-II)*std::log(x[0]);
  std::complex<long double> theta=-II*par.v*M_PIl+x[1]*M_PIl*(1.+2.*II*par.v);
  //std::complex<long double> PROVAT(std::complex<long double> N, std::complex<long double> b, void * par){
//  return (log_c_angle(N*N+b*b/4.,-M_PIl/2.));
//}
  y[0]=std::real(b/4.*(1.+2.*II*par.v)*par.bslope*(1.-II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  y[1]=std::imag(b/4.*(1.+2.*II*par.v)*par.bslope*(1.-II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  return 0;
}

std::complex<long double> integration::InverseBessel(int method, std::complex<long double> (Func)(std::complex<long double> , void*),
						     long double xp, long double bc, long double v, long double slope, void *pp){
  double prec=1e-8;
  int fail; double *err, *prob;
  InverBessel IBC;
  IBC.BB=Func;
  IBC.v=v;
  IBC.bslope=slope;
  IBC.xp=xp;
  IBC.param=pp;
  IBC.bc=bc;
  double *ris1,*ris2,*ris3;
  ris1=new double[2];
  ris2=new double[2];
  err=new double[2];
  prob=new double[2];
  CUBA(method,IB1,2,2,prec,&ris1,&fail,&err,&prob,&IBC);
  CUBA(method,IB2,2,2,prec,&ris2,&fail,&err,&prob,&IBC);
  CUBA(method,IB3,2,2,prec,&ris3,&fail,&err,&prob,&IBC);
  std::complex<long double> result;
  result=ris1[0]+ris2[0]+ris3[0]+II*(ris1[1]+ris2[1]+ris3[1]);
  return(result);  
}


int IT1(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverTotal par=*(InverTotal *)p;
  std::complex<long double> N=par.N0+par.Nslope*(1.-II)*std::log(x[0]);
  std::complex<long double> b=par.bc*x[1];
  
  y[0]=std::imag(-par.Nslope*(1.-II)/x[0]*par.bc/M_PIl*std::exp(-N*std::log(par.tau))*b/2.
      *(std::complex<long double>)(sp_bessel::besselJ(0.,(std::complex<double>)(b*std::sqrt(par.xp))))*par.TT(N,b,par.param));
  return 0;
}

int IT2(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverTotal par=*(InverTotal *)p;
  std::complex<long double> N=par.N0+par.Nslope*(1.-II)*std::log(x[0]);
  std::complex<long double> b=par.bc-par.bslope*(1.+II)*std::log(x[1]);
  std::complex<long double> theta=-II*par.v*M_PIl+x[2]*M_PIl*(-1.+2.*II*par.v);
  
  y[0]=std::imag(par.Nslope*(1.-II)/x[0]/M_PIl*std::exp(-N*std::log(par.tau))*b/4.*(-1.+2.*II*par.v)*par.bslope*(1.+II)/x[1]
      *std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.TT(N,b,par.param));
  return 0;
}

int IT3(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverTotal par=*(InverTotal *)p;
  std::complex<long double> N=par.N0+par.Nslope*(1.+II)*std::log(x[0]);
  std::complex<long double> b=par.bc-par.bslope*(1.-II)*std::log(x[1]);
  std::complex<long double> theta=-II*par.v*M_PIl+x[2]*M_PIl*(1.+2.*II*par.v);
  
  y[0]=std::imag(par.Nslope*(1.+II)/x[0]/M_PIl*std::exp(-N*std::log(par.tau))*b/4.*(1.+2.*II*par.v)*par.bslope*(1.-II)/x[1]
      *std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.TT(N,b,par.param));
  return 0;
}

double integration::InverseMellinBessel(int method, std::complex<long double> (Func)(std::complex<long double> , std::complex<long double>, void*),
						     long double xp, long double tau, void *pp){
  double prec=1e-8;
  int fail; double *err=NULL, *prob=NULL;
  InverTotal ITC1, ITC2;
  long double Nslope=1.5;
  long double bslope=2.5;
  long double v=15.;
  long double bc=5.;
  ITC1.TT=Func;
  ITC2.TT=Func;
  ITC1.xp=xp;
  ITC2.xp=xp;
  ITC1.tau=tau;
  ITC2.tau=tau;
  ITC1.Nslope=Nslope;
  ITC2.Nslope=Nslope;
  ITC1.v=v;
  ITC2.v=v;
  ITC1.param=pp;
  ITC2.param=pp;
  ITC1.bc=bc;
  ITC2.bc=bc;
  
  ITC1.N0=3.;
  ITC2.N0=2.;
  
  double *ris1=NULL,*ris2=NULL,*ris3=NULL;
  CUBA(method,IT1,3,1,prec,&ris1,&fail,&err,&prob,&ITC1);
  CUBA(method,IT2,3,1,prec,&ris2,&fail,&err,&prob,&ITC1);
  CUBA(method,IT3,3,1,prec,&ris3,&fail,&err,&prob,&ITC2);
  
  //double result= ris1[0];
  double result= ris1[0]+ris2[0]+ris3[0];
  return result;

}


















