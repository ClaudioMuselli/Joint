#include "integration.h"
#include <boost/concept_check.hpp>
 
 
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


int integration::CUBA( int method, int (Func)(int*, double *, int *, double *, void *), int ndim, int ncomp, 
		       double prec, double **res,int *fail, double **error, double **prob, void * par, int print){
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
      if (print == 0){
      std::cout << "Call to VEGAS" << std::endl;
      }
      Vegas(ndim,ncomp,func,par,nvec,prec, epsabs,verbose,seed,mineval,maxeval,nstart,
	    nincrease,nbatch,gridno,statefile,NULL,&neval, fail,ris,err,proba);
      break;
    case 1: 
      if (print == 0){
      std::cout << "Call to SUAVE" << std::endl;
      }
      Suave(ndim, ncomp, func, par,nvec, prec, epsabs, verbose, seed, mineval, maxeval, nnew, 5,
	    flatness,statefile,NULL, &nregions, &neval, fail, ris, err, proba);
      break;
    case 2:
      if (print == 0){
      std::cout << "Call to DIVONNE" << std::endl;
      }
      Divonne(ndim, ncomp, func, par,nvec, prec, epsabs, verbose, seed, mineval, maxeval, key1, 
	      key2, key3, maxpass, border, maxchisq, mindeviation, ngiven, ldxgiven, NULL, 
	      nextra, NULL,statefile,NULL, &nregions, &neval, fail, ris, err, proba);
      break;
    case 3: 
      if (print == 0){
      std::cout << "Call to CUHRE" << std::endl;
      }
      Cuhre(ndim, ncomp, func, par,nvec, prec, epsabs, verbose | last, mineval, maxeval, key, statefile, NULL,
	    &nregions, &neval, fail, ris, err, proba);
      break;
  }
  *error=err;
  *prob=proba;
  *res=ris;
  return 0;
}

int integration::PolynomialFit(int ncoeff, int ndat, double *xpoint, double *ypoint, 
			       double *errypoint, double **coeff, double **errcoeff, double *chisq){
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;
  X= gsl_matrix_alloc(ndat, ncoeff);
  y=gsl_vector_alloc(ndat);
  w=gsl_vector_alloc(ndat);
  
  c=gsl_vector_alloc(ncoeff);
  cov=gsl_matrix_alloc(ncoeff,ncoeff);
  
  double *errco=NULL, *co=NULL;
  
  errco=new double[ncoeff];
  co=new double[ncoeff];
  //std::cout << "CIAO" << std::endl;
  
  for( int i=0;i < ndat; i++){
    for (int j=0;j<ncoeff;j++){
      gsl_matrix_set(X,i,j, std::pow(xpoint[i],j));
    }
    gsl_vector_set(y,i,ypoint[i]);
    gsl_vector_set(w,i,(1.0/(errypoint[i]*errypoint[i])));
  }
  //std::cout << "CIAO" << std::endl;
  //Fitting with gsl
  gsl_multifit_linear_workspace *work =gsl_multifit_linear_alloc(ndat,ncoeff);
  gsl_multifit_wlinear(X,w,y,c,cov,chisq,work);
  gsl_multifit_linear_free(work);
  //std::cout << "CIAO" << std::endl;
  for( int i=0;i< ncoeff;i++){
    co[i]=gsl_vector_get(c,i);
    errco[i]=std::sqrt(gsl_matrix_get(cov,i,i));
  }
  *coeff=co;
  *errcoeff=errco;
  *chisq/=((double)(ndat-ncoeff));
  return 0;
}

int integration::RichardsonExtrapolation(int ncoeff, double hstart, double* ypoint, double* ris, double* errris) {
  double res=0.,errres=0.,res2=0.;
  std::vector<std::vector<double>> A1,A2;
  std::vector<double> yv;
  if (ncoeff==1){
    *ris=ypoint[0];
    *errris=hstart;
    return 0;
  }
  for (int i=0;i<ncoeff;i++){
    yv.push_back(ypoint[i]);
  }
  A1.push_back(yv);
  for (int i=1;i<ncoeff;i++){
    std::vector<double> mpoint;
    for (int j=0;j<ncoeff-i;j++){
      mpoint.push_back((std::pow(2.,i)*A1[i-1][j]-A1[i-1][j+1])/(std::pow(2.,i)-1.));
    }
    A1.push_back(mpoint);
  }
  res=A1[ncoeff-1.][0];
  for (int i=1;i<ncoeff;i++){
    yv.push_back(ypoint[i]);
  }
  A2.push_back(yv);
  for (int i=1;i<ncoeff-1;i++){
    std::vector<double> mpoint;
    for (int j=0;j<ncoeff-i-1;j++){
      mpoint.push_back((std::pow(2,i)*A2[i-1][j]-A2[i-1][j+1])/(std::pow(2.,i)-1.));
    }
    A2.push_back(mpoint);
  }
  res2=A2[ncoeff-2][0];
  *ris=res;
  double errmax=std::pow(hstart*std::pow(2.,ncoeff-1.),ncoeff);
  if (hstart*std::pow(2.,ncoeff-1.)>1.){
    std::cout << "WARNING: error may be not attendible; decrease hstart" << std::endl;
  }
  *errris=std::max(errmax,std::abs(res-res2));
  return 0;
  
  
  
}





struct InverseTotal{
   long double xp;
   long double tau;
   std::function<std::complex<long double> ( std::complex<long double>, std::complex<long double>, void *)> TT;
   void *param;
   long double N0;
   long double bc;
   long double epsilon;
};

int IT(int* ndim, double * x, int* ncomp, double* y, void *p){
  InverseTotal par=*(InverseTotal *)p;
  //long double z =(long double) (x[2]/(1.-x[2]));
  long double br=(long double) (x[1]/(1.-x[1]));
  long double yN=(long double) (x[0]/(1.-x[0]));
  std::complex<long double> N=par.N0+II*yN;
  std::complex<long double> b=br+par.bc*II;  
  std::complex<long double> y0= std::exp(-N*std::log(par.tau*std::pow(std::sqrt(1.+par.xp)-std::sqrt(par.xp),2)))*b/2.*LBesselJ(0.,b*std::sqrt(par.xp))
  *(std::exp(-par.epsilon*std::pow(br*br+yN*yN,1)))*par.TT(N,b,par.param)/std::pow((1.-x[0])*(1.-x[1]),2);
  y[0]= (double) (std::real(y0))/(std::pow(M_PI,1.));
}

int integration::InverseMellinBessel(int method, std::complex<long double> (Func)(std::complex<long double> , std::complex<long double>, void*),
						     long double xp, long double tau, void *pp, double *ris, double *error, double *chi,bool verb,
						     long double N0, int nump, double precint, double epstart, double eprate, bool fit){
  InverseTotal ITC;
  ITC.xp=xp;
  ITC.tau=tau;
  ITC.TT=Func;
  ITC.N0=N0;
  ITC.bc=0.;
  ITC.param=pp;
  int nn=nump;
  bool verbose=verb;
  double *epsilon=NULL,**res=NULL, prec=1e-14, *yy=NULL, *erryy=NULL;
  int fail; double **err=NULL, **prob=NULL;
  res=new double*[nn];
  err=new double*[nn];
  prob=new double*[nn];
  for (int i=0;i<nn; i++){
    res[i]=new double[1];
    err[i]=new double[1];
    prob[i]=new double[1];
  }
  epsilon=new double[nn];
  double epsilonstart=epstart,rate=eprate;
  epsilon[0]=epsilonstart;
  if (fit){
  for (int i=0;i<nn;i++){
    ITC.epsilon=epsilon[i];
    CUBA(method,IT,2,1,prec,&res[i],&fail,&err[i],&prob[i],&ITC);
    if (std::abs(err[i][0]/res[i][0]) > precint){
      i--;
      epsilon[i+1]+=rate;
    }
   else  epsilon[i+1]=epsilon[i]+rate;
  }
  }
  else{
    for (int i=0;i<nn;i++){
      ITC.epsilon=epsilon[i];
      CUBA(method,IT,2,1,prec,&res[i],&fail,&err[i],&prob[i],&ITC);
      if (std::abs(err[i][0]) > precint){
	std::cout << "WARNING: error out of precision for ep= " << epsilon[i] << std::endl;
      }
      epsilon[i+1]=epsilon[i]*2.;
    }
  }
  //print per prova poi cancellare
  yy=new double[nn];
  erryy=new double[nn];
  for (int i=0;i<nn;i++){
    yy[i]=res[i][0];
    erryy[i]=err[i][0];
    if (verbose){
    std::cout << "epsilon= " << epsilon[i] << " integral= " << yy[i] << " +- "<< erryy[i] << std::endl;
    }
  }
  if (fit){
  double *co=NULL,*errco=NULL,chisq;
  integration::PolynomialFit(3,nn,epsilon,yy,erryy, &co, &errco,&chisq);
  if (verbose){
  std::cout << co[0] << " + " << co[1] << " x + " << co[2] << " x^2 " << std::endl;
  std::cout << errco[0] << " " << errco[1] << " " << errco[2] << std::endl;
  std::cout << chisq << std::endl;
  }
  *ris=co[0];
  *error=errco[0];
  *chi=chisq;
  return 0; 
  }
  else {
    double res=0., errres=0.;
    integration::RichardsonExtrapolation(nump,epsilon[0],yy,&res,&errres);
    *ris=res;
    *error=errres;
    *chi=0.;
    return 0;
  }
}





























//Use For INTEGRATION
struct InverMellin{
  std::function<std::complex<long double> (std::complex<long double>, void*)> PP;
  long double N0;
  long double epsilon;
  long double tau;
  void *param;
};

struct InverBessel{
  std::function<std::complex<long double> (std::complex<long double>, void*)> BB;
  long double bc;
  long double bepsilon;
  long double xp;
  long double v;
  void *param;
};

//MELLIN INVERSE
int IM1(int* ndim, double * x, int* ncomp, double* y, void *p){
  InverMellin par= *(InverMellin *)p;
  std::complex<long double> N=par.N0+par.epsilon*(1.+II)*std::log(x[0]);
  
  y[0]=std::imag(par.epsilon*(1.+II)*std::exp(-N*std::log(par.tau))/M_PI/x[0]*par.PP(N,par.param));
  return 0;
}

int IMep(int* ndim, double* x, int* ncomp, double* y,void *p){
  InverMellin par= *(InverMellin *)p;
  std::complex<long double> N=par.N0-II/(1.-x[0]);
  y[0]=std::real(std::exp(-par.epsilon/std::pow(1.-x[0],2))*std::exp(-N*std::log(par.tau))/M_PI/std::pow(1.-x[0],2)*par.PP(N,par.param));
  std::complex<long double> N2=par.N0-II*x[1];
  y[1]=std::real(std::exp(-par.epsilon*x[1])*std::exp(-N2*std::log(par.tau))/M_PI*par.PP(N2,par.param));
  return 0;
}
  

double integration::InverseMellin_path(int method, std::complex<long double> (Func)(std::complex<long double> , void* ), 
				  long double tau, long double N0, long double slope, void *pp, long double *error) {
  double *res=NULL, prec=1e-8;
  int fail; double *err=NULL, *prob=NULL;
  res=new double[1];
  err=new double[1];
  prob=new double[1];
  InverMellin IMC;
  IMC.PP=Func;
  IMC.N0=N0;
  IMC.epsilon=slope;
  IMC.tau=tau;
  IMC.param=pp;
  CUBA(method,IM1,2,1,prec,&res,&fail,&err,&prob,&IMC);
  *error=err[0];
  return res[0];
}

double integration::InverseMellin_epsilon(int method,std::complex<long double> (Func)(std::complex<long double> , void* ), 
				  long double tau, long double N0,void *pp){
  int nn=10;
  double *epsilon=NULL,**res=NULL, prec=1e-8, *yy=NULL, *erryy=NULL;
  int fail; double **err=NULL, **prob=NULL;
  res=new double*[nn];
  err=new double*[nn];
  prob=new double*[nn];
  for (int i=0;i<nn; i++){
    res[i]=new double[2];
    err[i]=new double[2];
    prob[i]=new double[2];
  }
  InverMellin IMC;
  IMC.PP=Func;
  IMC.N0=N0;
  IMC.tau=tau;
  IMC.param=pp;
  epsilon=new double[nn];
  double epsilonstart=1e-5,rate=1e-5;
  epsilon[0]=epsilonstart;
  for (int i=0;i<nn;i++){
    IMC.epsilon=epsilon[i];
    CUBA(method,IMep,2,2,prec,&res[i],&fail,&err[i],&prob[i],&IMC);
    if ((err[i][0]+err[i][1]) >1e-4){
      i--;
      epsilon[i+1]+=rate;
    }
    else epsilon[i+1]=epsilon[i]+rate;
  }
  //print per prova poi cancellare
  yy=new double[nn];
  erryy=new double[nn];
  for (int i=0;i<nn;i++){
    yy[i]=res[i][0]+res[i][1];
    erryy[i]=err[i][0]+err[i][1];
    std::cout << "epsilon= " << epsilon[i] << " integral= " << yy[i] << " +- "<< erryy[i] << std::endl;
    
  }
  double *co=NULL,*errco=NULL,chisq;
  integration::PolynomialFit(3,nn,epsilon,yy,erryy, &co, &errco,&chisq);
  std::cout << co[0] << " + " << co[1] << " x + " << co[2] << " x^2 " << std::endl;
  std::cout << errco[0] << " " << errco[1] << " " << errco[2] << std::endl;
  std::cout << chisq << std::endl;
  return co[0];
}

int IB3(int* ndim, double *x, int* ncomp, double* y,void* p){
  InverBessel par=*(InverBessel *)p;
  std::complex<long double> b=par.bc*x[0];
  y[0]=std::real(par.bc*b/2.*(std::complex<long double>)(sp_bessel::besselJ(0.,(std::complex<double>)(b*std::sqrt(par.xp))))*par.BB(b,par.param));
  y[1]=std::imag(par.bc*b/2.*(std::complex<long double>)(sp_bessel::besselJ(0.,(std::complex<double>)(b*std::sqrt(par.xp))))*par.BB(b,par.param));
}

int IB1(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverBessel par=*(InverBessel *) p;
  std::complex<long double> b=par.bc-par.bepsilon*(1.+II)*std::log(x[0]);
  std::complex<long double> theta=-II*par.v*M_PIl+x[1]*M_PIl*(-1.+2.*II*par.v);
  
  y[0]=std::real(-b/4.*(-1.+2.*II*par.v)*par.bepsilon*(1.+II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  y[1]=std::imag(-b/4.*(-1.+2.*II*par.v)*par.bepsilon*(1.+II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  return 0;
  
}

int IB2(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverBessel par=*(InverBessel *) p;
  std::complex<long double> b=par.bc-par.bepsilon*(1.-II)*std::log(x[0]);
  std::complex<long double> theta=-II*par.v*M_PIl+x[1]*M_PIl*(1.+2.*II*par.v);
  
  y[0]=std::real(b/4.*(1.+2.*II*par.v)*par.bepsilon*(1.-II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  y[1]=std::imag(b/4.*(1.+2.*II*par.v)*par.bepsilon*(1.-II)/x[0]*std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.BB(b,par.param));
  return 0;
}

std::complex<long double> integration::InverseBessel_path(int method, std::complex<long double> (Func)(std::complex<long double> , void*),
						     long double xp, long double bc, long double v, long double slope, void *pp, std::complex<long double> *err){
  double prec=1e-12;
  int fail; double *err1,*err2,*err3, *prob;
  InverBessel IBC;
  IBC.BB=Func;
  IBC.v=v;
  IBC.bepsilon=slope;
  IBC.xp=xp;
  IBC.param=pp;
  IBC.bc=bc;
  double *ris1,*ris2,*ris3;
  ris1=new double[2];
  ris2=new double[2];
  err1=new double[2];
  err2=new double[2];
  err3=new double[2];
  prob=new double[2];
  CUBA(method,IB1,2,2,prec,&ris1,&fail,&err1,&prob,&IBC);
  CUBA(method,IB2,2,2,prec,&ris2,&fail,&err2,&prob,&IBC);
  CUBA(method,IB3,2,2,prec,&ris3,&fail,&err3,&prob,&IBC);
  std::complex<long double> result;
  *err=err1[0]+ err2[0]+err3[0]+II*(err1[1]+err2[1]+err3[1]);
  result=ris1[0]+ris2[0]+ris3[0]+II*(ris1[1]+ris2[1]+ris3[1]);
  return(result);  
}
 
 
int IBep(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverBessel par=*(InverBessel *)p;
  std::complex<long double> b=1./(1.-x[0])+par.bc*II;
  y[0]=std::real(std::exp(-par.bepsilon/std::pow(1.-x[0],2))*b/2.*LBesselJ(0.,b*std::sqrt(par.xp))*par.BB(b,par.param)/std::pow(1.-x[0],2));
  y[1]=std::imag(std::exp(-par.bepsilon/std::pow(1.-x[0],2))*b/2.*LBesselJ(0.,b*std::sqrt(par.xp))*par.BB(b,par.param)/std::pow(1.-x[0],2));
  std::complex<long double> b2= x[1]*(1.+par.bc*II);
  y[2]=std::real(std::exp(-par.bepsilon*x[1]*x[1])*b2/2.*LBesselJ(0.,b2*std::sqrt(par.xp))*par.BB(b2,par.param)*(1.+par.bc*II));
  y[3]=std::imag(std::exp(-par.bepsilon*x[1]*x[1])*b2/2.*LBesselJ(0.,b2*std::sqrt(par.xp))*par.BB(b2,par.param)*(1.+par.bc*II));
  
}

std::complex<long double> integration::InverseBessel_epsilon(int method,std::complex<long double> (Func)(std::complex<long double> , void* ), 
				  long double xp, long double bc,void *pp){
  int nn=10;
  double *epsilon=NULL,**res=NULL, prec=1e-12;
  double *yyreal=NULL, *erryyreal=NULL, *yyimag=NULL, *erryyimag=NULL;
  int fail; double **err=NULL, **prob=NULL;
  res=new double*[nn];
  err=new double*[nn];
  prob=new double*[nn];
  for (int i=0;i<nn; i++){
    res[i]=new double[4];
    err[i]=new double[4];
    prob[i]=new double[4];
  }
  InverBessel IBC;
  IBC.BB=Func;
  IBC.bc=bc;
  IBC.xp=xp;
  IBC.param=pp;
  epsilon=new double[nn];
  double epsilonstart=1e-7,rate=1e-7;
  epsilon[0]=epsilonstart;
  for (int i=0;i<nn;i++){
    IBC.bepsilon=epsilon[i];
    CUBA(method,IBep,2,4,prec,&res[i],&fail,&err[i],&prob[i],&IBC);
    if (((err[i][0]+err[i][2]) >1e-6)||((err[i][1]+err[i][3]) >1e-6)){
      i--;
      epsilon[i+1]+=rate;
    }
    else epsilon[i+1]=epsilon[i]+rate;
  }
  //print per prova poi cancellare
  yyreal=new double[nn];
  erryyreal=new double[nn];
  yyimag=new double[nn];
  erryyimag=new double[nn];
  for (int i=0;i<nn;i++){
    yyreal[i]=res[i][0]+res[i][2];
    erryyreal[i]=err[i][0]+err[i][2];
    yyimag[i]=res[i][1]+res[i][3];
    erryyimag[i]=err[i][1]+err[i][3];
    std::cout << "epsilon= " << epsilon[i] << " integral= " << yyreal[i]+II*yyimag[i] << " +- "<< erryyreal[i]+II*erryyimag[i] << std::endl;
    
  }
  double *coreal=NULL,*errcoreal=NULL,chisqreal, *coimag=NULL, *errcoimag=NULL, chisqimag;
  integration::PolynomialFit(3,nn,epsilon,yyreal,erryyreal, &coreal, &errcoreal,&chisqreal);
  integration::PolynomialFit(3,nn,epsilon,yyimag,erryyimag, &coimag, &errcoimag,&chisqimag);
  std::cout << coreal[0]+II*coimag[0] << " + " << coreal[1]+II*coimag[1] << " x + " << coreal[2]+II*coimag[2] << " x^2 " << std::endl;
  std::cout << errcoreal[0]+II*errcoimag[0] << " " <<errcoreal[1]+II*errcoimag[1] << " " << errcoreal[2]+II*errcoimag[2]<< std::endl;
  std::cout << chisqreal << " " << chisqimag << std::endl;
  std::complex<long double> ris= coreal[0]+II*coimag[0];
  return (ris);
}

/*

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
  //std::cout << par.v << " " << theta << std::endl;
  
  y[0]=std::imag(par.Nslope*(1.-II)/x[0]/M_PIl*std::exp(-N*std::log(par.tau))*b/4.*(-1.+2.*II*par.v)*par.bslope*(1.+II)/x[1]
      *std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.TT(N,b,par.param));
  //std::cout << y[0] << std::endl;
  return 0;
}

int IT3(int* ndim, double *x, int* ncomp, double* y, void* p){
  InverTotal par=*(InverTotal *)p;
  std::complex<long double> N=par.N0+par.Nslope*(1.+II)*std::log(x[0]);
  std::complex<long double> b=par.bc-par.bslope*(1.-II)*std::log(x[1]);
  std::complex<long double> theta=-II*par.v*M_PIl+x[2]*M_PIl*(1.+2.*II*par.v);
  
  y[0]=std::imag(par.Nslope*(1.+II)/x[0]/M_PIl*std::exp(-N*std::log(par.tau))*b/4.*(1.+2.*II*par.v)*par.bslope*(1.-II)/x[1]
      *std::exp(-II*b*std::sqrt(par.xp)*std::sin(theta))*par.TT(N,b,par.param));
  //std::cout << y[0] << std::endl;
  return 0;
}

double integration::InverseMellinBessel(int method, std::complex<long double> (Func)(std::complex<long double> , std::complex<long double>, void*),
						     long double xp, long double tau, void *pp){
  double prec=1e-8;
  int fail; double *err=NULL, *prob=NULL;
  InverTotal ITC1, ITC2;
  long double Nslope=1.5;
  long double bslope=2.5;
  long double v=5.;
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
  //CUBA(method,IT1,3,1,prec,&ris1,&fail,&err,&prob,&ITC1);
  CUBA(method,IT2,3,1,prec,&ris2,&fail,&err,&prob,&ITC1);
  //CUBA(method,IT3,3,1,prec,&ris3,&fail,&err,&prob,&ITC2);
  
  double result= ris2[0];
  //double result= ris1[0]+ris2[0]+ris3[0];
  return result;

}

*/
















