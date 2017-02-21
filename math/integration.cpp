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


double integration::CUBA( int method, int (Func)(int*, double *, int *, double *, void *), int ndim, double prec, int *fail, double error, double prob, void * par){
  static const double epsabs=1e-6;
  static const int verbose=2;
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
  static const int ncomp=1;
  static const int nvec=1;
  static const int spin=-1;
  
  int neval, nregions;
  double *ris=NULL, *err=NULL, *proba=NULL;
  ris=new double[ncomp];
  err=new double[ncomp];
  proba=new double[ncomp];
  integrand_t func = (integrand_t) Func;
  double res;
  
  switch (method){
    case 0:
      std::cout << "Call to VEGAS" << std::endl;
      Vegas(ndim,ncomp,func,par,nvec,prec, epsabs,verbose,seed,mineval,maxeval,nstart,
	    nincrease,nbatch,gridno,statefile,NULL,&neval, fail,ris,err,proba);
      res=ris[0];
      if(error) error=err[0];
      if(prob) prob=proba[0];
      break;
    case 1: 
      std::cout << "Call to SUAVE" << std::endl;
      Suave(ndim, ncomp, func, par,nvec, prec, epsabs, verbose, seed, mineval, maxeval, nnew, 5,
	    flatness,statefile,NULL, &nregions, &neval, fail, ris, err, proba);
      res=ris[0];
      if(error) error=err[0];
      if(prob) prob=proba[0];
      break;
    case 2:
      std::cout << "Call to DIVONNE" << std::endl;
      Divonne(ndim, ncomp, func, par,nvec, prec, epsabs, verbose, seed, mineval, maxeval, key1, 
	      key2, key3, maxpass, border, maxchisq, mindeviation, ngiven, ldxgiven, NULL, 
	      nextra, NULL,statefile,NULL, &nregions, &neval, fail, ris, err, proba);
      res=ris[0];
      if(error) error=err[0];
      if(prob) prob=proba[0]; 
      break;
    case 3: 
      std::cout << "Call to Cuhre" << std::endl;
      Cuhre(ndim, ncomp, func, par,nvec, prec, epsabs, verbose | last, mineval, maxeval, key, statefile, NULL,
	    &nregions, &neval, fail, ris, err, proba);
      res=ris[0];
      if(error) error=err[0];
      if(prob) prob=proba[0];
      break;
  }
  delete[] ris;
  delete[] err;
  delete[] proba;
  return res;
}