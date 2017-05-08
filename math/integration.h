#ifndef INTEGRATION_H
#define INTEGRATION_H


#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <cuba.h>
#include <vector>
#include <complex>
#include "complex_def.h"
#include <functional>
#include <complex_bessel.h>


using namespace sp_bessel;

namespace integration {
 
  //Integration GSL one dimension
  extern double GSL(double (func)(double , void *),double min, double max, double prec, double *err, void *par);

  //Integration CUBA one or more dimension VEGAS 0, SUAVE 1, DIVONNE 2, CUHRE 3 methods.
  //REMEMBER INTEGRATION LIMIT IS 0.0, 1.0 CHANGE VARIABLE
  extern int CUBA(int method, int (Func)(int*, double *, int*, double*, void*),
		     int ndim,int ncomp, double prec,double **res, int *fail, double **error, double **prob,void *par, int print=0);
  
  extern double InverseMellin_path(int method, std::complex<long double> (Func) ( std::complex<long double> , void *), long double tau , long double N0, long double slope, void *pp, long double *err);
  extern double InverseMellin_epsilon(int method, std::complex<long double> (Func) ( std::complex<long double> , void *), long double tau , long double N0, void *pp);
  
  extern std::complex<long double> InverseBessel_path(int method, std::complex<long double> (Func) (std::complex<long double> , void *), long double xp, long double bc, long double v, long double slope, void *pp);
  extern std::complex<long double> InverseBessel_epsilon(int method, std::complex<long double> (Func) (std::complex<long double> , void *), long double xp, long double bc, void *pp);

  extern int InverseMellinBessel(int method, std::complex<long double> (Func)(std::complex<long double> , std::complex<long double>, void*),
						     long double xp, long double tau, void *pp, double *ris, double *error, double *chi,bool verb,
				 long double N0=2.,int nump=10, double precint=1e-2, double epstart=5e-4, double eprate=5e-4, bool fit=true);
  
  extern int PolynomialFit(int ncoeff, int ndat, double *xpoint, double *ypoint, double *errypoint, double **coeff, double **errcoeff, double *chisq);
  //coeff and errcoeff are allocated inside
  
  extern int RichardsonExtrapolation(int ncoeff, double hstart, double *ypoint, double *ris, double *errris);
  //We apply a RichardsonExtrapolation of order ncoeff, evaluating the error as difference 
  //between the extrapolation of order ncoeff and of order ncoeff-1.
}





























#endif 