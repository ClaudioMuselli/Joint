#ifndef INTEGRATION_H
#define INTEGRATION_H


#include <iostream>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>
#include <cuba.h>
#include <vector>
#include <complex>
#include "complex_def.h"
#include <functional>

namespace integration {
 
  //Integration GSL one dimension
  extern double GSL(double (func)(double , void *),double min, double max, double prec, double *err, void *par);

  //Integration CUBA one or more dimension VEGAS 0, SUAVE 1, DIVONNE 2, CUHRE 3 methods.
  //REMEMBER INTEGRATION LIMIT IS 0.0, 1.0 CHANGE VARIABLE
  extern int CUBA(int method, int (Func)(int*, double *, int*, double*, void*),
		     int ndim,int ncomp, double prec,double **res, int *fail, double **error, double **prob,void *par);
  
  extern double InverseMellin(int method, std::complex<long double> (Func) ( std::complex<long double> , void *), long double tau , long double N0, long double slope, bool sign, void *pp);
  
  extern std::complex<long double> InverseBessel(int method, std::complex<long double> (Func) (std::complex<long double> , void *), long double xp, long double v, long double slope, void *pp); 

}


#endif 