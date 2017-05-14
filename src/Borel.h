#ifndef __BOREL_H__
#define __BOREL_H__

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_math.h>
#include "../math/integration.h"

using namespace std;


namespace Borel{
 int BorelJointC(int method, long double x, long double xp, long double asbar, std::complex<long double> (Func)(std::complex<long double>, std::complex<long double>, void*),
		 void *params,long double *ris, long double *err,long double Ccutoff, long double N0=3.,long double slope=1.5);
 int BorelJointCp(int method, long double x, long double xp, long double asbar, std::complex<long double> (Func)(std::complex<long double>, std::complex<long double>, void*),
		void *params,long double *ris, long double *err,long double Ccutoff, long double N0=3.,long double slope=1.5);
 int BorelJointpt(int method, long double x, long double xp, long double asbar, std::complex<long double> (Func)(std::complex<long double>, std::complex<long double>, void*),
		void *params,long double *ris, long double *err,long double Ccutoff, long double N0=3.,long double slope=1.5);
  
}
































































#endif