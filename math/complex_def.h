#ifndef __COMPLEXDEF_H__
#define __COMPLEXDEF_H__

#include <iostream>
#include <complex>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_dilog.h>
#include <gsl/gsl_sf_result.h>
#include <complex_bessel.h>

using namespace sp_bessel;



std::complex<long double> operator*(const std::complex<long double>& a, const double& b);

 std::complex<long double> operator*(const double& a,const std::complex<long double>& b);

 std::complex<long double> operator/(const std::complex<long double>& a, const double& b);

 std::complex<long double> operator/(const double& a,const std::complex<long double>& b);

 std::complex<long double> operator-(const std::complex<long double>& a, const double& b);

 std::complex<long double> operator-(const double& a,const std::complex<long double>& b);

 std::complex<long double> operator+(const std::complex<long double>& a, const double& b);

 std::complex<long double> operator+(const double& a,const std::complex<long double>& b);

 bool operator == (const std::complex<long double> &z,const double &a);

 bool operator != (const std::complex<long double> &z,const double &a);

 bool operator == (const double &a,const std::complex<long double> &z);

 bool operator != (const double &a,const std::complex<long double> &z);

//operation natural numbers
 std::complex<long double> operator + (const std::complex<long double> &z,const int n);

 std::complex<long double> operator - (const std::complex<long double> &z,const int n);

 std::complex<long double> operator * (const std::complex<long double> &z,const int n);

 std::complex<long double> operator / (const std::complex<long double> &z,const int n);

 std::complex<long double> operator + (const int n,const std::complex<long double> &z);

 std::complex<long double> operator - (const int n,const std::complex<long double> &z);

 std::complex<long double> operator * (const int n,const std::complex<long double> &z);

 std::complex<long double> operator / (const int n,const std::complex<long double> &z);

 std::complex<long double> operator + (const std::complex<long double> &z,const unsigned int n);

 std::complex<long double> operator - (const std::complex<long double> &z,const unsigned int n);

 std::complex<long double> operator * (const std::complex<long double> &z,const unsigned int n);

 std::complex<long double> operator / (const std::complex<long double> &z,const unsigned int n);

 std::complex<long double> operator + (const unsigned int n,const std::complex<long double> &z);

 std::complex<long double> operator - (const unsigned int n,const std::complex<long double> &z);

 std::complex<long double> operator * (const unsigned int n,const std::complex<long double> &z);

 std::complex<long double> operator / (const unsigned int n,const std::complex<long double> &z);

 bool operator == (const std::complex<long double> &z,const int n);

 bool operator != (const std::complex<long double> &z,const int n);

 bool operator == (const int n,const std::complex<long double> &z);

 bool operator != (const int n,const std::complex<long double> &z);

 bool operator == (const std::complex<long double> &z,const unsigned int n);

 bool operator != (const std::complex<long double> &z,const unsigned int n);

 bool operator == (const unsigned int n,const std::complex<long double> &z);

 bool operator != (const unsigned int n,const std::complex<long double> &z);

std::complex<long double> pow(const long double& a, const std::complex<long double> &b);
std::complex<long double> pow(const std::complex<long double>& a, const long double &b);
std::complex<long double> pow(const std::complex<long double>& a, const std::complex<long double> &b);

 long double inf_norm (const std::complex<long double> &z);
 bool isfinite (const std::complex<long double> &z);

//Useful Math functions
long double log_r(long double x);
long double dilog_r(long double x);

const std::complex<long double> II(0.0,1.0);
std::complex<long double> log_c(std::complex<long double> z);
std::complex<long double> dilog_c(std::complex<long double> z);

std::complex<long double> log_c_angle(std::complex<long double> z, long double angle);
std::complex<long double> expm1 (const std::complex<long double> &z);
std::complex<long double> log1p (const std::complex<long double> &z);

std::complex<long double> LogGamma(std::complex<long double> z);
std::complex<long double> Gamma_inv (const std::complex<long double> &z);


std::complex<long double> LBesselJ(long double k, std::complex<long double> z);
std::complex<long double> LBesselK(long double k, std::complex<long double> z);
std::complex<long double> CBesselK(std::complex<long double> nu, std::complex<long double> z, bool warn=true);

//Works only if abs(z)<1
std::complex<long double> Hyp2F1 (std::complex<long double> a,std::complex<long double> b, std::complex<long double> c,std::complex<long double> z);

#endif
