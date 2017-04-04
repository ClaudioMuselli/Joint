#include "complex_def.h"

std::complex<long double> operator*(const std::complex<long double>& a, const double& b){
	return a*static_cast<long double>(b);
}

std::complex<long double> operator*(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)*b;
}

std::complex<long double> operator/(const std::complex<long double>& a, const double& b){
	return a/static_cast<long double>(b);
}

std::complex<long double> operator/(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)/b;
}

std::complex<long double> operator-(const std::complex<long double>& a, const double& b){
	return a-static_cast<long double>(b);
}

std::complex<long double> operator-(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)-b;
}

std::complex<long double> operator+(const std::complex<long double>& a, const double& b){
	return a+static_cast<long double>(b);
}

std::complex<long double> operator+(const double& a,const std::complex<long double>& b){
	return static_cast<long double>(a)+b;
}

std::complex<long double> pow(const long double& a, const std::complex<long double> &b){
    std::complex<long double> A(a,0.);
    return(std::exp(b*std::log(A)));
  }
  std::complex<long double> pow(const std::complex<long double>& a, const long double &b){
    std::complex<long double> B(b,0.);
    return(std::exp(B*std::log(a)));
  }
  std::complex<long double> pow(const std::complex<long double>& a, const std::complex<long double> &b){
    return(std::exp(b*std::log(a)));
  }

//Useful Math functions
long double log_r(long double x){
//	std::cout << "log_r(" << x << ")" << std::endl;
	return gsl_sf_log(x);
}
long double dilog_r(long double x){
//	std::cout << "dilog_r(" << x << ")" << std::endl;
	return gsl_sf_dilog(x);
}

std::complex<long double> log_c(std::complex<long double> z){
	gsl_sf_result lnr;
	gsl_sf_result theta;
	gsl_sf_complex_log_e(std::real(z),std::imag(z),&lnr, &theta);
	std::complex<long double> ris(lnr.val,theta.val);
	return ris;
}
std::complex<long double> dilog_c(std::complex<long double> z){
	gsl_sf_result ris_r;
	gsl_sf_result ris_i;
	gsl_sf_complex_dilog_e(std::abs(z),std::arg(z),&ris_r, &ris_i);
	std::complex<long double> ris(ris_r.val,ris_i.val);
	return ris;
}

std::complex<long double> log_c_angle(std::complex<long double> z, long double angle){
  long double lnr=std::abs(z);
  long double theta=std::arg(z);
  long double epsilon=1e-20;
  if ((theta> (-M_PIl-epsilon))&&(theta < (angle+epsilon))){ 
    theta+=2.*M_PIl+angle;
  }
  else theta+=angle;
  std::complex<long double> ris(std::log(lnr),theta);
  return ris;
}

