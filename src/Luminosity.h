#ifndef LUM_H
#define LUM_H

#include <LHAPDF/LHAPDF.h>
#include <cmath>
#include <iostream>
#include <gsl/gsl_chebyshev.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_zeta.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_dilog.h>
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <fstream>
#include <string>

class Luminosity {
	public:
		typedef std::complex<long double> dcomplex;
		typedef LHAPDF::PDF* PDF_ptr;
	
		Luminosity(PDF_ptr thePDF, double MUF, unsigned short NF = 5, std::size_t order=20);
		Luminosity(const Luminosity& Lumi);
		virtual ~Luminosity();
		
		void Cheb_Lum(double );
		void c_large();

		dcomplex Lum_N(dcomplex N, unsigned channel) const;
		long double Lum_x(long double z, unsigned channel) const;

		dcomplex Lum_gg_N(dcomplex N) const {	return Lum_N(N,0);	}
		long double Lum_gg_x(long double z) const { return Lum_x(z,0);  }
		
		dcomplex Lum_gq_N(dcomplex N) const {	return Lum_N(N,1);	}
		long double Lum_gq_x(long double z) const { return Lum_x(z,1);  }
		
		dcomplex Lum_qqb_N(dcomplex N) const {	return Lum_N(N,2);	}
		long double Lum_qqb_x(long double z) const { return Lum_x(z,2);  }
		
		dcomplex Lum_qq_N(dcomplex N) const {	return Lum_N(N,3);	}
		long double Lum_qq_x(long double z) const { return Lum_x(z,3);  }
		
		dcomplex Lum_qQb_N(dcomplex N) const {	return Lum_N(N,4);	}
		long double Lum_qQb_x(long double z) const {	return Lum_x(z,4);	}
		
		dcomplex Lum_qQ_N(dcomplex N) const {	return Lum_N(N,5);	}
		long double Lum_qQ_x(long double z) const { return Lum_x(z,5);  }
		

		double xfxQ(int channel, double x){
		  return (_PDF->xfxQ(channel,x,_Muf)/x);
		}
		
		inline PDF_ptr get_PDF(){
		  return _PDF;
		}
		inline double get_muf(){
		  return _Muf;
		}
		inline size_t get_n(){
		  return _n;
		}
		
		struct par_struct{
			double Muf;
			PDF_ptr PDF;
			double t;
			unsigned int Nf;
		};
		
		void setNf(unsigned short NF){
		  if ( (NF < 1) || (NF > 6) ){
		    std::cout << "ERROR: INCONGRUENT NUMBER OF FLAVOUR (1-6 ACCEPTED)" << std::endl;
		    NF=5.;
		  }
		  _Nf=NF;
		}
		
		//Test Function
		void TestLum(std::string nameout);

	private:
		size_t _n;	//Order of Chebishev approximation
		std::array< gsl_cheb_series*,6>  Chebseries;
		std::array< std::vector<double> ,6>  C_larges;

		std::vector< std::vector< double > > _T;
		LHAPDF::PDF* _PDF;
		
		double _Muf;
		unsigned short _Nf;

};

#endif
