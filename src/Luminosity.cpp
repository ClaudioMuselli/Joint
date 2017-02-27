#include "Luminosity.h"
#include "complex_def.h"
#include <functional>


Luminosity::Luminosity(PDF_ptr thePDF, double MUF,  unsigned short NF, std::size_t n):_n(n),_Nf(NF),_Muf(MUF),_PDF(thePDF)
{	
	for (auto& it : Chebseries)
		it = gsl_cheb_alloc(_n);
	for (auto& it : C_larges)	
		it = std::vector < double >( _n+1 , 0 );

	for(std::size_t i = 0; i < (_n + 1) ; ++i)
		_T.push_back( std::vector < double >(_n+1,0) );

	_T[0][0]=1.0;
	_T[1][0]=0.0;
	_T[1][1]=1.0;
	for(std::size_t i = 2 ; i < (_n+1) ; i++){
		_T[i][0] = - _T[i-2][0];
		for(std::size_t j = 1 ; j < (i-1) ; j++){
			_T[i][j] = 2*_T[i-1][j-1] - _T[i-2][j];
		}
		_T[i][i-1] = 2*_T[i-1][i-2];
		_T[i][i] = 2*_T[i-1][i-1];
	}
	Cheb_Lum(MUF);

}
Luminosity::~Luminosity(){
 	for (auto& it : Chebseries)
		gsl_cheb_free(it);
}

//CORE FUNCTION (FOR INTERPOLATION)

//gg channel (0)

double f_gg(double x1,void *p){
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf; 
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
	return ( PDF-> xfxQ(0,x1,Muf)/x1 * PDF -> xfxQ(0,x2,Muf)/x2 ) / x1;
}


double _Lum_gg(double tau,void* p){
	Luminosity::par_struct* par = (Luminosity::par_struct*) p;
	par -> t = tau;
	double precision = 1e-4;
	double ris = 0.0, error = 0.0;
	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	Integrand.function = f_gg;
	Integrand.params = par;
	gsl_integration_qags (&Integrand, tau, 1.0, 0, precision, 1000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	return ris;
}

double _Lgg(double x, void *p){
	double z = std::exp(x);
	return (z*_Lum_gg(z,p));
}

//gq channel (1)

double f_gq(double x1,void *p){
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf;
	unsigned int Nf=p_struct-> Nf;
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
	double fS=0.;
	for (int i=0;i<Nf;i++){
	  fS+=(PDF->xfxQ(0,x1,Muf)/x1*(PDF->xfxQ(i+1,x2,Muf))/x2+PDF->xfxQ(0,x1,Muf)/x1*PDF-> xfxQ(-(i+1),x2,Muf)/x2);
	}
	return 2.*fS/x1;
}


double _Lum_gq(double tau,void* p){
	Luminosity::par_struct* par = (Luminosity::par_struct*) p;
	par -> t = tau;
	double precision = 1e-4;
	double ris = 0.0, error = 0.0;
	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	Integrand.function = f_gq;
	Integrand.params = par;
	gsl_integration_qags (&Integrand, tau, 1.0, 0, precision, 1000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	return ris;
}

double _Lgq(double x, void *p){
	double z = std::exp(x);
	return (z*_Lum_gq(z,p));
}

//qqbar channel (2)

double f_qqb(double x1,void *p){
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf; 
	unsigned int Nf=p_struct->Nf;
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
	double fqqbar=0.;
	for (int i=0;i<Nf;i++){
	  fqqbar+= PDF-> xfxQ((i+1),x1,Muf)/x1 * PDF -> xfxQ(-(i+1),x2,Muf)/x2;
	}
	return 2.*(fqqbar)/ x1;
}



double _Lum_qqb(double tau,void* p){
	Luminosity::par_struct* par = (Luminosity::par_struct*) p;
	par -> t = tau;
	double precision = 1e-4;
	double ris = 0.0, error = 0.0;
	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	Integrand.function = f_qqb;
	Integrand.params = par;
	gsl_integration_qags (&Integrand, tau, 1.0, 0, precision, 1000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	return ris;
}

double _Lqqb(double x, void *p){
	double z = std::exp(x);
	return (z*_Lum_qqb(z,p));
}

//qq channel (3)

double f_qq(double x1,void *p){
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf; 
	unsigned int Nf=p_struct->Nf;
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
	double fqq=0.;
	for (int i=0;i<Nf;i++){
	  fqq+= (PDF-> xfxQ((i+1),x1,Muf)/x1 * PDF -> xfxQ((i+1),x2,Muf)/x2)
	  +(PDF-> xfxQ(-(i+1),x1,Muf)/x1 * PDF -> xfxQ(-(i+1),x2,Muf)/x2);
	}
	return (fqq)/ x1;
}



double _Lum_qq(double tau,void* p){
	Luminosity::par_struct* par = (Luminosity::par_struct*) p;
	par -> t = tau;
	double precision = 1e-4;
	double ris = 0.0, error = 0.0;
	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	Integrand.function = f_qq;
	Integrand.params = par;
	gsl_integration_qags (&Integrand, tau, 1.0, 0, precision, 1000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	return ris;
}

double _Lqq(double x, void *p){
	double z = std::exp(x);
	return (z*_Lum_qq(z,p));
}

//qqbar prime channel (4)

double f_qQb(double x1,void *p){
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf; 
	unsigned int Nf=p_struct->Nf;
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
	double fqqbarprime=0.;
	for (int i=0;i<Nf;i++){
	  for (int j=0;j<Nf;j++){
	    if (i != j){
		fqqbarprime+= (PDF-> xfxQ(-(i+1),x1,Muf)/x1 * PDF -> xfxQ((j+1),x2,Muf)/x2)
		+(PDF-> xfxQ((i+1),x1,Muf)/x1 * PDF -> xfxQ(-(j+1),x2,Muf)/x2);
	    }
	  }
	}
	return (fqqbarprime)/ x1;
}



double _Lum_qQb(double tau,void* p){
	Luminosity::par_struct* par = (Luminosity::par_struct*) p;
	par -> t = tau;
	double precision = 1e-4;
	double ris = 0.0, error = 0.0;
	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	Integrand.function = f_qQb;
	Integrand.params = par;
	gsl_integration_qags (&Integrand, tau, 1.0, 0, precision, 1000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	return ris;
}

double _LqQb(double x, void *p){
	double z = std::exp(x);
	return (z*_Lum_qQb(z,p));
}

//qq prime channel (5)

double f_qQ(double x1,void *p){
	Luminosity::par_struct* p_struct = (Luminosity::par_struct *) p;
	double x2 = p_struct -> t /x1;
	double Muf = p_struct -> Muf; 
	unsigned int Nf=p_struct->Nf;
	Luminosity::PDF_ptr PDF = p_struct -> PDF;
	double fqqprime=0.;
	for (int i=0;i<Nf;i++){
	  for (int j=0;j<Nf;j++){
	    if (i != j){
		fqqprime+= (PDF-> xfxQ((i+1),x1,Muf)/x1 * PDF -> xfxQ((j+1),x2,Muf)/x2)
		+(PDF-> xfxQ((i+1),x1,Muf)/x1 * PDF -> xfxQ((j+1),x2,Muf)/x2);
	    }
	  }
	}
	return (fqqprime)/ x1;
}



double _Lum_qQ(double tau,void* p){
	Luminosity::par_struct* par = (Luminosity::par_struct*) p;
	par -> t = tau;
	double precision = 1e-4;
	double ris = 0.0, error = 0.0;
	gsl_function Integrand;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	Integrand.function = f_qQ;
	Integrand.params = par;
	gsl_integration_qags (&Integrand, tau, 1.0, 0, precision, 1000,	w, &ris, &error);
	gsl_integration_workspace_free (w);
	return ris;
}

double _LqQ(double x, void *p){
	double z = std::exp(x);
	return (z*_Lum_qQ(z,p));
}


//_____________________________________________________________________________________________________________________________________

//Function of Interpolation and final results


void Luminosity::Cheb_Lum(double muf){
	std::cout << "Computing parton luminosities at " << muf << " GeV..." << std::endl;
	_Muf = muf;
	
	double umin = -static_cast<double>(_n);

	for (auto& it : C_larges)	
		it = std::vector < double >( _n+1 , 0 );
	
	for (auto& it : Chebseries){
		gsl_cheb_free(it);
		it = gsl_cheb_alloc(_n);
	}
	gsl_function Large;
	par_struct par;
	par.PDF = _PDF;
	par.Muf = _Muf;
	par.Nf = _Nf;
	
	Large.params = &par;
	Large.function = &_Lgg;
	gsl_cheb_init(Chebseries[0],&Large,umin,0.0);
	Large.function = &_Lgq;
	gsl_cheb_init(Chebseries[1],&Large,umin,0.0);
	Large.function = &_Lqqb;
	gsl_cheb_init(Chebseries[2],&Large,umin,0.0);
	Large.function = &_Lqq;
	gsl_cheb_init(Chebseries[3],&Large,umin,0.0);
	Large.function = &_LqQb;
	gsl_cheb_init(Chebseries[4],&Large,umin,0.0);
	Large.function = &_LqQ;
	gsl_cheb_init(Chebseries[5],&Large,umin,0.0);
	
	c_large();
 }

void Luminosity::c_large(){

	double umin = -static_cast<double>(_n);
//	std::cout << " clarge n=" << _n << " umin =" << umin << std::endl;

	std::vector<double*> co_channels;
	for (auto it : Chebseries)
		co_channels.push_back(gsl_cheb_coeffs(it));
		
	std::size_t ser_length = Chebseries.size();
	
	for(unsigned k = 0; k < ser_length ; ++k){
		C_larges[k][0] = - co_channels[k][0]/2.;
		for(unsigned l = 0 ; l < (_n+1) ; l++ )
		 	for (unsigned i = l ; i < (_n+1) ; i++ )
				for(unsigned j = i ; j < (_n+1) ; j++ )
					C_larges[k][l] += pow(2./umin,l)*gsl_sf_fact(i)/gsl_sf_fact(i-l)*co_channels[k][j]*_T[j][i];
	}
}

Luminosity::dcomplex Luminosity::Lum_N(dcomplex N, unsigned channel) const {
	dcomplex ris(0.,0.);
	for (std::size_t i = 0 ; i < _n + 1 ; ++i)
		ris += C_larges[channel][i]/pow(N-1.,i+1.);
	return ris;
}
long double Luminosity::Lum_x(long double z, unsigned channel) const {
	long double ris=0.0;
	for (std::size_t i = 0 ; i < _n + 1 ; ++i)
		ris += C_larges[channel][i]/gsl_sf_fact(i)*std::pow(-1.,i)/z*std::pow(std::log(z),i);
	return ris;
}

//Test Function
void Luminosity::TestLum(std::string nameout){
  std::ofstream OUT;
  OUT.open(nameout);
  std::cout << "Test Luminosity written on file " << nameout << std::endl;
  par_struct pp;
  pp.PDF=_PDF;
  pp.Muf=_Muf;
  pp.Nf=_Nf;
  long double *x=NULL;
  int nt=100;
  long double xmin=1e-6; long double xmax=0.1;
  x=new long double[nt];
  x[0]=xmin;
  for (int i=1;i<nt;i++){
    x[i]=x[i-1]*1.14;
  }
  OUT << "x\t\trealPDF\t\tapproxPDF" << std::endl;
  OUT << "gg channel" << std::endl;
  for (int i=0;i < nt;i++){
    OUT << x[i] << "\t" << _Lum_gg(x[i],&pp)<< "\t" << Lum_gg_x(x[i]) << std::endl;
  }
  OUT << std::endl;
  OUT << "gq channel" << std::endl;
  for (int i=0;i < nt;i++){
    OUT << x[i] << "\t" << _Lum_gq(x[i],&pp)<< "\t" << Lum_gq_x(x[i]) << std::endl;
  }
  OUT << std::endl;
  OUT << "qqb channel" << std::endl;
  for (int i=0;i < nt;i++){
    OUT << x[i] << "\t" << _Lum_qqb(x[i],&pp)<< "\t" << Lum_qqb_x(x[i]) << std::endl;
  }
  OUT << std::endl;
  OUT << "qq channel" << std::endl;
  for (int i=0;i < nt;i++){
    OUT << x[i] << "\t" << _Lum_qq(x[i],&pp)<< "\t" << Lum_qq_x(x[i]) << std::endl;
  }
  OUT << std::endl;
   OUT << "qQb channel" << std::endl;
  for (int i=0;i < nt;i++){
    OUT << x[i] << "\t" << _Lum_qQb(x[i],&pp)<< "\t" << Lum_qQb_x(x[i]) << std::endl;
  }
  OUT << std::endl;
  OUT << "qQ channel" << std::endl;
  for (int i=0;i < nt;i++){
    pp.t=x[i];
    OUT << x[i] << "\t" << _Lum_qQ(x[i],&pp)<< "\t" << Lum_qQ_x(x[i]) << std::endl;
  }
  OUT << std::endl;
  OUT.close();
  return;
  
}



