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
  
  
long double inf_norm (const std::complex<long double> &z)
{
  return std::max (std::abs (std::real (z)),std::abs (std::imag (z)));
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

std::complex<long double> expm1 (const std::complex<long double> &z)
{
  const long double x = std::real (z),y = std::imag (z);

  if ((std::abs (x) >= 1.0) || (std::abs (y) >= 1.0)) return (std::exp (z) - 1.0);

  const long double expm1_x = std::expm1 (x),exp_x = 1.0 + expm1_x,sin_y_over_two = std::sin (0.5*y),sin_y = std::sin (y);

  return std::complex<long double> (expm1_x - 2.0*exp_x*sin_y_over_two*sin_y_over_two,exp_x*sin_y);
}

std::complex<long double> log1p (const std::complex<long double> &z)
{
  const long double x = std::real (z),y = std::imag (z); 

  const long double xp1 = 1.0 + x,abs_x = std::abs (x), abs_y = std::abs (y);

  if ((abs_x >= 1.0) || (abs_y >= 1.0)) return std::log (1.0 + z);
 
  const long double y_over_xp1 = y/xp1;

  return std::complex<long double> (std::log1p (x) + 0.5*std::log1p (y_over_xp1*y_over_xp1),std::atan2 (y,xp1));
}


std::complex<long double> LogGamma(std::complex<long double> z){
 int g=7;
 long double p[9];
 p[0]=0.99999999999980993;
 p[1]=676.5203681218851;
 p[2]=-1259.1392167224028;
 p[3]=771.32342877765313;
 p[4]=-176.61502916214059;
 p[5]=12.507343278686905;
 p[6]=-0.13857109526572012;
 p[7]=9.9843695780195716e-6;
 p[8]=1.5056327351493116e-7;
 std::complex<long double> sum,ris;
 z -=1;
 sum=p[0];
 for(int i=1;i<(g+2);i++)
   sum += p[i]/(z+i);
 ris=0.5*std::log(2*M_PIl)+(z+0.5)*std::log(z+g+0.5)-(z+g+0.5)+std::log(sum);
 return ris;
}

std::complex<long double> Gamma_inv (const std::complex<long double> &z)
{
  if ((!std::isfinite (std::real(z)))||(!std::isfinite(std::imag(z)))) std::cout<<"z is not finite in Gamma_inv."<<std::endl, abort ();

  const long double x = real (z);

  if (x >= 0.5)
  {
    const long double log_sqrt_2Pi = 0.91893853320467274177,g = 4.7421875;
    const std::complex<long double> z_m_0p5 = z - 0.5, z_pg_m0p5 = z_m_0p5 + g, zm1 = z - 1.0;
    const long double c[15] = {0.99999999999999709182,
			  57.156235665862923517,
			  -59.597960355475491248,
			  14.136097974741747174,
			  -0.49191381609762019978,
			  0.33994649984811888699E-4,
			  0.46523628927048575665E-4,
			  -0.98374475304879564677E-4,
			  0.15808870322491248884E-3,
			  -0.21026444172410488319E-3,
			  0.21743961811521264320E-3,
			 -0.16431810653676389022E-3,
			  0.84418223983852743293E-4,
			  -0.26190838401581408670E-4,
			  0.36899182659531622704E-5};
      
    std::complex<long double> sum = c[0];
    for (int i = 1 ; i < 15 ; i++) sum += c[i]/(zm1 + static_cast<long double>(i));

    const std::complex<long double> Gamma_inv_z = std::exp (z_pg_m0p5 - z_m_0p5*std::log (z_pg_m0p5) - log_sqrt_2Pi)/sum;

    return Gamma_inv_z;
  }
  else
  {
    const int n = static_cast<int> (std::rint (x));

    const std::complex<long double> eps = z - (static_cast<long double>(n));

    if (n%2 == 0)
      return (std::sin (M_PI*eps)*M_1_PI)/Gamma_inv (1.0 - z);
    else
      return (-std::sin (M_PI*eps)*M_1_PI)/Gamma_inv (1.0 - z);
  }
}

std::complex<long double> Gamma_ratio_diff_small_eps (const std::complex<long double> &z,const std::complex<long double> &eps)
{
  const std::complex<long double> zero(0.,0.);
  const long double g = 4.7421875;
  if (inf_norm (eps) > 0.1) std::cout<<"One must have |eps|oo < 0.1 in Gamma_ratio_diff_small_eps.", abort ();

  const std::complex<long double> eps_pz = z + eps,z_m_0p5 = z - 0.5,z_pg_m0p5 = z_m_0p5 + g,eps_pz_pg_m0p5 = z_pg_m0p5 + eps, zm1 = z - 1.0,zm1_p_eps = zm1 + eps;
  const long double x = std::real (z),eps_px = std::real (eps_pz);
  const int n = static_cast<int> (std::rint (x)),m = static_cast<int> (std::rint (eps_px));
  long double nl=static_cast<long double> (n);
  long double ml=static_cast<long double> (m);
  if ((std::real(z) == static_cast<long double>(n))&& (std::imag(z)<1e-20) && (n <= 0)) std::cout<<"z is negative integer in Gamma_ratio_diff_small_eps.", abort ();
  if ((std::real(eps_pz) == static_cast<long double> (m))&&(std::imag(z)<1e-20) && (m <= 0)) std::cout<<"z+eps is negative integer in Gamma_ratio_diff_small_eps.", abort ();

  const long double c[15] = {0.99999999999999709182,
			57.156235665862923517,
			-59.597960355475491248,
			14.136097974741747174,
			-0.49191381609762019978,
			0.33994649984811888699E-4,
			0.46523628927048575665E-4,
			-0.98374475304879564677E-4,
			0.15808870322491248884E-3,
			-0.21026444172410488319E-3,
			0.21743961811521264320E-3,
			-0.16431810653676389022E-3,
			0.84418223983852743293E-4,
			-0.26190838401581408670E-4,
			0.36899182659531622704E-5};

  if ((x >= 0.5) || (eps_px >= 0.5))
  {
    std::complex<long double> sum_num = 0.0, sum_den = c[0];

    for (int i = 1 ; i < 15 ; i++) 
    {
      long double il=static_cast<long double> (i);
      const std::complex<long double> ci_zm1_pi_inv = c[i]/(zm1 + il);
      sum_num += ci_zm1_pi_inv/(zm1_p_eps + il), sum_den += ci_zm1_pi_inv;
    }

    if (eps != zero)
      return expm1 (z_m_0p5*log1p (eps/z_pg_m0p5) + eps*std::log (eps_pz_pg_m0p5) - eps + log1p (-eps*sum_num/sum_den))/eps;
    else
      return (z_m_0p5/z_pg_m0p5 + std::log (eps_pz_pg_m0p5) - 1.0 - sum_num/sum_den);
  }
  else
  {
    if (eps != zero)
    {
      
      const std::complex<long double> Pi_eps = M_PIl*eps,term = std::sin (Pi_eps)/std::tan (M_PIl*(z-nl)),T1_eps_z = (cos (Pi_eps) + term)*Gamma_ratio_diff_small_eps (1.0 - z,-eps);
      const std::complex<long double> sin_Pi_2_eps = sin (M_PI_2*eps),T2_eps_z = (2.0*sin_Pi_2_eps*sin_Pi_2_eps - term)/eps, T_eps_z = T1_eps_z + T2_eps_z;
      
      return (T_eps_z/(1.0 - eps*T_eps_z));
    }
    else
      return (Gamma_ratio_diff_small_eps (1.0 - z,0.0) - M_PIl/tan (M_PIl*(z-nl)));
  }
}

std::complex<long double> Gamma_inv_diff_eps (const std::complex<long double> &z,const std::complex<long double> &eps)
{
 
  const std::complex<long double> eps_pz = z + eps;
  const long double x = std::real (z),eps_px = std::real (eps_pz);
  const int n = static_cast<int> (std::rint (x)), m = static_cast<int> (std::rint (eps_px));
  
  const bool is_z_negative_integer = ((std::real(z) == static_cast<long double>(n))&&(std::imag(z)<1e-20)) && (n <= 0);
  const bool is_eps_pz_negative_integer = ((std::real(eps_pz) == static_cast<long double> (m))&&(std::imag(eps_pz)<1e-20)) && (m <= 0);

  if (inf_norm (eps) > 0.1)
    return (Gamma_inv (z) - Gamma_inv (eps_pz))/eps;
  else if (eps_pz != z)
  { 
    if (is_z_negative_integer)
      return (-Gamma_inv (eps_pz)/eps);
    else if (is_eps_pz_negative_integer)
      return (Gamma_inv (z)/eps);
    else
    {
      const long double z_neg_int_distance = inf_norm (z + std::abs (x)),eps_pz_neg_int_distance = inf_norm (eps_pz + std::abs (eps_px));

      if (z_neg_int_distance < eps_pz_neg_int_distance)
	return Gamma_ratio_diff_small_eps (z,eps)*Gamma_inv (eps_pz);
      else
	return Gamma_ratio_diff_small_eps (eps_pz,-eps)*Gamma_inv (z);
    }
  }
  else if (is_z_negative_integer && is_eps_pz_negative_integer)
  {
    long double fact = -1.0;
    for (int k = -1 ; k >= n ; k--) fact *= static_cast<long double>(k);
    
    return fact;
  }
  else
    return Gamma_ratio_diff_small_eps (z,eps)*Gamma_inv (eps_pz);
}



std::complex<long double> LBesselJ(long double k, std::complex<long double> z){
  const long double zreal=static_cast<double>(std::real(z)),zimag=static_cast<double>(std::imag(z));
  double kk= static_cast<double> (k);
  std::complex<double> zz(zreal,zimag);
  std::complex<double> IId(0.,1.);
  std::complex<double> ris;
  if (std::abs(zz) < 25.){
  ris=sp_bessel::besselJ(kk,zz);
  }
  else {
    ris=1./std::sqrt(2.*M_PI*zz)*(std::exp(IId*(1./4.*(2.*kk+1.)*M_PI-zz))*(1.+3675./32768./std::pow(zz,4.)-3229./6144.*kk*kk/std::pow(zz,4.)
    +329./1024.*std::pow(kk/zz,4.)-7./128.*std::pow(kk,6.)/std::pow(zz,4.)+std::pow(kk,8.)/384./std::pow(zz,4.)
    -75.*IId/1024./std::pow(zz,3.)+259.*IId*kk*kk/768./std::pow(zz,3)-35.*IId*std::pow(kk,4.)/192./std::pow(zz,3.)+IId*std::pow(kk,6.)/48.*std::pow(zz,3.)
    -9./128./(zz*zz)+5./16.*kk*kk/zz/zz-std::pow(kk,4.)/(8.*zz*zz)+IId/8./zz-IId *kk*kk/2./zz)+std::exp(-IId*(1./4.*(2.*kk+1.)*M_PI-zz))*
    (1.+3675./32768./std::pow(zz,4.)-3229./6144.*kk*kk/std::pow(zz,4.)
    +329./1024.*std::pow(kk/zz,4.)-7./128.*std::pow(kk,6.)/std::pow(zz,4.)+std::pow(kk,8.)/384./std::pow(zz,4.)
    +75.*IId/1024./std::pow(zz,3.)-259.*IId*kk*kk/768./std::pow(zz,3)+35.*IId*std::pow(kk,4.)/192./std::pow(zz,3.)-IId*std::pow(kk,6.)/48.*std::pow(zz,3.)
    -9./128./(zz*zz)+5./16.*kk*kk/zz/zz-std::pow(kk,4.)/(8.*zz*zz)-IId/8./zz+IId *kk*kk/2./zz));
    
  }
  const long double risreal=static_cast<long double>(std::real(ris)),risimag=static_cast<long double>(std::imag(ris));
  std::complex<long double> res(risreal,risimag);
  return res;
}

std::complex<long double> LBesselK(long double k, std::complex<long double> z){
  const long double zreal=static_cast<double>(std::real(z)),zimag=static_cast<double>(std::imag(z));
  double kk= static_cast<double> (k);
  std::complex<double> zz(zreal,zimag);
  std::complex<double> IId(0.,1.);
  std::complex<double> ris;
  if (std::abs(zz) < 25.){
  ris=sp_bessel::besselK(kk,zz);
  }
  else {
    ris=std::sqrt(M_PI/(2.*zz))*std::exp(-zz)*(1.-59535./262144./std::pow(zz,5.)+352407.*kk*kk/327680./std::pow(zz,5.)
    -17281./24576.*std::pow(kk,4.)/std::pow(zz,5.)+1463.*std::pow(kk,6.)/10240./std::pow(zz,5.)
    -11./1024.*std::pow(kk,8.)/std::pow(zz,5.)+std::pow(kk,10.)/3840./std::pow(zz,5.)+3675./32768./std::pow(zz,4.)-3229./6144.
    *kk*kk/std::pow(zz,4.)+329./1024.*std::pow(kk/zz,4.)-7./128.*std::pow(kk,6.)/std::pow(zz,4.)+std::pow(kk,8.)/384./std::pow(zz,4.)
    -75./1024./std::pow(zz,3.)+259./768.*kk*kk/std::pow(zz,3.)-35./192.*std::pow(kk,4.)/std::pow(zz,3.)-35./192.*std::pow(kk,4.)/std::pow(zz,3.)
    +std::pow(kk,6.)/48./std::pow(zz,3.)+9./128./zz/zz-5./16.*kk*kk/zz/zz+std::pow(kk,4.)/(8.*zz*zz)-1./(8.*zz)+kk*kk/(2.*zz));
    
  }
  const long double risreal=static_cast<long double>(std::real(ris)),risimag=static_cast<long double>(std::imag(ris));
  std::complex<long double> res(risreal,risimag);
  return res;
}

//Useful Function For the evaluation of the Hypergeometric Function 2F1 //VEDERE AEAE per maggiori informazioni

std::complex<long double> A_sum_init (const int m,const std::complex<long double> &eps,const std::complex<long double> &Gamma_inv_one_meps)
{
  long double ml=static_cast<long double> (m);
  const std::complex<long double> one_meps = 1.0 - eps;
  
  if (one_meps-ml != 1.-ml)
  {
    std::complex<long double> Gamma_inv_one_meps_mm = Gamma_inv_one_meps;
    for (int i = 1 ; i <= ml ; i++) Gamma_inv_one_meps_mm *= one_meps - i;
    return Gamma_inv_one_meps_mm/eps;
  }
  else
  {
    long double fact = 1.0;
    for (int n = 2 ; n < m ; n++) fact *= n;

    return (m%2 == 0) ? (fact) : (-fact);
  }
}

std::complex<long double> log_A_sum_init (const int m,const std::complex<long double> &eps)
{
  long double ml=static_cast<long double> (m);
  const std::complex<long double> one_meps_mm = 1.0 - ml - eps;

  if (one_meps_mm != 1.-ml)
    return (-LogGamma(one_meps_mm) - std::log (eps));
  else
  {
    const std::complex<long double> i_Pi(0,M_PIl);

    long double log_fact = 0.0;
    for (int n = 2 ; n < m ; n++) log_fact += std::log (static_cast<long double>(n));

    return (m%2 == 0) ? (log_fact) : (log_fact + i_Pi);
  }
}

std::complex<long double> B_sum_init_PS_one (const std::complex<long double> &a,const std::complex<long double> &b,const std::complex<long double> &c,
				   const std::complex<long double> &Gamma_c,const std::complex<long double> &Gamma_inv_one_meps,
				   const std::complex<long double> &Gamma_inv_eps_pa_pm,const std::complex<long double> &Gamma_inv_eps_pb_pm,
				   const std::complex<long double> &one_minus_z,const int m,const std::complex<long double> &eps)
{
  long double ml=static_cast<long double> (m);
  const long double inf_norm_eps = inf_norm (eps),phase = (m%2 == 0) ? (1.) : (-1.);
  const std::complex<long double> a_pm = a + ml,b_pm = b + ml,one_meps = 1.0 - eps,Pi_eps = M_PIl*eps,Pi_eps_pm = M_PIl*(eps + ml);

  std::complex<long double> Gamma_inv_one_meps_mm = Gamma_inv_one_meps;
  for (int i = 1 ; i <= m ; i++) Gamma_inv_one_meps_mm *= one_meps - (static_cast<long double>(i));
  if (inf_norm_eps > 0.1)
  {
    const std::complex<long double> Gamma_inv_eps_pm_p1 = phase*std::sin(Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
    std::complex<long double> prod1 = Gamma_inv_one_meps*Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm;
    for (int n = 0 ; n < m ; n++) {
      long double nl=static_cast<long double>(n);
      prod1 *= (a + nl)*(b + nl)/(nl + 1.0);
    }
    const std::complex<long double> prod2 = Gamma_inv (a)*Gamma_inv (b)*Gamma_inv_eps_pm_p1*pow (one_minus_z,eps),res = Gamma_c*(prod1 - prod2)/eps;
    return res;
  }
  else
  {
    long double Gamma_inv_mp1 = 1.0;
    std::complex<long double> prod_ab = 1.0;
    for (int n = 0 ; n < m ; n++) {
      long double nl=static_cast<long double>(n);
      Gamma_inv_mp1 /= nl + 1.0; 
      prod_ab *= (a + nl)*(b + nl);
    }

    const bool is_eps_non_zero = (one_meps-ml != 1.-ml);
    const std::complex<long double> Gamma_inv_eps_pm_p1 = (is_eps_non_zero) ? (phase*std::sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm)) : (Gamma_inv_mp1);
    const std::complex<long double> Gamma_inv_a_pm = Gamma_inv (a_pm),Gamma_inv_b_pm = Gamma_inv (b_pm);
    const std::complex<long double> z_term = (is_eps_non_zero) ? (expm1 (eps*std::log (one_minus_z))/eps) : (std::log (one_minus_z));

    const std::complex<long double> prod1 = Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm*(Gamma_inv_mp1*Gamma_inv_diff_eps (1.0,-eps) + Gamma_inv_diff_eps (m+1,eps));
    const std::complex<long double> prod2 = Gamma_inv_eps_pm_p1*(Gamma_inv_eps_pb_pm*Gamma_inv_diff_eps (a_pm,eps) + Gamma_inv_a_pm*Gamma_inv_diff_eps (b_pm,eps));
    const std::complex<long double> prod3 = Gamma_inv_a_pm*Gamma_inv_b_pm*Gamma_inv_eps_pm_p1*z_term;

    const std::complex<long double> res = Gamma_c*prod_ab*(prod1 - prod2 - prod3); 
    if ((std::isfinite (std::real(res)))&&(std::isfinite(std::imag(res))))
      return res;
    else
    {
      const std::complex<long double> Gamma_inv_eps_pm_p1 = phase*std::sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
      std::complex<long double> prod1 = Gamma_inv_one_meps*Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm;
      for (int n = 0 ; n < m ; n++) {
	long double nl=static_cast<long double> (nl);
	prod1 *= (a + nl)*(b + nl)/(nl + 1.0);
      }
      const std::complex<long double> prod2 = Gamma_inv (a)*Gamma_inv (b)*Gamma_inv_eps_pm_p1*pow (one_minus_z,eps),res_default = Gamma_c*(prod1 - prod2)/eps;
      return res_default;
    }
  }
}

std::complex<long double> B_sum_init_PS_infinity (const std::complex<long double> &a,const std::complex<long double> &c,
					const std::complex<long double> &Gamma_c,const std::complex<long double> &Gamma_inv_cma,
					const std::complex<long double> &Gamma_inv_one_meps,const std::complex<long double> &Gamma_inv_eps_pa_pm,
					const std::complex<long double> &z,const int m,const std::complex<long double> &eps)
{
  long double ml=static_cast<long double> (m);
  const long double inf_norm_eps = inf_norm(eps),phase = (m%2 == 0) ? (1.) : (-1.);
  const std::complex<long double> cma = c - a,a_mc_p1 = 1.0 - c + a,a_mc_p1_pm = a_mc_p1 + ml,cma_meps = cma - eps,eps_pa_mc_p1 = eps + a_mc_p1,a_pm = a + m;
  const std::complex<long double> Gamma_inv_cma_meps = Gamma_inv (cma_meps),one_meps = 1.0 - eps,Pi_eps = M_PIl*eps,Pi_eps_pm = M_PIl*(eps + ml);

  std::complex<long double> Gamma_inv_one_meps_mm = Gamma_inv_one_meps;
  for (int i = 1 ; i <= m ; i++) {
    long double il=static_cast<long double> (i);
    Gamma_inv_one_meps_mm *= one_meps - il;
  } 
  if (inf_norm_eps > 0.1)
  { 
    const std::complex<long double> Gamma_inv_eps_pm_p1 = phase*std::sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
    std::complex<long double> prod1 = Gamma_inv_cma*Gamma_inv_eps_pa_pm*Gamma_inv_one_meps,prod2 = Gamma_inv (a)*Gamma_inv_cma_meps*Gamma_inv_eps_pm_p1*std::pow (-z,-eps);

    for (int n = 0 ; n < m ; n++)
    {
      long double nl=static_cast<long double> (n);
      const long double n_p1 = nl + 1.0;
      const std::complex<long double> a_pn = a + nl,a_mc_p1_pn = a_mc_p1 + nl,eps_pa_mc_p1_pn = eps_pa_mc_p1 + nl;
      prod1 *= a_pn*a_mc_p1_pn/n_p1, prod2 *= eps_pa_mc_p1_pn;
    }

    const std::complex<long double> res = Gamma_c*(prod1 - prod2)/eps;
    return res;
  }
  else
  {
    const int n0 = -static_cast<int> (std::rint (std::real (a_mc_p1)));
    const bool is_eps_non_zero = (one_meps-ml != 1.-ml),is_n0_here = (n0 >= 0) && (n0 < m);

    long double Gamma_inv_mp1 = 1.0;
    std::complex<long double> prod_a = 1.0,prod_a_mc_p1 = 1.0,prod_eps_pa_mc_p1_n0 = (is_n0_here) ? (1.0) : (0.0),prod_eps_pa_mc_p1 = 1.0,sum = 0.0;
    for (int n = 0 ; n < m ; n++)
    {
      long double nl=static_cast<long double> (n);
      const std::complex<long double> a_pn = a + nl,a_mc_p1_pn = a_mc_p1 + nl,eps_pa_mc_p1_pn = eps_pa_mc_p1 + nl;
      prod_a *= a_pn, prod_a_mc_p1 *= a_mc_p1_pn, prod_eps_pa_mc_p1 *= eps_pa_mc_p1_pn, Gamma_inv_mp1 /= nl+1.0;

      if (n != n0) 
      {
	if (is_n0_here) prod_eps_pa_mc_p1_n0 *= eps_pa_mc_p1_pn;
	sum += (is_eps_non_zero) ? (log1p (eps/a_mc_p1_pn)) : (1.0/a_mc_p1_pn);
      }
    }

    const std::complex<long double> Gamma_inv_eps_pm_p1 = (is_eps_non_zero) ? (phase*std::sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm)) : (Gamma_inv_mp1);
    const std::complex<long double> sum_term = (is_eps_non_zero) ? (expm1 (sum)/eps) : (sum),prod_diff_eps = prod_eps_pa_mc_p1_n0 + prod_a_mc_p1*sum_term;

    const std::complex<long double> z_term = (is_eps_non_zero) ? (expm1 (-eps*std::log (-z))/eps) : (-std::log (-z));
    const std::complex<long double> Gamma_inv_a_pm = Gamma_inv (a_pm),Gamma_prod1 = Gamma_inv_cma*Gamma_inv_eps_pa_pm;
    const std::complex<long double> prod1 = Gamma_prod1*Gamma_inv_mp1*(Gamma_inv_diff_eps (1.0,-eps)*prod_eps_pa_mc_p1 - Gamma_inv_one_meps*prod_diff_eps);
    const std::complex<long double> prod_2a = Gamma_prod1*Gamma_inv_diff_eps (m+1,eps),prod_2b = Gamma_inv_cma*Gamma_inv_eps_pm_p1*Gamma_inv_diff_eps (a_pm,eps);
    const std::complex<long double> prod_2c = Gamma_inv_eps_pm_p1*Gamma_inv_a_pm*(Gamma_inv_diff_eps (cma,-eps) + Gamma_inv_cma_meps*z_term);
    const std::complex<long double> prod2 = prod_eps_pa_mc_p1*(prod_2a - prod_2b - prod_2c),res = Gamma_c*prod_a*(prod1 + prod2);

    if ((std::isfinite (std::real(res)))&&(std::isfinite(std::imag(res))))
      return res;
    else
    { 
      const std::complex<long double> Gamma_inv_eps_pm_p1 = phase*std::sin (Pi_eps)/(Pi_eps_pm*Gamma_inv_one_meps_mm);
      std::complex<long double> prod1 = Gamma_inv_cma*Gamma_inv_eps_pa_pm*Gamma_inv_one_meps,prod2 = Gamma_inv (a)*Gamma_inv_cma_meps*Gamma_inv_eps_pm_p1*std::pow (-z,-eps);
      
      for (int n = 0 ; n < m ; n++)
      {
	long double nl=static_cast<long double> (n);
	const double n_p1 = nl + 1.0;
	const std::complex<long double> a_pn = a + nl,a_mc_p1_pn = a_mc_p1 + nl,eps_pa_mc_p1_pn = eps_pa_mc_p1 + nl;
	prod1 *= a_pn*a_mc_p1_pn/n_p1, prod2 *= eps_pa_mc_p1_pn;
      }

      const std::complex<long double> res_default = Gamma_c*(prod1 - prod2)/eps;
      return res_default;
    }
  }
}


void cv_poly_der_tab_calc (const std::complex<long double> &a,const std::complex<long double> &b,const std::complex<long double> &c,const std::complex<long double> &z,long double cv_poly_der_tab[])
{
  const long double mod_a2 = std::norm (a),mod_b2 = std::norm (b),mod_c2 = std::norm (c),mod_z2 = std::norm (z);
  const long double Re_a = std::real (a),Re_b = std::real (b),Re_c = std::real (c);

  cv_poly_der_tab[0] = 2.0*((Re_a*mod_b2 + Re_b*mod_a2)*mod_z2 - Re_c - mod_c2);
  cv_poly_der_tab[1] = 2.0*((mod_a2 + mod_b2 + 4.0*Re_a*Re_b)*mod_z2 - 1.0 - 4.0*Re_c - mod_c2);
  cv_poly_der_tab[2] = 6.0*((Re_a + Re_b)*mod_z2 - Re_c - 1.0);
  cv_poly_der_tab[3] = 4.0*(mod_z2 - 1.0);
}

long double cv_poly_der_calc (const long double cv_poly_der_tab[],const long double x)
{
  const long double Px = cv_poly_der_tab[0] + x*(cv_poly_der_tab[1] + x*(cv_poly_der_tab[2] + x*cv_poly_der_tab[3]));
  return Px;
}

int min_n_calc (const long double cv_poly_der_tab[])
{
  const long double C1 = cv_poly_der_tab[1],C2 = cv_poly_der_tab[2],three_C3 = 3.0*cv_poly_der_tab[3];
  const long double Delta = C2*C2 - three_C3*C1;

  if (Delta <= 0.0) 
    return 0;
  else
  {
    const double largest_root = -(C2 + std::sqrt (Delta))/three_C3;
    return std::max (static_cast<int> (std::ceil (largest_root)),0);
  }
}

// Calculation of the 2F1 power series converging for |z| < 1
std::complex<long double> hyp_PS_zero (const std::complex<long double> &a,const std::complex<long double> &b,const std::complex<long double> &c,const std::complex<long double> &z)
{  
  std::complex<long double> term = 1.0,sum = 1.0;
  
  const int na = std::abs (static_cast<int> (std::rint (std::real (a)))),nb = std::abs (static_cast<int> (std::rint (std::real (b))));

  if ((static_cast<int>(std::real(a)) == -na)&&(std::imag(a)<1e-20))
    for (int n = 0 ; n < na ; n++)
    {
      long double nl=static_cast<long double> (n);
      term *= z*(a+nl)*(b+nl)/((nl+1.0)*(c+nl));
      sum += term;
    }
  else if ((static_cast<int>(std::real(b)) == -nb)&&(std::imag(b)<1e-20))
    for (int n = 0 ; n < nb ; n++)
    {
      long double nl=static_cast<long double> (n);
      term *= z*(a+nl)*(b+nl)/((nl+1.0)*(c+nl));
      sum += term;
    }
  else
  {
    long double cv_poly_der_tab[4];
    cv_poly_der_tab_calc (a,b,c,z,cv_poly_der_tab);

    const int min_n = min_n_calc (cv_poly_der_tab);
    bool possible_false_cv = true;
    int n = 0;
    long double nl=static_cast<long double> (n);
    while (possible_false_cv || (inf_norm (term) > 1E-15))
    {
      term *= z*(a+nl)*(b+nl)/((nl+1.0)*(c+nl));
      sum += term;
      if (possible_false_cv && (n >= min_n)) possible_false_cv = (cv_poly_der_calc (cv_poly_der_tab,nl) > 0);
      n++;
      nl=static_cast<long double> (n);
    }
  }

  return sum;
}


std::complex<long double> hyp_PS_one(const std::complex<long double> &a,const std::complex<long double> &b, const std::complex<long double> &c,const std::complex<long double> &one_minus_z)
{
  const std::complex<long double> zero(0.,0.);
  const int m = static_cast<int> (std::rint (std::real (c - a - b))),phase = (m%2 == 0) ? (1) : (-1),m_m1 = m - 1, m_p1 = m + 1;
  const long double ml=(static_cast<long double>(m));
  const std::complex<long double> eps = c - a - b - ml,eps_pm = eps + ml,eps_pm_p1 = eps_pm + 1.0,a_pm = a + ml,b_pm = b + ml,one_meps = 1.0 - eps,one_meps_mm = one_meps - ml;
  const std::complex<long double> eps_pa = a + eps,eps_pb = b + eps,eps_pa_pm = eps_pa + ml,eps_pb_pm = eps_pb + ml,Pi_eps = M_PI*eps,Gamma_c = 1.0/Gamma_inv (c);
  const std::complex<long double> Gamma_inv_eps_pa_pm = Gamma_inv (eps_pa_pm),Gamma_inv_eps_pb_pm = Gamma_inv (eps_pb_pm);
  const std::complex<long double> Gamma_prod = Gamma_c*Gamma_inv_eps_pa_pm*Gamma_inv_eps_pb_pm,Gamma_inv_one_meps = Gamma_inv (one_meps);
  
  const std::complex<long double> A_first_term = (m > 0) ? (Gamma_prod*A_sum_init (m,eps,Gamma_inv_one_meps)) : (0.0);
  std::complex<long double> A_sum = A_first_term,A_term = A_first_term;

  if ((!std::isfinite (std::real(A_first_term)))||(!std::isfinite(std::imag(A_first_term)))) 
  {
    A_term = std::exp (LogGamma(c) - LogGamma(eps_pa_pm) - LogGamma(eps_pb_pm) + log_A_sum_init (m,eps));
    A_sum=A_term;
    if ((std::imag (a) == 0.0) && (std::imag (b) == 0.0) && (std::imag (c) == 0.0)) {
      A_term = std::real (A_term);
      A_sum=A_term;
    }
  }

  const std::complex<long double> pow_mzp1_m = std::pow (one_minus_z,m);
  const std::complex<long double> B_first_term = B_sum_init_PS_one (a,b,c,Gamma_c,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm,Gamma_inv_eps_pb_pm,one_minus_z,m,eps)*pow_mzp1_m;

  std::complex<long double> prod_B = pow_mzp1_m;
  for (int n = 0 ; n < m_m1 ; n++) 
  {
    long double nl=static_cast<long double> (nl);
    const std::complex<long double> ratio = (a + nl)*(b + nl)/(nl + 1.0);
    A_term *= one_minus_z*ratio/(nl + one_meps_mm), A_sum += A_term, prod_B *= ratio;
  }

  if (m > 0) prod_B *= (a + ml - 1.0)*(b + ml - 1.0)/ml;

  std::complex<long double> B_extra_term = prod_B*Gamma_prod*Gamma_inv_one_meps,B_term = B_first_term,B_sum = B_first_term;
  const long double B_prec = 1E-15*inf_norm (B_first_term);

  long double cv_poly1_der_tab[4],cv_poly2_der_tab[4];
  cv_poly_der_tab_calc (a,b,one_meps_mm,one_minus_z,cv_poly1_der_tab), cv_poly_der_tab_calc (eps_pb_pm,eps_pa_pm,eps_pm_p1,one_minus_z,cv_poly2_der_tab);

  const int min_n = std::max (min_n_calc (cv_poly1_der_tab),min_n_calc (cv_poly2_der_tab));
  bool possible_false_cv = true;

  int n = 0;
  while (possible_false_cv || (inf_norm (B_term) > B_prec))
  {
    long double nl=static_cast<long double> (n);
    const int n_pm_p1 = n + m_p1,n_p1 = n + 1;
    long double n_pm_p1l=static_cast<long double> (n_pm_p1);
    long double n_p1l=static_cast<long double> (n_p1);
    const std::complex<long double> a_pm_pn = a_pm + nl,b_pm_pn = b_pm + nl,eps_pm_p1_pn = eps_pm_p1 + nl,n_p1_meps = one_meps + nl;
    const std::complex<long double> eps_pa_pm_pn = eps_pa_pm + nl,eps_pb_pm_pn = eps_pb_pm + nl,eps_pm_pn = eps_pm + nl;
    const std::complex<long double> prod1 = eps_pa_pm_pn*eps_pb_pm_pn,prod2 = eps_pm_p1_pn*n_p1l,prod3 = a_pm_pn*b_pm_pn;

    B_term = one_minus_z*(B_term*prod1/prod2 + B_extra_term*(prod3/n_pm_p1l - a_pm_pn - b_pm_pn - eps + prod1/n_p1l)/(eps_pm_p1_pn*n_p1_meps));
    B_sum += B_term;
    B_extra_term *= one_minus_z*prod3/(n_pm_p1l*n_p1_meps);
    if (possible_false_cv && (n >= min_n)) possible_false_cv = (cv_poly_der_calc (cv_poly1_der_tab,nl) > 0) || (cv_poly_der_calc (cv_poly2_der_tab,nl) > 0);
    n++;
  }

  return (eps == zero) ? (phase*(A_sum + B_sum)) : ((A_sum + B_sum)*phase*Pi_eps/std::sin (Pi_eps));
} 

std::complex<long double> hyp_PS_infinity (const std::complex<long double> &a,const std::complex<long double> &b,const std::complex<long double> &c,const std::complex<long double> &z)
{
  const std::complex<long double> zero(0.,0.);
  const int m = static_cast<int> (std::rint (std::real (b - a))),phase = (m%2 == 0) ? (1) : (-1), m_m1 = m - 1,m_p1 = m + 1;
  long double ml=static_cast<long double> (m);
  long double m_m1l=static_cast<long double> (m_m1);
  long double m_p1l=static_cast<long double> (m_p1);
  const std::complex<long double> eps = b - a - ml,a_mc_p1 = 1.0 - c + a,one_meps = 1.0 - eps,one_meps_mm = one_meps - ml,a_pm = a + ml,a_mc_p1_pm = a_mc_p1 + ml,cma = c - a;
  const std::complex<long double> eps_pa = eps + a, eps_pm_p1 =  eps + ml + 1,eps_pa_mc_p1_pm = eps + a_mc_p1_pm,Pi_eps = M_PI*eps;
  const std::complex<long double> eps_pa_pm = eps_pa + ml,eps_pm = eps + ml,Gamma_c = 1.0/Gamma_inv (c),Gamma_inv_eps_pa_pm = Gamma_inv (eps_pa_pm);
  const std::complex<long double> Gamma_inv_cma = Gamma_inv (cma),z_inv = 1.0/z,pow_mz_ma = pow (-z,-a),Gamma_inv_one_meps = Gamma_inv (one_meps);
  const std::complex<long double> Gamma_prod = Gamma_c*Gamma_inv_cma*Gamma_inv_eps_pa_pm;

  const std::complex<long double> A_first_term = (m > 0) ? (Gamma_prod*A_sum_init (m,eps,Gamma_inv_one_meps)) : (0.0);
  std::complex<long double> A_sum = A_first_term,A_term = A_first_term;

  if ((!std::isfinite (std::real(A_first_term)))||(!std::isfinite(std::imag(A_first_term)))) 
  {
    A_term = std::exp (LogGamma (c) - LogGamma (cma) - LogGamma (b) + log_A_sum_init (m,eps));
    A_sum=A_term;
    if ((std::imag (a) == 0.0) && (std::imag (b) == 0.0) && (std::imag (c) == 0.0)) {
      A_term = std::real (A_term);
      A_sum=A_term;
    }
  }

  const std::complex<long double> pow_z_inv_m = std::pow (z_inv,m);
  const std::complex<long double> B_first_term = B_sum_init_PS_infinity (a,c,Gamma_c,Gamma_inv_cma,Gamma_inv_one_meps,Gamma_inv_eps_pa_pm,z,m,eps)*pow_z_inv_m;

  std::complex<long double> prod_B = pow_z_inv_m;
  for (int n = 0 ; n < m_m1 ; n++) 
  {
    long double nl=static_cast<long double> (n);
    const std::complex<long double> ratio = (a + nl)*(a_mc_p1 + nl)/(nl + 1.0);
    A_term *= z_inv*ratio/(nl + one_meps_mm), A_sum += A_term, prod_B *= ratio;
  }

  if (m > 0) prod_B *= (a + ml - 1.0)*(a_mc_p1 + ml - 1.0)/ml;

  std::complex<long double> B_extra_term = prod_B*Gamma_prod*Gamma_inv_one_meps,B_term = B_first_term,B_sum = B_first_term;
  const long double B_prec = 1E-15*inf_norm (B_first_term);
  
  long double cv_poly1_der_tab[4],cv_poly2_der_tab[4];
  cv_poly_der_tab_calc (a,a_mc_p1,one_meps_mm,z_inv,cv_poly1_der_tab), cv_poly_der_tab_calc (b,eps_pa_mc_p1_pm,eps_pm_p1,z_inv,cv_poly2_der_tab);

  const int min_n = std::max (min_n_calc (cv_poly1_der_tab),min_n_calc (cv_poly2_der_tab));
  bool possible_false_cv = true;

  int n = 0;
  while (possible_false_cv || (inf_norm (B_term) > B_prec))
  {
    long double nl=static_cast<long double> (n);
    const int n_pm_p1 = n + m_p1,n_p1 = n + 1;
    long double n_pm_p1l=static_cast<long double> (n_pm_p1);
    long double n_p1l=static_cast<long double> (n_p1);
    const std::complex<long double> a_pm_pn = a_pm + nl,a_mc_p1_pm_pn = a_mc_p1_pm + nl,eps_pm_p1_pn = eps_pm_p1 + nl,n_p1_meps = one_meps + nl;
    const std::complex<long double> eps_pa_pm_pn = eps_pa_pm + nl,eps_pa_mc_p1_pm_pn = eps_pa_mc_p1_pm + nl,eps_pm_pn = eps_pm + nl;
    const std::complex<long double> prod1 = eps_pa_pm_pn*eps_pa_mc_p1_pm_pn,prod2 = eps_pm_p1_pn*n_p1,prod3 = a_pm_pn*a_mc_p1_pm_pn;

    B_term = z_inv*(B_term*prod1/prod2 + B_extra_term*(prod3/n_pm_p1l - a_pm_pn - a_mc_p1_pm_pn - eps + prod1/n_p1l)/(eps_pm_p1_pn*n_p1_meps));
    B_sum += B_term;
    B_extra_term *= z_inv*prod3/(n_pm_p1l*n_p1_meps);
    if (possible_false_cv && (n >= min_n)) possible_false_cv = (cv_poly_der_calc (cv_poly1_der_tab,nl) > 0) || (cv_poly_der_calc (cv_poly2_der_tab,nl) > 0);
    n++;
  }

  return (eps == zero) ? (phase*pow_mz_ma*(A_sum + B_sum)) : ((A_sum + B_sum)*phase*pow_mz_ma*Pi_eps/std::sin (Pi_eps));
}

std::complex<long double> hyp_PS_complex_plane_rest (const std::complex<long double> &a,const std::complex<long double> &b,const std::complex<long double> &c,const std::complex<long double> &z)
{
  const long double abs_z = std::abs (z);
  const bool is_abs_z_small = abs_z < 1.0;
  
  const std::complex<long double> z0 = (is_abs_z_small) ? (0.9*z/abs_z) : (1.1*z/abs_z),zc = z - z0;
  const std::complex<long double> zc_z0_ratio = zc/(z0*(1.0 - z0)),z0_term1 = 2.0*z0 - 1.0,z0_term2 = c - (a + b + 1.0)*z0;
  
  const std::complex<long double> hyp_PS_z0 = (is_abs_z_small) ? (hyp_PS_zero (a,b,c,z0)) : (hyp_PS_infinity (a,b,c,z0));
  const std::complex<long double> dhyp_PS_z0 = (is_abs_z_small) ? (hyp_PS_zero (a + 1.0,b + 1.0,c + 1.0,z0)*a*b/c) : (hyp_PS_infinity (a + 1.0,b + 1.0,c + 1.0,z0)*a*b/c);

  int n = 0;
  std::complex<long double> an = hyp_PS_z0, anp1 = zc*dhyp_PS_z0, sum = an + anp1;

  const double prec = 1E-15*(inf_norm (an) + inf_norm (anp1));

  while (inf_norm (an) + inf_norm (anp1) > prec)
  {
    long double nl=static_cast<long double> (n);
    const std::complex<long double> anp2 = zc_z0_ratio*(anp1*(nl*z0_term1 - z0_term2) + an*zc*(a + nl)*(b + nl)/(nl + 1.))/(nl + 2.);
    sum += anp2;
    n++;
    an = anp1;
    anp1 = anp2;
  }

  return sum;
}




std::complex<long double> Hyp2F1(std::complex<long double> a,std::complex<long double> b, std::complex<long double> c, std::complex<long double> z)
{
  const double Re_a = std::real (a), Re_b = std::real (b), Re_c = std::real (c);

  const int na = static_cast<int> (std::rint (Re_a)),nb = static_cast<int> (std::rint (Re_b)),nc = static_cast<int> (std::rint (Re_c));
  const bool is_a_neg_int = ((std::real(a) == static_cast<long double> (na))&&(std::imag(a)<1e-20)) && (na <= 0),is_b_neg_int = ((std::real(b) == static_cast<long double> (nb))&&(std::imag(b)<1e-20)) && (nb <= 0),is_c_neg_int = ((std::real(c) == static_cast<long double> (nc))&&(std::imag(c)<1e-20)) && (nc <= 0);
  const std::complex<long double> zm1 = z-1.0;
  const std::complex<long double> uno(1.0,0.0);

  if (is_c_neg_int)
  {
    if (is_a_neg_int && (nc < na))
    {
      const std::complex<long double> z_over_zm1 = z/zm1;
      return ((z == uno) || (std::abs (z) < std::abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (std::pow (-zm1,-a)*hyp_PS_zero (a,c-b,c,z_over_zm1));
    }
    else if (is_b_neg_int && (nc < nb))
    {
      const std::complex<long double> z_over_zm1 = z/zm1;
      return ((z == uno) || (std::abs (z) < std::abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (std::pow (-zm1,-b)*hyp_PS_zero (b,c-a,c,z_over_zm1));
    }
    else 
      std::cout<<"2F1 undefined"<<std::endl, abort ();
  }

  if (is_a_neg_int)
  {
    const std::complex<long double> z_over_zm1 = z/zm1;
    return ((z == uno) || (std::abs (z) < std::abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (std::pow (-zm1,-a)*hyp_PS_zero (a,c-b,c,z_over_zm1));
  }
  else if (is_b_neg_int)
  {
    const std::complex<long double> z_over_zm1 = z/zm1;
    return ((z == uno) || (std::abs (z) < std::abs (z_over_zm1))) ? (hyp_PS_zero (a,b,c,z)) : (std::pow (-zm1,-b)*hyp_PS_zero (b,c-a,c,z_over_zm1));
  }

  const std::complex<long double> z_shift(std::real (z),-1E-307);
  if ((std::real (z) >= 1.0) && (std::imag (z) == 0.0)) return Hyp2F1 (a,b,c,z_shift);

  const bool ab_condition = (Re_b >= Re_a), cab_condition = (Re_c >= Re_a + Re_b);
  if (!ab_condition || !cab_condition)
  {
    if (!ab_condition && cab_condition) return Hyp2F1 (b,a,c,z);
    else if (!cab_condition && ab_condition) return std::pow (-zm1,c-a-b)*Hyp2F1 (c-b,c-a,c,z);
    else return std::pow (-zm1,c-a-b)*Hyp2F1 (c-a,c-b,c,z);
  }

  const long double abs_zm1 = std::abs (zm1);
  if (abs_zm1 < 1E-5) return hyp_PS_one (a,b,c,-zm1);

  const long double abs_z = std::abs (z),abs_z_inv = 1.0/abs_z,abs_z_over_zm1 = abs_z/abs_zm1,abs_zm1_inv = 1.0/abs_zm1,abs_zm1_over_z = 1.0/abs_z_over_zm1;
  const bool are_ac_small = (inf_norm (a) < 5.0) && (inf_norm (c) < 5.0),is_cmb_small = (inf_norm (c-b) < 5.0);
  const bool are_abc_small = are_ac_small && (inf_norm (b) < 5.0),are_a_cmb_c_small = are_ac_small && is_cmb_small;
  
  const double R_tab[5] = {0.5,0.6,0.7,0.8,0.9};
  
  for (unsigned int i = 0 ; i < 5 ; i++)
  {
    const double R = R_tab[i];

    if (abs_z <= R) return hyp_PS_zero (a,b,c,z);
    if (is_cmb_small && (abs_z_over_zm1 <= R)) return std::pow (-zm1,-a)*hyp_PS_zero (a,c-b,c,z/zm1);
  }
  
  for (unsigned int i = 0 ; i < 5 ; i++)
  {
    const double R = R_tab[i];

    if (abs_z_inv <= R) return hyp_PS_infinity (a,b,c,z);
    if (is_cmb_small && (abs_zm1_over_z <= R)) return std::pow (-zm1,-a)*hyp_PS_infinity (a,c-b,c,z/zm1);

    if (are_abc_small && (abs_zm1 <= R)) return hyp_PS_one (a,b,c,-zm1);
    if (are_a_cmb_c_small && (abs_zm1_inv <= R)) return std::pow (-zm1,-a)*hyp_PS_one (a,c-b,c,-1.0/zm1);
  }

  return hyp_PS_complex_plane_rest (a,b,c,z);  
}