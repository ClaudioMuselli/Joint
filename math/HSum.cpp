#include "HSum.h"



//Table of interpolated function




HSum::HSum(bool verbose,bool testinterfun, bool testharmsums):
_verbose(verbose),_testinterpolatedfunction(testinterfun),_testharmonicsums(testharmsums){
  
  (_verbose) ? std::cout<< "HSUM is a class to evaluate Harmonic Sums up to four order" << endl:std::cout<< "";
  (_verbose) ? std::cout<< "Inizializing coefficients" << std::endl:std::cout<< "";
  InizializeConst();
  (_verbose) ? std::cout << "DONE" << endl:std::cout<< "";
  
  
  
  if (_testinterpolatedfunction==true){
    std::cout << "Testing Semi-Analiting Mellin Trasformed function" << endl;
    TestInterpolatedFunction();
    std::cout << "DONE" << endl;
  }
  if (_testharmonicsums==true){
    std::cout << "Testing Harmonic Sums on integers value" << endl;
    std::cout << "Set the positive integer value you want to test:" << endl;
    std::cout << "Positive Integer Number: ";
    int n=0;
    std::cin >> n;
    if (n<=0) {
      std::cout << "ERROR: negative Test number - Test will be evaluated with n=7" << endl;
      n=7;
    }
    TestHarmonicSums(n);
    std::cout << "DONE " << endl;
  }
  (_verbose) ? std::cout << "Ready to use" << endl:std::cout<< "";
  
		    
}

std::complex<long double> HSum::HS(int i, std::complex<long double> N){
  switch (i){
    case(-4):
      return H_m4(N);
    case(-3):
      return H_m3(N);
    case(-2):
      return H_m2(N);
    case(-1):
      return H_m1(N);
    case(1):
      return H_1(N);
    case(2):
      return H_2(N);
    case(3):
      return H_3(N);
    case(4):
      return H_4(N);
  }
  cout << "ERROR: not valid Harmonic Sums index (Implemented Weight up to four yet)"<< endl;
  return (0.,0.);
      
}

std::complex<long double> HSum::HS(int i,int j, std::complex<long double> N){
  switch (i){
    case(-3):{
      switch (j){
	case(-1):
	  return H_m3_m1(N);
	case(1):
	  return H_m3_1(N);
      }
    }
    case(-2):{
      switch (j){
	case(-2):
	  return H_m2_m2(N);
	case(-1):
	  return H_m2_m1(N);
	case(1):
	  return H_m2_1(N);
	case(2):
	  return H_m2_2(N);
      }
    }
    case(-1):{
      switch (j){
	case(-3):
	  return H_m1_m3(N);
	case(-2):
	  return H_m1_m2(N);
	case(-1):
	  return H_m1_m1(N);
	case(1):
	  return H_m1_1(N);
	case(2):
	  return H_m1_2(N);
	case(3):
	  return H_m1_3(N);
      }
    }
    case(1):{
      switch (j){
	case(-3):
	  return H_1_m3(N);
	case(-2):
	  return H_1_m2(N);
	case(-1):
	  return H_1_m1(N);
	case(1):
	  return H_1_1(N);
	case(2):
	  return H_1_2(N);
	case(3):
	  return H_1_3(N);
      }
    }
    case(2):{
      switch (j){
	case(-2):
	  return H_2_m2(N);
	case(-1):
	  return H_2_m1(N);
	case(1):
	  return H_2_1(N);
	case(2):
	  return H_2_2(N);
      }
    }
    case(3):{
      switch (j){
	case(-1):
	  return H_3_m1(N);
	case(1):
	  return H_3_1(N);
      }
    }
  }
  cout << "ERROR: not valid Harmonic Sums index (Implemented Weight up to four yet)"<< endl;
  return (0.,0.);
      
}

std::complex<long double> HSum::HS(int i, int j, int k, std::complex<long double> N){
  switch (i){
    case(-2):{
      switch(j){
	case(-1):{
	  switch(k){
	    case(-1):
	      return H_m2_m1_m1(N);
	    case(1):
	      return H_m2_m1_1(N);
	  }
	}
	case(1):{
	  switch(k){
	    case(-1):
	      return H_m2_1_m1(N);
	    case(1):
	      return H_m2_1_1(N);
	  }
	}
      }
    }
    case(-1):{
      switch (j){
	case(-2):{
	  switch(k){
	    case(-1):
	      return H_m1_m2_m1(N);
	    case(1):
	      return H_m1_m2_1(N);
	  }
	}
	case(-1):{
	  switch(k){
	    case(-2):
	      return H_m1_m1_m2(N);
	    case(-1):
	      return H_m1_m1_m1(N);
	    case(1):
	      return H_m1_m1_1(N);
	    case(2):
	      return H_m1_m1_2(N);
	  }
	}
	case(1):{
	  switch(k){
	    case(-2):
	      return H_m1_1_m2(N);
	    case(-1):
	      return H_m1_1_m1(N);
	    case(1):
	      return H_m1_1_1(N);
	    case(2):
	      return H_m1_1_2(N);
	  }
	}
	case(2):{
	  switch(k){
	    case(-1):
	      return H_m1_2_m1(N);
	    case(1):
	      return H_m1_2_1(N);
	  }
	}
      }
    }
    case(1):{
      switch (j){
	case(-2):{
	  switch(k){
	    case(-1):
	      return H_1_m2_m1(N);
	    case(1):
	      return H_1_m2_1(N);
	  }
	}
	case(-1):{
	  switch(k){
	    case(-2):
	      return H_1_m1_m2(N);
	    case(-1):
	      return H_1_m1_m1(N);
	    case(1):
	      return H_1_m1_1(N);
	    case(2):
	      return H_1_m1_2(N);
	  }
	}
	case(1):{
	  switch(k){
	    case(-2):
	      return H_1_1_m2(N);
	    case(-1):
	      return H_1_1_m1(N);
	    case(1):
	      return H_1_1_1(N);
	    case(2):
	      return H_1_1_2(N);
	  }
	}
	case(2):{
	  switch(k){
	    case(-1):
	      return H_1_2_m1(N);
	    case(1):
	      return H_1_2_1(N);
	  }
	}
      }
    }
    case(2):{
      switch(j){
	case(-1):{
	  switch(k){
	    case(-1):
	      return H_2_m1_m1(N);
	    case(1):
	      return H_2_m1_1(N);
	  }
	}
	case(1):{
	  switch(k){
	    case(-1):
	      return H_2_1_m1(N);
	    case(1):
	      return H_2_1_1(N);
	  }
	}
      }
    }  
  }
  cout << "ERROR: not valid Harmonic Sums index (Implemented Weight up to four yet)"<< endl;
  return (0.,0.);
      
}

std::complex<long double> HSum::HS(int i,int j, int k, int m, std::complex<long double> N){
  switch (i){
    case(-1):{
      switch (j){
	case(-1):{
	  switch(k){
	    case(-1):{
	      switch(m){
		case(-1):
		  return H_m1_m1_m1_m1(N);
		case(1):
		  return H_m1_m1_m1_1(N);
	      }
	    }
	    case(1):{
	      switch(m){
		case(-1):
		  return H_m1_m1_1_m1(N);
		case(1):
		  return H_m1_m1_1_1(N);
	      }
	    }
	  }
	}
	case(1):{
	  switch(k){
	    case(-1):{
	      switch(m){
		case(-1):
		  return H_m1_1_m1_m1(N);
		case(1):
		  return H_m1_1_m1_1(N);
	      }
	    }
	    case(1):{
	      switch(m){
		case(-1):
		  return H_m1_1_1_m1(N);
		case(1):
		  return H_m1_1_1_1(N);
	      }
	    }
	  }
	}
      }
    }
    case(1):{
      switch (j){
	case(-1):{
	  switch(k){
	    case(-1):{
	      switch(m){
		case(-1):
		  return H_1_m1_m1_m1(N);
		case(1):
		  return H_1_m1_m1_1(N);
	      }
	    }
	    case(1):{
	      switch(m){
		case(-1):
		  return H_1_m1_1_m1(N);
		case(1):
		  return H_1_m1_1_1(N);
	      }
	    }
	  }
	}
	case(1):{
	  switch(k){
	    case(-1):{
	      switch(m){
		case(-1):
		  return H_1_1_m1_m1(N);
		case(1):
		  return H_1_1_m1_1(N);
	      }
	    }
	    case(1):{
	      switch(m){
		case(-1):
		  return H_1_1_1_m1(N);
		case(1):
		  return H_1_1_1_1(N);
	      }
	    }
	  }
	}
      }
    }
	      
  }
  cout << "ERROR: not valid Harmonic Sums index (Implemented Weight up to four yet)"<< endl;
  return (0.,0.);
      
}

HSum::~HSum(){
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//CROSSCHECK FUNCTIONS

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  






std::complex<long double> HSum::HS_int(int i,int n){
  long double ris=0.;
  long double neg=1.;
  if (i < 0){
    i=std::abs(i);
    neg=-1.;
  }
  long double I=((long double) i);
  for (int j=1;j<=n;j++){
    long double J=((long double) j);
    ris+=pow(neg,J)/(pow(J,I));
  }
  std::complex<long double> RIS(ris,0.);
  return RIS;
}

std::complex<long double> HSum::HS_int(int i,int j,int n){
  long double ris=0.,ris2=0.;
  long double neg1=1.,neg2=1.;
  if (i < 0){
    i=std::abs(i);
    neg1=-1.;
  }
  if (j < 0){
    j=std::abs(j);
    neg2=-1.;
  }
  long double I=((long double) i);
  long double J=((long double) j);
  for (int jj=1;jj<=n;jj++){
    ris2=0.;
    long double JJ=((long double) jj);
    for (int kk=1;kk<=jj;kk++){
      long double KK=((long double) kk);
      ris2+=pow(neg2,KK)/pow(KK,J);
    }
    ris+=pow(neg1,JJ)*ris2/(pow(JJ,I));
  }
  std::complex<long double> RIS(ris,0.);
  return RIS;
}

std::complex<long double> HSum::HS_int(int i,int j,int k,int n){
  long double ris=0.,ris2=0.,ris3=0.;
  long double neg1=1.,neg2=1.,neg3=1.;
  if (i < 0){
    i=std::abs(i);
    neg1=-1.;
  }
  if (j < 0){
    j=std::abs(j);
    neg2=-1.;
  }
  if (k < 0){
    k=std::abs(k);
    neg3=-1.;
  }
  long double I=((long double) i);
  long double J=((long double) j);
  long double K=((long double) k);
  for (int jj=1;jj<=n;jj++){
    ris2=0.;
    long double JJ=((long double) jj);
    for (int kk=1;kk<=jj;kk++){
      ris3=0.;
      long double KK=((long double) kk);
      for (int mm=1;mm<=kk;mm++){
	long double MM=((long double) mm);
	ris3+=pow(neg3,MM)/pow(MM,K);
      }
      ris2+=pow(neg2,KK)*ris3/pow(KK,J);
    }
    ris+=pow(neg1,JJ)*ris2/(pow(JJ,I));
  }
  std::complex<long double> RIS(ris,0.);
  return RIS;
}

std::complex<long double> HSum::HS_int(int i,int j,int k,int m,int n){
  long double ris=0.,ris2=0.,ris3=0., ris4=0.;
  long double neg1=1.,neg2=1.,neg3=1., neg4=1.;
  if (i < 0){
    i=std::abs(i);
    neg1=-1.;
  }
  if (j < 0){
    j=std::abs(j);
    neg2=-1.;
  }
  if (k < 0){
    k=std::abs(k);
    neg3=-1.;
  }
  if (m < 0){
    m=std::abs(m);
    neg4=-1.;
  }
  long double I=((long double) i);
  long double J=((long double) j);
  long double K=((long double) k);
  long double M=((long double) m);
  for (int jj=1;jj<=n;jj++){
    ris2=0.;
    long double JJ=((long double) jj);
    for (int kk=1;kk<=jj;kk++){
      ris3=0.;
      long double KK=((long double) kk);
      for (int mm=1;mm<=kk;mm++){
	ris4=0.;
	long double MM=((long double) mm);
	for (int ll=1;ll<=mm;ll++){
	  long double LL=((long double) ll);
	  ris4+=pow(neg4,LL)/pow(LL,M);
	}
	ris3+=pow(neg3,MM)*ris4/pow(MM,K);
      }
      ris2+=pow(neg2,KK)*ris3/pow(KK,J);
    }
    ris+=pow(neg1,JJ)*ris2/(pow(JJ,I));
  }
  std::complex<long double> RIS(ris,0.);
  return RIS;
}


void HSum::TestInterpolatedFunction(){
  //Values of CROSSCHECK
  std::vector<std::complex<long double>> g_v;
  std::complex<long double> NTest=2.+II;
  std::complex<long double> g_v_data[]={{0.0959473930360097,-0.028037163902461058},
  {0.055864908089439794,-0.013368158396651264},{0.17989000024656945,-0.04506116321742105},
  {-0.11083533556687188,+0.03156645059393223},{-0.039484446737595116,+0.021639724835958503},
  {0.14916093025637875,-0.03976343536944421},{-0.11969459327246698,+0.03361342326614195},
  {0.07221345274187332,-0.011765100701642627},{0.01695153326475704,-0.0039424452598304585},
  {-0.06575363790727289,+0.01310217472374898},{-0.3928282582240154,+0.05522661213910878},
  {0.22776185633616544,-0.03627358096229505}, {-0.06475111922019319,+0.015157058455765633},
  {0.20111362850218645,-0.06414491621631654}, {0.17127361976121863,-0.04077048445788456},
  {-0.1049810481804515,+0.02835961498705548},{0.1125205320368638,-0.027364003573713103},
  {0.9588133851510667,-0.24106914649348782}, {-0.21475568434004566,0.0738638239530835},
  {0.44028201016202234,-0.13732029537036497}, {1.7352370584367085,-0.30488004079023684},
  {0.3624728911057466,-0.09263683350215834}, {0-0.2511865365413727,+0.08517250266884471},
  {-0.6279919034848717,+0.14707979613696398}, {0.06804238802027272,-0.021226000425085434},
  {0.812027569319294,-0.08362066947521188}, {-0.31799014934197334,+0.062087083834581655},
  {0.15916325750396956,-0.056264163173798044},{-2.7155627365652166,+0.1469351851157078},
  {0.033305643746930014,-0.00676025817399897},{0.19263263117248935,-0.05662868509088942},
  {-0.004572972623993649,+0.0031181197404396054},{0.02006930899203301,-0.014493838276617852},
  {-0.1521324495518403,+0.051397028812878255},{0.001001015438086814,-0.0011574868971451393},
  {0.01016382069860128,-0.0056565347910271895},{-0.12177019428170494,+0.01702805482268213},
  {0.29997819614117216,-0.054266573935924864},{-1.3056857877637884,+0.17204278167513112}};
  g_v.insert(g_v.end(),&g_v_data[0],&g_v_data[sizeof(g_v_data)/(sizeof(g_v_data[0]))]);
  std::vector<std::complex<long double> >g;
  
  g.push_back(g1(NTest));g.push_back(g2(NTest));g.push_back(g3(NTest));g.push_back(g4(NTest));g.push_back(g5(NTest));
  g.push_back(g6(NTest));g.push_back(g7(NTest));g.push_back(g8(NTest));g.push_back(g9(NTest));g.push_back(g10(NTest));
  g.push_back(g11(NTest));g.push_back(g12(NTest));g.push_back(g13(NTest));g.push_back(g14(NTest));g.push_back(g15(NTest));
  g.push_back(g16(NTest));g.push_back(g17(NTest));g.push_back(g18(NTest));g.push_back(g19(NTest));g.push_back(g20(NTest));
  g.push_back(g21(NTest));g.push_back(g22(NTest));g.push_back(g23(NTest));g.push_back(g24(NTest));g.push_back(g25(NTest));
  g.push_back(g26(NTest));g.push_back(g27(NTest));g.push_back(g28(NTest));g.push_back(g29(NTest));g.push_back(g30(NTest));
  g.push_back(g31(NTest));g.push_back(g32(NTest));g.push_back(g33(NTest));g.push_back(g34(NTest));g.push_back(g35(NTest));
  g.push_back(g36(NTest));g.push_back(g37(NTest));g.push_back(g38(NTest));g.push_back(g39(NTest));
  bool test=true;
  std::vector<int >error;
  for (int i=0;i<g.size();i++){
    std::complex<long double> DIFF=g[i]-g_v[i];
    if ((std::abs(std::real(DIFF))<2e-8)&&(std::abs(std::imag(DIFF))<2e-8)){
    }
    else{
      test=false;
      error.push_back(i+1.);
    }
  }
  if (test==false){
    std::cout << "ERROR in functions:" << endl;
    for (int i=0;i<error.size();i++){
      std::cout << "g" << error[i] << " ";
    }
  }
  else
    std::cout << "All Interpolated Function are OK" << endl;
  return;
}


void HSum::TestHarmonicSums(int n){
  bool test1=true,test2=true,test3=true,test4=true;
  std::vector<int> error;
  std::complex<long double> N((long double) n,0.);
  for (int i=-4;i<=4;i++){
    if ((i!=0)&&(std::abs(i)<5)){
      std::complex<long double> DIFF=HS(i,N)-HS_int(i,n);
      if ((std::abs(std::real(DIFF))<2e-8)&&(std::abs(std::imag(DIFF))<2e-8)){
      }
      else{
	test1=false;
	error.push_back(i);
	error.push_back(0);
      }
    }
  }
  for (int i=-3;i<=3;i++){
    for (int j=-3;j<=3;j++){
      if ((i!=0)&&(j!=0)&&(std::abs(j)+std::abs(i)<5)){
	std::complex<long double> DIFF=HS(i,j,N)-HS_int(i,j,n);
	if ((std::abs(std::real(DIFF))<2e-8)&&(std::abs(std::imag(DIFF))<2e-8)){
	}
	else{
	  test2=false;
	  error.push_back(i);
	  error.push_back(j);
	  error.push_back(0);
	}
      }
    }
  }
  for (int i=-2;i<=2;i++){
    for (int j=-2;j<=2;j++){
      for(int k=-2;k<=2;k++){
	if ((i!=0)&&(j!=0)&&(k!=0)&&(std::abs(k)+std::abs(j)+std::abs(i)<5)){
	  std::complex<long double> DIFF=HS(i,j,k,N)-HS_int(i,j,k,n);
	  if ((std::abs(std::real(DIFF))<2e-8)&&(std::abs(std::imag(DIFF))<2e-8)){
	  }
	  else{
	    test3=false;
	    error.push_back(i);
	    error.push_back(j);
	    error.push_back(k);
	    error.push_back(0);
	  }
	}
      }
    }
  }
  for (int i=-1;i<=1;i++){
    for( int j=-1;j<=1;j++){
      for(int k=-1;k<=1;k++){
	for(int m=-1;m<=1;m++){
	  if ((i!=0)&&(j!=0)&&(k!=0)&&(m!=0)&&(std::abs(i)+std::abs(j)+std::abs(k)+std::abs(m)<5)){
	    std::complex<long double> DIFF=HS(i,j,k,m,N)-HS_int(i,j,k,m,n);
	    if ((std::abs(std::real(DIFF))<2e-8)&&(std::abs(std::imag(DIFF))<2e-8)){
	      }
	    else{
	      test4=false;
	      error.push_back(i);
	      error.push_back(j);
	      error.push_back(k);
	      error.push_back(m);
	      error.push_back(0);
	    }
	  }
	}
      }
    }
  }
  if ((test1==false)||(test2==false)||(test3==false)||(test4==false)){
    cout << "ERROR in Harmonic Sums with indices" << endl;
    for( int i=0;i<error.size();i++){
      if (error[i] !=0){
	cout << error[i];
      }
      else 
	cout << endl;
    }
  }
  else
    cout << "All Harmonic Sums OK" << endl;
    
}
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//CORE FUNCTIONS

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Inizialize Constants
void HSum::InizializeConst(){
  //Inizialize constant
  zeta2=gsl_sf_zeta(2.);
  zeta3=gsl_sf_zeta(3.);
  Li4=0.5174790616738993863;
  log2=std::log(2);
  log2q=log2*log2;
  log2c=log2*log2*log2;
  zeta2q=zeta2*zeta2;
  EulerGamma=0.577215664901532860606512090082;  
  
  
  //Definition of the expansion coefficients
  long double a1_data[]={0.999999974532238,-0.499995525889840,0.333203435557262,
		    -0.248529457782640,0.191451164719161,-0.137466222728331,
		    0.0792107412244877,-0.0301109656912626,0.00538406208663153};
  a1.insert(a1.end(),&a1_data[0],&a1_data[sizeof(a1_data)/(sizeof(a1_data[0]))]);


  long double a2_data[]={0.,0.999999980543793,-0.999995797779624,0.916516447393493,
		    -0.831229921350708,0.745873737923571,-0.634523908078600,
		    0.467104011423750,-0.261348046799178,0.0936814286867420,
		    -0.0156249375012462};
  a2.insert(a2.end(),&a2_data[0],&a2_data[sizeof(a2_data)/(sizeof(a2_data[0]))]);
		    
  long double a3_data[]={0.,0.,0.999999989322696,-1.49999722020708,1.74988008499745,
		    -1.87296689068405,1.91539974617231,-1.85963744001295,
		    1.62987195424434,-1.17982353224299,0.628710122994999,
		    -0.211307487211713,0.0328953352932140};
  a3.insert(a3.end(),&a3_data[0],&a3_data[sizeof(a3_data)/(sizeof(a3_data[0]))]);
		    
  long double b1_data[]={0.693147166991375,-0.306850436868254,0.193078041088284,
		    -0.139403892894644,0.105269615988049,-0.0746801353858524,
		    0.0427339135378207,-0.0161809049989783,0.00288664611077007};
  b1.insert(b1.end(),&b1_data[0],&b1_data[sizeof(b1_data)/(sizeof(b1_data[0]))]);		    
  
  long double b2_data[]={0.480453024731510,0.480450679641120,-0.519463586324817,
		    0.479285947990175,-0.427765744446172,0.360855321373065,
		    -0.263827078164263,0.146927719341510,-0.0525105367350968,
		    0.00874144396622167};
  b2.insert(b2.end(),&b2_data[0],&b2_data[sizeof(b2_data)/(sizeof(b2_data[0]))]);		    
  
  long double b3_data[]={0.33302465198526926,0.33302465294458333,0.33302458698245874,
		    -0.6669733788751905,0.8329914797601391,-0.9166248356766676,
		    0.9555225313767898,-0.9626278376733954,0.9315917979228248,
		    -0.8411435378489305,0.6749616810627938,-0.4530769450120006,
		    0.23754927911496704,-0.08942075739064419,0.021221995497457305,
		    -0.00236584329575064};
  b3.insert(b3.end(),&b3_data[0],&b3_data[sizeof(b3_data)/(sizeof(b3_data[0]))]);
		    
  long double c1_data[]={2.2012182965269744e-8,2.833327652357064,-1.8330909624101532,
		    0.7181879191200942,-0.0280403220046588,-0.181869786537805,
		    0.532318519269331,-1.07281686995035,1.38194913357518,
		    -1.11100841298484,0.506649587198046,-0.100672390783659};
  c1.insert(c1.end(),&c1_data[0],&c1_data[sizeof(c1_data)/(sizeof(c1_data[0]))]);
		    
  long double P21_data[]={11./6.,-3.,3./2.,-1./3.};
   P21.insert(P21.end(),&P21_data[0],&P21_data[sizeof(P21_data)/(sizeof(P21_data[0]))]);
  
  long double c2_data[]={1.989197140339627e-8,-6.050453690953361e-6,2.12530461606213,
		    -1.0523034829446278,0.160180000661971,-0.351982379713689,
		    1.41033369447519,-3.53344124579927,5.93934899678262,
		    -6.60019998525006,4.66330491799074,-1.89825521858848,
		    0.339773000152805};
  c2.insert(c2.end(),&c2_data[0],&c2_data[sizeof(c2_data)/(sizeof(c2_data[0]))]);
		    
  long double P22_data[]={-1.,5./2.,-2.,1./2.};
  P22.insert(P22.end(),&P22_data[0],&P22_data[sizeof(P22_data)/(sizeof(P22_data[0]))]);  

  long double c3_data[]={-2.878388705163104e-8,1.423616247405256,-0.08001203559240111,
		    -0.39875367195395994,0.339241791547134,-0.0522116678353452,
		    -0.0648354706049337,0.0644165053822532,-0.0394927322542075,
		    0.0100879370657869};
  c3.insert(c3.end(),&c3_data[0],&c3_data[sizeof(c3_data)/(sizeof(c3_data[0]))]);
		    
  long double P23_data[]={205./144.,-25./12.,23./24.,-13./36.,1./16.};
  P23.insert(P23.end(),&P23_data[0],&P23_data[sizeof(P23_data)/(sizeof(P23_data[0]))]);  
  long double P33_data[]={-25./24.,2.,-3./2.,2./3.,-1./8.};
  P33.insert(P33.end(),&P33_data[0],&P33_data[sizeof(P33_data)/(sizeof(P33_data[0]))]);  
  
  long double c4_data[]={-1.84461401708802e-8,2.2150086978693073,-0.9133677154535804,
		    3.4783104357500143,-2.823955592989266,0.992890266001707,
		    -1.30026190226546,3.41870577921103,-5.76763902370864,
		    6.45554138192407,-4.59405622046138,1.88510809558304,
		    -0.340476080290674};
  c4.insert(c4.end(),&c4_data[0],&c4_data[sizeof(c4_data)/(sizeof(c4_data[0]))]);
		    
  long double P24_data[]={-167./36.+25./6.*zeta2,235./18.-8.*zeta2,-40./3.+6.*zeta2,109./18.-8./3.*zeta2,-41./36.+1./2.*zeta2};
  P24.insert(P24.end(),&P24_data[0],&P24_data[sizeof(P24_data)/(sizeof(P24_data[0]))]);
  long double P34_data[]={35./12.,-26./3.,19./2.,-14./3.,11./12.};
  P34.insert(P34.end(),&P34_data[0],&P34_data[sizeof(P34_data)/(sizeof(P34_data[0]))]);
  
  long double c5_data[]={-0.822467033400776,0.0887664705657325,-0.0241549406045162,
		    0.00965074750946139,-0.00470587487919749,0.00246014308378549,
		    -0.00116431121874067,0.000395705193848026,-0.0000664699010014505};
  c5.insert(c5.end(),&c5_data[0],&c5_data[sizeof(c5_data)/(sizeof(c5_data[0]))]);
		    
  long double d5_data[]={-0.822467033400776,0.999999974532241,-0.249997762945014,
		    0.111067811851394,-0.0621323644338330,0.0382902328987004,
		    -0.0229110370338977,0.0113158200819689,-0.00376387065979726,
		    0.000598229109013054};
  d5.insert(d5.end(),&d5_data[0],&d5_data[sizeof(d5_data)/(sizeof(d5_data[0]))]);
		    
  long double q1_data[]={-0.09475300423010691,0.1454735027234842,-0.08098683369213346,
		    0.05023079316164129,-0.03411803198282981,0.02472443734108192,
		    -0.01877010647266778,0.014708747453969606,-0.011671217826519803,
		    0.009069802023937482,-0.006547951759267011,0.0041015946864745925,
		    -0.0020613218541571006,0.00075748485509974,-0.00017760906428930553,
		    0.000019714636373319215};
  q1.insert(q1.end(),&q1_data[0],&q1_data[sizeof(q1_data)/(sizeof(q1_data[0]))]);
		    
  long double q2_data[]={0.537213193604428,-1.119453718988143,1.174906969846667,
		    -1.1890779043586297,1.1945050676798088,-1.1967467424151426,
		    1.1955129514326055,-1.1835779969286169,1.1403589239212737,
		    -1.0314757738318858,0.8310004277895442,-0.5601432110623664,
		    0.29472391510409235,-0.11124910312554286,0.026456840745986296,
		    -0.002953839416623553};
  q2.insert(q2.end(),&q2_data[0],&q2_data[sizeof(q2_data)/(sizeof(q2_data[0]))]);
		    
  long double q3_data[]={-0.5372131936080373,0.045027332856608177,-0.010425994173158348,
		    0.0037472572760287814,-0.0017183961106026168,0.0009183269608058933,
		    -0.0005444461658565019,0.0003469566797451847,-0.00023098812103251377,
		    0.00015458246246608968,-0.00009826699232598133,0.00005519314304989385,
		    -0.000025223902811693763,8.518715418948123e-6,-1.8505304740610155e-6,
		    1.9151017720947663e-7};
  q3.insert(q3.end(),&q3_data[0],&q3_data[sizeof(q3_data)/(sizeof(q3_data[0]))]);
		    
  long double q4_data[]={0.0947530042271157,-0.3349795102848428,0.5614397835582168,
		    -0.6926554895103129,0.7769724000297276,-0.8354872769601905,
		    0.8767455156397378,-0.8995864697499111,0.8896397543874993,
		    -0.8197443755832042,0.6688529284198806,-0.45468132147833146,
		    0.24058512756020864,-0.09115874520857946,0.021735979874238424,
		    -0.0024313049233671805};
  q4.insert(q4.end(),&q4_data[0],&q4_data[sizeof(q4_data)/(sizeof(q4_data[0]))]);
		    
  long double q5_data[]={1.4092730215810628e-11,0.5822405222318227,-1.5665076685110302,
		    2.260578785582297,-2.7781329026065738,3.184698704229937,
		    -3.509076133374793,3.7417122180120175,-3.8132268023357767,
		    3.5937383966892074,-2.97986443052705,2.0483750064405015,
		    -1.0921389915653676,0.41599509143779606,-0.09955775776575748,
		    0.011165962058518938};
  q5.insert(q5.end(),&q5_data[0],&q5_data[sizeof(q5_data)/(sizeof(q5_data[0]))]);
		    
  long double q6_data[]={-2.3184816518117378e-11,6.968818347250131e-9,0.9999995028321177,
		    -1.999984562032115,2.916402389377294,-3.7471734064835496,
		    4.490720489784361,-5.106704825230779,5.465791341869803,
		    -5.340171712347202,4.541827792247,-3.1768555249894264,
		    1.7139895199002004,-0.6582025319003859,0.15843088157139373,
		    -0.01784285460130432};
  q6.insert(q6.end(),&q6_data[0],&q6_data[sizeof(q6_data)/(sizeof(q6_data[0]))]);
		    
  long double q7_data[]={-0.5822405264650113,0.11090665409384426,-0.042519755533205965,
		    0.02186263484057763,-0.013183846224852545,0.00877828110191841,
		    -0.006247889444605592,0.00465184457099542,-0.003542248910792736,
		    0.002663393028672393,-0.0018730231186815261,0.0011486239470476912,
		    -0.0005670390849428907,0.00020510083224129793,-0.000047391008407603993,
		    5.187375208790952e-6};
  q7.insert(q7.end(),&q7_data[0],&q7_data[sizeof(q7_data)/(sizeof(q7_data[0]))]);
			    

}



//Table of Interpolated functions

std::complex <long double> HSum::g1(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a2.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a2[i]/(N+k+1.);
}
return (0.5*(log2q-N*sum));
}
std::complex <long double> HSum::g2(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a3.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a3[i]/(N+k+1.);
}
return (1./3.*(log2c-N*sum));
  
}
std::complex <long double> HSum::g3(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a1.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a1[i]*(N/(N+k+1.)*zeta2+(k+1.)/(N+k+1.)/(N+k+1.)*H_1(N+k+1.));
}
return (zeta2*log2-sum);
}

std::complex <long double> HSum::g4(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a1.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a1[i]*(N/(N+k+1.)*zeta2/2.+(k+1.)/(N+k+1.)/(N+k+1.)*(log2-B(0,N+k+2.)));
}
return (-0.5*zeta2*log2+sum);
}

std::complex <long double> HSum::g5(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a1.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a1[i]*((k+1)/(N+k+1.)/(N+k+1.)*(zeta2+PolyGamma(1,N+k+2.)-2.*H_1(N+k+1.)/(N+k+1.)));
}
return (-sum);
}

std::complex <long double> HSum::g6(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a1.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a1[i]*(N/(N+k+1.)*zeta3+(k+1.)/(N+k+1.)/(N+k+1.)*(zeta2-H_1(N+k+1.)/(N+k+1.)));
}
return (zeta3*log2-sum);
}

std::complex <long double> HSum::g7(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a1.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a1[i]*(N/(N+k+1.)*zeta3*3./4.+(k+1.)/(N+k+1.)/(N+k+1.)*0.5*zeta2-(k+1.)/pow(N+k+1.,3.)*(log2-B(0,N+k+2.)));
}
return (-3./4.*zeta3*log2+sum);
}

std::complex <long double> HSum::g8(std::complex<long double > N){
std::complex<long double> sum;
for (int i=0;i<a1.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum+=a1[i]*(N/(N+k+1.)*zeta3+(k+1.)/(N+k+1.)/(N+k+1.)*0.5*(H_1(N+k+1.)*H_1(N+k+1.)+H_2(N+k+1.)));
}
return (zeta3*log2-sum);
}

std::complex <long double> HSum::g9(std::complex<long double > N){
std::complex<long double> sum,sum3;
for (int i=0;i<a1.size();i++){
  std::complex<long double> sum2;
  std::complex<long double> k((long double) i,0.);
  for (int j=0;j<a2.size();j++){
    std::complex<long double> l((long double) j,0.);
    sum2+=a2[j]/(N+k+l+2.);
  }
  sum+=a1[i]*N/(N+k+1.)*(zeta3/8.-0.5*sum2);
  
}
for (int i=0;i<a3.size();i++){
  std::complex<long double> k((long double) i,0.);
  sum3+=0.5*a3[i]/(N+k+1.);
}
return (1/8.*zeta3*log2-sum-sum3);
}

std::complex <long double> HSum::g10(std::complex<long double > N){
std::complex<long double> sum,sum3;
for (int i=0;i<a1.size();i++){
  std::complex<long double> sum2;
  std::complex<long double> k((long double) i,0.);
  for (int j=0;j<a1.size();j++){
    std::complex<long double> l((long double) j,0.);
    sum2+=a1[j]/(N+k+l+2.)*H_1(N+k+l+2.);
  }
  sum+=a1[i]*N/(N+k+1.)*(5.*zeta3/8.-sum2);
}
for (int h=0;h<a2.size();h++){
  std::complex<long double> m((long double) h,0.);
  sum3+=a2[h]/(N+m+1.)*H_1(N+m+1.);
}
return (-5./8.*zeta3*log2+sum+sum3);
}

std::complex <long double> HSum::g11(std::complex<long double > N){
  std::complex<long double> sum1;
  for (int i=0;i<c1.size();i++){
    std::complex<long double> k((long double) i,0.);
    std::complex<long double> sum2;
    for (int j=0;j<a1.size();j++){
      std::complex<long double> l((long double) j,0.);
      sum2+=a1[j]*((l+1.)/(N+k+l+1.)*H_1(N+k+l+1.)-H_1(l+1.));
    }
    sum1+=c1[i]*(0.5*(log2q-zeta2)-sum2);
  }
  std::complex<long double> sum3;
  for(int i=0;i<P21.size();i++){
    std::complex<long double> k((long double) i,0.);
    std::complex<long double> sum4;
    for (int j=0;j<a1.size();j++){
      std::complex<long double> l((long double) j,0.);
      sum4+=a1[j]*((l+1.)/(N+l+k+1.)*(H_1(N+k+l+1.)*H_1(N+k+l+1.)+H_2(N+l+k+1.))-H_1(l+1.)*H_1(l+1.)-H_2(l+1.));
    }
    sum3+=P21[i]*(7./4.*zeta3-zeta2*log2+1./3.*log2c+sum4);
  }
  return (sum1+sum3);
}

std::complex <long double> HSum::g12(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double)i,0.);
    std::complex<long double> sum2;
    for(int j=0;j<a2.size();j++){
      std::complex<long double> l((long double)j,0.);
      sum2+=a2[j]*(N+k+1.)/(N+k+l+2.);
    }
    sum+=a1[i]/(k+1.)*(0.5*(log2q-sum2)-B(1,N+k+2.)+B(0,N+k+2.)*(H_1(N+k+1.)-log2));
  }
  return sum;
}

std::complex <long double> HSum::g13(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a3.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=a3[i]/(N+k+1.);
  }
  for(int i=0;i<a2.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=a2[i]*N/(N+k+1.)*(0.5*zeta2-(log2-B(0,N+k+2.))/(N+k+1.));
  }
  return(-1./4.*zeta2*log2q+0.5*sum);
}

std::complex <long double> HSum::g14(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<b2.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=b2[i]/(N+k+1.);
  }
  return sum;
}

std::complex <long double> HSum::g15(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<b1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=b1[i]/(N+k+1.)*(zeta2-H_1(N+k+1.)/(N+k+1.));
  }
  return sum;
}

std::complex <long double> HSum::g16(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<b1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=b1[i]/(N+k+1.)*(-0.5*zeta2+(log2-B(0,N+k+2.))/(N+k+1.));
  }
  return sum;
}

std::complex <long double> HSum::g17(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<b2.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=b2[i]/(N+k+1.)/(N+k+1.);
  }
  return (-sum+log2q*PolyGamma(1,N+1.));
}

std::complex <long double> HSum::g18(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<c1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=c1[i]*N/(N+k)*H_1(N+k);
  }
  for (int i=0;i<P21.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=-P21[i]*N/(N+k)*(H_1(N+k)*H_1(N+k)+H_2(N+k));
  }
  return (1./N*(H_1(N)*H_1(N)+H_2(N))-zeta2*H_1(N)+sum);
}

std::complex <long double> HSum::g19(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=a1[i]/(k+1.)*H_1(N+k+1.);
  }
  return (0.5*zeta2*H_1(N)-sum);
}
std::complex <long double> HSum::g20(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<c2.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=c2[i]*N/(N+k)*H_1(N+k);
  }
  for (int i=0;i<P22.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=-P22[i]*N/(N+k)*(H_1(N+k)*H_1(N+k)+H_2(N+k));
  }
  for (int i=0;i<c4.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=-0.5*c4[i]*N/(N+k);
  }
  for (int i=0;i<P24.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=N*0.5*P24[i]*H_1(N+k)/(N+k);
  }
  for (int i=0;i<P34.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=-N*0.5*P34[i]*(H_1(N+k)*H_1(N+k)+H_2(N+k))/(N+k);
  }
  return (0.5*zeta2q-zeta3*H_1(N)+sum);
}

std::complex <long double> HSum::g21(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<c3.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=c3[i]*N/(N+k)*H_1(N+k);
  }
  for (int i=0;i<P33.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=N/(N+k)*(P33[i]*(H_1(N+k)*H_1(N+k)*H_1(N+k)+3.*H_1(N+k)*H_2(N+k)+2.*H_3(N+k))-P23[i]*(H_1(N+k)*H_1(N+k)+H_2(N+k)));
  }
  return (-zeta3*H_1(N)+1./(2.*N)*(H_1(N)*H_1(N)*H_1(N)+3.*H_1(N)*H_2(N)+2.*H_3(N))+sum);
}

std::complex <long double> HSum::g22(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<c1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=c1[i]*PolyGamma(1,N+k+1.);
  }
  for (int i=0;i<P21.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=-P21[i]*(H_1(N+k)*PolyGamma(1,N+k+1.)-0.5*PolyGamma(2,N+k+1.));
  }
  return (sum);
}

std::complex <long double> HSum::g23(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=a1[i]*N/(N+k+1.)*(zeta3-zeta2/(N+k+1.)-zeta2/(k+1.)+H_1(N+k+1.)*(1./(N+k+1.)/(N+k+1.)+1./(k+1.)/(k+1.)+1./((k+1.)*(N+k+1.))));
  }
  return (-0.5*zeta2*zeta2+zeta3*log2+3./4.*zeta3*H_1(N)-g6(N)-sum);
}

std::complex <long double> HSum::g24(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=a1[i]*N/(N+k+1.)*(2.*zeta3-(H_1(N+k+1.)*H_1(N+k+1.)+H_2(N+k+1.))/(N+k+1.));
  }
  for (int i=0;i<c5.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=c5[i]*N/(N+k+1.)*H_1(N+k+1.);
  }
  for (int i=0;i<d5.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=-d5[i]*N/(N+k)*(H_1(N+k)*H_1(N+k)+H_2(N+k));
  }
  return (-2.*zeta3*log2+2.*g8(N)+5./8.*zeta3*H_1(N)+sum);
}

std::complex <long double> HSum::g25(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a2.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=0.5*a2[i]/(k+1.)*N/(N+k+1.)*H_1(N+k+1.);
  }
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    std::complex<long double> sum2;
    for (int j=0;j<a1.size();j++){
      std::complex<long double> l((long double) j,0.);
      sum2+=a1[j]*H_1(N+k+l+2.)/(N+k+l+2.);
    }
    sum+=0.5*a1[i]*N/(N+k+1.)*(-5./8.*zeta3+sum2);
  }
  return (5./16.*zeta3*log2+0.5*g10(N)-1./8.*zeta3*H_1(N)+sum);
}

std::complex <long double> HSum::g26(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=a1[i]*((k+1.)/(N+k+1.)*(H_1(N+k+1.)*H_1(N+k+1.)+H_2(N+k+1.))-(H_1(k+1.)*H_1(k+1.)+H_2(k+1.)));
  }
  return( sum+7./4.*zeta3-zeta2*log2+1./3.*log2c);
}

std::complex <long double> HSum::g27(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum-=a1[i]*((k+1.)/(N+k+1.)*H_1(N+k+1.)-H_1(k+1.));
  }
  return( 0.5*(log2q-zeta2)+sum);
}

std::complex <long double> HSum::g28(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<b1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=b1[i]/(N+k+1.);
  }
  return sum;
}

std::complex <long double> HSum::g29(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum-=a1[i]*((k+1.)/(N+k+1.)*(H_1(N+k+1.)*H_1(N+k+1.)*H_1(N+k+1.)+3.*H_1(N+k+1.)*H_2(N+k+1.)+2.*H_3(N+k+1.))
    -(H_1(k+1.)*H_1(k+1.)*H_1(k+1.)+3.*H_1(k+1.)*H_2(k+1.)+2.*H_3(k+1.)));
  }
  return (sum-6.*Li4);
}

std::complex <long double> HSum::g30(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<a1.size();i++){
    std::complex<long double> k((long double) i,0.);
    std::complex<long double> sum2;
    for (int j=0;j<a3.size();j++){
      std::complex<long double> l((long double) j,0.);
      sum2+=a3[j]/(N+l+k+2.);
    }
    sum-=a1[i]*(N*sum2+3.*g2(N+k+1.));
  }
  return (log2q*log2q+sum);
}

std::complex <long double> HSum::g31(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<b3.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=b3[i]/(N+k+1.);
  }
  return sum;
}
    
std::complex <long double> HSum::g32(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<q1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=q1[i]/(N+k+1.);
  }
  return sum;
}

std::complex <long double> HSum::g33(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<q2.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=q2[i]/(N+k+1.);
  }
  return sum;
}
    
std::complex <long double> HSum::g34(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<q3.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=q3[i]/(N+k+1.);
  }
  return sum;
}    
    
std::complex <long double> HSum::g35(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<q4.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=q4[i]/(N+k+1.);
  }
  return sum;
}    
    
std::complex <long double> HSum::g36(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<q5.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum+=q5[i]/(N+k+1.);
  }
  return sum;
}

std::complex <long double> HSum::g37(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<q6.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum-=q6[i]/(N+k+1.)*H_1(N+k+1.);
  }
  return sum;
}

std::complex <long double> HSum::g38(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<q7.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum-=q7[i]/(N+k+1.)*H_1(N+k+1.);
  }
  return sum;
} 

std::complex <long double> HSum::g39(std::complex<long double > N){
  std::complex<long double> sum;
  for (int i=0;i<b1.size();i++){
    std::complex<long double> k((long double) i,0.);
    sum-=b1[i]/(N+k+1.)*(H_1(N+k+1.)*H_1(N+k+1.)+H_2(N+k+1.)+H_1(N+k+1.)*log2);
  }
  return (-g38(N)+sum);
}






//Table of Harmonic Sums
//Weight 1

std::complex<long double> HSum::H_1(std::complex<long double> N){
  return(PolyGamma(0,N+1.)+EulerGamma);
}
std::complex<long double> HSum::H_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)*B(0,N+1.)-log2);
}

//Weight 2

std::complex<long double> HSum::H_2(std::complex<long double> N){
  return(-PolyGamma(1,N+1.)+zeta2);
}
std::complex<long double> HSum::H_m2(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)*B(1,N+1.)-0.5*zeta2);
}

std::complex<long double> HSum::H_m1_m1(std::complex<long double> N){
  return(0.5*(H_m1(N)*H_m1(N)+H_2(N)));
}

std::complex<long double> HSum::H_m1_1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)*g1(N)+H_1(N)*H_m1(N)+H_m2(N)+(H_1(N)-H_m1(N))*log2-0.5*log2q);
}

std::complex<long double> HSum::H_1_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)*g1(N)-(H_1(N)-H_m1(N))*log2+0.5*log2q);
}

std::complex<long double> HSum::H_1_1(std::complex<long double> N){
  return(0.5*(H_1(N)*H_1(N)+H_2(N)));
}

//Weight 3

std::complex<long double> HSum::H_3(std::complex<long double> N){
  return(1./2.*PolyGamma(2,N+1.)+zeta3);
}

std::complex<long double> HSum::H_m3(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)/2.*B(2,N+1.)-3./4.*zeta3);
}

std::complex<long double> HSum::H_m2_m1(std::complex<long double> N){
  return(-g19(N)+log2*(H_2(N)-H_m2(N))-5./8.*zeta3);
}

std::complex<long double> HSum::H_m2_1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)*g3(N)+zeta2*H_m1(N)-5./8.*zeta3+zeta2*log2);
}

std::complex<long double> HSum::H_2_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)*g4(N)-log2*(H_2(N)-H_m2(N))-0.5*zeta2*H_m1(N)+1./4.*zeta3-0.5*zeta2*log2);
}

std::complex<long double> HSum::H_2_1(std::complex<long double> N){
  return(-g18(N)+2.*zeta3);
}

std::complex<long double> HSum::H_m1_m2(std::complex<long double> N){
  return(-H_m2_m1(N)+H_m2(N)*H_m1(N)+H_3(N));
}

std::complex<long double> HSum::H_m1_2(std::complex<long double> N){
  return(-H_2_m1(N)+H_2(N)*H_m1(N)+H_m3(N));
}

std::complex<long double> HSum::H_1_2(std::complex<long double> N){
  return(-H_2_1(N)+H_2(N)*H_1(N)+H_3(N));
}

std::complex<long double> HSum::H_1_m2(std::complex<long double> N){
  return(-H_m2_1(N)+H_m2(N)*H_1(N)+H_m3(N));
}

std::complex<long double> HSum::H_m1_m1_m1(std::complex<long double> N){
  return(1./6.*(2.*H_m3(N)+std::pow(H_m1(N),3.)+3.*H_m1(N)*H_2(N)));
}

std::complex<long double> HSum::H_m1_1_m1(std::complex<long double> N){
  return(0.5*g14(N)+log2*(H_m1_m1(N)-H_m1_1(N))+0.5*log2q*H_m1(N)+1./8.*zeta3-0.5*zeta2*log2+1./3.*log2c);
}

std::complex<long double> HSum::H_1_1_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(0.5*pow(-1.,nn)*g2(N)+log2*(H_1_m1(N)-H_1_1(N))+0.5*log2q*(H_1(N)-H_m1(N))-1./6.*log2c);
}

std::complex<long double> HSum::H_1_m1_m1(std::complex<long double> N){
  return(0.5*(-H_m1_1_m1(N)+H_m1(N)*H_1_m1(N)-H_2_1(N)-H_m1_m2(N)+H_1(N)*H_2(N)+H_m1(N)*H_m2(N)+2.*H_3(N)));
}

std::complex<long double> HSum::H_m1_m1_1(std::complex<long double> N){
  return(0.5*(-H_m1_1_m1(N)+H_m1(N)*H_m1_1(N)+H_m1_m2(N)+H_2_1(N)));
}

std::complex<long double> HSum::H_m1_1_1(std::complex<long double> N){
  return(H_1_1_m1(N)-H_1(N)*H_1_m1(N)+H_m2_1(N)+H_m1_2(N)+0.5*(H_1(N)*H_1(N)*H_m1(N)-H_m1(N)*H_2(N))-H_m3(N));
}

std::complex<long double> HSum::H_1_m1_1(std::complex<long double> N){
  return(-2.*H_1_1_m1(N)+H_1(N)*H_1_m1(N)-H_m2_1(N)-H_m1_2(N)+H_1(N)*H_m2(N)+H_m1(N)*H_2(N)+2.*H_m3(N));
}

std::complex<long double> HSum::H_1_1_1(std::complex<long double> N){
  return(1./6.*(H_1(N)*H_1(N)*H_1(N)+3.*H_1(N)*H_2(N)+2.*H_3(N)));
}

//Weight 4

std::complex<long double> HSum::H_4(std::complex<long double> N){
  return(-1./6.*PolyGamma(3,N+1.)+2./5.*zeta2q);
}

std::complex<long double> HSum::H_m4(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)/6.*B(3,N+1.)-7./20.*zeta2q);
}

std::complex<long double> HSum::H_m2_m2(std::complex<long double> N){
  return(0.5*(H_m2(N)*H_m2(N)+H_4(N)));
}

std::complex<long double> HSum::H_m3_1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)*g6(N)+zeta2*H_m2(N)-zeta3*H_m1(N)-3./5.*zeta2q+2.*Li4+3./4.*zeta3*log2-0.5*zeta2*log2q+1./12.*log2c*log2);
}

std::complex<long double> HSum::H_m2_2(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)*g5(N)-2.*H_m3_1(N)+2.*zeta2*H_m2(N)+3./40.*zeta2q);
}

std::complex<long double> HSum::H_2_m2(std::complex<long double> N){
  return(-H_m2_2(N)+H_m4(N)+H_m2(N)*H_2(N));
}

std::complex<long double> HSum::H_2_2(std::complex<long double> N){
  return(0.5*(H_2(N)*H_2(N)+H_4(N)));
}

std::complex<long double> HSum::H_3_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)*g7(N)-log2*(H_3(N)-H_m3(N))-0.5*zeta2*H_m2(N)+3./4.*zeta3*H_m1(N)-1./8.*zeta2q+3./4.*zeta3*log2);
}

std::complex<long double> HSum::H_m3_m1(std::complex<long double> N){
  return(g23(N)+log2*(H_3(N)-H_m3(N))-0.5*zeta2*H_2(N)-2.*Li4+11./10.*zeta2q-7./4.*zeta3*log2+0.5*zeta2*log2q-1./12.*log2c*log2);
}

std::complex<long double> HSum::H_3_1(std::complex<long double> N){
  return(0.5*g22(N)-1./4.*H_4(N)-1./4.*H_2(N)*H_2(N)+zeta2*H_2(N)-3./20.*zeta2q);
}

std::complex<long double> HSum::H_m1_m3(std::complex<long double> N){
  return(-H_m3_m1(N)+H_m3(N)*H_m1(N)+H_4(N));
}

std::complex<long double> HSum::H_m1_3(std::complex<long double> N){
  return(-H_3_m1(N)+H_3(N)*H_m1(N)+H_m4(N));
}

std::complex<long double> HSum::H_1_m3(std::complex<long double> N){
  return(-H_m3_1(N)+H_m3(N)*H_1(N)+H_m4(N));
}

std::complex<long double> HSum::H_1_3(std::complex<long double> N){
  return(-H_3_1(N)+H_3(N)*H_1(N)+H_4(N));
}

std::complex<long double> HSum::H_m2_1_1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)*g8(N)+zeta3*H_m1(N)-Li4+1./8.*zeta2q+1./8.*zeta3*log2+1./4.*zeta2*log2q-1./24.*log2c*log2);
}

std::complex<long double> HSum::H_m2_1_m1(std::complex<long double> N){
  return(-g25(N)-log2*(H_m2_1(N)-H_m2_m1(N))+0.5*log2q*(H_m2(N)-H_2(N))+3./40.*zeta2q);
}

std::complex<long double> HSum::H_m2_m1_m1(std::complex<long double> N){
  return(0.5*(-H_m1_m2_m1(N)+H_m2_2(N)+H_3_m1(N)+H_m1(N)*H_m2_m1(N)));
}


std::complex<long double> HSum::H_m2_m1_1(std::complex<long double> N){
  return(zeta2q/4.-g24(N)-H_2_m1_m1(N)+(-H_2_m1(N)+H_2_1(N))*log2-2.*Li4
	 +1./12.*(-H_m2(N)*(6.*zeta2-6.*log2q)+H_2(N)*(6*zeta2-6.*log2q)
	 -log2*(-6.*zeta2*log2+log2c+21.*zeta3)));
}

std::complex<long double> HSum::H_2_1_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)*g9(N)-log2*(H_2_1(N)-H_2_m1(N))+0.5*log2q*(H_2(N)-H_m2(N))
	+1./8.*zeta3*H_m1(N)+3.*Li4-6./5.*zeta2q+11./4.*zeta3*log2
	-3./4.*zeta2*log2q+1./8.*log2q*log2q);
}

std::complex<long double> HSum::H_1_2_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)*g13(N)-1./40.*(48.*zeta2q+20.*zeta2*H_1_m1(N)+80.*H_2_1_m1(N)
  -40.*H_1_m2(N)*log2+40.*H_1_2(N)*log2-80.*H_2_m1(N)*log2
  +80.*H_2_1(N)*log2-5.*log2q*log2q-10.*zeta2*H_m1(N)*2.*log2+10.*zeta2*H_1(N)*2.*log2+20.*zeta2*log2q
  +20.*H_m2(N)*2.*log2q-20.*H_2(N)*2.*log2q-120.*Li4-10.*H_1(N)*zeta3-105.*log2*zeta3));
}

std::complex<long double> HSum::H_m1_2_m1(std::complex<long double> N){
  return(g16(N)-2.*H_m2_1_m1(N)-log2*(2.*(H_m2_1(N)-H_m2_m1(N))+H_m1_2(N)-H_m1_m2(N))-0.5*zeta2*H_m1_m1(N)
  +log2q*(H_m2(N)-H_2(N))+1./4.*zeta3*H_1(N)+(0.25*zeta3-0.5*zeta2*log2)*(H_m1(N)-H_1(N))
  -log2*(H_m2_m1(N)-log2*(H_2(N)-H_m2(N))+0.5*zeta2*H_1(N))+33./20.*zeta2q-4.*Li4-13./4.*zeta3*log2+3./4.*zeta2*log2q-1./6.*log2q*log2q);
}

std::complex< long double > HSum::H_m1_2_1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(g11(N)*pow(-1.,nn)+1./1440.*(1044.*zeta2q+1440.*zeta2*H_m1_1(N)-2880.*H_m2_1_1(N)+360.*zeta2*log2q-180.*(log2q*log2q+24.*Li4)+2880.*H_m1(N)*zeta3));
}

std::complex<long double> HSum::H_m1_1_m2(std::complex<long double> N){
  return(-g17(N)/2.+(7./4.*zeta2q-zeta2*H_m1_1(N)-2.*H_m2_1_m1(N)-2.*H_m1_2_m1(N)
  +5./2.*zeta2*log2q-H_2(N)*log2q-log2q*log2q/6.+(H_3(N)-H_m2_1(N)-H_m1_2(N))*2.*log2+H_m2(N)*(log2q+H_m1(N)*2.*log2)-4.*Li4+1./4.*(H_m1(N)-21.*log2)*zeta3)/2.);
}

std::complex<long double> HSum::H_2_m1_m1(std::complex<long double> N){
  return(0.5*(-H_m1_2_m1(N)+H_m3_m1(N)+H_m1(N)*H_2_m1(N)+H_2_2(N)));
}

std::complex<long double> HSum::H_m1_m1_2(std::complex<long double> N){
  return(0.5*(-H_m1_2_m1(N)+H_m1_m3(N)+H_m1(N)*H_m1_2(N)+H_2_2(N)));
}

std::complex<long double> HSum::H_2_1_1(std::complex<long double> N){
  return(-g21(N)+6./5.*zeta2q);
}

std::complex<long double> HSum::H_1_2_1(std::complex<long double> N){
  return(-2.*H_2_1_1(N)+H_3_1(N)+H_1(N)*H_2_1(N)+H_2_2(N));
}

std::complex<long double> HSum::H_1_1_2(std::complex<long double> N){
  return(H_2_1_1(N)+0.5*(H_1(N)*(H_1_2(N)-H_2_1(N))+H_1_3(N)-H_3_1(N)));
}

std::complex<long double> HSum::H_1_m2_1(std::complex<long double> N){
  return(-2.*H_m2_1_1(N)+H_m3_1(N)+H_1(N)*H_m2_1(N)+H_m2_2(N));
}

std::complex<long double> HSum::H_1_1_m2(std::complex<long double> N){
  return(H_m2_1_1(N)+H_m2(N)*H_2(N)-H_m2_2(N)-H_m2(N)*H_1_1(N)+H_1(N)*H_1_m2(N)+H_1_m3(N)-H_1(N)*H_m3(N));
}

std::complex<long double> HSum::H_m1_m2_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(2.*pow(-1.,nn)*g12(N)-1./12.*(12.*zeta2*H_m1_1(N)+18.*zeta2*log2q
  -4.*log2q*log2q+12.*H_2(N)*(zeta2-log2q)+12.*H_m2(N)*(-zeta2+log2q+H_m1(N)*2.*log2)
  +12.*(3.*zeta2q+H_m2_2(N)+H_3_m1(N)+2.*H_2_m1_1(N)+H_3(N)*2.*log2-(H_m2_1(N)+H_m1_2(N))*2.*log2-8.*Li4)
  -63.*log2*zeta3+3.*H_m1(N)*(4.*H_m2_m1(N)+5.*zeta3)));
}

std::complex<long double> HSum::H_m1_m2_1(std::complex<long double> N){
  return(g15(N)-(19.*zeta2q/40.+H_m2_m1_1(N)+H_2_m1_m1(N)-Li4+1./24.*(-12.*zeta2*H_m1(N)*H_m1(N)-24.*zeta2*H_2(N)+24.*H_2_m1(N)*log2
  -6.*zeta2*log2q+12.*H_2(N)*log2q-log2q*log2q+2.*H_m2(N)*(6.*zeta2-6.*log2q)-6.*log2*zeta3+H_m1(N)*(-24.*zeta2*log2+15.*zeta3))));
}

std::complex<long double> HSum::H_2_m1_1(std::complex<long double> N){
  return(-H_1_2_m1(N)-H_2_1_m1(N)+H_3_m1(N)+H_2(N)*(H_m1_1(N)+H_1_m1(N))+H_2_m2(N)-H_2(N)*H_m2(N)-H_1(N)*H_m1_2(N)+H_1(N)*H_m3(N));
}

std::complex<long double> HSum::H_m1_m1_m2(std::complex<long double> N){
  return(0.5*(-H_m1_m2_m1(N)+H_2_m2(N)+H_m1(N)*H_m1_m2(N)+H_m1_3(N)));
}

std::complex<long double> HSum::H_1_m1_m2(std::complex<long double> N){
  return(-H_m1_m2_1(N)-H_m1_1_m2(N)+H_m1_m3(N)+H_m2_m2(N)+H_1(N)*H_m1_m2(N));
}

std::complex<long double> HSum::H_1_m2_m1(std::complex<long double> N){
  return(H_m1_m2_1(N)-H_m1(N)*H_m2_1(N)-H_m1_m3(N)+H_m1(N)*H_m3(N)+H_1(N)*H_m2_m1(N)+H_1_3(N)-H_1(N)*H_3(N));
}

std::complex<long double> HSum::H_1_m1_2(std::complex<long double> N){
  return(-H_1_2_m1(N)-H_2_1_m1(N)+H_1_m3(N)+H_3_m1(N)+H_2(N)*H_1_m1(N));
}

std::complex<long double> HSum::H_m1_1_2(std::complex<long double> N){
  return(H_1_2_m1(N)+H_2_1_m1(N)-H_m1_2_1(N)+H_m2_2(N)+H_m3_1(N)-H_m4(N)+H_m1_3(N)-H_3_m1(N)-H_2(N)*H_1_m1(N)+H_1(N)*H_m1_2(N)-H_1(N)*H_m3(N));
}

std::complex<long double> HSum::H_m1_m1_m1_m1(std::complex<long double> N){
  return (1./4.*H_4(N)+1./8.*H_2(N)*H_2(N)+1./3.*H_m3(N)*H_m1(N)+1./4.*H_2(N)*H_m1(N)*H_m1(N)+1./24.*H_m1(N)*H_m1(N)*H_m1(N)*H_m1(N));
}

std::complex<long double> HSum::H_m1_1_1_1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)/6.*g29(N)-Li4);
}

std::complex<long double> HSum::H_m1_1_1_m1(std::complex<long double> N){
  return(-1./6.*g31(N)+1./6.*(12.*zeta2q/5.-6.*H_m1_1_1(N)*log2+1./2.*6.*zeta2*log2q+H_m1_1_m1(N)*6.*log2
  -6.*Li4-log2*(log2q*(H_m1(N)+log2)+(H_m1_m1(N)-H_m1_1(N))*3.*log2+6.*zeta3)));
}

std::complex<long double> HSum::H_1_m1_1_m1(std::complex<long double> N){
  return(-1./20.*zeta2q-g32(N)-H_1_m1_1(N)*log2-0.5*H_1_m1(N)*log2q-1./6.*log2c*H_1(N)-1./24.*log2q*log2q+1./8.*zeta3*H_1(N)+1./8.*log2*zeta3);
}

std::complex<long double> HSum::H_m1_1_m1_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(pow(-1.,nn)*g33(N)+1./288.*(36.*zeta2q-288.*H_m1_1_m1(N)*log2-12.*6.*zeta2*log2q
  -48.*H_m1(N)*log2c-12.*log2q*log2q+24.*6.*H_m1_1(N)*(zeta2-log2q)
  +12.*6.*zeta2*H_m1(N)*2.*log2+12.*6.*zeta2*2.*log2q-252.*H_m1(N)*zeta3-252.*log2*zeta3));
}

std::complex<long double> HSum::H_m1_m1_1_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)*g35(N)-1./120.*(4.*36.*zeta2q+120.*H_m1_m1_1(N)*log2+15.*6.*zeta2*log2q
  +60.*H_m1_m1(N)*log2q+20.*H_m1(N)*log2c-10.*log2q*log2q-360.*Li4
  -15.*H_m1(N)*zeta3-330.*log2*zeta3));
}

std::complex<long double> HSum::H_1_1_m1_m1(std::complex<long double> N){
  return(g34(N)-1./24.*(24.*H_1_1_m1(N)*log2+4.*H_1(N)*log2c-12.*H_1_1(N)*(zeta2-log2q)
  -12.*zeta2*H_1(N)*log2-24.*Li4+21.*H_1(N)*zeta3));
}

std::complex<long double> HSum::H_1_1_1_m1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)/6.*g30(N)+log2*(H_1_1_m1(N)-H_1_1_1(N))+0.5*log2q*(H_1_1(N)-H_1_m1(N))-1./6.*log2c*(H_1(N)-H_m1(N))+1./24.*log2q*log2q);
}

std::complex<long double> HSum::H_1_1_1_1(std::complex<long double> N){
  return (1./4.*H_4(N)+1./8.*H_2(N)*H_2(N)+1./3.*H_3(N)*H_1(N)+1./4.*H_2(N)*H_1(N)*H_1(N)+1./24.*H_1(N)*H_1(N)*H_1(N)*H_1(N));
}

std::complex<long double> HSum::H_1_m1_m1_m1(std::complex<long double> N){
  long double frac,n=0.;
  frac=modf(std::real(N),&n);
  return(pow(-1.,n)*g36(N)-1./120.*(4.*36.*zeta2q+240.*H_m1_m1_1_m1(N)+10.*6.*zeta2*H_m1(N)*log2-10.*6.*zeta2*log2*H_1(N)
  -120.*H_m1_m1_m1(N)*log2+240.*H_m1_m1_1(N)*log2+120.*H_1_m1_m1(N)*log2+20.*6.*zeta2*log2q-20.*H_m1(N)*log2c+20.*H_1(N)*log2c
    -20.*log2q*log2q-10.*H_1_m1(N)*(6.*zeta2-6.*log2q)-360.*Li4-30.*H_m1(N)*zeta3+30.*H_1(N)*zeta3-345.*log2*zeta3));
}

std::complex<long double> HSum::H_1_1_m1_1(std::complex<long double> N){
  long double frac,nn=0.;
  frac=modf(std::real(N),&nn);
  return(-pow(-1.,nn)/2.*g37(N)-1./360.*(4.*36.*zeta2q+360.*H_m1_m1_1_m1(N)+360.*H_1_m1_m1_m1(N)
  +30*6.*zeta2*log2*(H_m1(N)-H_1(N))-360.*log2*(H_m1_m1_m1(N)-H_m1_m1_1(N)-H_1_m1_m1(N)+H_1_m1_1(N))
    +30.*6.*zeta2*log2q-180.*log2q*(H_m1_m1(N)-H_m1_1(N))-120.*log2c*(H_m1(N)-H_1(N))-60.*log2q*log2q
    -30.*6.*H_1_m1(N)*(zeta2-log2q)+30.*6.*H_1_1(N)*(zeta2-log2q)-360.*Li4-45.*H_m1(N)*zeta3+45.*H_1(N)*zeta3-360.*log2*zeta3)
  );
}

std::complex<long double> HSum::H_1_m1_1_1(std::complex<long double> N){
  return(1./6.*(H_m1(N)*H_1(N)*H_1(N)*H_1(N)+3.*H_1(N)*H_1(N)*H_m2(N)+3.*H_m1(N)*H_1(N)*H_2(N)
  +3.*H_m2(N)*H_2(N)+6.*H_1(N)*H_m3(N)+2.*H_m1(N)*H_3(N)+6.*H_m4(N))-H_m1_1_1_1(N)-H_1_1_m1_1(N)-H_1_1_1_m1(N)
  );
}

std::complex<long double> HSum::H_m1_m1_m1_1(std::complex<long double> N){
  return(1./6.*(H_1(N)*H_m1(N)*H_m1(N)*H_m1(N)+3.*H_m1(N)*H_m1(N)*H_m2(N)+3.*H_1(N)*H_m1(N)*H_2(N)
  +3.*H_m2(N)*H_2(N)+6.*H_m1(N)*H_3(N)+2.*H_1(N)*H_m3(N)+6.*H_m4(N))-H_1_m1_m1_m1(N)-H_m1_1_m1_m1(N)-H_m1_m1_1_m1(N)
  );
}

std::complex<long double> HSum::H_1_m1_m1_1(std::complex<long double> N){
  return(g38(N)-(2.*H_1_1_m1_m1(N)+1./12.*H_1_m1(N)*6.*(zeta2-log2q)+H_1_1(N)*(-zeta2+log2q)+H_1_1_m1(N)*log2*2.+Li4) );
}

std::complex<long double> HSum::H_m1_m1_1_1(std::complex<long double> N){
  return(0.5*g39(N)-1./48.*(24.*H_1_m1_m1_1(N)+24.*log2*(H_m1_m1_1(N)+H_1_m1_m1(N))+2.*6.*zeta2*log2q-4.*log2c*(H_m1(N)-H_1(N))
  -4.*log2q*log2q+2.*6.*H_1_m1(N)*(zeta2+log2q)-6.*2.*log2*zeta2*H_m1(N)+6.*zeta2*H_1(N)*2.*log2-72.*Li4
  +42.*H_m1(N)*zeta3-42.*H_1(N)*zeta3-21.*log2*zeta3)
  );
}

std::complex<long double> HSum::H_m1_1_m1_1(std::complex<long double> N){
  return(1./4.*(H_1(N)*H_1(N)*H_m1(N)*H_m1(N)+H_1(N)*H_1(N)*H_2(N)+H_m1(N)*H_m1(N)*H_2(N)+4.*H_m1(N)*H_1(N)*H_m2(N)
  +H_2(N)*H_2(N)+2.*H_m2(N)*H_m2(N)+4.*H_1(N)*H_3(N)+4.*H_m1(N)*H_m3(N)+6.*H_4(N))-H_m1_m1_1_1(N)-H_1_m1_m1_1(N)-H_1_1_m1_m1(N)
  -H_1_m1_1_m1(N)-H_m1_1_1_m1(N)
  );
}
















    