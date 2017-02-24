#include "AP.h"





GammaAP::GammaAP(const long double  Nf, const long double Nc): HAP(false,false,false){
  NF=Nf;
  CA=Nc;
  CF=(Nc*Nc-1.)/(2.*Nc); 
}
GammaAP::~GammaAP(){
}

void GammaAP::ComputeGamma(std::complex<long double> N, int order){
  if(order>=0){
    gg0=gammagg0(N);
    gS0=gammagS0(N);
    Sg0=gammaSg0(N);
    WW0=gammansplus0(N);
    SS0=WW0;
    TT0=WW0;
    VV0=WW0;
    plus0=0.5*(gg0+SS0+std::sqrt((gg0-SS0)*(gg0-SS0)+4.*Sg0*gS0));
    minus0=0.5*(gg0+SS0-std::sqrt((gg0-SS0)*(gg0-SS0)+4.*Sg0*gS0));
    
  }
  if(order>=1){
    sums(N);
    DEF1(N);
    gg1=gammagg1(N);
    gS1=gammagS1(N);
    Sg1=gammaSg1(N);
    WW1=gammansplus1(N);
    SS1=WW1+gammaps1(N);
    TT1=gammansminus1(N);
    VV1=TT1;
  }
  if(order>=2){
    DEF2(N);
    gg2=gammagg2(N);
    gS2=gammagS2(N);
    Sg2=gammaSg2(N);
    WW2=gammansplus2(N);
    SS2=WW2+gammaps2(N); 
    TT2=gammansminus2(N);
    VV2=gammansval2(N);
  }
  return;
}

// Analytic continuations of the occuring sums
void GammaAP::sums(std::complex<long double> N) {
  //
  g1s.S1  = HAP.HS(1,N);
  g1s.S2  = HAP.HS(2,N);
  g1s.S3 = HAP.HS(3,N);
  g1s.S4 = HAP.HS(4,N);
  //
  g1s.NI = 1./N;
  g1s.NI2 = g1s.NI*g1s.NI;
  g1s.NI3 = g1s.NI*g1s.NI2;
  g1s.NS  = N * N;
  g1s.NT  = g1s.NS * N;
  g1s.NFO = g1s.NT * N;
  g1s.NFI = g1s.NFO * N;
  g1s.NSI = g1s.NFI * N;
  g1s.NSE = g1s.NSI * N;
  g1s.NE  = g1s.NSE * N;
  g1s.NN  = g1s.NE * N;
  // 
  g1s.NM  = N - 1.;
  g1s.NMI = 1./g1s.NM;
  g1s.NMI2 = g1s.NMI*g1s.NMI;
  g1s.N1  = N + 1.;
  g1s.N1I = 1./g1s.N1;
  g1s.N1I2 = g1s.N1I*g1s.N1I;
  g1s.N1I3 = g1s.N1I*g1s.N1I2;  
  g1s.N2  = N + 2.;
  g1s.N2I = 1./g1s.N2;
  g1s.NMS = g1s.NM * g1s.NM;
  g1s.N1S = g1s.N1 * g1s.N1;
  g1s.N1T = g1s.N1S * g1s.N1;
  g1s.N2S = g1s.N2 * g1s.N2;
  g1s.N2T = g1s.N2S * g1s.N2;
  //
  g1s.N3  = N + 3.;
  g1s.N4  = N + 4.;
  g1s.N5  = N + 5.;
  g1s.N6  = N + 6.;
  g1s.S1M = g1s.S1 - g1s.NI;
  g1s.S2M = g1s.S2 - g1s.NI2;
  g1s.S11 = g1s.S1  + 1./g1s.N1;
  g1s.S12 = g1s.S11 + 1./g1s.N2;
  g1s.S13 = g1s.S12 + 1./g1s.N3;
  g1s.S14 = g1s.S13 + 1./g1s.N4;
  g1s.S15 = g1s.S14 + 1./g1s.N5;
  g1s.S16 = g1s.S15 + 1./g1s.N6;
  g1s.S21 = g1s.S2 + g1s.N1I2;
  g1s.S31 = g1s.S3 + g1s.N1I3;
  g1s.SPMOM = ( 1.0000 * (gsl_sf_zeta(2) - g1s.S1 / N ) / N  -
		0.9992 * (gsl_sf_zeta(2) - g1s.S11/ g1s.N1) / g1s.N1 +
		0.9851 * (gsl_sf_zeta(2) - g1s.S12/ g1s.N2) / g1s.N2 -
		0.9005 * (gsl_sf_zeta(2) - g1s.S13/ g1s.N3) / g1s.N3 +
		0.6621 * (gsl_sf_zeta(2) - g1s.S14/ g1s.N4) / g1s.N4 -
		0.3174 * (gsl_sf_zeta(2) - g1s.S15/ g1s.N5) / g1s.N5 +
		0.0699 * (gsl_sf_zeta(2) - g1s.S16/ g1s.N6) / g1s.N6 );
  //
  g1s.SLC = - 5./8. * gsl_sf_zeta(3);
  g1s.SLV = - gsl_sf_zeta(2)/2.* (HAP.HS(1,g1s.N1/2.-1.) - HAP.HS(1,N/2.-1.)) + g1s.S1/g1s.NS + g1s.SPMOM;
  g1s.SSCHLM = g1s.SLC - g1s.SLV;
  g1s.SSTR2M = HAP.HS(2,g1s.N1/2.-1.);
  g1s.SSTR3M = HAP.HS(3,g1s.N1/2.-1.);
  //
  g1s.SSCHLP = g1s.SLC + g1s.SLV;
  g1s.SSTR2P = HAP.HS(2,g1s.N2/2.-1.);
  g1s.SSTR3P = HAP.HS(3,g1s.N2/2.-1.);
  //
  return;
}

//LO anomalous dimensions
std::complex<long double> GammaAP::gammagg0(std::complex<long double> N){
  return (-1./(4.*M_PIl)*(CA*(4.*(HAP.HS(1,N-2.)-2.*HAP.HS(1,N-1.)-2.*HAP.HS(1,N+1.)+HAP.HS(1,N+2.)+3.*HAP.HS(1,N))-11./3.)+2./3.*NF));
  
}
std::complex<long double> GammaAP::gammagS0(std::complex<long double> N){
  return -1./(4.*M_PIl)*(2.*CF*(2.*HAP.HS(1,N-2.)-4.*HAP.HS(1,N-1.)-HAP.HS(1,N+1.)+3.*HAP.HS(1,N)));
}
std::complex<long double> GammaAP::gammaSg0(std::complex<long double> N){
  return -1./(4.*M_PIl)*(2.*NF*(HAP.HS(1,N-1.)+4.*HAP.HS(1,N+1.)-2.*HAP.HS(1,N+2)-3.*HAP.HS(1,N)));
}
std::complex<long double> GammaAP::gammansplus0(std::complex<long double> N){
  return -1./(4.*M_PIl)*CF*(2.*(HAP.HS(1,N-1.)+HAP.HS(1,N+1.))-3.);
}

//NLO important definition
void GammaAP::DEF1(std::complex<long double> N) {
  PNPA = ( 16.* g1s.S1 * (2.* N + 1.) / (g1s.NS * g1s.N1S) +
	   16.* (2.* g1s.S1 - 1./(N * g1s.N1)) * ( g1s.S2 - g1s.SSTR2P ) +
	   64.* g1s.SSCHLP + 24.* g1s.S2 - 3. - 8.* g1s.SSTR3P -
	   8.* (3.* g1s.NT + g1s.NS -1.) / (g1s.NT * g1s.N1T) -
	   16.* (2.* g1s.NS + 2.* N +1.)/(g1s.NT * g1s.N1T) ) * (-0.5);
  PNMA = ( 16.* g1s.S1 * (2.* N + 1.) / (g1s.NS * g1s.N1S) +
	   16.* (2.* g1s.S1 - 1./(N * g1s.N1)) * ( g1s.S2 - g1s.SSTR2M ) +
	   64.* g1s.SSCHLM + 24.* g1s.S2 - 3. - 8.* g1s.SSTR3M -
	   8.* (3.* g1s.NT + g1s.NS -1.) / (g1s.NT * g1s.N1T) +
	   16.* (2.* g1s.NS + 2.* N +1.) / (g1s.NT * g1s.N1T) ) * (-0.5);
  PNSB = ( g1s.S1 * (536./9. + 8.* (2.* N + 1.) / (g1s.NS * g1s.N1S)) -
	   (16.* g1s.S1 + 52./3.- 8./(N * g1s.N1)) * g1s.S2 - 43./6. -
	   (151.* g1s.NFO + 263.* g1s.NT + 97.* g1s.NS + 3.* N + 9.) *
	   4./ (9.* g1s.NT * g1s.N1T) ) * (-0.5);
  PNSC = ( -160./9.* g1s.S1 + 32./3.* g1s.S2 + 4./3. +
	   16.*(11.*g1s.NS+5.*N-3.)/(9.* g1s.NS * g1s.N1S))*(-0.5);

	   
  PPSA = (5.* g1s.NFI + 32.* g1s.NFO + 49.* g1s.NT+38.* g1s.NS + 28.* N + 8.) 
    / (g1s.NM * g1s.NT * g1s.N1T * g1s.N2S) * 2.;
    
  PQGA = (-2.* g1s.S1 * g1s.S1 + 2.* g1s.S2 - 2.* g1s.SSTR2P) 
    * (g1s.NS + N + 2.) / (N * g1s.N1 * g1s.N2) 
    + (8.* g1s.S1 * (2.* N + 3.)) / (g1s.N1S * g1s.N2S)
    + 2.* (g1s.NN + 6.* g1s.NE + 15.* g1s.NSE + 25.* g1s.NSI + 36.* g1s.NFI
	   + 85.* g1s.NFO + 128.* g1s.NT + 104.* g1s.NS + 64.* N + 16.)
    / (g1s.NM * g1s.NT * g1s.N1T * g1s.N2T);
  PQGB = (2.* g1s.S1 * g1s.S1 - 2.* g1s.S2 + 5.) * (g1s.NS + N + 2.)
    / (N * g1s.N1 * g1s.N2) - 4.* g1s.S1 / g1s.NS
    + (11.* g1s.NFO + 26.* g1s.NT + 15.* g1s.NS + 8.* N + 4.)
    / (g1s.NT * g1s.N1T * g1s.N2) ; 
    
  PGQA = (- g1s.S1 * g1s.S1 + 5.* g1s.S1 - g1s.S2) * (g1s.NS + N + 2.) 
    / (g1s.NM * N * g1s.N1)  -  2.* g1s.S1 / g1s.N1S
    - (12.* g1s.NSI + 30.* g1s.NFI + 43.* g1s.NFO + 28.* g1s.NT - g1s.NS
       - 12.* N - 4.) / (2.* g1s.NM * g1s.NT * g1s.N1T) ;
  PGQB = (g1s.S1*g1s.S1 + g1s.S2 - g1s.SSTR2P) * (g1s.NS + N + 2.) / (g1s.NM * N * g1s.N1)
    - g1s.S1 * (17.* g1s.NFO + 41.* g1s.NS - 22.* N - 12.) 
    / (3.* g1s.NMS * g1s.NS * g1s.N1)
    + (109.* g1s.NN + 621.* g1s.NE + 1400.* g1s.NSE + 1678.* g1s.NSI
       + 695.* g1s.NFI - 1031.* g1s.NFO - 1304.* g1s.NT - 152.* g1s.NS
       + 432.* N + 144.) / (9.* g1s.NMS * g1s.NT * g1s.N1T * g1s.N2S);
  PGQC = (g1s.S1 - 8./3.) * (g1s.NS + N + 2.) / (g1s.NM * N * g1s.N1) + 1./ g1s.N1S;
  PGQC *= 4./3.;
  
  PGGA = - (2.* g1s.NFI + 5.* g1s.NFO + 8.* g1s.NT + 7.* g1s.NS- 2.* N - 2.)
    * 8.* g1s.S1 / (g1s.NMS * g1s.NS * g1s.N1S * g1s.N2S) -  67./9.* g1s.S1 + 8./3.
    - 4.* g1s.SSTR2P * (g1s.NS + N + 1.) / (g1s.NM * N * g1s.N1 * g1s.N2)
    + 2.* g1s.S1 * g1s.SSTR2P - 4.* g1s.SSCHLP + 0.5 * g1s.SSTR3P
    + (457.* g1s.NN + 2742.* g1s.NE + 6040.* g1s.NSE + 6098.* g1s.NSI
       + 1567.* g1s.NFI - 2344.* g1s.NFO - 1632.* g1s.NT + 560.* g1s.NS
       + 1488.* N + 576.) / (18.* g1s.NMS * g1s.NT * g1s.N1T * g1s.N2T);
  PGGB = (38.* g1s.NFO + 76.* g1s.NT + 94.* g1s.NS + 56.* N + 12.) *(-2.)
    / (9.* g1s.NM * g1s.NS * g1s.N1S * g1s.N2)  +  20./9.* g1s.S1  -  4./3.;
  PGGC= (2.*g1s.NSI + 4. * g1s.NFI + g1s.NFO - 10.* g1s.NT - 5. * g1s.NS - 4.* N
	 - 4.) * (-2.) / (g1s.NM * g1s.NT * g1s.N1T * g1s.N2) -1. ;  
}


//NLO anomalous dimensions
std::complex<long double> GammaAP::gammagg1(std::complex<long double> N){
  return (CA*CA*PGGA + 0.5*NF*(CA*PGGB+CF*PGGC))*4. /pow(4*M_PIl,2);
}
std::complex<long double> GammaAP::gammagS1(std::complex<long double> N){
  return (CF*CF*PGQA + CF*CA*PGQB+0.5*NF*CF*PGQC)*4. /pow(4*M_PIl,2);
}
std::complex<long double> GammaAP::gammaSg1(std::complex<long double> N){
  return 0.5*NF * (CA * PQGA + CF * PQGB)*4. /pow(4*M_PIl,2);
}
std::complex<long double> GammaAP::gammaps1(std::complex<long double> N){
  return  (0.5*NF*CF*PPSA*4. /pow(4*M_PIl,2));
}
std::complex<long double> GammaAP::gammansplus1(std::complex<long double> N){
  return CF *((CF-CA/2.)* PNPA + CA* PNSB + 0.5*(NF)* PNSC) /pow(4*M_PIl,2);
}
std::complex<long double> GammaAP::gammansminus1(std::complex<long double> N){
  return CF *((CF-CA/2.)* PNMA + CA* PNSB + 0.5*NF* PNSC) /(pow(4*M_PIl,2));
}

//NNLO important definitions
void GammaAP::DEF2(std::complex<long double> N){
  std::complex<long double> A0  = - g1s.S1M;
  std::complex<long double> B1  = - g1s.S1 * g1s.NI;
  std::complex<long double> B1M;
  if (( abs(imag(N)) < (1.-5.))&& (abs(real(N)-1.)<(1.-5.)))
         B1M = - gsl_sf_zeta(2);
  else
     B1M = - g1s.S1M * g1s.NMI;
  std::complex<long double> B11 = - g1s.S11 * g1s.N1I;
  std::complex<long double> B12 = - g1s.S12 * g1s.N2I;
  std::complex<long double> B2  = (pow(g1s.S1,2) + g1s.S2) * g1s.NI;
  std::complex<long double> B2M = (pow(g1s.S1M,2) + g1s.S2M) * g1s.NMI;
  std::complex<long double> B21 = (pow(g1s.S11,2) + g1s.S21) * g1s.N1I;
  std::complex<long double> B3  = - (pow(g1s.S1,3) + 3.*g1s.S1*g1s.S2 + 2.*g1s.S3) * g1s.NI;
  std::complex<long double> B31 = - (pow(g1s.S11,3) + 3.*g1s.S11*g1s.S21 + 2.*g1s.S31) * g1s.N1I;
  std::complex<long double> B4  = (pow(g1s.S1,4) + 6.*pow(g1s.S1,2)*g1s.S2 + 8.*g1s.S1*g1s.S3 + 3.*pow(g1s.S2,2) + 6.*g1s.S4) * g1s.NI;
  std::complex<long double> C0 = g1s.NI;
  std::complex<long double> CM = g1s.NMI;
  std::complex<long double> C1 = g1s.N1I;
  std::complex<long double> C2 = g1s.N2I;
  std::complex<long double> C3 = 1./(N+3.);
  std::complex<long double> C4 = 1./(N+4.);
  std::complex<long double> D1  = - g1s.NI2;
  std::complex<long double> D1M = - g1s.NMI2;
  std::complex<long double> D11 = - g1s.N1I2;
  std::complex<long double> D2  = 2.* g1s.NI3;
  std::complex<long double> D21 = 2.* g1s.N1I3;
  std::complex<long double> D3  = - 6.* g1s.NI2*g1s.NI2;
  std::complex<long double> D31 = - 6.* g1s.N1I2*g1s.N1I2;
  std::complex<long double> D4  = 24.* g1s.NI2*g1s.NI3;
  std::complex<long double> D41 = 24.* g1s.N1I2*g1s.N1I3;
  std::complex<long double> E1  = g1s.S1*g1s.NI2 + (g1s.S2-gsl_sf_zeta(2))*g1s.NI;
  std::complex<long double> E11 = g1s.S11*g1s.N1I2 + (g1s.S21-gsl_sf_zeta(2))*g1s.N1I;
  std::complex<long double> E2  = 2.* ( - g1s.S1*g1s.NI3 + (gsl_sf_zeta(2)-g1s.S2)*g1s.NI2 - (g1s.S3-gsl_sf_zeta(3))*g1s.NI );
  std::complex<long double> PP2 = + 1174.898 * A0 + 1295.384 + 714.1 * B1 - 522.1 * C3
            + 243.6 * C2 - 3135.* C1 + 1641.1 * C0 + 1258.* D1
            + 294.9 * D2 + 800/27. * D3 + 128/81. * D4
            + 563.9 * E1 + 256.8 * E2
          + NF * ( - 183.187 * A0 - 173.924 - 5120/81. * B1
            + 44.79 * C3 + 72.94 * C2 + 381.1 * C1 - 197.0 * C0
            - 152.6 * D1 - 2608./81. * D2 - 192./81. * D3
           - 56.66 * E1 - 1.497 * D31 );
  std::complex<long double> PM2 = + 1174.898 * A0 + 1295.470 + 714.1 * B1 - 433.2 * C3
            + 297.0 * C2 - 3505.*C1 + 1860.2 * C0 + 1465.2 * D1
            + 399.2 * D2 + 320./9. * D3 + 116./81. * D4
            + 684.0 * E1 + 251.2 * E2
          + NF * ( - 183.187 * A0 - 173.933 - 5120/81. * B1
            + 34.76 * C3 + 77.89 * C2 + 406.5 * C1 - 216.62 * C0
            - 172.69 * D1 - 3216./81 * D2 - 256./81. * D3
            - 65.43 * E1 - 1.136 * D31 );
   std::complex<long double> PSS2 = - 163.9 * (B1M-B1)-7.208 * (B11-B12) + 4.82 * (C3-C4)
            - 43.12 * (C2-C3) + 44.51 * (C1-C2) + 151.49 * (C0-C1)
            + 178.04 * D1 + 6.892 * D2 - 40./27. * (2.*D3 - D4)
            - 173.1 * E1 + 46.18 * E2;
   std::complex<long double> PF2 = - ( 17./72. - 2./27. * g1s.S1 - 10./27. * g1s.S2
       + 2./9. * g1s.S3 - (12.* pow(N,4) + 2.* pow(N,3) - 12.* pow(N,2)
     - 2.* N + 3.)/(27.* pow(N,3) * pow(g1s.N1,3)) ) * 32./3.;
   P2PLSN = PP2 + pow(NF,2) * PF2;
   P2MINN = PM2 + pow(NF,2)* PF2;
   P2VALN = P2MINN + NF * PSS2;
   std::complex<long double> PS1 = - 3584./27.* (D1M-D1) - 506.* (CM-C0) + 160./27.* (D4-D41)
            - 400./9.* (D3-D31) + 131.4 * (D2-D21) - 661.6 * (D1-D11)
            - 5.926 * (B3-B31) - 9.751 * (B2-B21) - 72.11 * (B1-B11)
            + 177.4 * (C0-C1) + 392.9 * (C1-C2) - 101.4 * (C2-C3)
            - 57.04 * (E1-E11);
  std::complex<long double> PS2 = 256./81.* (CM-C0) + 32./27.* (D3-D31) + 17.89 * (D2-D21)
            + 61.75 * (D1-D11) + 1.778 * (B2-B21) + 5.944 * (B1-B11)
            + 100.1 * (C0-C1) - 125.2 * (C1-C2) + 49.26 * (C2-C3)
            - 12.59 * (C3-C4) - 1.889 * (E1-E11);
  std::complex<long double> QG1 = - 896./3.* D1M - 1268.3 * CM + 536./27.* D4 - 44./3.* D3 
            + 881.5 * D2 + 424.9 * D1 + 100./27.* B4 - 70./9.* B3 
            - 120.5 * B2 + 104.42 * B1 + 2522.* C0 - 3316. * C1 
            + 2126. * C2 + 1823. * E1 - 25.22 * E2 - 252.5 * D31;
  std::complex<long double> QG2 =   1112./243.* CM - 16./9.* D4 - 376./27.* D3 - 90.8 * D2
            - 254.0 * D1 + 20./27.* B3 + 200./27.* B2 - 5.496 * B1
            - 252.0 * C0 + 158.0 * C1 + 145.4 * C2 - 139.28 * C3
            - 53.09 * E1 - 80.616 * E2 - 98.07 * D21 + 11.70 * D31;
  std::complex<long double> GQ0 = 1189.3 * D1M + 6163.1 * CM - 4288./81. * D4 + 1568./9.* D3
            - 1794.* D2 + 4033.* D1 + 400./81.* B4 + 2200./27.* B3
            + 606.3 * B2 + 2193.* B1 - 4307.* C0 + 489.3 * C1 
            + 1452.* C2 + 146.* C3 - 447.3 * E2 - 972.9 * D21; 
  std::complex<long double> GQ1 =   71.082 * D1M - 46.41 * CM + 128./27.* D4 + 704/81.* D3
            + 20.39 * D2 + 174.8 * D1 - 400./81.* B3 - 68.069 * B2 
            - 296.7 * B1 - 183.8 * C0 + 33.35 * C1 - 277.9 * C2 
            + 108.6 * D21 - 49.68 * E1;
  std::complex<long double> GQ2 = ( 64.* (- CM + C0 + 2.* C1) + 320.* (B1M - B1 + 0.8 * B11)
            + 96.* (B2M - B2 + 0.5 * B21) ) / 27. ;
  std::complex<long double> GG0 = 2675.8 * D1M + 14214.* CM - 144. * D4 + 72. * D3
            - 7471.* D2 + 274.4 * D1 - 20852.* C0 + 3968.* C1 
            - 3363.* C2 + 4848.* C3 + 7305.* E1 + 8757.* E2 
            + 3589.* B1 + 4425.894 + 2643.521 * A0; 
  std::complex<long double> GG1 =   157.27 * D1M + 182.96 * CM + 512./27. * D4 
            + 832./9. * D3 + 491.3 * D2 + 1541.* D1 - 350.2 * C0
            + 755.7 * C1 - 713.8 * C2 + 559.3 * C3 + 26.15 * E1
            - 808.7 * E2 - 320. * B1 - 528.723 - 412.172 * A0;
  std::complex<long double> GG2 = - 680./243. * CM - 32./27. * D3 + 9.680 * D2 
            - 3.422 * D1 - 13.878 * C0 + 153.4 * C1 - 187.7 * C2
            + 52.75 * C3 - 115.6 * E1 + 85.25 * E11 - 63.23 * E2
            + 6.4630 - 16./9. * A0;
  P2PSN =       NF * ( PS1 + NF * PS2 );
  P2QGN = NF * ( QG1 + NF * QG2 );
  P2GQN = GQ0 + NF * ( GQ1 + NF * GQ2 );
  P2GGN = GG0 + NF * ( GG1 + NF * GG2 );
}

//NNLO anomalous dimensions
std::complex<long double> GammaAP::gammagg2(std::complex<long double> N){
  return P2GGN/pow(4*M_PIl,3);
}
std::complex<long double> GammaAP::gammagS2(std::complex<long double> N){
  return P2GQN/pow(4*M_PIl,3);
}
std::complex<long double> GammaAP::gammaSg2(std::complex<long double> N){
  return P2QGN/pow(4*M_PIl,3);
}
std::complex<long double> GammaAP::gammaps2(std::complex<long double> N){
  return P2PSN/pow(4*M_PIl,3);
}
std::complex<long double> GammaAP::gammansplus2(std::complex<long double> N){
  return P2PLSN/pow(4*M_PIl,3);
}
std::complex<long double> GammaAP::gammansminus2(std::complex<long double> N){
  return P2MINN/pow(4*M_PIl,3);
}
std::complex<long double> GammaAP::gammansval2(std::complex<long double> N){
  return P2VALN/pow(4*M_PIl,3);
}











