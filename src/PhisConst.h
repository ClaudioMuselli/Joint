#ifndef __CONST_H__
#define __CONST_H__


class ConstResum {
public:
  ConstResum(long double MUF, long double MUR, long double as, long double Q, long double sigma0Higgs);
  ConstResum(ConstResum & CC);
  long double _MUF;
  long double _MUR;
  long double _as;
  long double _Q;
  long double _sigma0Higgs;
private:
  long double Gf = 0.00001166364;// Fermi constant in GeV^-2
  long double GeVtopb = 389379304.;// GeV^-2 to pb conversion factor == (hc)^2 
};

#endif