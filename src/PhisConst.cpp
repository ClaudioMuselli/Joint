#include "PhisConst.h"

ConstResum::ConstResum(long double MUF, long double MUR, long double as, long double Q, long double sigma0Higgs){
  _MUF=MUF;_MUR=MUR;_as=as;_Q=Q;_sigma0Higgs=sigma0Higgs;
}
ConstResum::ConstResum(ConstResum &CC){
  _MUF=CC._MUF;_MUR=CC._MUR;_as=CC._as;_Q=CC._Q;_sigma0Higgs=CC._sigma0Higgs;
}