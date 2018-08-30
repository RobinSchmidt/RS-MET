#pragma once

template<class T>
class PoleZeroMapperAnalog
{

public:

  typedef std::complex<T> Cmplx;

  // N: order
  // pk,pp,pz: prototype gain/poles/zeros 
  // k,p,z: mapped gain7poles/zeros

  void protoToLowpass( size_t N, T pk, Cmplx* pp, Cmplx* pz, T* k, Cmplx* p, Cmplx* z, T wc);
  void protoToHighpass(size_t N, T pk, Cmplx* pp, Cmplx* pz, T* k, Cmplx* p, Cmplx* z, T wc);
  void protoToBandpass(size_t N, T pk, Cmplx* pp, Cmplx* pz, T* k, Cmplx* p, Cmplx* z, T wl, T wu);
  void protoToBandstop(size_t N, T pk, Cmplx* pp, Cmplx* pz, T* k, Cmplx* p, Cmplx* z, T wl, T wu);

protected:

  //int mode = LOWPASS;
  //T wl = 1;
  //T wu = 2;

};

