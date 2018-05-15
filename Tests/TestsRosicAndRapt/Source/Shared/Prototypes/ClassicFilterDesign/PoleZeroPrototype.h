#pragma once

/** A class for computing poles and zeros (and gain factors) for analog unit cutoff lowpass or 
low-shelving prototype filters. */

template<class T>
class PoleZeroPrototype
{

public:

  typedef std::complex<T> Complex;  // for convenience

  //-----------------------------------------------------------------------------------------------
  // \name static gain/pole/zero computation functions

  // lowpass prototypes (for allpole filters, zeros will be set to inf):
  // N: order, k: gain, p: poles, z: zeros, ripple and rejection is given in dB
  static void butterworth(size_t N, T* k, Complex* p, Complex* z);
  static void bessel(     size_t N, T* k, Complex* p, Complex* z);
  static void papoulis(   size_t N, T* k, Complex* p, Complex* z);
  static void halpern(    size_t N, T* k, Complex* p, Complex* z);
  static void gaussian(   size_t N, T* k, Complex* p, Complex* z);
  static void chebychev(  size_t N, T* k, Complex* p, Complex* z, T ripple);
  static void chebychev2( size_t N, T* k, Complex* p, Complex* z, T rejection);
  static void elliptic(   size_t N, T* k, Complex* p, Complex* z, T ripple, T rejection);

  // low-shelving prototypes:
  // G0: reference gain, G: shelving-gain
  static void butterworth(size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  static void bessel(     size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  static void papoulis(   size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  static void halpern(    size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  static void gaussian(   size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  static void chebychev(  size_t N, T G0, T G, T* k, Complex* p, Complex* z, T ripple);
  static void chebychev2( size_t N, T G0, T G, T* k, Complex* p, Complex* z, T ripple);
  static void elliptic(   size_t N, T G0, T G, T* k, Complex* p, Complex* z, T ripple);

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setOrder(size_t newOrder);
  void setShelvingGain(T newGain);
  void setPassbandRipple(T newRipple);
  void setStopbandRejection(T newRejection);

  //-----------------------------------------------------------------------------------------------
  // \name Retrieving Poles/Zeros/Gain

  void getPolesZerosAndGain(Complex* poles, Complex* zeros, T* gain);
    // should be a simple dispatcher function


protected:

  //int method = BUTTERWORTH;
  T G0 = T(0);
  T G  = T(1);

};

// maybe the order of the functions should be be bessel, gaussian, butterworth, papoulis, halpern, 
// chebychev2, chebychev, elliptic - from time-domain to frequency-domain superiority (roughly)

// class should return only upper-left quarter-plane poles and zeros, their number is (order+1)/2
// using integer division