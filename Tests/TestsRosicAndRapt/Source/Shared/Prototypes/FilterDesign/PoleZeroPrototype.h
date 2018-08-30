#pragma once

/** This class determines the locations of poles and zeros in the s-plane for a continuous time, 
unit cutoff lowpass or low-shelving prototype filter. It supports Butterworth, Chebychev, inverse 
Chebychev, elliptic, Bessel, Papoulis, Halpern and Gaussian designs. The low-shelving design is a 
generalization of a unit cutoff lowpass filter and and its magnitude response can be seen as the 
lowpass-response raised to a pedestal. The height of this pedestal is called the reference gain 
(which is zero in the lowpass-case). The gain the 'passband' (represented by the member variable G) 
is either a boost (when G > 1) or an attenuation/cut (when G < 1) - for standard lowpass designs, 
this gain would be unity.

References:
-(1) Sophocles J. Orfanidis: Lecture Notes on Elliptic Filter Design
-(2) Sophocles J. Orfanidis: High-Order Elliptical Equalizer Design
-(3) Larry D.Paarmann: Design and Analysis of Analog Filters */

template<class T>
class PoleZeroPrototype
{

public:

  /** This is an enumeration of the available approximation methods. */
  enum approximationMethods
  {
    //BUTTERWORTH = 1,   ///< maximally flat at DC
    BUTTERWORTH = 0,   ///< maximally flat at DC
    CHEBYCHEV,         ///< equiripple in passband, monotonic in stopband
    CHEBYCHEV_INV,     ///< equiripple in stopband, monotonic in passband
    ELLIPTIC,          ///< equiripple in passband and stopband, maximally steep transition
    BESSEL,            ///< approximates linear phase
    PAPOULIS,          ///< maximizes steepness at cutoff (selectivity) under constraint of monotonicity

    HALPERN,           ///< minimizes ratio of bandwidths at specified magnitudes (shaping factor) under constraint of monotonicity
                       ///< ...less steep at cutoff but steeper in stopband than Papoulis
    GAUSSIAN,          ///< smallest timelength*bandwidth product, good time response (no overshoot?)

    NUM_APPROXIMATION_METHODS
  };

  typedef std::complex<T> Complex;  // for convenience

  //-----------------------------------------------------------------------------------------------
  // \name static gain/pole/zero computation functions

  // lowpass prototypes (for allpole filters, zeros will be set to inf):
  // N: order, k: gain, p: poles, z: zeros, ripple and rejection is given in dB
  static void butterworth(size_t N, T* k, Complex* p, Complex* z);
  //static void bessel(     size_t N, T* k, Complex* p, Complex* z);
  //static void papoulis(   size_t N, T* k, Complex* p, Complex* z);
  //static void halpern(    size_t N, T* k, Complex* p, Complex* z);
  //static void gaussian(   size_t N, T* k, Complex* p, Complex* z);
  //static void chebychev(  size_t N, T* k, Complex* p, Complex* z, T ripple);
  //static void chebychev2( size_t N, T* k, Complex* p, Complex* z, T rejection);
  //static void elliptic(   size_t N, T* k, Complex* p, Complex* z, T ripple, T rejection);

  // low-shelving prototypes:
  // G0: reference gain, G: shelving-gain
  static void butterworth(size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  //static void bessel(     size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  //static void papoulis(   size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  //static void halpern(    size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  //static void gaussian(   size_t N, T G0, T G, T* k, Complex* p, Complex* z);
  //static void chebychev(  size_t N, T G0, T G, T* k, Complex* p, Complex* z, T ripple);
  //static void chebychev2( size_t N, T G0, T G, T* k, Complex* p, Complex* z, T ripple);
  //static void elliptic(   size_t N, T G0, T G, T* k, Complex* p, Complex* z, T ripple, T reject);

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  void setOrder(int newOrder) { order = newOrder; }
  void setApproximationMethod(int newMethod) { method = newMethod; }
  //void setShelvingGain(T newGain);
  //void setPassbandRipple(T newRipple);
  //void setStopbandRejection(T newRejection);

  //-----------------------------------------------------------------------------------------------
  // \name Retrieving Poles/Zeros/Gain

  void getPolesZerosAndGain(Complex* poles, Complex* zeros, T* gain);
    // should be a simple dispatcher function


protected:

  int method = BUTTERWORTH;
  int order  = 1;
  T ripple = 1;    // use Ap as member later (avoid recomputations)
  T reject = 60;   // use As as member later
  T G0 = T(0);
  T G  = T(1);     
  //T GB = RS_SQRT2_INV; // maybe later for more generality (may also avoid re-compuations of GB)
  // maybe rename to g0 (gain at 0), gi (gain at infinity), g1 (gain at 1)
  // or gZero, gOne, gInf


};

