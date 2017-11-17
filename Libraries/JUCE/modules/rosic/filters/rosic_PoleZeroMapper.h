#ifndef rosic_PoleZeroMapper_h
#define rosic_PoleZeroMapper_h

namespace rosic
{

/** This class contains static functions for converting a set poles and zeros of some prototype 
filter in the complex s-plane (or z-plane) to another set of poles and zeros of the corresponding 
unnormalized filter. Additionaly, it also implements a mapping between the s-plane and z-plane via 
the bilinear transform. */

class rsPoleZeroMapper
{

public:

  /** Transforms an analog lowpass prototype filter represented by its zeros, poles and gain (in 
  "z", "p", "k") into a corresponding low-shelving filter with reference gain G0 and low-frequency 
  gain G. */
  static void sLowpassToLowshelf(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, 
    double *kNew, int N, double G0, double G);

  /** Transforms an analog lowpass prototype filter represented by its zeros, poles and gain (in 
  "z", "p", "k") with nominal radian cutoff frequency of unity to an analog lowpass with nominal 
  radian cutoff frequency of wc. The new zeros, poles and gain factor are returned in "zNew", 
  "pNew" and "kNew". */
  static void sLowpassToLowpass(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, 
    double *kNew, int N, double wc);

  /** Transforms an analog lowpass prototype filter to an analog highpass with nominal radian 
  cutoff frequency of wc.  @see sLowpassToLowpass */
  static void sLowpassToHighpass(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, 
    double *kNew, int N, double wc);

  /** Transforms an analog lowpass prototype filter to an analog bandpass with lower and upper 
  nominal radian bandedge frequencies of wl and wu. @see sLowpassToLowpass */
  static void sLowpassToBandpass(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, 
    double *kNew, int N, double wl, double wu);

  /** Transforms an analog lowpass prototype filter to an analog bandreject with lower and upper 
  nominal radian bandedge frequencies of wl and wu. @see sLowpassToLowpass */
  static void sLowpassToBandreject(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, 
    double *kNew, int N, double wl, double wu);

  //static void zLowpassToLowpass(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, 
  //  double *kNew, int N, double wc);

  /** Performs an s-plane lowpass-to-lowpass transform. This transform is done in place which 
  means, the result is stored in the same arrays as the input. These arrays must be of length 
  numPoles and numZeros respectively. The target cutoff frequency is assumed to be expressed as 
  radian frequency (Omega = 2*pi*frequency). ...to be deprecated. */
  static void prototypeToAnalogLowpass(Complex *poles, int numPoles, Complex *zeros, int numZeros, 
    double *gain, double targetCutoff);

  /** \todo comment this function */
  static void sPlanePrototypeToLowpass(Complex *prototypePoles, Complex *prototypeZeros, 
    Complex *targetPoles, Complex *targetZeros, int prototypeOrder, double targetCutoff);

  /** Performs an s-plane lowpass-to-highpass transform. This transform is done in place which 
  means, the result is stored in the same arrays as the input. These arrays must BOTH be of length 
  max(numPoles, numZeros) because we get additional finite zeros at s=0 for zeros at s=inf in the 
  prototype. ...to be deprecated. */
  static void prototypeToAnalogHighpass(Complex *poles, int numPoles, Complex *zeros, int numZeros, 
    double *gain, double targetCutoff);

  /** \todo comment this function */
  static void sPlanePrototypeToHighpass(Complex *prototypePoles, Complex *prototypeZeros, 
    Complex *targetPoles, Complex *targetZeros, int prototypeOrder, double targetCutoff);

  /** Same as prototypeToAnalogLowpass but additionally exchanges the roles of poles and zeros. */
  static void prototypeToAnalogHighShelv(Complex *poles, int numPoles, Complex *zeros, 
    int numZeros, double *gain, double targetCutoff);

  static void sPlanePrototypeToHighShelv(Complex *prototypePoles, Complex *prototypeZeros, 
    Complex *targetPoles, Complex *targetZeros, int prototypeOrder, double targetCutoff);

  /** Performs an s-plane lowpass-to-bandpass transform. This transform is done in place which 
  means, the result is stored in the same arrays as the input. These arrays must BOTH be of length
  2*max(numPoles, numZeros) because each finite pole and zero maps to a pair of finite poles and 
  zeros and we get additional finite zeros for zeros at s=inf in the prototype. */
  static void prototypeToAnalogBandpass(Complex *poles, int numPoles, Complex *zeros, int numZeros, 
    double *gain, double targetLowCutoff, double targetHighCutoff);

  static void sPlanePrototypeToBandpass(Complex *prototypePoles, Complex *prototypeZeros, 
    Complex *targetPoles, Complex *targetZeros, int prototypeOrder, double targetLowerCutoff, 
    double targetUpperCutoff);

  /** Performs an s-plane lowpass-to-bandpass transform. This transform is done in place which 
  means, the result is stored in the same arrays as the input. These arrays must BOTH be of length 
  2*max(numPoles, numZeros) because each finite pole and zero maps to a pair of finite poles and 
  zeros and we get additional finite zeros for zeros at s=inf in the prototype. */
  static void prototypeToAnalogBandstop(Complex *poles, int numPoles, Complex *zeros, int numZeros, 
    double *gain, double targetLowCutoff, double targetHighCutoff);

  static void sPlanePrototypeToBandreject(Complex *prototypePoles, Complex *prototypeZeros, 
    Complex *targetPoles, Complex *targetZeros, int prototypeOrder, double targetLowerCutoff, 
    double targetUpperCutoff);

  /** Maps zeros, poles and gain factor of a digital prototype lowpass with normalized radian 
  cutoff frequency "wc" to the new zeros, poles and gain factor of a digital targer lowpass with 
  normalized radian frequency "wt". */
  static void zLowpassToLowpass(Complex *z, Complex *p, double *k, Complex *zNew, Complex *pNew, 
    double *kNew, int N, double wc, double wt);

  /*
  void analogToDigitalLowpassToLowpass(Complex *poles, int numPoles,
  Complex *zeros, int numZeros, double *gain, double targetCutoff);
  */

  /** Performs a bilinear transform from the s-plane to the z-plane. */
  static void bilinearAnalogToDigital(Complex *poles, int numPoles, Complex *zeros, int numZeros, 
    double sampleRate, double *gain);

  /** Performs a bilinear transform from the z-plane to the s-plane. */
  //void bilinearDigitalToAnalog(Complex *poles, int numPoles, Complex *zeros, int numZeros, 
  //  double sampleRate, double *gain);

  // make function sLowpassToLowshelf
  // include a gain factor k into the mapping functions -> include it also into the 
  // PrototypeDesigner update function names sLowpassToLowpass etc.
  // include different s-to-z mappings (matched z-tansform (aka impulse invariant), 
  // magnitude matched, phase matched etc.)


protected:

  /** Transforms an array of roots (poles or zeros) of an analog lowpass prototype filter to 
  corresponding roots of an analog highpass with nominal radian cutoff frequency of wc. */
  static void sLowpassToHighpass(Complex *r, Complex *rNew, int N, double wc);

  /** Transforms an array of roots (poles or zeros) of an analog lowpass prototype filter to 
  corresponding roots of an analog bandpass with nominal radian center frequency of wc and 
  bandwidth bw. */
  static void sLowpassToBandpass(Complex *r, Complex *rNew, int N, double wc, double bw);

  /** Transforms an array of roots (poles or zeros) of an analog lowpass prototype filter to 
  corresponding roots of an analog bandreject with nominal radian center frequency of wc and 
  bandwidth bw. */
  static void sLowpassToBandreject(Complex *r, Complex *rNew, int N, double wc, double bw);

};

}

#endif
