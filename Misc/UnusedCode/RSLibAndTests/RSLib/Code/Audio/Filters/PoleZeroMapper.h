#ifndef RS_POLEZEROMAPPER_H
#define RS_POLEZEROMAPPER_H

namespace RSLib
{

  /**

  This class contains static functions for converting a set poles and zeros of some prototype 
  filter in the complex s-plane (or z-plane) to another set of poles and zeros of the corresponding
  unnormalized filter. Additionaly, it also implements a mapping between the s-plane and z-plane 
  via the bilinear transform.

  \todo make coding style consistent, comment uncommented functions, make function 
  sLowpassToLowshelf, include a gain factor k into the mapping functions (include it also into the 
  PrototypeDesigner), update function names sLowpassToLowpass etc., include different s-to-z 
  mappings (matched z-tansform (aka impulse invariant), magnitude matched, phase matched etc.), 
  include z-plane LP->LP, LP->HP, LP->BP, LP->BR mappings (Constantinides formulas)

  // use convention: *z, Nz, *p, Np, k for inputs and *zm, *pm, *km for outputs, for functionNames:
  // sLP2LP, bilinearS2Z
  // adopt a specific convention for the ordering of the poles and zeros (this should be adhered
  // to in the rsPrototypeDesigner as well


  http://en.wikipedia.org/wiki/Prototype_filter describes also a lowpass->multiband transform -
  maybe, we can implement this

  */

  class RSLib_API rsPoleZeroMapper
  {

  public:

    /** Transforms an analog lowpass prototype filter represented by its zeros, poles and gain (in 
    "z", "p", "k") into a corresponding low-shelving filter with reference gain G0 and 
    low-frequency gain G. */
    static void sLowpassToLowshelf(rsComplexDbl *z, rsComplexDbl *p, double *k, rsComplexDbl *zNew,
      rsComplexDbl *pNew, double *kNew, int N, double G0, double G);

    /** Transforms an analog lowpass prototype filter represented by its zeros, poles and gain (in
    "z", "p", "k") with nominal radian cutoff frequency of unity to an analog lowpass with nominal
    radian cutoff frequency of wc. The new zeros, poles and gain factor are returned in "zNew", 
    "pNew" and "kNew". */
    static void sLowpassToLowpass(rsComplexDbl *z, rsComplexDbl *p, double *k, rsComplexDbl *zNew, 
      rsComplexDbl *pNew, double *kNew, int N, double wc);

    /** Transforms an analog lowpass prototype filter to an analog highpass with nominal radian 
    cutoff frequency of wc. @see sLowpassToLowpass */
    static void sLowpassToHighpass(rsComplexDbl *z, rsComplexDbl *p, double *k, rsComplexDbl *zNew,
      rsComplexDbl *pNew, double *kNew, int N, double wc);

    /** Transforms an analog lowpass prototype filter to an analog bandpass with lower and upper 
    nominal radian bandedge frequencies of wl and wu. @see sLowpassToLowpass */
    static void sLowpassToBandpass(rsComplexDbl *z, rsComplexDbl *p, double *k, rsComplexDbl *zNew,
      rsComplexDbl *pNew, double *kNew, int N, double wl, double wu);

    /** Transforms an analog lowpass prototype filter to an analog bandreject with lower and upper 
    nominal radian bandedge frequencies of wl and wu. @see sLowpassToLowpass */
    static void sLowpassToBandreject(rsComplexDbl *z, rsComplexDbl *p, double *k, 
      rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double wl, double wu);

    //static void zLowpassToLowpass(rsComplexDbl *z, rsComplexDbl *p, double *k, 
    //  rsComplexDbl *zNew, rsComplexDbl *pNew, double *kNew, int N, double wc);

    /** Performs an s-plane lowpass-to-lowpass transform. This transform is done in place which 
    means, the result is stored in the same arrays as the input. These arrays must be of length 
    numPoles and numZeros respectively. The target cutoff frequency is assumed to be expressed as 
    radian frequency (Omega = 2*pi*frequency). ...to be deprecated. */
    static void prototypeToAnalogLowpass(rsComplexDbl *poles, int numPoles, rsComplexDbl *zeros, 
      int numZeros, double *gain, double targetCutoff);

    /** \todo comment this function */
    static void sPlanePrototypeToLowpass(rsComplexDbl *prototypePoles, 
      rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
      int prototypeOrder, double targetCutoff);

    /** Performs an s-plane lowpass-to-highpass transform. This transform is done in place which 
    means, the result is stored in the same arrays as the input. These arrays must BOTH be of 
    length max(numPoles, numZeros) because we get additional finite zeros at s=0 for zeros at 
    s=inf in the prototype. ...to be deprecated. */
    static void prototypeToAnalogHighpass(rsComplexDbl *poles, int numPoles, rsComplexDbl *zeros, 
      int numZeros, double *gain, double targetCutoff);

    /** \todo comment this function */
    static void sPlanePrototypeToHighpass(rsComplexDbl *prototypePoles, 
      rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
      int prototypeOrder, double targetCutoff);

    /** Same as prototypeToAnalogLowpass but additionally exchanges the roles of poles and 
    zeros. */
    static void prototypeToAnalogHighShelv(rsComplexDbl *poles, int numPoles, rsComplexDbl *zeros, 
      int numZeros, double *gain, double targetCutoff);

    /** \todo comment this function */
    static void sPlanePrototypeToHighShelv(rsComplexDbl *prototypePoles, 
      rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
      int prototypeOrder, double targetCutoff);

    /** Performs an s-plane lowpass-to-bandpass transform. This transform is done in place which 
    means, the result is stored in the same arrays as the input. These arrays must BOTH be of 
    length 2*max(numPoles, numZeros) because each finite pole and zero maps to a pair of finite 
    poles and zeros and we get additional finite zeros for zeros at s=inf in the prototype. */
    static void prototypeToAnalogBandpass(rsComplexDbl *poles, int numPoles, rsComplexDbl *zeros, 
      int numZeros, double *gain, double targetLowCutoff, double targetHighCutoff);

    /** \todo comment this function */
    static void sPlanePrototypeToBandpass(rsComplexDbl *prototypePoles, 
      rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
      int prototypeOrder, double targetLowerCutoff, double targetUpperCutoff);

    /** Performs an s-plane lowpass-to-bandpass transform. This transform is done in place which 
    means, the result is stored in the same arrays as the input. These arrays must BOTH be of 
    length 2*max(numPoles, numZeros) because each finite pole and zero maps to a pair of finite 
    poles and zeros and we get additional finite zeros for zeros at s=inf in the prototype. */
    static void prototypeToAnalogBandstop(rsComplexDbl *poles, int numPoles, rsComplexDbl *zeros, 
      int numZeros, double *gain, double targetLowCutoff, double targetHighCutoff);

    /** \todo comment this function */
    static void sPlanePrototypeToBandreject(rsComplexDbl *prototypePoles, 
      rsComplexDbl *prototypeZeros, rsComplexDbl *targetPoles, rsComplexDbl *targetZeros, 
      int prototypeOrder, double targetLowerCutoff, double targetUpperCutoff);

    /** Maps zeros, poles and gain factor of a digital prototype lowpass with normalized radian 
    cutoff frequency "wc" to the new zeros, poles and gain factor of a digital targer lowpass with 
    normalized radian frequency "wt". */
    static void zLowpassToLowpass(rsComplexDbl *z, rsComplexDbl *p, double *k, rsComplexDbl *zNew, 
      rsComplexDbl *pNew, double *kNew, int N, double wc, double wt);

    //void analogToDigitalLowpassToLowpass(rsComplexDbl *poles, int numPoles, 
    //  rsComplexDbl *zeros, int numZeros, double *gain, double targetCutoff);

    /** Performs a bilinear transform from the s-plane to the z-plane. */
    static void bilinearAnalogToDigital(rsComplexDbl *poles, int numPoles, rsComplexDbl *zeros, 
      int numZeros, double sampleRate, double *gain);

    //void bilinearDigitalToAnalog(rsComplexDbl *poles, int numPoles, rsComplexDbl *zeros, 
    //  int numZeros, double sampleRate, double *gain);
    /**< Performs a bilinear transform from the z-plane to the s-plane. */


  protected:

    // why are these not public?

    /** Transforms an array of roots (poles or zeros) of an analog lowpass prototype filter to 
    corresponding roots of an analog highpass with nominal radian cutoff frequency of wc. */
    static void sLowpassToHighpass(rsComplexDbl *r, rsComplexDbl *rNew, int N, double wc);

    /** Transforms an array of roots (poles or zeros) of an analog lowpass prototype filter to 
    corresponding roots of an analog bandpass with nominal radian center frequency of wc and 
    bandwidth bw. */
    static void sLowpassToBandpass(rsComplexDbl *r, rsComplexDbl *rNew, int N, double wc, 
      double bw);

    /** Transforms an array of roots (poles or zeros) of an analog lowpass prototype filter to 
    corresponding roots of an analog bandreject with nominal radian center frequency of wc and 
    bandwidth bw. */
    static void sLowpassToBandreject(rsComplexDbl *r, rsComplexDbl *rNew, int N, double wc, 
      double bw);

  };

}

#endif
