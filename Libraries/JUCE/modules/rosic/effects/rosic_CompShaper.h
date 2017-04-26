#ifndef rosic_CompShaper_h
#define rosic_CompShaper_h

// rosic-indcludes:
#include "../math/rosic_Interpolation.h"

namespace rosic
{

  /**

  This class implements a waveshaper with a (softknee) compressor alike transfer function. As there 
  is no dynamic  behaviour here, there are no attack/release parameters. Threshold, ratio, knee and 
  the boolean limiter flag work entirely similar as in class SoftKneeCompressor, only that here the 
  transfer function here acts in the linear amplitude domain (as opposed to the dB-domain).

  \todo: implement oversampling

  */

  class CompShaper
  {

  public:

    enum kneeShapes
    {
      CUBIC,
      QUARTIC
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    CompShaper();

    /** Destructor */
    ~CompShaper();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the input signal gain / drive in dB. */
    void setDrive(double newDrive) { driveFactor = dB2amp(newDrive); }

    /** Sets the amount of DC (direct current) to applied to the iput signal (after the drive, before 
    the waveshaper. */
    void setDC(double newDC) { dc = newDC; }

    /** Sets the threshold above which the transfer curve snaps off. */
    void setThreshold(double newThreshold) { threshold = dB2amp(newThreshold); calculateCoefficients(); }

    /** Sets the ratio - that is the reciprocal of the transfer function's slope above the the
    threshold. */
    void setRatio(double newRatio) 
    { ratio = newRatio; oneOverRatio = 1.0/ratio; calculateCoefficients(); }

    /** Switches the waveshaper into clipper mode - levels above the threshold will be mapped to
    the threshold itself - this corresponds to an infinite ratio. */
    void setToClipperMode(bool shouldClip) { clipperMode = shouldClip; calculateCoefficients(); }

    /** Sets the transition width between the two slopes as a percentage of the threshold.
    
    
    (the value is specified in decibels despite
    the fact that the transfer curve acts on linear amplitude). */
    void setKneeWidth(double newKneeWidth);

    /** Sets the shape of the transition between threshold-knee/2 and  threshold+knee/2. 
    @see kneeShapes */
    void setKneeShape(int newKneeShape);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //=============================================================================================

    /** Returns the value of the transfer curve at a given x - x is assumed to be the current
    amplitude (envelope) of some signal, expressed in dB. */
    INLINE double transferCurveAt(double x);

  protected:

    /** Computes the polynomial coefficients. */
    void calculateCoefficients();

    /** Calculates the polynomial coefficients for the quartic transition. */
    void calculateQuarticCoefficients();

    /** Calculates the polynomial coefficients for the cubic transition. */
    void calculateCubicCoefficients();

    double a0, a1, a2, a3, a4; // polynomial coefficients
    double driveFactor, dc, outputGainFactor, amount, threshold, ratio, oneOverRatio, kneeWidth;
    int    kneeShape;
    bool   clipperMode;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE double CompShaper::transferCurveAt(double x)
  {
    double xs = sign(x);
    double xa = fabs(x);

    if( xa < threshold-0.5*kneeWidth*0.01*threshold )
      return x;
    else if( xa > threshold+0.5*kneeWidth*0.01*threshold )
    {
      if( clipperMode )
        return xs * threshold;
      else
        return xs * (threshold + oneOverRatio * (xa-threshold)); 
    }
    else
      return xs * evaluateQuartic(xa, a0, a1, a2, a3, a4);
  }

  INLINE void CompShaper::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double tmpL = driveFactor * (*inOutL) + dc;
    double tmpR = driveFactor * (*inOutR) + dc;  

    tmpL = transferCurveAt(tmpL);
    tmpR = transferCurveAt(tmpR);

    *inOutL  = tmpL; 
    *inOutR  = tmpR; 
  }

} // end namespace rosic

#endif // #ifndef rosic_CompShaper_h
