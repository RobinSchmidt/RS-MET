#ifndef rosic_StereoPan_h
#define rosic_StereoPan_h

//// rosic-indcludes:
//#include "../basics/rosic_ChannelMatrix2x2.h"
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements a gain based stereo-pan with different pan-laws (including laws that
  involve cross-mixing). The pan-laws are normalized such that in center position the total gain
  of each of the two channels is unchanged (when there's no cross-mixing).

  \todo: make also a delay-based pan class and one that combines both approaches in a consistent
  way

  */

  class StereoPan : public ChannelMatrix2x2
  {

  public:

    enum panLaws
    {
      LINEAR = 0,
      SINCOS,
      SQUARE_ROOT,
      LINEAR_CLIPPED,
      LINEAR_SQUARE_NORMALIZED,
      LINEAR_CLIPPED_SQUARE_NORMALIZED,
      LINEAR_CROSSMIX_SQUARE_NORMALIZED,
      SINCOS_CROSSMIX_SQUARE_NORMALIZED,
      SQRT_CROSSMIX_SQUARE_NORMALIZED,
    };

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    StereoPan();

    /** Destructor. */
    ~StereoPan();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the panorama position as value in the range -1...+1 where -1 is hard-left and +1 is
    hard-right. */
    void setPanoramaPosition(double newPosition)
    { p = RAPT::rsClip(newPosition, -1.0, 1.0); calculateGainFactors(); }

    /** Sets a global gain in decibels */
    void setGain(double newGain)
    { g = RAPT::rsDbToAmp(newGain); calculateGainFactors(); }

    /** Selects one of the pan-laws. @see panLaws */
    void setPanLaw(int newLaw) { panLaw = newLaw; calculateGainFactors(); }

    /** Selects whether the channel gain-factors should be normalized to unity in the center - if
    false, they will be normalized at the extreme sides. */
    void setNormalizeAtCenter(bool shouldNormalizeAtCenter)
    { normalizeAtCenter = shouldNormalizeAtCenter; calculateGainFactors(); }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the panorama position as value in the range -1...+1 where -1 is hard-left and +1 is
    hard-right. */
    double getPanoramaPosition() const { return p; }

    /** Returns the index of the currently selected pan-law. @see panLaws */
    int getPanLaw() const { return panLaw; }

    /** Returns true when the selected pan-law applies some channel cross-mixing, false
    otherwise. */
    bool doesPanLawApplyCrossMix() const { return panLaw >= 6; }

    //---------------------------------------------------------------------------------------------
    // static member functions for the different pan-laws:

    static void normalizeSum(double &x1, double &x2, double &x3, double &x4, double targetSum)
    {
      double s = targetSum / (x1+x2+x3+x4);  // scaler to be applied
      x1 *= s; x2 *= s; x3 *= s; x4 *= s;
    }

    static void normalizeSumOfSquares(double &x1, double &x2, double &x3, double &x4,
      double targetSum)
    {
      double s = sqrt( targetSum / (x1*x1+x2*x2+x3*x3+x4*x4) );  // scaler to be applied
      x1 *= s; x2 *= s; x3 *= s; x4 *= s;
    }

    /** Computes the gain factors for a linear channel weighting without cross-mixing. */
    static void panLawLinear(double p, double &gLL, double &gRL, double &gLR, double &gRR,
      bool centerNormalized = false)
    {
      double pn = 0.5*p+0.5;  // normalized to 0...1
      gLL = 1.0-pn;
      gRL = 0.0;
      gLR = 0.0;
      gRR = pn;
      if( centerNormalized == true )
      {
        gLL *= 2.0;
        gRR *= 2.0;
      }
    }

    /** Computes the gain factors for a sin/cos based channel weighting (giving a constant sum of
    the squares of the weights) without cross-mixing. */
    static void panLawSinCos(double p, double &gLL, double &gRL, double &gLR, double &gRR,
      bool centerNormalized = false)
    {
      equalPowerGainFactors(p, &gLL, &gRR, -1.0, 1.0);
      gRL  = 0.0;
      gLR  = 0.0;
      if( centerNormalized == true )
      {
        gLL *= SQRT2;
        gRR *= SQRT2;
      }
    }

    /** Computes the gain factors for a square-root based channel weighting (giving a constant sum
    of the squares of the weights) without cross-mixing. */
    static void panLawSquareRoot(double p, double &gLL, double &gRL, double &gLR,
      double &gRR, bool centerNormalized = false)
    {
      double pn = 0.5*p+0.5;  // normalized to 0...1
      gLL = sqrt(1.0-pn);
      gRL = 0.0;
      gLR = 0.0;
      gRR = sqrt(pn);
      if( centerNormalized == true )
      {
        gLL *= SQRT2;
        gRR *= SQRT2;
      }
    }

    /** Computes the gain factors for a channel weighting without cross-mixing that increases
    linearly until it reaches 1 at the center position. */
    static void panLawLinearClipped(double p, double &gLL, double &gRL, double &gLR, double &gRR,
      bool /*centerNormalized = false*/)
    {
      double pn = 0.5*p+0.5;  // normalized to 0...1
      gLL = rmin(2.0*(1.0-pn), 1.0);
      gRL = 0.0;
      gLR = 0.0;
      gRR = rmin(2.0*pn, 1.0);
    }

    /** Computes the gain factors for a linear channel weighting without cross-mixing and
    re-normalization of the sum-of-squares. */
    static void panLawLinearSquareNormalized(double p, double &gLL, double &gRL, double &gLR,
      double &gRR, bool centerNormalized = false)
    {
      double pn = 0.5*p+0.5;  // normalized to 0...1
      gLL = 1.0-pn;
      gRL = 0.0;
      gLR = 0.0;
      gRR = pn;
      if( centerNormalized == true )
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 2.0);
      else
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 1.0);
    }

    /** Computes the gain factors for a channel weighting without cross-mixing that increases
    linearly until it reaches 1 at the center position and then re-normalizes of the
    sum-of-squares. */
    static void panLawLinearClippedSquareNormalized(double p, double &gLL, double &gRL,
      double &gLR, double &gRR, bool centerNormalized = false)
    {
      double pn = 0.5*p+0.5;  // normalized to 0...1
      gLL = rmin(2.0*(1.0-pn), 1.0);
      gRL = 0.0;
      gLR = 0.0;
      gRR = rmin(2.0*pn, 1.0);
      if( centerNormalized == true )
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 2.0);
      else
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 1.0);
    }

    /** Calculates the four gain factors for a pan-law that first applies a linear channel
    weighting, then cross-mixes in the opposite channel whenever the pan-postion is off center and
    finally renormalizes the sum-of-the-squares of the gain factors. */
    static void panLawLinearCrossMixSquareNormalized(double p, double &gLL, double &gRL,
      double &gLR, double &gRR, bool centerNormalized = false)
    {
      // calculate preliminary gain factors:
      if( p >= 0.0 )
      {
        gRL = 0.0;
        gLR = p;
        gRR = (0.5*p)+0.5;
        gLL = 1.0-gRR;
      }
      else
      {
        gRL = -p;
        gLR = 0.0;
        gLL = (-0.5*p)+0.5;
        gRR = 1.0-gLL;
      }
      if( centerNormalized == true )
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 2.0);
      else
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 1.0);

    }

    /** Calculates the four gain factors for a pan-law that first applies a sin/cos channel
    weighting, then cross-mixes in the opposite channel whenever the pan-postion is off center and
    finally renormalizes the sum-of-the-squares of the gain factors. */
    static void panLawSinCosCrossMixSquareNormalized(double p, double &gLL, double &gRL,
      double &gLR, double &gRR, bool centerNormalized = false)
    {
      // calculate preliminary gain factors:
      equalPowerGainFactors(p, &gLL, &gRR, -1.0, 1.0);
      if( p >= 0.0 )
      {
        gRL = 0.0;
        gLR = sin(0.5*PI*p);
      }
      else
      {
        gRL = sin(-0.5*PI*p);
        gLR = 0.0;
      }
      if( centerNormalized == true )
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 2.0);
      else
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 1.0);
    }

    /** Calculates the four gain factors for a pan-law that first applies a sqrt-based channel
    weighting, then cross-mixes in the opposite channel whenever the pan-postion is off center and
    finally renormalizes the sum-of-the-squares of the gain factors. */
    static void panLawSqrtCrossMixSquareNormalized(double p, double &gLL, double &gRL,
      double &gLR, double &gRR, bool centerNormalized = false)
    {
      // calculate preliminary gain factors:
      double pn = 0.5*p+0.5;  // normalized to 0...1
      gLL = sqrt(1.0-pn);
      gRR = sqrt(pn);
      if( p >= 0.0 )
      {
        gRL = 0.0;
        gLR = sqrt(p);
      }
      else
      {
        gRL = sqrt(-p);
        gLR = 0.0;
      }
      if( centerNormalized == true )
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 2.0);
      else
        normalizeSumOfSquares(gLL, gRL, gLR, gRR, 1.0);
    }

    //=============================================================================================

  protected:

    /** Calculates the 4 gain-factors gLL, gRL, gLR, and gRR according to the desired panorama
    position and the selected pan-law. */
    void calculateGainFactors();

    double g, p;
    int    panLaw;
    bool   normalizeAtCenter;

  };

} // end namespace rosic

#endif // rosic_StereoPan_h
