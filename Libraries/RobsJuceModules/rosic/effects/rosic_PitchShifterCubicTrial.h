#ifndef rosic_PitchShifter_h
#define rosic_PitchShifter_h

//// rosic-indcludes:
////#include "../math/rosic_RealFunctions.h"
//#include "../filters/rosic_EllipticSubbandFilter.h"
//#include "../basics/rosic_Interpolator.h"


namespace rosic
{

  /**

  This is a simple grain-based pitch shifter. It uses a delayline with two output tap pointers 
  which are spaced apart by half of the delayline length. The individual tap-outputs are multiplied 
  by a cos^2 shaped window-function (or grain-envelope) which has its maximum when the pointer 
  passes through the point of half the delayline length and its minima at the two 
  (forward/backward) wraparound points. These two output tap signals are then summed. For upward 
  pitch-shifting, the input can optionally be lowpass-filtered before entering the delayline in 
  order to avoid aliasing.

  */

  class PitchShifter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    PitchShifter();   
    ///< Constructor.

    ~PitchShifter();  
    ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate. */

    void setDetuneCoarse(double newDetuneCoarse);
    /**< Sets the coarse detuning in semitones. */

    void setDetuneFine(double newDetuneFine);
    /**< Sets the coarse detuning in cents. */

    void setGrainLength(double newGrainLength);
    /**< Sets the length of the grains - this is the maximum distance which the output pointers fall
    behind the input pointer (in milliseconds) - it can not exceed the length of the delayline. */

    void setFeedback(double newFeedback);
    /**< Sets the amount by which the pitch shifted output is fed back to the input 
    (in percent). */

    void setDryWet(double newDryWet);
    /**< Sets the ratio between dry and wet signal (in percent wet). */

    void setAntiAliasing(bool shouldAntiAlias);
    /**< Switches an anti-aliasing filter on or off. */ 

    //---------------------------------------------------------------------------------------------
    // inquiry:

    double getDetuneCoarse();
    /**< Returns the coarse detuning in semitones. */

    double getDetuneFine();
    /**< Returns the coarse detuning in cents. */

    double getGrainLength();
    /**< Returns the size of the grains - this is the maximum distance which the output pointers fall
    behind the input pointer. */

    double getFeedback();
    /**< Returns the amount by which the pitch shifted output is fed back to the input (in 
    percent). */

    double getDryWet();
    /**< Returns the ratio between dry and wet signal (in percent wet). */

    bool getAntiAliasing();
    /**< Informs whether the anti-alias filter is on or off. */ 

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double* inL, double* inR, double* outL, double* outR);
    /**< Calculates a stereo-ouput frame. */

    //---------------------------------------------------------------------------------------------
    // others:

    void reset();
    /**< Resets the content of the delaylines to all zeros. */

    //=============================================================================================

  protected:

    void updateIncrement();
    /**< Updates the increment for the two tapOut pointers according to the desired 
    pitch-shift. */

    void initTapDistances();
    /**< Initializes the distances of the output taps with respect to the input tap. */

    void updateDistanceVariables();
    /**< Updates the member variables which have to do with the distances of between the write
    and read pointers. */


    int    delayLineLength;
    double *delayLineL, *delayLineR;
    int    tapIn;

    double distance1, distance2, maxDistance, maxDistanceRec, maxDistanceHalf;
    /**< Distances which the two output taps are behind the input tap (in samples) and the maximum 
    value which they assume ...and the reciprocal of the maximum value. */

    double distanceIncrement;
    /**< Distance increment (or decrement) for the ouput pointers. */

    double detuneCoarse, detuneFine;

    double sampleRate;
    double grainLength;
    double feedbackFactor, dry, wet;

    bool antiAliasingIsOn;

    EllipticSubBandFilter antiAliasFilterL, antiAliasFilterR;
  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void PitchShifter::getSampleFrameStereo(double* inL,  double* inR,
                                                 double* outL, double* outR)
  { 
    double fracPart, tmpL1, tmpR1, tmpL, tmpR, weight;
    int    intPart, tapOut;

    // get the ouput form the first output-tap:

    intPart  = (int) floor(distance1);
    fracPart = distance1 - intPart;
    tapOut   = tapIn - intPart;
    if( tapOut < 2 )
      tapOut += delayLineLength-2;
    if( tapOut > delayLineLength-2 )
      tapOut = 2;

    // get the sample from the delayline with linear interpolation:
    //tmpL1 = Interpolator::getSampleLinear( 1.0-fracPart, &(delayLineL[tapOut-1]) );
    // tmpL1 = Interpolator::getSampleHermite4p3o( 1.0-fracPart, &(delayLineL[tapOut-1]) );
    tmpL1 = delayLineL[tapOut] + fracPart * ( delayLineL[tapOut-1] - delayLineL[tapOut] );
    tmpR1 = delayLineR[tapOut] + fracPart * ( delayLineR[tapOut-1] - delayLineR[tapOut] );




    // apply the cos^2 weighting function and store the value:
    //weight = cosSquaredApprox( maxDistanceRec*(2.0*distance1-maxDistance) );
    weight = cos(0.5*PI* maxDistanceRec*(2.0*distance1-maxDistance) ); weight *= weight;
    tmpL   = tmpL1 * weight;
    tmpR   = tmpR1 * weight;

    // add the ouput form the second output-tap:
    if( distance1 >= maxDistanceHalf )
      distance2 = distance1 - maxDistanceHalf;
    else 
      distance2 = distance1 + maxDistanceHalf;

    intPart  = (int) floor(distance2);
    fracPart = distance2 - intPart;
    tapOut   = tapIn - intPart;
    if( tapOut < 2 )
      tapOut += delayLineLength-2;
    if( tapOut > delayLineLength-2 )
      tapOut = 2;

    // get the sample from the delayline with linear interpolation:
    //tmpL1 = Interpolator::getSampleLinear( 1.0-fracPart, &(delayLineL[tapOut-1]) );
    //tmpL1 = Interpolator::getSampleHermite4p3o( 1.0-fracPart, &(delayLineL[tapOut-1]) );
    tmpL1 = delayLineL[tapOut] + fracPart * ( delayLineL[tapOut-1] - delayLineL[tapOut] );
    tmpR1 = delayLineR[tapOut] + fracPart * ( delayLineR[tapOut-1] - delayLineR[tapOut] );

    // apply the cos^2 weighting function and store the value:
    //weight  = cosSquaredApprox( maxDistanceRec*(2.0*distance2-maxDistance) );
    weight = cos(0.5*PI* maxDistanceRec*(2.0*distance2-maxDistance) ); weight *= weight;
    tmpL   += tmpL1 * weight;
    tmpR   += tmpR1 * weight;

    // write the incoming samples plus the feedback signal into the delaylines:
    if( antiAliasingIsOn )
    {
      delayLineL[tapIn] = antiAliasFilterL.getSampleDirect1(*inL + feedbackFactor * tmpL);
      delayLineR[tapIn] = antiAliasFilterR.getSampleDirect1(*inR + feedbackFactor * tmpR);
    }
    else
    {
      delayLineL[tapIn] = *inL + feedbackFactor * tmpL;
      delayLineR[tapIn] = *inR + feedbackFactor * tmpR;
    }

    // repeat the sample in the last (additional) cell at index zero  for the linear interpolator:
    if( tapIn == delayLineLength-1 )
    {
      delayLineL[1] = delayLineL[tapIn];
      delayLineR[1] = delayLineR[tapIn];
    }
    if( tapIn == delayLineLength-2 )
    {
      delayLineL[0] = delayLineL[tapIn];
      delayLineR[0] = delayLineR[tapIn];
    }


    // decrement distance and wrap around if necesarry:
    distance1 += distanceIncrement;
    if( distance1 <= 0.0 )
      distance1 += maxDistance;
    if( distance1 > maxDistance )
      distance1 -= maxDistance;

    // increment tapIn-pointer and wrap around if necesarry:
    tapIn++;
    if( tapIn > delayLineLength-2 )
      tapIn = 2;

    // mix dry wet and store the result:
    *outL = dry * (*inL) + wet * tmpL;
    *outR = dry * (*inR) + wet * tmpR;
  }

} // end namespace rosic

#endif // rosic_PitchShifter_h
