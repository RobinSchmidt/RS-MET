#ifndef rosic_PitchShifter_h
#define rosic_PitchShifter_h

//// rosic-indcludes:
//#include "../filters/rosic_EllipticSubBandFilterDirectForm.h"
////#include "../basics/rosic_Interpolator.h"

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

    /** Constructor. */
    PitchShifter();

    /** Destructor. */
    ~PitchShifter();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the coarse detuning in semitones. */
    void setDetuneCoarse(double newDetuneCoarse);

    /** Sets the coarse detuning in cents. */
    void setDetuneFine(double newDetuneFine);

    /** Sets the length of the grains - this is the maximum distance which the output pointers fall
    behind the input pointer (in milliseconds) - it can not exceed the length of the delayline. */
    void setGrainLength(double newGrainLength);

    /** Sets the amount by which the pitch shifted output is fed back to the input
    (in percent). */
    void setFeedback(double newFeedback);

    /** Sets the ratio between dry and wet signal (in percent wet). */
    void setDryWet(double newDryWet);

    /** Sets the playback direction to reverse. */
    void setReversePlayback(bool shouldPlayReverse);

    /** Inverts the polarity of the wet signal. */
    void setNegativePolarity(bool shouldBeNegative);

    /** Switches an anti-aliasing filter on or off. */
    void setAntiAliasing(bool shouldAntiAlias);

    /** Switches between stereo and mono mode. */
    //void setMonoMode(bool newMonoMode) { mono = newMonoMode; }

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the coarse detuning in semitones. */
    double getDetuneCoarse() const { return detuneCoarse; }

    /** Returns the coarse detuning in cents. */
    double getDetuneFine() const { return detuneFine; }

    /** Returns the size of the grains - this is the maximum distance which the output pointers fall
    behind the input pointer. */
    double getGrainLength() const { return grainLength; }

    /** Returns the amount by which the pitch shifted output is fed back to the input (in
    percent). */
    double getFeedback() const { return 100.0 * feedbackFactor; }

    /** Returns the ratio between dry and wet signal (in percent wet). */
    double getDryWet() const { return 100.0 * wet; }

    /** Informs whether the playback direction is reverse (true) or normal (false). */
    bool isPlaybackReverse() const { return reversePlayback; }

    /** Informs whether the polarity of the wet signal is inverted (true) or normal (false). */
    bool isPolarityNegative() const { return (wetPolarity == -1.0); }

    /** Informs whether the anti-alias filter is on or off. */
    bool getAntiAliasing() const { return antiAliasingIsOn; }

    /** Informs whether the object is in mono processing mode. */
    //bool isInMonoMode() const { return mono; }

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates a stereo-ouput frame. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the content of the delaylines to all zeros. */
    void reset();

    //=============================================================================================

  protected:

    /** Updates the increment for the two tapOut pointers according to the desired
    pitch-shift. */
    void updateIncrement();

    /** Initializes the distances of the output taps with respect to the input tap. */
    void initTapDistances();

    /** Updates the member variables which have to do with the distances of between the write
    and read pointers. */
    void updateDistanceVariables();

    int    delayLineLength;
    double *delayLineL, *delayLineR;
    int    tapIn;

    double distance1, distance2, maxDistance, maxDistanceRec, maxDistanceHalf;
      // Distances which the two output taps are behind the input tap (in samples) and the maximum
      // value which they assume ...and the reciprocal of the maximum value. */

    double distanceIncrement;
      // Distance increment (or decrement) for the ouput pointers. */

    //double incrementSign;
      // A sign factor to be applied to the distance incrment (for reverse playback). */

    double detuneCoarse, detuneFine;

    double sampleRate;
    double grainLength;
    double feedbackFactor, dry, wet;
    double wetPolarity;

    bool   reversePlayback;
    bool   antiAliasingIsOn;
    //bool   mono;


    //EllipticSubBandFilterDirectForm antiAliasFilterL, antiAliasFilterR;
    rsEllipticSubBandFilter antiAliasFilterL, antiAliasFilterR;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void PitchShifter::getSampleFrameStereo(double *inOutL,  double *inOutR)
  {
    double fracPart, tmpL1, tmpR1, tmpL, tmpR, weight;
    int    intPart, tapOut;

    // get the ouput form the first output-tap:

    //intPart  = (int) floor(distance1);
    intPart  = floorInt(distance1);
    fracPart = distance1 - (double) intPart;
    tapOut   = tapIn - intPart;
    if( tapOut < 1 )
      tapOut += delayLineLength-1;
    if( tapOut > delayLineLength-1 )
      tapOut = 1;

    // get the sample from the delayline with linear interpolation:
    //tmpL1 = Interpolator::getSampleLinear( 1.0-fracPart, &(delayLineL[tapOut-1]) );
    tmpL1 = delayLineL[tapOut] + fracPart * ( delayLineL[tapOut-1] - delayLineL[tapOut] );
    tmpR1 = delayLineR[tapOut] + fracPart * ( delayLineR[tapOut-1] - delayLineR[tapOut] );

    // apply the cos^2 weighting function and store the value:
    weight = RAPT::rsCosSquaredApprox( maxDistanceRec*(2.0*distance1-maxDistance) );
    //weight = cos(0.5*PI* maxDistanceRec*(2.0*distance1-maxDistance) ); weight *= weight;
    tmpL   = tmpL1 * weight;
    tmpR   = tmpR1 * weight;

    // add the ouput from the second output-tap:
    if( distance1 >= maxDistanceHalf )
      distance2 = distance1 - maxDistanceHalf;
    else
      distance2 = distance1 + maxDistanceHalf;

    //intPart  = (int) floor(distance2);
    intPart  = floorInt(distance2);
    fracPart = distance2 - (double) intPart;
    tapOut   = tapIn - intPart;
    if( tapOut < 1 )
      tapOut += delayLineLength-1;
    if( tapOut > delayLineLength-1 )
      tapOut = 1;

    // get the sample from the delayline with linear interpolation:
    //tmpL1 = Interpolator::getSampleLinear( 1.0-fracPart, &(delayLineL[tapOut-1]) );
    tmpL1 = delayLineL[tapOut] + fracPart * ( delayLineL[tapOut-1] - delayLineL[tapOut] );
    tmpR1 = delayLineR[tapOut] + fracPart * ( delayLineR[tapOut-1] - delayLineR[tapOut] );

    // apply the cos^2 weighting function and store the value:
    weight  = RAPT::rsCosSquaredApprox( maxDistanceRec*(2.0*distance2-maxDistance) );
    //weight = cos(0.5*PI* maxDistanceRec*(2.0*distance2-maxDistance) ); weight *= weight;
    tmpL   += tmpL1 * weight;
    tmpR   += tmpR1 * weight;

    // write the incoming samples plus the feedback signal into the delaylines:
    if( antiAliasingIsOn )
    {
      delayLineL[tapIn] = antiAliasFilterL.getSample(*inOutL + feedbackFactor * tmpL);
      delayLineR[tapIn] = antiAliasFilterR.getSample(*inOutR + feedbackFactor * tmpR);
    }
    else
    {
      delayLineL[tapIn] = *inOutL + feedbackFactor * tmpL;
      delayLineR[tapIn] = *inOutR + feedbackFactor * tmpR;
    }

    // repeat the sample in the last (additional) cell at index zero  for the linear interpolator:
    if( tapIn == delayLineLength-1 )
    {
      delayLineL[0] = delayLineL[tapIn];
      delayLineR[0] = delayLineR[tapIn];
    }

    // increment or decrement distance and wrap around if necesarry:
    if( reversePlayback == true )
      distance1 += (2.0 - distanceIncrement);
    else
      distance1 += distanceIncrement;

    if( distance1 <= 0.0 )
      distance1 += maxDistance;
    if( distance1 > maxDistance )
      distance1 -= maxDistance;

    // increment tapIn-pointer and wrap around if necesarry:
    tapIn++;
    if( tapIn >= delayLineLength )
      tapIn = 1;

    // mix dry wet and store the result:
    *inOutL = dry * (*inOutL) + wet * wetPolarity * tmpL;
    *inOutR = dry * (*inOutR) + wet * wetPolarity * tmpR;
  }

} // end namespace rosic

#endif // rosic_PitchShifter_h
