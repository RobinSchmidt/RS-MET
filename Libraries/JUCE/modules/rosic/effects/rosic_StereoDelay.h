#ifndef rosic_StereoDelay_h
#define rosic_StereoDelay_h

// rosic-indcludes:
#include "../math/rosic_ElementaryFunctionsReal.h"
#include "../delaylines/rosic_AllpassDiffusor.h"
#include "../filters/rosic_LowpassHighpass.h"

namespace rosic
{

  /**

  This StereoDelay class consists of two delay-lines with subsequent diffusor alpasses and
  lowpass/highpass filters. The delaylines have the option for feedback and crossfeedback. The 
  delay-times are tuned in beats with respect to some tempo which is assumed to be given in 
  BPM.

  */

  class StereoDelay
  {

  public:

    enum delaylines
    {
      LEFT = 0,
      RIGHT
    };
    /**< Enumeration of the delaylines. */

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    StereoDelay();   
    ///< Constructor.

    ~StereoDelay();  
    ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    void setSampleRate(double newSampleRate);
    /**< Sets the sample-rate. */

    void setBeatsPerMinute(double newBpm);
    /**< Sets the tempo in beats per minute. */

    void setDelayInBeats(double newDelay, int delayline);
    /**< Sets the delay time in beats for one of the two delay lines. */

    void setDelayScale(double newFactor, int delayline);
    /**< Sets a multiplicator for the delay time for one of the two delay lines. */

    void setInjection(double newPercentage, int sourceChannel, int delayline);
    /**< Adjusts, how much of the signal from the 'sourceChannel' (LEFT or RIGHT) goes into the 
    delay line  specified by 'delayline' (LEFT or RIGHT) - in percent, can be negative. */

    void setDiffusorTimeInMilliseconds(double newDiffusorTime, int delayline);
    /**< Sets the delay-time (in seconds milliseconds) for the diffusor-allpass for one of the 
    delay lines. */

    void setDiffusorAmount(double newDiffusorAmount, int delayline);
    /**< Sets the amount of diffusion (in percent) for one of the delay lines. */

    void setLowpassCutoff(double newCutoff, int delayline);
    /**< Sets the cutoff frequency of the lowpass-filter (in Hz) for one of the delay lines. */

    void setHighpassCutoff(double newCutoff, int delayline);
    /**< Sets the cutoff frequency of the highpass-filter (in Hz) for one of the delay lines. */

    void setCutoffScale(double newCutoffScale);
    /**< Sets a scaling factor for the cutoff-frequencies of all the lowpass/highpass filters. */

    void setFeedback(double newFeedback, int source, int target);
    /**< Sets the amount of feedback or crossfeedback from one delay-line to itself or to the other 
    one (in percent, negative values allowed). */

    void setOutputMix(double newPercentage, int delayline, int targetChannel);
    /**< Adjusts, how much of the signal from the delay line  specified by 'delayline' 
    (LEFT or RIGHT) goes out to the 'targetChannel' (LEFT or RIGHT) - in percent, can be 
    negative. */

    void setWetDelayInBeats(double newDelay, int channel);
    /**< Sets the 'pre'-delay time for the wet signal in beats for one of the two delay lines. */

    void setDryWet(double newDryWet);
    /**< Sets the mix between the dry and the wet signal (in percent wet). */

    //---------------------------------------------------------------------------------------------
    // inquiry:

    double getDelayInBeats(int delayline);
    /**< Returns the delay time in whole notes for one of the two delay lines. */

    double getDelayScale(int delayline);
    /**< Returns a multiplicator for the delay time for one of the two delay lines. */

    double getInjection(int sourceChannel, int delayline);
    /**< Returns, how much of the signal from the 'sourceChannel' (LEFT or RIGHT) goes into the 
    delay line  specified by 'delayline' (LEFT or RIGHT) - in percent, can be negative. */

    double getDiffusorTimeInMilliseconds(int delayline);
    /**< Returns the delay-time (in milliseconds) for the diffusor-allpass for one of the delay 
    lines. */

    double getDiffusorAmount(int delayline);
    /**< Returns the amount of diffusion (in percent) for one of the delay lines. */

    double getLowpassCutoff(int delayline);
    /**< Returns the cutoff frequency of the lowpass-filter (in Hz) for one of the delay lines. */

    double getHighpassCutoff(int delayline);
    /**< Returns the cutoff frequency of the highpass-filter (in Hz) for one of the delay lines. */

    double getCutoffScale();
    /**< Returns a scaling factor for the cutoff-frequencies of all the lowpass/highpass 
    filters. */

    double getFeedback(int source, int target);
    /**< Returns the amount of feedback or crossfeedback from one delay-line to itself or to the 
    other one (in percent, negative values allowed). */

    double getOutputMix(int delayline, int targetChannel);
    /**< Adjusts, how much of the signal from the delay line  specified by 'delayline' 
    (LEFT or RIGHT) goes out to the 'targetChannel' (LEFT or RIGHT) - in percent, can be 
    negative. */

    double getWetDelayInBeats(int channel);
    /**< Returns the pre-delay time for the output signal in beats for one of the two delay 
    lines. */

    double getDryWet();
    /**< Returns the mix between the dry and the wet signal (in percent wet). */

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double* inOutL,  double* inOutR);
    /**< Calculates a stereo-ouput frame. */

    //---------------------------------------------------------------------------------------------
    // others:

    void reset();
    /**< Resets the internal buffers. */

    //=============================================================================================

  protected:

    void updateDelayTimes();
    /**< Updates the delay-times of the delay-lines according to the sampleRate, the tempo, the 
    delay-length in whole notes and the delaytime fineTuning-factor. */

    void updateFilterFrequencies();
    /**< Updates the filters frequencies acording to their nominal settings and the cutoffScale 
    factor. */

    double inL2L, inL2R, inR2L, inR2R;
    /**< Input scale factors (injection), first L or R specifies input delayline, second L/R 
    specifies delayline. */

    double fbL2L, fbL2R, fbR2L, fbR2R;
    /**< Feedback and crossfeedback factors: left->left, left->right, right->left and 
    right->right. */

    double outL2L, outL2R, outR2L, outR2R;
    /**< Output scale factors (injection), first L or R specifies delayline, second L/R 
    specifies output delayline. */

    double dryGain, wetGain;
    /**< Amplitude factors for the dry and wet signals. */

    double wetDelayInBeatsL, wetDelayInBeatsR;
    /**< 'Pre'-delay times for left and right wet signal in beats. */

    IntegerDelayLine delayLineL, delayLineR;
    /**< The embedded delaylines for left and right delayline. */

    AllpassDiffusor diffusorL, diffusorR;    
    /**< The embedded diffusors for left and right delayline. */

    LowpassHighpass filterL, filterR;   
    /**< The embedded bandpasses for left and right delayline. */

    IntegerDelayLine wetDelayLineL, wetDelayLineR;
    /**< The embedded delaylines for left and right wet signal. */

    double sampleRate;
    double bpm;
    double delayInBeatsL, delayInBeatsR;
    double delayScaleL, delayScaleR;
    double lowpassCutoffL, highpassCutoffL, lowpassCutoffR, highpassCutoffR;
    double cutoffScale;
    double dryWet; // dry/wet ratio in % wet
  };

  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  INLINE void StereoDelay::getSampleFrameStereo(double* inOutL,  double* inOutR)
  {
    // establish the input signals to the delaylines by weighted sums of the inputs:
    double tmpL = inL2L * (*inOutL)  +  inR2L * (*inOutR);
    double tmpR = inL2R * (*inOutL)  +  inR2R * (*inOutR);

    // apply the current input samples to the delaylines but do not increment their tap-pointers
    // yet, because we must still add the feedback signals to the input taps (which are not yet
    // available at this time):
    tmpL = delayLineL.getSampleSuppressTapIncrements(tmpL);
    tmpR = delayLineR.getSampleSuppressTapIncrements(tmpR);

    // apply the diffusor allpass and bandpass filtering:
    tmpL = filterL.getSample(diffusorL.getSample(tmpL));
    tmpR = filterR.getSample(diffusorR.getSample(tmpR));

    // apply feedback and crossfeedback around the delaylines:
    delayLineL.addToInput(fbL2L*tmpL);
    delayLineL.addToInput(fbR2L*tmpR);    
    delayLineR.addToInput(fbR2R*tmpR);
    delayLineR.addToInput(fbL2R*tmpL);

    // increment the tap-pointers in the delaylines:
    delayLineL.incrementTapPointers();
    delayLineR.incrementTapPointers();

    // establish the output signals by weighted sums of the delayline-outputs:
    tmpL = outL2L*tmpL + outR2L*tmpR;
    tmpR = outL2R*tmpL + outR2R*tmpR;

    // apply pre-delay to the output-signals:
    tmpL = wetDelayLineL.getSample(tmpL);
    tmpR = wetDelayLineR.getSample(tmpR);

    // apply dry/wet control and store the calculated output samples in the output slots:
    *inOutL = dryGain * (*inOutL) + wetGain * tmpL;
    *inOutR = dryGain * (*inOutR) + wetGain * tmpR;
  }

} // end namespace rosic

#endif // rosic_StereoDelay_h
