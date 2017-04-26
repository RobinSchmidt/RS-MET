#ifndef rosic_LadderFilter_h
#define rosic_LadderFilter_h

// standard-library includes:
#include <stdlib.h>          // for the NULL macro

// includes for the STL:
#include <vector>

// rosic-indcludes:
#include "../math/rosic_ComplexFunctions.h"
#include "rosic_OnePoleFilterStereo.h"

namespace rosic
{

  /**

  This class defines the user parameters for the LadderFilter class.

  \todo: feedbackByCutoff via a laguerre type curve with parameter 'a', 'a' is determined via
  a = -(cutoff-b)/c with (for example b = 20 and c = 20000 -> these represent frequencies, at 20 Hz
  'a' will be zero and at 20 kHz 'a' will be close to 1, maybe dabble a bit with higher values for
  'c' to make it less close to 1 but don't use lower values)
  then: y = (x-a*x) / (1-2*a*x+a) where x is the feedback (normalized to 0...1) and y is the
  frequency-scaled feedback
  -or: scale the resonance with the logarithm of frequency, 500 Hz serving as point where the 
  resonance is scaled down most
  -morph between various slopes for lowpass and highpass
  -squelch parameter s (idea): modifies input signal according to: x' = (1-s)*x + s*x^3


  */

  class LadderFilterParameters
  {

  public:

    enum modes
    {
      LOWPASS = 0,
      HIGHPASS,
      BANDPASS_6_12,
      BANDPASS_6_18,
      BANDPASS_12_12,
      BANDPASS_18_6,
      BANDPASS_12_6,
      HIGH_BAND_LOW_MORPH,

      NUM_MODES
    };

    double resonanceRaw;     // resonance parameter
    double resonanceSkewed;  // the mapped/skewed resonance for better adjustment
    double allpassFreq;      // characteristic frequency of the allpass before the Moog filter
    double driveFactor;      // a drive factor for distortion
    //double makeUpByBoost;    // the amount of low-frequency compensation by a shelver
    double makeUp;           // the amount of remaining gain compensation in percent
    double dcOffset;         // a dc-offset for odd harmonics
    double sampleRate;       // the sample-rate
    double morph;            // for morphing between highpass through bandpass to lowpass
    int    outputStage;      // filter stage from which the output is taken
    int    mode;             // switches between different modes

    LadderFilterParameters()
    {
      resonanceRaw    = 0.0;
      resonanceSkewed = 0.0;
      allpassFreq     = 20000.0;
      //makeUpByBoost   = 0.0;
      makeUp          = 0.0;
      driveFactor     = 1.0;
      //dcOffset        = 0.05;
      dcOffset        = 0.0;
      sampleRate      = 44100.0;
      morph           = 0.0;       // 4th order lowpass
      outputStage     = 4;
      mode            = LOWPASS;
    }

  };

  /**

  This class is a Moog style filter similar to the MoogFilter class. This class here, however
  mnakes several tweaks to make it more CPU-friendly and customize the characteristic sound.

  \todo factor out the aspect of polyphony, provide different modes instead of just the selection of the slope

  */

  class LadderFilter
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    LadderFilter();

    /** Destructor. */
    ~LadderFilter();

    //---------------------------------------------------------------------------------------------
    // parameter settings:

    /** Sets the sample-rate for this filter. */
    void setSampleRate(double newSampleRate);

    /** Sets the cutoff frequency for this filter - the actual coefficient calculation may be 
    supressed by passing 'false' as second parameter, in this case, it should be triggered
    manually later by calling calculateCoefficients. */
    INLINE void setCutoff(double newCutoff, bool updateCoefficients = true);

    void setCutoff(double newCutoff) { setCutoff(newCutoff, true); } 
      // overloaded to use as target for callbacks - todo make an extra function setCutoffWithoutUpdate

    /** Sets the amount of inverted feedback from the output to the input for the resonant peak. */
    INLINE void setResonance(double newResonance, bool updateCoefficients = true);

    /** Sets the amount of resonance in percent. */
    void setResonanceInPercent(double newResonanceInPercent) { setResonance(0.01*newResonanceInPercent, true); }

    /** Sets the characteristic frequnecy of an allpass filter which is applied before the actual
    Moog filter model in order to pre-shape the input waveform - this affects mainly the sound of
    the distortion due to the nonlinearity. */
    void setAllpassFreq(double newAllpassFreq);

    /** Sets the input drive in decibels. */
    void setDrive(double newDrive);

    /** Sets the input dc-offset. */
    void setDcOffset(double newDrive);

    /** Selects the filter stage from which the output is taken (1...4). */
    void setOutputStage(int newOutputStage);

    /** Selects the mode of the filter. @see LadderFilterParameters::modes */
    void setMode(int newMode);

    /** Sets the morphing parameter to smoothly morph between highpass through bnadpass to 
    lowpass (mode HIGH_BAND_LOW_MORPH must be chosen). */
    void setMorph(double newMorph);

    /** The makeUp parameter controls the amount by which the low frequency loss due to high
    resonance will be compesated by an overall gain factor. A value of 0.0 means no
    compensation, 100.0 is full compensation, higher values lead to overcompensation and
    negative values will even emphasize the low frequency loss (the shelver will apply a
    negative gain then). */
    void setMakeUp(double newMakeUp, bool updateCoefficients = true);

    void setMakeUp(double newMakeUp) { setCutoff(newMakeUp, true); } 

    //---------------------------------------------------------------------------------------------
    // inquiry:

    /** Returns the cutoff frequency of this filter. */
    double getCutoff();

    /** Returns the resonance parameter of this filter. */
    double getResonance();

    /** Returns the characteristic frequency of the allpass. */
    double getAllpassFreq();

    /** Returns the MakeUp gain of this filter. */
    double getMakeUp();

    /** Returns the drive parameter of this filter. */
    double getDrive();

    /** Returns the filter stage from which the output is taken (1...4). */
    int getOutputStage();

    /** Returns the value of the filters transfer function (including resonance) at the complex
    value z (the filter is viewed as a linear filter for that matter, which it is not). The last
    parameter selects the filter stage from which the output will be taken. */
    Complex getTransferFunctionAt(Complex z, bool withFeedback = true, bool withMakeUpBoost = true,
      bool withMakeUpGain = true, int stage = 4);

    /** Returns the value of the magnitude response of the filter at the given frequency. */
    double getMagnitudeAt(double frequency, bool withFeedback = true, bool withMakeUpBoost = true,
      bool withMakeUpGain = true, int stage = 4);

    /** Calculates the magnitudes of the frequency-response at the frequencies given in the array
    "frequencies" (in Hz) and stores them in the array "magnitudes". Both arrays are assumed to be
    "numBins" long. "inDecibels" indicates, if the frequency response should be returned in
    decibels. If "accumulate" is true, the magnitude response of this biquad-cascade will be
    multiplied with (or added to, when "inDecibels" is true) to the magnitudes which are already
    there in the "magnitudes"-array. This is useful for calculating the magnitude response
    of several biquad-cascades in series. */
    void getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins,
      bool inDecibels = false, bool accumulate = false);

    //---------------------------------------------------------------------------------------------
    // audio processing:


    /** under construction for test */
    double getSampleTest(double in);

    /** Calculates one output sample at a time. */
    INLINE double getSample(double in);

    /** Calculates one stereo output sample frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // master/slave configuration:

    /** Adds a  slave instance to the vector of slaves - this will also cause the 'isMaster' flag
    of the slave to be set to false, redirect the slaves parameters-pointer to the one of this
    instance and delete the old (now unused) parameters pointer of the slave. */
    void addSlave(LadderFilter* newSlave);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Causes the filter to re-calculate the coeffiecients - is not called automatically by
    setCutoff(). */
    INLINE void calculateCoefficients();

    //INLINE void calcuateCoefficientsApprox();




    /** Resets the internal state variables. */
    void reset();

    //=============================================================================================

  protected:

    double b0, a1, a1Old; // coefficients for the first order sections

    double y1L, y1R, y2L, y2R, y3L, y3R, y4L, y4R, yOldL, yOldR;
      // Output signals of the 4 filter stages and the overall output signal from the previous
      // iteration for left and right channel

    double c0, c1, c2, c3, c4;
      // weights for the different output taps

    double k, kOld;
      // Feedback factor in the loop

    double makeupGain; //, makeupGainRec, makeupGainSq;
      // The factor resulting from the makeup (compensation for low frequency loss)

    double cutoff;
      // Cutoff frequency of this filter instance - this is not to be shared among voices

    OnePoleFilterStereo allpass;
      // This is a first order allpass filter which is apllied before the actual ladder filter

    LadderFilterParameters* parameters;
      // A pointer to the parameters which are potentially shared by among instances

    std::vector<LadderFilter*> slaves;
      // A vector of pointers to other instances of this class which shall be kept in sync to this
      // instance with regard to their parameters. */

    bool isMaster;
      // A flag which indicates whether or not this instance is a master which controls other
      // instances of this class - this will also determine whether or not this objects will delete 
      // the pointer to the parameter set on destruction. By default, instances will be constructed 
      // as master, later they can be re-configured as slaves by adding them as slaves via addSlave 
      // to another instance

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed 
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void LadderFilter::setCutoff(double newCutoff, bool updateCoefficients)
  {
    if( newCutoff != cutoff )
    {
      if( newCutoff < 20.0 )
        cutoff = 20.0;
      else if( newCutoff > 20000.0 )
        cutoff = 20000.0;
      else
        cutoff = newCutoff;

      if( updateCoefficients == true )
        calculateCoefficients();
    }
  }

  INLINE void LadderFilter::setResonance(double newResonance, bool updateCoefficients)
  {
    if( newResonance != parameters->resonanceRaw )
    {
      parameters->resonanceRaw    = newResonance;
      parameters->resonanceSkewed = (1.0-exp(-3.0*newResonance)) / (1.0-exp(-3.0));
      if( updateCoefficients == true )
      {
        calculateCoefficients();
        for(unsigned int s=0; s<slaves.size(); s++)
          slaves[s]->calculateCoefficients();
      }
    }
  }

  INLINE void LadderFilter::calculateCoefficients()
  {
    // calculate intermediate variables:
    double wc = 2 * PI * cutoff / parameters->sampleRate;
    double s, c;
    sinCos(wc, &s, &c);             // c = cos(wc); s = sin(wc);
    double t  = tan(0.25*(wc-PI));
    double r  = parameters->resonanceSkewed;

    // calculate filter a1-coefficient tuned such the resonance frequency is just right:
    double a1_fullRes = t / (s-c*t);

    // calculate filter a1-coefficient as if there were no resonance:
    double x        = exp(-wc);
    double a1_noRes = -x;

    // use a weighted sum between the resonance-tuned and no-resonance coefficient:
    a1 = r*a1_fullRes + (1.0-r)*a1_noRes;

    // calculate the b0-coefficient from the condition that each stage should be a leaky
    //integrator:
    b0 = 1.0+a1;

    // calculate feedback factor by dividing the resonance parameter by the magnitude at the
    // resonant frequency:
    double gsq = b0*b0 / (1.0 + a1*a1 + 2.0*a1*c);
    k          = r / (gsq*gsq);

    // evaluate the magnitude response at DC:
    double b0_4   = b0*b0*b0*b0;
    //double dr     = ( (((a1+4.0)*a1+6.0)*a1+4.0)*a1 + k*b0_4 + 1.0);
    double dcGain = b0_4 / ( (((a1+4.0)*a1+6.0)*a1+4.0)*a1 + k*b0_4 + 1.0);

    // derive the required makeup gain from the dc-gain with optimized calculations for the two
    // extreme settings of the makeUp parameter:
    if( parameters->makeUp == 100.0 )
      makeupGain = 1.0 / dcGain;
    else if( parameters->makeUp == 0.0 )
      makeupGain = 1.0;
    else
      makeupGain = pow(dcGain, -0.01*parameters->makeUp);

    /*
    double dcGainDb = amp2dB(dcGain);


    // calculate the required gain, its reciprocal and its square (the latter two are needed
    // because we apply the inverse gain befor the nonlinearity and its square after it in order to
    // decouple the makeup-parameter from the amount of distortion):
    double boostDb = -(0.01*parameters->makeUp) * dcGainDb;
    makeupGain     = dB2amp(boostDb);
    */
  }

  INLINE double LadderFilter::getSample(double in)
  {
    // apply drive, feedback and DC-offset:
    double y0L = parameters->driveFactor*in - k*yOldL + parameters->dcOffset;

    // cascade of 4 1st order sections:
    y1L = b0*y0L - a1*y1L;
    y2L = b0*y1L - a1*y2L;
    y3L = b0*y2L - a1*y3L;
    y4L = b0*y3L - a1*y4L;

    // apply the (first) nonlinearity:
    yOldL = y4L / (1.0 + 0.25*y4L*y4L);

    // boost the lows, apply makeup gain, apply the saturation and return the sample:
    //return lowBooster.getSample(yOldL);
    //return makeupGainSq * tanhApprox(makeupGainRec*yOldL + 0.05) - 0.05;
    //return makeupGain * tanhApprox(yOldL + 0.05) - 0.05;
    return 0.0; // to do: copy and edit code form the getSampleFarmeStereo
  }

  INLINE void LadderFilter::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // apply the allpass-filter:
    double tmpL, tmpR;
    allpass.getSampleFrameStereo(inOutL, inOutR, &tmpL, &tmpR);


    // \todo try, if time-varying behavior may be improved by scaling yOld by some function of 
    // k/kOld where k is the feedback now and kOld is the feedback one sample ago. the idea is to 
    // counteract the change of stored energy in y4 when the feedback factor changes
    // maybe: y4 *= kOld/k or y4 *= k/kOld - or maybe use (k/kOld)^2 - experiment


    // apply drive, feedback and DC-offset:
    double y0L = parameters->driveFactor*tmpL - k*yOldL + parameters->dcOffset;
    double y0R = parameters->driveFactor*tmpR - k*yOldR + parameters->dcOffset;

    // cascade of 4 1st order sections:
    y1L = b0*y0L - a1*y1L;
    y1R = b0*y0R - a1*y1R;
    y2L = b0*y1L - a1*y2L;
    y2R = b0*y1R - a1*y2R;
    y3L = b0*y2L - a1*y3L;
    y3R = b0*y2R - a1*y3R;
    y4L = b0*y3L - a1*y4L;
    y4R = b0*y3R - a1*y4R;

    // apply the (first) nonlinearity:
    yOldL = y4L / (1.0 + 0.25*y4L*y4L);
    yOldR = y4R / (1.0 + 0.25*y4R*y4R);
    //yOldL = 1.25*y4L / (1.0 + 0.25*y4L*y4L);
    //yOldR = 1.25*y4R / (1.0 + 0.25*y4R*y4R);

    // select the signal of the output stage:
    //double tmpL, tmpR;

    switch( parameters->mode )
    {
    case LadderFilterParameters::LOWPASS:
      {
        switch( parameters->outputStage )
        {
        case 0: { tmpL = y0L;  tmpR = y0R; } break;
        case 1: { tmpL = y1L;  tmpR = y1R; } break;
        case 2: { tmpL = y2L;  tmpR = y2R; } break;
        case 3: { tmpL = y3L;  tmpR = y3R; } break;
        case 4: { tmpL = y4L;  tmpR = y4R; } break;
        }
      }
      break; // end case LOWPASS
    case LadderFilterParameters::HIGHPASS:
      {
        switch( parameters->outputStage )
        {
        case 0: { tmpL = y0L;                       tmpR = y0R;                       } break;
        case 1: { tmpL = y0L-  y1L;                 tmpR = y0R-  y1R;                 } break;
        case 2: { tmpL = y0L-2*y1L+  y2L;           tmpR = y0R-2*y1R+  y2R;           } break;
        case 3: { tmpL = y0L-3*y1L+3*y2L-  y3L;     tmpR = y0R-3*y1R+3*y2R-  y3R;     } break;
        case 4: { tmpL = y0L-4*y1L+6*y2L-4*y3L+y4L; tmpR = y0R-4*y1R+6*y2R-4*y3R+y4R; } break;
        }
      }
      break; // end case HIGHPASS

    case LadderFilterParameters::BANDPASS_6_12:
      {
        tmpL = y2L-y3L; 
        tmpR = y2R-y3R; 
      }
      break;
    case LadderFilterParameters::BANDPASS_6_18:
      {
        tmpL = y3L-y4L; 
        tmpR = y3R-y4R; 
      }
      break;
    case LadderFilterParameters::BANDPASS_12_12:
      {
        tmpL = y2L-2*y3L+y4L; 
        tmpR = y2R-2*y3R+y4R; 
      }
      break;
    case LadderFilterParameters::BANDPASS_18_6:
      {
        tmpL = y1L-3*y2L+3*y3L-y4L;
        tmpR = y1R-3*y2R+3*y3R-y4R;
      }
      break;
    case LadderFilterParameters::BANDPASS_12_6:
      {
        tmpL = y1L-2*y2L+y3L; 
        tmpR = y1R-2*y2R+y3R; 
      }
      break;

    case LadderFilterParameters::HIGH_BAND_LOW_MORPH:
      {
        tmpL = c0*y0L + c1*y1L + c2*y2L + c3*y3L + c4*y4L; 
        tmpR = c0*y0R + c1*y1R + c2*y2R + c3*y3R + c4*y4R; 
      }
      break;

    }

    /*
    //apply makeup gain, apply the saturation and return the samples:
    double dc       = parameters->dcOffset;
    double scaledDc = makeupGainRec * dc;
    *outL           = tanhApprox(makeupGainRec*tmpL + scaledDc) - scaledDc;
    *outR           = tanhApprox(makeupGainRec*tmpR - scaledDc) + scaledDc;
    *outL          *= makeupGainSq;
    *outR          *= makeupGainSq;
    */

    double dc = parameters->dcOffset;
    *inOutL   = tanhApprox(tmpL + dc) - dc;
    *inOutR   = tanhApprox(tmpR - dc) + dc;
    *inOutL  *= makeupGain;
    *inOutR  *= makeupGain;
  }

}

#endif // rosic_LadderFilter_h
