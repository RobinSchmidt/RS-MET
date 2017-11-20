#ifndef RS_FOURPOLEFILTER_H
#define RS_FOURPOLEFILTER_H

namespace RSLib
{

  /**

  This class defines the user parameters for the rsFourPoleFilter class.

  \todo: remove the polyphony-handling - this aspect should be handled somewhere else

  */

  class RSLib_API rsFourPoleFilterParameters
  {

  public:

    /** This is an enumeration of the available filter modes. */
    enum modes
    {
      BYPASS = 0,

      LOWPASS_6,
      LOWPASS_12,
      LOWPASS_18,
      LOWPASS_24,

      HIGHPASS_6,
      HIGHPASS_12,
      HIGHPASS_18,
      HIGHPASS_24,

      BANDPASS_RBJ,
      BANDREJECT_RBJ,
      PEAK_OR_DIP_RBJ,
      LOW_SHELV_1ST,
      LOW_SHELV_RBJ,
      HIGH_SHELV_1ST,
      HIGH_SHELV_RBJ,
      ALLPASS_1ST,
      ALLPASS_RBJ,
      LOWPASS_HIGHPASS_RBJ,   // first stage highpass, second stage lowpass

      MORPH_LP_BP_HP,         // lowpass, bandpass, highpass
      MORPH_LP_PK_HP,         // lowpass, peak, highpass

      MORPH_LP_RES_HP,        // lowpass, resonator, highpass
      MORPH_LP_HP,            // lowpass, resonator, highpass
      MORPH_LP_N_HP,          // lowpass, notch, highpass
      MORPH_BP_BN_N,          // bandpass, bi-notch, notch

      RESON,                  // two-pole resonator (Steiglitz)
      RESONZ,                 // two-pole-one-zero resonator (Steiglitz)

      RC_LPF_HPF_SERIES,      // 1-pole-LPF and -HPF in series
      RC_2LPFS_SERIES,        // 2 1-pole LPFs in series
      RC_2HPFS_SERIES,        // 2 1-pole HPFs in series
      RC_LPF_HPF_PARALLEL,    // 1-pole-LPF and -HPF in parallel

      BUTTER_LPF,             // butterworth-LPF of 2nd order
      BUTTER_HPF,             // butterworth-HPF of 2nd order
      BUTTER_BPF,             // butterworth-BPF of 2nd order
      BUTTER_BRF,             // butterworth-BRF of 2nd order

      ONE_POLE_APF_SERIES,    // 2 1-pole APFs in series

      FOF_FILTER,           // realizes a FOF-filter specified by the
      // FOF-parameters fofA, fofAlpha, fofOmega, fofPhi
    };

    double sampleRate;          // the sample-rate
    double sampleRateRec;       // reciprocal of the sample-rate
    double resonance;
    double q;
    double gainDb;              // gain in dB
    double gainFactor;          // gain in dB as raw amplitude-factor
    double morph;               ///< \todo: do not share this among voices
    double transition;          // morph transition parameter
    double freq2Scale;          // scale factor of the second frequency
    double freq2Offset;         // additive offset of the second frequency
    double q2Scale;             // scale factor of the second q-factor
    double gain2Scale;          // scale factor of the second gain (in dB)

    int    mode;                // mode of the filter
    //int    order;             // order of the filter
    bool   twoStageSwitch;      // indicator flag for the second stage

    rsFourPoleFilterParameters()
    {
      sampleRate          =  44100.0;
      sampleRateRec       =  1.0 / sampleRate;
      resonance           =  0.0;
      q                   =  rsSqrt(0.5);
      gainDb              =  0.0;
      gainFactor          =  1.0;
      morph               = -1.0;
      transition          =  0.0;
      freq2Scale          =  2.0;
      freq2Offset         =  0.0;
      q2Scale             =  1.0;
      gain2Scale          =  1.0;
      mode                =  0;
      twoStageSwitch      = false;
    }

  };


  /**

  This class implements a cascade of two biquad filters, where the second stage is optional. The
  filter can be set up to the desired specifications and can then calculate its filter
  coefficients. The (re) calculation of the coefficients, however, does not happen automatically
  each time a parameter changes - instead it must be triggered from outside by means of a call to
  updateFilterCoefficients(). This was done so, in order to avoid multiple calculations of the same
  coefficient-set when more than one parameter changes at a given time. The two biquad stages are
  each implemented in direct form 1, because i found that this form responds best to rapidly
  varying parameters.

  \todo: mode Notch-Peak-Notch - one pole-pair near the circle surrounded by two zero pairs on the 
    circle, the othe two poles are used to normalize the gain at DC and pi

  */

  class RSLib_API rsFourPoleFilter
  {

  public:

    /** \name Construction/Destruction */

    /** Constructor. */
    rsFourPoleFilter();

    /** Destructor. */
    ~rsFourPoleFilter();


    /** \name Setup */

    /** Sets up the sample-rate. */
    void setSampleRate(double newSampleRate);

    /** Sets the mode of the filter. @see modes */
    void setMode(int newMode);

    /** Decides, whether only the first stage or both should be used. */
    void useTwoStages(bool shouldUseTwoStages);

    /** Sets the order of the filter (1...4). This will behave differently depending on the chosen
    filter mode. For example, if the mode is LOWPASS and the order is 3, a cookbook resonant
    lowpass and a first order lowpass will be used in series connection. */
    //void setOrder(double newOrder);

    // parameters that are supposed to be updated at sample-rate (inlined
    // access-functions):

    /** Sets the cuttoff/center frequency of the filter. */
    RS_INLINE void setFrequency(double newFrequency);

    /** Sets Q of the filter. */
    RS_INLINE void setQ(double newQ);

    /** Sets resonance of the filter. */
    RS_INLINE void setResonance(double newResonance);

    /** Sets the gain value (in dB) for shelving and peaking modes. */
    RS_INLINE void setGain(double newGain);

    /** Sets the morph parameter for the morphable modes. */
    RS_INLINE void setMorph(double newMorph);

    /** Sets the morph-transition parameter for the morphable modes. */
    RS_INLINE void setTransition(double newTransition);

    /** Sets a scale factor for the characteristic frequency of the second filter stage with
    respect to the first. */
    RS_INLINE void setSecondFrequencyScaleFactor(double newFactor);

    /** Sets an offset (in Hz) for the characteristic frequency of the second filter stage with
    respect to the first. */
    RS_INLINE void setSecondFrequencyOffset(double newOffset);

    /** Sets a scale factor for the quality-factor 'Q' of the second filter stage with
    respect to the first. */
    RS_INLINE void setSecondQScaleFactor(double newFactor);

    /** Sets a scale factor for the gain (in dB) of the second filter stage with respect to the
    first. */
    RS_INLINE void setSecondGainScaleFactor(double newFactor);

    /** Calculates filter coefficients from filter parameters. Has to be called each time Freq, Q
    or Gain changes -> it is NOT called automatically by setFreq, setQ or setGain to avoid multiple
    calculations of the coeffs when more than one of these parameters is changed at the same time
    (for example by an envelope).
    costs (AMD Athlon 64 3200+): 240-375 CPU-cycles (depending on the mode) */
    RS_INLINE void updateFilterCoefficients();


    /** \name Inquiry */

    /** Returns the mode of the filter. @see modes */
    int getMode() { return parameters->mode; }

    /** Informs, whether only the first stage or both are used. */
    bool usesTwoStages() { return parameters->twoStageSwitch; }

    /** Returns the cuttoff/center frequency of the filter. */
    double getFrequency() { return freq; }

    /** Returns the Q of the filter. */
    double getQ() { return parameters->q; }

    /** Returns the gain value (in dB) for shelving and peaking modes. */
    double getGain() { return parameters->gainDb; }

    /** Returns the mroph parameter for morphable modes. */
    double getMorph() { return parameters->morph; }

    /** Returns the scale factor for the characteristic frequency of the second filter stage with
    respect to the first. */
    double getSecondFrequencyScaleFactor() { return parameters->freq2Scale; }

    /** Returns the offset (in Hz) for the characteristic frequency of the second filter stage with
    respect to the first. */
    double getSecondFrequencyOffset() { return parameters->freq2Offset; }

    /** Returns the scale factor for the quality factor 'Q' of the second filter stage with
    respect to the first. */
    double getSecondQScaleFactor() { return parameters->q2Scale; }

    /** Returns the scale factor for the gain (in dB) of the second filter stage with
    respect to the first. */
    double getSecondGainScaleFactor() { return parameters->gain2Scale; }

    /** Returns the value of the magnitude response of the filter at the given frequency. */
    double getMagnitudeAt(double frequency);

    /** Calculates the magnitudes of the frequency-response at the frequencies given in the array
    "frequencies" (in Hz) and stores them in the array "magnitudes". Both arrays are assumed to be
    "numBins" long. "inDecibels" indicates, if the frequency response should be returned in
    decibels. If "accumulate" is true, the magnitude response of this biquad-cascade will be
    multiplied with (or added to, when "inDecibels" is true) to the magnitudes which are already
    there in the "magnitudes"-array. This is useful for calculating the magnitude response
    of several biquad-cascades in series. */
    void getMagnitudeResponse(double *frequencies, double *magnitudes, int numBins,
      bool inDecibels = false, bool accumulate = false);


    //static RS_INLINE void calculateCookbookPeakFilterCoeffs(double* b0, double* b1, double* b2,
    //  double* a1, double* a2, double centerFrequency, double q, double gain);


    /** \name Audio Processing */

    /** Calculates a single filtered output-sample via Direct Form 1.
    costs (AMD Athlon 64 3200+): ~22.5 cycles per stage. */
    RS_INLINE double getSample(double in);

    /** Calculates one output stereo sample-frame at a time. */
    RS_INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);


    /** \name Polyphony handling */

    /** Adds a  slave instance to the vector of slaves - this will also cause the 'isMaster' flag
    of the slave to be set to false, redirect the slaves parameters-pointer to the one of this
    instance and delete the old (now unused) parameters pointer of the slave. */
    void addSlave(rsFourPoleFilter* newSlave);


    /** \name Misc */

    /** Sets the buffers for the previous input and output samples of both biquad stages to
    zero. */
    void reset();

    /** Sets the buffers for the previous input and output samples of biquad stage 2 to zero. */
    void resetBuffersStage2();

    /** Copies the coefficients form the first biquad stage into the ones for the seconds stage
    also in order to use two identical coefficient sets for both biquad stages. */
    RS_INLINE void copyFirstStageCoeffsIntoSecondStageCoeffs();


  protected:

    /** Maps the morph parameter in order to make the behaviour more musical. */
    RS_INLINE double mapMorphParameter(double m);

    // direct form coefficients for stage 1 and 2:
    doubleA a1_s1, a2_s1, b0_s1, b1_s1, b2_s1,
            a1_s2, a2_s2, b0_s2, b1_s2, b2_s2;

    // buffer variables for the first biquad stage:
    doubleA xL_s1_d1,   // left input, 1st stage, 1 sample delay
            xR_s1_d1,   // right input, 1st stage, 1 sample delay
            xL_s1_d2,   // left input, 1st stage, 2 samples delay
            xR_s1_d2,   // right input, 1st stage, 2 samples delay
            yL_s1_d1,   // left output, 1st stage, 1 sample delay
            yR_s1_d1,   // right output, 1st stage, 1 sample delay
            yL_s1_d2,   // left output, 1st stage, 2 samples delay
            yR_s1_d2;   // right output, 1st stage, 2 samples delay

    // buffer variables for the second biquad stage:
    doubleA xL_s2_d1,   // left input, 2nd stage, 1 sample delay
            xR_s2_d1,   // right input, 2nd stage, 1 sample delay
            xL_s2_d2,   // left input, 2nd stage, 2 samples delay
            xR_s2_d2,   // right input, 2nd stage, 2 samples delay
            yL_s2_d1,   // left output, 2nd stage, 1 sample delay
            yR_s2_d1,   // right output, 2nd stage, 1 sample delay
            yL_s2_d2,   // left output, 2nd stage, 2 samples delay
            yR_s2_d2;   // right output, 2nd stage, 2 samples delay

    doubleA preGain;  // overall gain applied to the input

    // filter parameters:
    doubleA freq;
    bool secondStageIsActive; // indicator flag for the second stage - redundant now (?)
    //doubleA resonance;
    //doubleA q;
    //doubleA gain, A; // gain in dB and as raw amplitude-factor
    //doubleA morph; ///< todo: do not share morph among voices


    rsFourPoleFilterParameters* parameters;
      // A pointer to the parameters which are potentially shared by among instances. */

    std::vector<rsFourPoleFilter*> slaves;
      // A vector of pointers to other instances of this class which shall be kept in sync to this
      // instance with regard to their parameters. */

    bool isMaster;
      // A flag which indicates whether or not this instance is a master which controls other
      // instances of this class - this will also determine whether or not this objects will delete 
      // the pointer to the parameter set on destruction. By default, instances will be constructed 
      // as master, later they can be re-configured as slaves by adding them as slaves via addSlave 
      // to another instance.

  };


  //-----------------------------------------------------------------------------------------------
  // inlined functions:

  RS_INLINE void rsFourPoleFilter::setFrequency(double newFrequency)
  {
    // restrict cutoff frequency to the range between 20 and 20000 Hz:
    if( newFrequency <= 20.0 )
      freq = 20.0;
    else if( newFrequency >= 20000.0 )
      freq = 20000.0;
    else
      freq = newFrequency;
    updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::setQ(double newQ)
  {
    if( newQ >= 0.1 )
    {
      parameters->q = newQ;
      updateFilterCoefficients();
      for(unsigned int s=0; s<slaves.size(); s++)
        slaves[s]->updateFilterCoefficients();
    }
  }

  RS_INLINE void rsFourPoleFilter::setResonance(double newResonance)
  {
    if( newResonance >= 0.0 )
      parameters->resonance = newResonance;
  }

  RS_INLINE void rsFourPoleFilter::setGain(double newGain)
  {
    parameters->gainDb     = newGain;
    parameters->gainFactor = pow(10, (0.025*parameters->gainDb) );
    updateFilterCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::setMorph(double newMorph)
  {
    parameters->morph = newMorph;
    updateFilterCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::setTransition(double newTransition)
  {
    parameters->transition = newTransition;
    updateFilterCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::setSecondFrequencyScaleFactor(double newFactor)
  {
    parameters->freq2Scale = newFactor;
    updateFilterCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::setSecondFrequencyOffset(double newOffset)
  {
    parameters->freq2Offset = newOffset;
    updateFilterCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::setSecondQScaleFactor(double newFactor)
  {
    parameters->q2Scale = newFactor;
    updateFilterCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::setSecondGainScaleFactor(double newFactor)
  {
    parameters->gain2Scale = newFactor;
    updateFilterCoefficients();
    for(unsigned int s=0; s<slaves.size(); s++)
      slaves[s]->updateFilterCoefficients();
  }

  RS_INLINE void rsFourPoleFilter::updateFilterCoefficients()
  {
    preGain = 1.0;

    double q          = parameters->q;
    double gainFactor = parameters->gainFactor;

    if(  parameters->mode == rsFourPoleFilterParameters::PEAK_OR_DIP_RBJ
      || parameters->mode == rsFourPoleFilterParameters::LOW_SHELV_1ST
      //|| parameters->mode == rsFourPoleFilterParameters::LOW_SHELV_RBJ
      || parameters->mode == rsFourPoleFilterParameters::HIGH_SHELV_1ST    )
      //|| parameters->mode == rsFourPoleFilterParameters::HIGH_SHELV_RBJ     )
    {
      if( secondStageIsActive )
        gainFactor = rsSqrt(gainFactor);
    }
    else if(  parameters->mode == rsFourPoleFilterParameters::LOW_SHELV_RBJ
           || parameters->mode == rsFourPoleFilterParameters::HIGH_SHELV_RBJ  )
    {
      if( secondStageIsActive )
      {
        gainFactor = rsSqrt(gainFactor);
        q = rsSqrt(q);
      }
    }
    else
    {
      if( secondStageIsActive )
        q = rsSqrt(q);
    }

    switch( parameters->mode )
    {
    case rsFourPoleFilterParameters::LOWPASS_6:
      {
        rsBiquadDesigner::calculateFirstOrderLowpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq);
        //BiquadDesigner::makeBypassBiquad(b0_s2, b1_s2, b2_s2, a1_s2, a2_s2);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::LOWPASS_12:
      {
        rsBiquadDesigner::calculateCookbookLowpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        //BiquadDesigner::makeBypassBiquad(b0_s2, b1_s2, b2_s2, a1_s2, a2_s2);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::LOWPASS_18:
      {
        rsBiquadDesigner::calculateCookbookLowpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        rsBiquadDesigner::calculateFirstOrderLowpassCoeffs(b0_s2, b1_s2, b2_s2, a1_s2, a2_s2,
          parameters->sampleRateRec, freq);
      }
      break;
    case rsFourPoleFilterParameters::LOWPASS_24:
      {
        rsBiquadDesigner::calculateCookbookLowpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;


    case rsFourPoleFilterParameters::HIGHPASS_6:
      {
        rsBiquadDesigner::calculateFirstOrderHighpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq);
        //BiquadDesigner::makeBypassBiquad(b0_s2, b1_s2, b2_s2, a1_s2, a2_s2);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::HIGHPASS_12:
      {
        rsBiquadDesigner::calculateCookbookHighpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        //BiquadDesigner::makeBypassBiquad(b0_s2, b1_s2, b2_s2, a1_s2, a2_s2);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::HIGHPASS_18:
      {
        rsBiquadDesigner::calculateCookbookHighpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        rsBiquadDesigner::calculateFirstOrderHighpassCoeffs(b0_s2, b1_s2, b2_s2, a1_s2, a2_s2,
          parameters->sampleRateRec, freq);
      }
      break;
    case rsFourPoleFilterParameters::HIGHPASS_24:
      {
        rsBiquadDesigner::calculateCookbookHighpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;

      /*
    case HIGHPASS_6:
      {
        BiquadDesigner::calculateFirstOrderHighpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case HIGHPASS_12:
      {
        BiquadDesigner::calculateCookbookHighpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
      */

    case rsFourPoleFilterParameters::BANDPASS_RBJ:
      {
        rsBiquadDesigner::calculateCookbookBandpassConstSkirtCoeffsViaQ(b0_s1, b1_s1, b2_s1, a1_s1,
          a2_s1, parameters->sampleRateRec, freq, q);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::BANDREJECT_RBJ:
      {
        rsBiquadDesigner::calculateCookbookBandrejectCoeffsViaQ(b0_s1, b1_s1, b2_s1, a1_s1,
          a2_s1, parameters->sampleRateRec, freq, q);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::PEAK_OR_DIP_RBJ:
      {
        rsBiquadDesigner::calculateCookbookPeakFilterCoeffsViaQ(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q, gainFactor);
        copyFirstStageCoeffsIntoSecondStageCoeffs();

      }
      break;
    case rsFourPoleFilterParameters::LOW_SHELV_1ST:
      {
        rsBiquadDesigner::calculateFirstOrderLowShelvCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, gainFactor);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::LOW_SHELV_RBJ:
      {
        rsBiquadDesigner::calculateCookbookLowShelvCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q, gainFactor);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::HIGH_SHELV_1ST:
      {
        rsBiquadDesigner::calculateFirstOrderHighShelvCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, gainFactor);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::HIGH_SHELV_RBJ:
      {
        rsBiquadDesigner::calculateCookbookHighShelvCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q, gainFactor);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::ALLPASS_1ST:
      {
        rsBiquadDesigner::calculateFirstOrderAllpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::ALLPASS_RBJ:
      {
        rsBiquadDesigner::calculateCookbookAllpassCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1,
          parameters->sampleRateRec, freq, q);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    case rsFourPoleFilterParameters::MORPH_LP_PK_HP:
      {
        rsBiquadDesigner::calculateLowPeakHighMorphCoeffs(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1, preGain,
          parameters->sampleRate, parameters->sampleRateRec, freq, q, parameters->morph,
          secondStageIsActive);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
      break;
    default:
      {
        rsBiquadDesigner::makeBypassBiquad(b0_s1, b1_s1, b2_s1, a1_s1, a2_s1);
        copyFirstStageCoeffsIntoSecondStageCoeffs();
      }
    } // end of switch(mode)
  }

  RS_INLINE double rsFourPoleFilter::getSample(double in)
  {
    doubleA out;

    // calculate output of the first stage:
    in *= preGain;
    out = b0_s1*in + (b1_s1*xL_s1_d1 + b2_s1*xL_s1_d2) + (a1_s1*yL_s1_d1 + a2_s1*yL_s1_d2);

    // update buffers of the first stage:
    xL_s1_d2 = xL_s1_d1;
    xL_s1_d1 = in;
    yL_s1_d2 = yL_s1_d1;
    yL_s1_d1 = out;

    // optionally apply the second stage:
    if( secondStageIsActive )
    {
      out = b0_s2*out + (b1_s2*xL_s2_d1 + b2_s2*xL_s2_d2) + (a1_s2*yL_s2_d1 + a2_s2*yL_s2_d2);

      // update buffers of the first stage:
      xL_s2_d2 = xL_s2_d1;
      xL_s2_d1 = yL_s1_d1;  // this contains still the sample which was input to this stage
      yL_s2_d2 = yL_s2_d1;
      yL_s2_d1 = out;
    }

    return out;
  }

  RS_INLINE void rsFourPoleFilter::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    // we need some temporary variables for the filtering because we cannot take for granted that
    // *inL is distinct from *outL (the same for right channel)
    double tmpL, tmpR, tmpL1, tmpR1;

    // calculate output of the first stage:
    tmpL1 = preGain * (*inOutL);
    tmpR1 = preGain * (*inOutR);
    tmpL  = b0_s1*tmpL1 + (b1_s1*xL_s1_d1 + b2_s1*xL_s1_d2) + (a1_s1*yL_s1_d1 + a2_s1*yL_s1_d2);
    tmpR  = b0_s1*tmpR1 + (b1_s1*xR_s1_d1 + b2_s1*xR_s1_d2) + (a1_s1*yR_s1_d1 + a2_s1*yR_s1_d2);

    // update buffers of the first stage:
    xL_s1_d2 = xL_s1_d1;
    xL_s1_d1 = tmpL1;
    yL_s1_d2 = yL_s1_d1;
    yL_s1_d1 = tmpL;

    xR_s1_d2 = xR_s1_d1;
    xR_s1_d1 = tmpR1;
    yR_s1_d2 = yR_s1_d1;
    yR_s1_d1 = tmpR;

    // optionally apply the second stage:
    if( secondStageIsActive )
    {
      tmpL = b0_s2*tmpL + (b1_s2*xL_s2_d1 + b2_s2*xL_s2_d2) + (a1_s2*yL_s2_d1 + a2_s2*yL_s2_d2);
      tmpR = b0_s2*tmpR + (b1_s2*xR_s2_d1 + b2_s2*xR_s2_d2) + (a1_s2*yR_s2_d1 + a2_s2*yR_s2_d2);

      // update buffers of the second stage:
      xL_s2_d2 = xL_s2_d1;
      xL_s2_d1 = yL_s1_d1;  // this contains still the sample which was input to this stage
      yL_s2_d2 = yL_s2_d1;
      yL_s2_d1 = tmpL;
      xR_s2_d2 = xR_s2_d1;
      xR_s2_d1 = yR_s1_d1;
      yR_s2_d2 = yR_s2_d1;
      yR_s2_d1 = tmpR;
    }

    *inOutL = tmpL;
    *inOutR = tmpR;
  }

  RS_INLINE double rsFourPoleFilter::mapMorphParameter(double m)
  {
    //double x1=x;
    //double c = cosSquaredApprox(x);
    //return 1.0-c;

    return fabs(m);
    //return m*m;

    /*
    m = m*m;
    double a  = 3.0;
    double y1 = 0.5*(tanh(a*(m-0.5))+1.0);
    double o  = 0.5*(tanh(a*(0.0-0.5))+1.0);
    double y2 = y1-o;
    double y  = (1.0/(1.0-2.0*o))*y2;
    return y;
    */
  }

  RS_INLINE void rsFourPoleFilter::copyFirstStageCoeffsIntoSecondStageCoeffs()
  {
    a1_s2 = a1_s1;
    a2_s2 = a2_s1;
    b0_s2 = b0_s1;
    b1_s2 = b1_s1;
    b2_s2 = b2_s1;
  }

}

#endif
