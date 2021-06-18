#ifndef RAPT_BIQUADCASCADE_H_INCLUDED
#define RAPT_BIQUADCASCADE_H_INCLUDED

/** This class implements a cascade of biquad-filter stages where each of the biquads implements 
the difference equation:

y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]

Each stage has its own set of coefficients which has to be set from outside this class - this class 
does not do the filter-design. The coefficients can be calculated by one of the "Designer" classes 
such as for example the BiquadDesigner class. */

template<class TSig, class TCoef>  // types for signal and coefficients
class rsBiquadCascade // rename to BiquadChain
{
  typedef const TSig&  CRSig;   // const reference to a signal value
  typedef const TCoef& CRCoef;  // const reference to a coefficient value

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. You may pass the maximum number of stages that this cascade will be able to 
  realize here (default is 12). */
  rsBiquadCascade(int newMaxNumStages = 12);

  /** Destructor. */
  ~rsBiquadCascade();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets a new sample-rate. */
  //void setSampleRate(double newSampleRate);

  /** Sets the number of biquad stages. */
  void setNumStages(int newNumStages);

  /** Sets the order of the filter. The number of stages will we either half this value (if the 
  order is even) or (order+1)/2 if the order is odd. */
  void setOrder(int newOrder);

  /** Sets up the global gain factor (by multiplying the feedforward coefficients of the first 
  stage by the factor). */
  void setGlobalGainFactor(CRCoef newGainFactor);

  /** Copies the settings (numStages and the coefficients) from another instance of this class. */
  void copySettingsFrom(rsBiquadCascade *other);

  /** Turns this biquad-cascade into an allpass filter that has the same poles as the original 
  filter. The zeros are moevd to positions that are reflections of the poles in the unit circle. */
  void turnIntoAllpass();

  /** Multiplies all feedforward ('b') coefficients by some factor. */
  //void multiplyFeedforwardCoeffsBy(double factor);

  /** Allows the user to set the filter coefficients for the individual biquad-stages. The 
  difference-equation of each of the biquad stages is: 
  \f[ y[n] = b_0 x[n] + b_1 x[n-1] + b_2 x[n-2] - a_1 y[n-1] - a_2 y[n-2] \f] */
  inline void setCoeffs(TCoef* newB0, TCoef* newB1, TCoef* newB2, TCoef* newA1, TCoef* newA2); 
  //, double newGain = 1.0);

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the current number of stages. */
  int getNumStages() const { return numStages; }

  /** Returns the memory-address of the b0 array. */
  TCoef* getAddressB0() const { return b0; }

  /** Returns the memory-address of the b1 array. */
  TCoef* getAddressB1() const { return b1; }

  /** Returns the memory-address of the b2 array. */
  TCoef* getAddressB2() const { return b2; }

  /** Returns the memory-address of the a1 array. */
  TCoef* getAddressA1() const { return a1; }

  /** Returns the memory-address of the a2 array. */
  TCoef* getAddressA2() const { return a2; }

  /** Returns the global gain factor. */
  //TCoef getGlobalGainFactor() const { return gain; }

  /** Writes the complex frequency-response of a biquad-cascade at the normalized radian 
  frequencies given in 'w' into the array 'H'. */
  void getFrequencyResponse(TCoef* w, std::complex<TCoef>* H, int numBins, 
    int accumulationMode = rsFilterAnalyzer<TCoef>::NO_ACCUMULATION) const;

  /** Writes the magnitdue-response of a biquad-cascade at the normalized radian frequencies given 
  in 'w' into the array 'magnitudes'. Both arrays are assumed to be "numBins" long. "inDecibels" 
  indicates, if the frequency response should be returned in decibels. If "accumulate" is true, the 
  magnitude response of this biquad-cascade will be multiplied with (or added to, when "inDecibels" 
  is true) to the magnitudes which are already there in the "magnitudes"-array. This is useful for 
  calculating the magnitude response of several biquad-cascades in series. */
  void getMagnitudeResponse(TCoef* w, TCoef* magnitudes, int numBins, bool inDecibels = false, 
    bool accumulate = false) const;

  /** Writes the magnitdue-response of a biquad-cascade at the physical frequencies given in 
  'frequencies' into the array 'magnitudes'. */
  void getMagnitudeResponse(TCoef* frequencies, CRCoef sampleRate, TCoef* magnitudes, 
    int numBins, bool inDecibels = false, bool accumulate = false) const;

  /** Returns an estimate of the time (in samples) it takes for the impulse response to ring 
  down to "threshold". The estimate is based on the pole which has the highest Q (i.e. is 
  closest to the unit circle). */
  TCoef getRingingTimeEstimate(CRCoef threshold) const;
  // todo: write a function getRingingFrequency - the normalized radian ringing frequency is
  // the angle of the pole that is closest to the unit circle

  //-----------------------------------------------------------------------------------------------
  /** \name Audio Processing */

  /** Calculates a single filtered output-sample via a cascade of biquads in Direct-Form 1 */
  inline TSig getSampleDirect1(CRSig in);

  /** Calculates a single filtered output-sample via a cascade of biquads in Direct-Form 2 */
  inline TSig getSampleDirect2(CRSig in);

  /** Calculates a single filtered output-sample via a cascade of biquads. */
  inline TSig getSample(CRSig in) { return getSampleDirect1(in); }

  //-----------------------------------------------------------------------------------------------
  /** \name Misc */

  /** Initializes the biquad coefficients to neutral values. */
  void initBiquadCoeffs();

  /** Sets the buffers for the previous input and output doubles of all biquad stages to zero. */
  void reset();


protected:

  TCoef *a1, *a2, *b0, *b1, *b2;  // filter coefficients
  TSig  *x1, *x2, *y1, *y2;       // buffering
  int numStages;                  // current number of biquad-stages
  int maxNumStages;               // maximum number of biquad-stages
};

//-------------------------------------------------------------------------------------------------
// inlined functions:

template<class TSig, class TCoef>
inline void rsBiquadCascade<TSig, TCoef>::setCoeffs(TCoef* newB0, TCoef* newB1, TCoef* newB2, 
  TCoef* newA1, TCoef* newA2)
{
  for(int i = 0; i < numStages; i++)
  {
    b0[i] = newB0[i];
    b1[i] = newB1[i];
    b2[i] = newB2[i];
    a1[i] = newA1[i];
    a2[i] = newA2[i];
  }
}

template<class TSig, class TCoef>
inline TSig rsBiquadCascade<TSig, TCoef>::getSampleDirect1(CRSig in)
{
  TSig tmp, tmp2;
  int  i;  // for the loop through the stages 

  tmp = in;

  // calculate current output-sample (y[n]) of all the BiQuad-stages (the output of one stage is 
  // the input for the next stage):
  for(i = 0; i < numStages; i++)
  {
    tmp2 = tmp; // for x1[i]

    // calculate current output-sample (y[n]) of BiQuad-stage i:
    tmp = b0[i]*tmp + (b1[i]*x1[i] + b2[i]*x2[i]) - (a1[i]*y1[i] + a2[i]*y2[i]);

    // set x[n-1], x[n-2], y[n-1] and y[n-2] for the next iteration:
    x2[i] = x1[i];
    x1[i] = tmp2;
    y2[i] = y1[i];
    y1[i] = tmp;
  }

  return tmp;
}

template<class TSig, class TCoef>
inline TSig rsBiquadCascade<TSig, TCoef>::getSampleDirect2(CRSig in)
{
  TSig x, y, g;
  y = x = in; // y = ... to make it work also with 0 stages as bypass

  // calculate current output-sample (y[n]) of all the BiQuad-stages (the output of one stage is 
  // the input for the next stage):
  for(int i = 0; i < numStages; i++)
  {
    // calculate current output-sample (y[n]) of BiQuad-stage i:
    //g = x - a1[i]*y1[i] - a2[i]*y2[i];
    g = x - (a1[i]*y1[i] + a2[i]*y2[i]);
    //g = x + a1[i]*y1[i] + a2[i]*y2[i]; // this is wrong (a has wrong sign) - for performance test
    y = b0[i]*g + b1[i]*y1[i] + b2[i]*y2[i];

    // set g[n-1], g[n-2] for the next iteration:
    y2[i] = y1[i];
    y1[i] = g;

    x = y; // output of one stage is input to the next
  }

  return y;
}

#endif