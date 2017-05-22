#ifndef rosic_CookbookFilter_h
#define rosic_CookbookFilter_h

//// rosic-indcludes:
//#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

 /**

 This class implements a cascade of up to 5 (= const short maxStages) simple
 BiQuad stages. The filter can calculate it's own coefficients according
 to its member variables mode, freq, q, gain, sampleRate which are set
 up via the corresponding setFoo()-functions. To do the coefficient
 calculation, the function calcCoeffs() has to be called (it is not called
 automatically when one of the parameters changes). This function uses the
 filter design formulas taken from "Cookbook formulae for audio EQ biquad
 filter coefficients" by Robert Bristow-Johnson. Alternatively, the
 coefficients can be set manually via setCoeffs(). Each stage will use the
 same set of filter-coefficients. Simply use the setNumStages()-function to
 specify the number of stages you need. The class has been designed to
 facilitate real-time (i.e. per sample) adjustments of the filters parameters.
 The inconvinience, that the user has to call the calcCoeffs()-function to
 trigger a coefficient recalculation after changing some filter-parameter
 avoids multiple calculations of the same coefficient set when more than one
 parameter changes at a time. The resulting biquad (-cascade) is implemented
 in various topologies, namely direct form 1 & 2 and in normalized ladder form
 1 & 2. If you want to use the ladder-form (which is recommendable for rapidly
 time-varying filters), make sure to always call convertDirectToLadder() after
 each coefficient update.

 */

 class CookbookFilter
 {
 public:

  /** This is an enumeration of the 9 available filter modes. */
  enum modes
  {
   BYPASS = 0,
   LOWPASS,
   HIGHPASS,
   BANDPASS_CONST_SKIRT,
   BANDPASS_CONST_PEAK,
   BANDREJECT,
   ALLPASS,
   PEAK,
   LOW_SHELF,
   HIGH_SHELF
  };

  //---------------------------------------------------------------------------
  // construction/destruction:

  CookbookFilter();   ///< Constructor.
  ~CookbookFilter();  ///< Destructor.

  //---------------------------------------------------------------------------
  // parameter settings:

  void setSampleRate(double newSampleRate);
  ///< Overrides the setSampleRate() method of the AudioModule base class.

  void setMode(int newMode);
  /**< Sets the mode of the filter. 9 modes are available: 1:lowpass,
       2:highpass, 3:bandpass(constant skirt gain), 4:bandpass(constant
       peak gain), 5:bandreject, 6:allpass, 7:peaking, 8:low-shelf,
       9:high-shelf. */

  void setNumStages(int newNumStages);
  ///< Lets the user set the number of biquad stages.

  // parameters that are supposed to be updated at sample-rate (inlined
  // access-functions):

  INLINE void setFreq(double newFreq);
  ///< Sets the cuttoff/center frequency of the filter.

  INLINE void setQ(double newQ);
  ///< Sets Q of the filter.

  INLINE void setGain(double newGain);
  //< Sets the gain value for shelving and peaking modes.

  INLINE void setCoeffs(double newA1, double newA2,
                        double newB0, double newB1,
                        double newB2);
  ///< Allows the user to set the filter coefficients manually.

  INLINE void calcCoeffs();
  /**< Calculates filter coefficients from filter parameters. Has to be called
       each time Freq, Q or Gain changes -> it is NOT called automatically
       by setFreq, setQ or setGain to avoid multiple calculations of the coeffs
       when more than one of these parameters is changed at the same time
       (for example by an envelope).

       costs (AMD Athlon 64 3200+): 240-375 CPU-cycles (depending on the mode) */

  INLINE void convertDirectToLadder();
  /**< This function converts the set of direct form coefficients
       (b0, b1, b2, a1, a2) into an alternative set of coefficients
       (k1, k2, c1, c2) which will be used in the normalized ladder structure.
       If you want to use the noramlized ladder-topology (which is
       recommendable for rapidly time-varying filters), call this function after
       a call to calcCoeffs() or setCoeffs(). */

  //---------------------------------------------------------------------------
  // audio processing:

  INLINE double getSample(double in);
  /**< Calculates a single filtered output-sample via the chosen implementaion strcuture.  */

  INLINE double getSampleDirect1(double in);
  /**< Calculates a single filtered output-sample via Direct Form 1.
       costs (AMD Athlon 64 3200+): ~22.5 cycles per stage */

  INLINE double getSampleDirect2(double in);
  /**< Calculates a single filtered output-sample via Direct Form 2. */

  INLINE double getSampleLadder1(double in);
  /**< Calculates a single filtered output-sample via Normalized Ladder Form
       in pole-before-zero configuration. */

  //INLINE double getSampleLadder2(double in);
  /**< Calculates a single filtered output-sample via Normalized Ladder Form
       in zero-before-pole configuration. */

  INLINE double getSampleAutoChoose(double in);
  /**< Calculates a single filtered output-sample. The most appropriate
       structure is chosen automatically according to the direction of the
       cutoff-change. */

  //---------------------------------------------------------------------------
  // others:

  void resetBuffers ();
  /**< Sets the buffers for the previous input and output samples of all biquad
       stages to zero. */

  //===========================================================================

 protected:

  static const int maxNumStages = 5;

  // direct form coefficients:
  doubleA a0, a1, a2, b0, b1, b2;

  // normalized ladder coefficients:
  doubleA k1, k2, c1, c2, ladderGain;

  // filter parameters:
  intA numStages;
  intA mode;
  doubleA freq;
  doubleA oldFreq; // to determine direction of freq-change in order to choose
                   // the most appropriate filter-topology
  doubleA q;
  doubleA gain, A; // gain in dB and as raw amplitude-factor

  doubleA sampleRate;
  doubleA sampleRateRec;  // reciprocal of the sample-rate

  // buffering:
  doubleA x1[maxNumStages];   //arrays of previous biquad input samples,
  doubleA x2[maxNumStages];   //array index indicates the biquad-stage
  doubleA y1[maxNumStages];   //arrays of previous output samples
  doubleA y2[maxNumStages];

  doubleA g1[maxNumStages]; // arrays of previous intermediate samples for DF2
  doubleA g2[maxNumStages];

  // buffers for the intermediate variables in the ladder, first number is an
  // index to the signal, second number indicates a time-delay:
  doubleA e2_0[maxNumStages];
  doubleA e2_1[maxNumStages];
  doubleA e1_0[maxNumStages];
  doubleA e1_1[maxNumStages];
  doubleA e0_0[maxNumStages];
  doubleA e0_1[maxNumStages];
  doubleA e0_2[maxNumStages];
  doubleA e1t_0[maxNumStages];
  doubleA e1t_1[maxNumStages];
  doubleA e2t_0[maxNumStages];
  doubleA e2t_1[maxNumStages];
 };

 //----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):
 INLINE void CookbookFilter::setFreq(double newFreq)
 {
  oldFreq = freq; // remember the previous freq-setting

  //restrict cutoff frequency to the range between 20 and 20000 Hz:
  if( newFreq <= 20.0 )
   freq = 20.0;
  else if( newFreq >= 20000.0 )
   freq = 20000.0;
  else
   freq = newFreq;
 }

 INLINE void CookbookFilter::setQ(double newQ)
 {
  if( newQ > 0 )
   q = newQ;
 }

 INLINE void CookbookFilter::setGain(double newGain)
 {
  gain = newGain;
  A    = pow(10, (0.025*gain) );
 }

 INLINE void CookbookFilter::setCoeffs(double newA1,
                                       double newA2,
                                       double newB0,
                                       double newB1,
                                       double newB2)
 {
  a1 = newA1;
  a2 = newA2;
  b0 = newB0;
  b1 = newB1;
  b2 = newB2;
  return;
 }

 INLINE void CookbookFilter::calcCoeffs()
 {
  static doubleA omega, sine, cosine, alpha, beta, a0Rec;

  if(mode == BYPASS)
  {
   a0 = 1.0;
   a1 = 0.0;
   a2 = 0.0;
   b0 = 1.0;
   b1 = 0.0;
   b2 = 0.0;
   return;
  }

  // calculate intermediate variables:
  omega  = (2.0*PI)*freq*sampleRateRec;
  //sine   = sin(omega);
  //cosine = cos(omega);
  sinCos(omega, &sine, &cosine); // cheaper than sin and cos
  alpha  = sine/(2.0*q);

  switch(mode)
  {
   //calculate LPF-coefficients:
   case LOWPASS:
   {
    a0 = 1+alpha;
    a1 = -2*cosine;
    a2 = 1-alpha;

    b1 = 1-cosine;
    b0 = b1/2;
    b2 = b0;
   }
   break;

   //calculate HPF-coefficients:
   case HIGHPASS:
   {
    a0 = 1+alpha;
    a1 = -2*cosine;
    a2 = 1-alpha;

    b1 = -(1+cosine);
    b0 = (-b1)/2;
    b2 = b0;
   }
   break;

   //calculate BPF (constant skirt) - coefficients:
   case BANDPASS_CONST_SKIRT:
   {
    a0 = 1+alpha;
    a1 = -2*cosine;
    a2 = 1-alpha;

    b1 = 0;
    b0 = q*alpha;
    b2 = -b0;
   }
   break;

   //calculate BPF (constant peak) - coefficients:
   case BANDPASS_CONST_PEAK:
   {
    a0 = 1+alpha;
    a1 = -2*cosine;
    a2 = 1-alpha;

    b0 = alpha;
    b1 = 0;
    b2 = -alpha;
   }
   break;

   //calculate BRF-coefficients:
   case BANDREJECT:
   {
    a0 = 1+alpha;
    a1 = -2*cosine;
    a2 = 1-alpha;

    b0 = 1;
    b1 = -2*cosine;
    b2 = 1;
   }
   break;

   //calculate APF-coefficients:
   case ALLPASS:
   {
    a0 = 1+alpha;
    a1 = -2*cosine;
    a2 = 1-alpha;

    b0 = 1-alpha;
    b1 = -2*cosine;
    b2 = 1+alpha;
   }
   break;

   //calculate peaking-coefficients:
   case PEAK:
   {
    a0 = 1+alpha/A;
    a1 = -2*cosine;
    a2 = 1-alpha/A;

    b0 = 1+alpha*A;
    b1 = -2*cosine;
    b2 = 1-alpha*A;
   }
   break;

   //calculate low-shelf-coefficients:
   case LOW_SHELF:
   {
    beta   = sqrt(A)/q;

    a0 =         (A+1) + (A-1)*cosine + beta*sine;
    a1 =  -2*  ( (A-1) + (A+1)*cosine              );
    a2 =         (A+1) + (A-1)*cosine - beta*sine;

    b0 =   A*  ( (A+1) - (A-1)*cosine + beta*sine  );
    b1 = 2*A*  ( (A-1) - (A+1)*cosine              );
    b2 =   A*  ( (A+1) - (A-1)*cosine - beta*sine  );
   }
   break;

   //calculate high-shelf-coefficients:
   case HIGH_SHELF:
   {
    beta   = sqrt(A)/q;

    a0 =          (A+1) - (A-1)*cosine + beta*sine;
    a1 =    2*  ( (A-1) - (A+1)*cosine              );
    a2 =          (A+1) - (A-1)*cosine - beta*sine;

    b0 =    A*  ( (A+1) + (A-1)*cosine + beta*sine  );
    b1 = -2*A*  ( (A-1) + (A+1)*cosine              );
    b2 =    A*  ( (A+1) + (A-1)*cosine - beta*sine  );
   }
   break;

   //if mode not in 1..9 calculate LPF-coefficients by default
   default:
   {
    a0 = 1+alpha;
    a1 = -2*cosine;
    a2 = 1-alpha;

    b1 = 1-cosine;
    b0 = b1/2;
    b2 = b0;
   }
   //the equations are taken from "Cookbook formulae for audio EQ biquad
   //filter coefficients" by Robert Bristow-Johnson

  } //end switch(mode)

  // scale all coefficients by (1/a0):
  a0Rec  = 1.0 / a0;
  a1    *= a0Rec;
  a2    *= a0Rec;
  b0    *= a0Rec;
  b1    *= a0Rec;
  b2    *= a0Rec;
 }

 INLINE void CookbookFilter::convertDirectToLadder()
 {
  k2         = a2;
  k1         = (a1-a2*a1) / (1-a2*a2);
  c2         = sqrt(1-k2*k2);
  c1         = sqrt(1-k1*k1);
  ladderGain = 1/(c1*c2);
 }

 INLINE double CookbookFilter::getSample(double in)
 {
   return getSampleDirect1(in);
 }

 INLINE double CookbookFilter::getSampleDirect1(double in)
 {
  doubleA tmp, tmp2;
  intA    i;  // for the loop through the stages

  tmp = in;

  // calculate current output-sample (y[n]) of all the BiQuad-stages
  // (the output of one stage is the input for the next stage):
  for (i=0; i<numStages; i++)
  {
   tmp2 = tmp; //for x_1[i]

   // calculate current output-sample (y[n]) of BiQuad-stage i:
   tmp = b0*(tmp) + b1*x1[i] + b2*x2[i]
                  - a1*y1[i] - a2*y2[i];

   // set x[n-1], x[n-2], y[n-1] and y[n-2] for the next call:
   x2[i] = x1[i];
   x1[i] = tmp2;
   y2[i] = y1[i];
   y1[i] = tmp;
  }

  return tmp;
 }

 INLINE double CookbookFilter::getSampleDirect2(double in)
 {
  static doubleA x, y, g;
  static intA    i;  //for the loop through the stages

  x = in;

  // calculate current output-sample (y[n]) of all the BiQuad-stages
  // (the output of one stage is the input for the next stage):
  for (i=0; i<numStages; i++)
  {
   // calculate intermediate signal g[n] and output-signal y[n] of
   // BiQuad-stage i:
   g = x - a1*g1[i] - a2*g2[i];
   y = b0*g + b1*g1[i] + b2*g2[i];

   // set g[n-1], g[n-2] for the next call:
   g2[i] = g1[i];
   g1[i] = g;

   x = y; // output of one stage is input to the next
  }

  return y;
 }

 INLINE double CookbookFilter::getSampleLadder1(double in)
 {
  static doubleA x, y;
  static intA    i;  // for the loop through the stages

  x = in;

  // calculate current output-sample (y[n]) of all the Ladder-stages
  // (the output of one stage is the input for the next stage):
  for (i=0; i<numStages; i++)
  {
   e2_0[i]  = x;
   e1_0[i]  = c2*e2_0[i]  - k2*e1t_1[i];
   e0_0[i]  = c1*e1_0[i]  - k1*e0_1[i];  // e0t equal to e0, no extra variable
   e1t_0[i] = c1*e0_1[i]  + k1*e1_0[i];
   e2t_0[i] = c2*e1t_1[i] + k2*e2_0[i];

   // calculate ouput of this stage (apply the zeros and the gain):
   y = ladderGain * (b0*e0_0[i] + b1*e0_1[i] + b2*e0_2[i]);

   // update the buffers:
   e2_1[i]  = e2_0[i];
   e2t_1[i] = e2t_0[i];
   e1_1[i]  = e1_0[i];
   e1t_1[i] = e1t_0[i];
   e0_2[i]  = e0_1[i];
   e0_1[i]  = e0_0[i];

   x = y; // output of one stage is input to the next
  }

  return y;
 }

 //INLINE double CookbookFilter::getSampleLadder2(double in)
 //{
 // static doubleA x, y;
 // static intA    i;  // for the loop through the stages

 // x = in;

 // // calculate current output-sample (y[n]) of all the Ladder-stages
 // // (the output of one stage is the input for the next stage):
 // for (i=0; i<numStages; i++)
 // {

 // }

 // return y;
 //}

 INLINE double CookbookFilter::getSampleAutoChoose(double in)
 {
  static doubleA x, y, g, tmp;
  static intA    i;  //for the loop through the stages

  x = in;

  // calculate current output-sample (y[n]) of all the BiQuad-stages
  // (the output of one stage is the input for the next stage):
  for (i=0; i<numStages; i++)
  {
   // calculate intermediate signal g[n] and output-signal y[n] of
   // BiQuad-stage i:
   g = x - a1*g1[i] - a2*g2[i];

   if( oldFreq >= freq ) // freq sweeps up - use DF1-topology
   {
    // calculate current output-sample (y[n]) of BiQuad-stage i:
    y = b0*(tmp) + b1*x1[i] + b2*x2[i]
                 - a1*y1[i] - a2*y2[i];
   }
   else  // freq sweeps down - use DF2-topology
   {
    y = b0*g + b1*g1[i] + b2*g2[i];
   }

   // set x[n-1], x[n-2], y[n-1], y[n-2] g[n-1], g[n-2] for the next call:
   x2[i] = x1[i];
   x1[i] = x;
   y2[i] = y1[i];
   y1[i] = y;
   g2[i] = g1[i];
   g1[i] = g;

   x = y; // output of one stage is input to the next
  }

  return y;
 }


} // end namespace rosic

#endif // rosic_CookbookFilter_h
