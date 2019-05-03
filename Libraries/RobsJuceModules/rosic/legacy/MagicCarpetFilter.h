#ifndef MagicCarpetFilter_h
#define MagicCarpetFilter_h

//#include "AudioModule.h"
#include "Definitions.h"
#include "MoreMath.h"
using namespace MoreMath;

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
parameter changes at a time. 

*/

class MagicCarpetFilter
{
public:

 /** This is an enumeration of the 9 available filter modes. */
 enum modes
 {
  BYPASS = 0,
  LOWPASS_6,
  LOWPASS_12,
  HIGHPASS_6,
  HIGHPASS_12,
  BANDPASS,
  BANDREJECT,
  PEAK,
  LOW_SHELF,
  HIGH_SHELF
 };

 //---------------------------------------------------------------------------
 // construction/destruction:

 MagicCarpetFilter();   ///< Constructor.     
 ~MagicCarpetFilter();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 void setMode(int newMode);      
 /**< Sets the mode of the filter. 9 modes are available: 1:lowpass, 
      2:highpass, 3:bandpass(constant skirt gain), 4:bandpass(constant
      peak gain), 5:bandreject, 6:allpass, 7:peaking, 8:low-shelf,
      9:high-shelf. */

 void setTwoStages(bool newTwoStagesSwitch);    
 ///< Lets the user switch between one or two biquad stages.

 // parameters that are supposed to be updated at sample-rate (inlined
 // access-functions):

 INLINE void setFreq(double newFreq);      
 ///< Sets the cuttoff/center frequency of the filter.

 INLINE void setQ(double newQ);         
 ///< Sets Q of the filter.

 INLINE void setGain(double newGain); 
 //< Sets the gain value for shelving and peaking modes.

 INLINE void calcCoeffs(); 
 /**< Calculates filter coefficients from filter parameters. Has to be called
      each time Freq, Q or Gain changes -> it is NOT called automatically
      by setFreq, setQ or setGain to avoid multiple calculations of the coeffs
      when more than one of these parameters is changed at the same time
      (for example by an envelope). 

      costs (AMD Athlon 64 3200+): 240-375 CPU-cycles (depending on the mode) */



 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR); 
 /**< Calculates a single filtered output-sample via Direct Form 1. Make sure
      that the input and output-slots are distinct, or it won't work.  */

 //---------------------------------------------------------------------------
 // others:

 void resetBuffers (); 
 /**< Sets the buffers for the previous input and output samples of all biquad
      stages to zero. */

 //===========================================================================

protected:

 // direct form coefficients:
 doubleA a1, a2, b0, b1, b2;

 // buffering:
 doubleA x_s1_d1_L;   // input to stage 1, delayed by 1, left channel
 doubleA x_s1_d2_L;   // input to stage 1, delayed by 2, left channel
 doubleA y_s1_d1_L;   // output of stage 1, delayed by 1, left channel
 doubleA y_s1_d2_L;   // output of stage 1, delayed by 2, left channel

 doubleA x_s1_d1_R;   // input to stage 1, delayed by 1, right channel
 doubleA x_s1_d2_R;   // input to stage 1, delayed by 2, right channel
 doubleA y_s1_d1_R;   // output of stage 1, delayed by 1, right channel
 doubleA y_s1_d2_R;   // output of stage 1, delayed by 2, right channel

 //doubleA x_s2_d1_L;   // input to stage 2, delayed by 1, left channel
 //doubleA x_s2_d2_L;   // input to stage 2, delayed by 2, left channel
 doubleA y_s2_d1_L;   // output of stage 2, delayed by 1, left channel
 doubleA y_s2_d2_L;   // output of stage 2, delayed by 2, left channel

 //doubleA x_s2_d1_R;   // input to stage 2, delayed by 1, right channel
 //doubleA x_s2_d2_R;   // input to stage 2, delayed by 2, right channel
 doubleA y_s2_d1_R;   // output of stage 2, delayed by 1, right channel
 doubleA y_s2_d2_R;   // output of stage 2, delayed by 2, right channel

 // filter parameters:
 doubleA freq; 
 doubleA q;
 doubleA gain, A; // gain in dB and as raw amplitude-factor
 intA    mode;  
 bool    twoStages;

 doubleA sampleRate;
 doubleA sampleRateRec;  // reciprocal of the sample-rate
};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):
INLINE void MagicCarpetFilter::setFreq(double newFreq)
{
 // restrict cutoff frequency to the range between 20 and 20000 Hz:
 if( newFreq <= 20.0 )
  freq = 20.0;
 else if( newFreq >= 20000.0 )
  freq = 20000.0;
 else 
  freq = newFreq;
}

INLINE void MagicCarpetFilter::setQ(double newQ)
{
 if( newQ > 0.01 )
  q = newQ;
}

INLINE void MagicCarpetFilter::setGain(double newGain)
{
 gain = newGain;
 A    = pow(10, (0.025*gain) );
}

INLINE void MagicCarpetFilter::calcCoeffs()
{
 static doubleA a0, omega, sine, cosine, alpha, beta, a0Rec;

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
 //MoreMath::sincos(omega, &sine, &cosine); // cheaper than sin and cos
 //alpha  = sine/(2.0*q);

 switch(mode)
 {
  //calculate LPF-coefficients:
  case LOWPASS_6:
  {
   a1 = -exp( -omega );
   b0 = 1.0 + a1;
   b1 = 0.0;
   b2 = 0.0;
   a2 = 0.0;
  }
  break;

  //calculate LPF-coefficients:
  case LOWPASS_12:
  {
   // calculate sine and cosine of omega::
   MoreMath::sinCos(omega, &sine, &cosine); // cheaper than sin and cos
   alpha  = sine/(2.0*q);

   // set up parameters
   a0 = 1+alpha;
   a1 = -2*cosine;
   a2 = 1-alpha;

   b1 = 1-cosine;
   b0 = b1/2;
   b2 = b0;
  }
  break;

  //calculate HPF-coefficients:
  case HIGHPASS_6:
  {
   a1 = -exp( -omega );
   b0 =  0.5*(1.0-a1);
   b1 = -0.5*(1.0-a1);
   b2 = 0.0;
   a2 = 0.0;
  }
  break;

  //calculate HPF-coefficients:
  case HIGHPASS_12:
  {
   MoreMath::sinCos(omega, &sine, &cosine);
   alpha  = sine/(2.0*q);

   a0 = 1+alpha;
   a1 = -2*cosine;
   a2 = 1-alpha;

   b1 = -(1+cosine);
   b0 = (-b1)/2;
   b2 = b0;
  }
  break;

  //calculate BPF (constant skirt) - coefficients:
  case BANDPASS:
  {
   MoreMath::sinCos(omega, &sine, &cosine);
   alpha  = sine/(2.0*q);

   a0 = 1+alpha;
   a1 = -2*cosine;
   a2 = 1-alpha;

   b1 = 0;
   b0 = q*alpha;
   b2 = -b0;
  }
  break;

  /*
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
  */

  //calculate BRF-coefficients:
  case BANDREJECT:
  {
   MoreMath::sinCos(omega, &sine, &cosine);
   alpha  = sine/(2.0*q);

   a0 = 1+alpha;
   a1 = -2*cosine;
   a2 = 1-alpha;

   b0 = 1;
   b1 = -2*cosine;
   b2 = 1;
  }
  break;

  /*
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
  */

  //calculate peaking-coefficients:
  case PEAK:
  {
   MoreMath::sinCos(omega, &sine, &cosine);
   alpha  = sine/(2.0*q);

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
   MoreMath::sinCos(omega, &sine, &cosine);
   alpha  = sine/(2.0*q);
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
   MoreMath::sinCos(omega, &sine, &cosine);
   alpha  = sine/(2.0*q);
   beta   = sqrt(A)/q;

   a0 =          (A+1) - (A-1)*cosine + beta*sine;
   a1 =    2*  ( (A-1) - (A+1)*cosine              );
   a2 =          (A+1) - (A-1)*cosine - beta*sine;

   b0 =    A*  ( (A+1) + (A-1)*cosine + beta*sine  );
   b1 = -2*A*  ( (A-1) + (A+1)*cosine              );
   b2 =    A*  ( (A+1) + (A-1)*cosine - beta*sine  );
  }
  break;

  // if mode is out of range, calculate BYPASS-coefficients by default
  default:
  {
  a0 = 1.0;
  a1 = 0.0;
  a2 = 0.0; 
  b0 = 1.0;
  b1 = 0.0;
  b2 = 0.0;  
  }
  //the equations are taken from "Cookbook formulae for audio EQ biquad
  //filter coefficients" by Robert Bristow-Johnson

 } //end switch(mode)

 // scale all coefficients by (1/a0):
 a0Rec  = 1.0 / a0;
 a1     = -a1*a0Rec;
 a2     = -a2*a0Rec;
 b0    *= a0Rec;
 b1    *= a0Rec;
 b2    *= a0Rec;
}



INLINE void MagicCarpetFilter::getSampleFrameStereo(double *inL, 
                                                    double *inR, 
                                                    double *outL, 
                                                    double *outR)
{
 static doubleA tmpL, tmpR;

 tmpL = *inL;
 tmpR = *inR;

 // apply first stage:
 tmpL = b0*tmpL + b1*x_s1_d1_L + b2*x_s1_d2_L
                + a1*y_s1_d1_L + a2*y_s1_d2_L
                + TINY;
 tmpR = b0*tmpR + b1*x_s1_d1_R + b2*x_s1_d2_R
                + a1*y_s1_d1_R + a2*y_s1_d2_R
                + TINY;

 // update first buffer:
 x_s1_d2_L = x_s1_d1_L;
 x_s1_d1_L = *inL;
 y_s1_d2_L = y_s1_d1_L;
 y_s1_d1_L = tmpL;

 x_s1_d2_R = x_s1_d1_R;
 x_s1_d1_R = *inR;
 y_s1_d2_R = y_s1_d1_R;
 y_s1_d1_R = tmpR;

 *outL = tmpL;
 *outR = tmpR;

 if( twoStages == true )
 {
  // apply second stage:
  *outL = b0*tmpL + b1*y_s1_d1_L + b2*y_s1_d2_L
                  + a1*y_s2_d1_L + a2*y_s2_d2_L;
  *outR = b0*tmpR + b1*y_s1_d1_R + b2*y_s1_d2_R
                  + a1*y_s2_d1_R + a2*y_s2_d2_R;

  // update second buffer:
  //x_s2_d2_L = x_s2_d1_L;
  //x_s2_d1_L = tmpL;
  y_s2_d2_L = y_s2_d1_L;
  y_s2_d1_L = *outL;

  //x_s2_d2_R = x_s2_d1_R;
  //x_s2_d1_R = tmpR;
  y_s2_d2_R = y_s2_d1_R;
  y_s2_d1_R = *outR;
 }
 else
 {
  *outL = tmpL;
  *outR = tmpR;
 }

}

#endif // MagicCarpetFilter_h
