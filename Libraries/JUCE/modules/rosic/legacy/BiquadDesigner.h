#ifndef BiquadDesigner_h
#define BiquadDesigner_h

//#include "AudioModule.h"
#include "MoreMath.h"
using namespace MoreMath;


/**

This class calculates the b0, b1, b2, a1, a2 coefficients for a digital
biquad-filter of the form:
y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2]
               - a1*y[n-1] - a2*y[n-2]
Several different design formulas are provided including the 
"Cookbook"-designs by Robert Bristow Johnson (RBJ).

*/

class BiquadDesigner
{
public:

 /**< This is an enumeration of the 9 available filter modes. */
 enum modes
 {
  BYPASS = 0,           // bypass

  RBJ_LPF,              // RBJ-cookbook lowpass
  RBJ_HPF,              // RBJ-cookbook highpass
  RBJ_BPF1,             // RBJ-cookbook bandpass (constant skirt width)
  RBJ_BPF2,             // RBJ-cookbook bandpass (constant peak width)
  RBJ_BRF,              // RBJ-cookbook bandreject
  RBJ_APF,              // RBJ-cookbook allpass
  RBJ_PEAK,             // RBJ-cookbook peak/bell
  RBJ_LOW_SHELV,        // RBJ-cookbook low-shelving
  RBJ_HIGH_SHELV,       // RBJ-cookbook high-shelving

  RESON,                // two-pole resonator (Steiglitz)
  RESONZ,               // two-pole-one-zero resonator (Steiglitz)

  ONE_POLE_LPF,         // one-pole lowpass ("RC-filter" from dspguide)
  ONE_POLE_HPF,         // one-pole highpass ("RC-filter" from dspguide)
  ONE_POLE_APF,         // one-pole allpass (from DAFX)
  ONE_POLE_LOW_SHELV,   // one-pole low-shelving (from DAFX)
  ONE_POLE_HIGH_SHELV,  // one-pole high-shelving (from DAFX)

  RC_LPF_HPF_SERIES,    // 1-pole-LPF and -HPF in series
  RC_2LPFS_SERIES,      // 2 1-pole LPFs in series
  RC_2HPFS_SERIES,      // 2 1-pole HPFs in series
  RC_LPF_HPF_PARALLEL,  // 1-pole-LPF and -HPF in parallel

  BUTTER_LPF,           // butterworth-LPF of 2nd order
  BUTTER_HPF,           // butterworth-HPF of 2nd order
  BUTTER_BPF,           // butterworth-BPF of 2nd order
  BUTTER_BRF,           // butterworth-BRF of 2nd order

  NOTCH,                // deep notch: with zeros on the unit circle and poles
                        // near the unit-circle to control bandwidth

  LPF_BPF_HPF_MORPH,    // 2-pole-resonator witzh two zeros placed either both 
                        // at z=1, one at z=1 and the other at z=-1 or both at 
                        // z=-1 - this gives a morph between LPF/BPF/HPF 
  LPF_RESON_HPF_MORPH,  // similar as above, but the zeros move together from 
                        // z=1 to z=-1 - this gives a morph between 
                        // LPF/RESON/HPF

  ONE_POLE_APFS,        // 2 1-pole APFs in series

  FOF_FILTER,           // realizes a FOF-filter specified by the 
                        // FOF-parameters fofA, fofAlpha, fofOmega, fofPhi

  FREE_IN_THE_PLANE     // lets the use specify the positions of the poles 
                        // and zeros in the z-plane

 };

 //---------------------------------------------------------------------------
 // construction/destruction:

          BiquadDesigner();  ///< Constructor.     
 virtual ~BiquadDesigner();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual void setMode(int newMode);      
 /**< Sets the mode of the filter. For the available modes see the 
      enumaration above. */

 // parameters that are supposed to be updated at double-rate (inlined
 // access-functions):

 virtual INLINE void setFreq(double newFreq);      
 ///< Sets the cuttoff/center frequency of the filter.

 virtual INLINE void setQ(double newQ);         
 ///< Sets Q of the filter.

 virtual INLINE void setGain(double newGain); 
 //< Sets the gain value for shelving and peaking modes.

 virtual INLINE void getCoeffs(double *b0, double *b1, double *b2,
                                           double *a1, double *a2);
 /**< Calculates biquad-coefficients according to the filter parameters and 
      stores them at the adresses which have been passed by the pointers. */

 //===========================================================================

protected:

 doubleA sampleRate;        // the sample-rate
 doubleA sampleRateRec;     // reciprocal of the sample-rate

 // filter parameters:
 intA mode;  
 doubleA freq;
 doubleA q;
 doubleA gain;

};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):
INLINE void BiquadDesigner::setFreq(double newFreq)
{
 //restrict cutoff frequency to the range between 20 and 20000 Hz:
 if( newFreq <= 20.0 )
  freq = 20.0;
 else if( newFreq >= 20000.0 )
  freq = 20000.0;
 else 
  freq = newFreq;
}

INLINE void BiquadDesigner::setQ(double newQ)
{
 if( newQ > 0.001 )
  q = newQ;
}

INLINE void BiquadDesigner::setGain(double newGain)
{
 gain = newGain;
}


INLINE void BiquadDesigner::getCoeffs(double *b0, double *b1, double *b2,
                                                  double *a1, double *a2)
{
 static doubleA A, omega, sine, cosine, alpha, beta, a0, a0Rec;

 switch(mode)
 {
  // "calculate" Bypass-coefficients:
  case BYPASS:
  {
   *b0 = 1.0;
   *b1 = 0.0;
   *b2 = 0.0;

   a0  = 1.0;
   *a1 = 0.0;
   *a2 = 0.0; 
  }
  break;

  //calculate LPF-coefficients:
  case RBJ_LPF:
  {
   omega  = 2*PI*freq*sampleRateRec;
   cosine = cos(omega);
   sine   = sin(omega);
   alpha  = sine/(2*q);

   a0  = 1+alpha;
   *a1 = -2*cosine;
   *a2 = 1-alpha;

   *b1 = 1-cosine;
   *b0 = *b1/2;
   *b2 = *b0;
  }
  break;

  //calculate HPF-coefficients:
  case RBJ_HPF:
  {
   omega  = 2*PI*freq*sampleRateRec;
   cosine = cos(omega);
   sine   = sin(omega);
   alpha  = sine/(2*q);

   a0  = 1+alpha;
   *a1 = -2*cosine;
   *a2 = 1-alpha;

   *b1 = -(1+cosine);
   *b0 = (-*b1)/2;
   *b2 = *b0;
  }
  break;

  //calculate BPF (constant skirt) - coefficients:
  case RBJ_BPF1:
  {
   omega  = 2*PI*freq*sampleRateRec;
   cosine = cos(omega);
   sine   = sin(omega);
   alpha  = sine/(2*q);

   a0  = 1+alpha;
   *a1 = -2*cosine;
   *a2 = 1-alpha;

   *b1 = 0;
   *b0 = q*alpha;
   *b2 = -*b0;
  }
  break;

  //calculate BPF (constant peak) - coefficients:
  case RBJ_BPF2:
  {
   omega  = 2*PI*freq*sampleRateRec;
   cosine = cos(omega);
   sine   = sin(omega);
   alpha  = sine/(2*q);

   a0  = 1+alpha;
   *a1 = -2*cosine;
   *a2 = 1-alpha;

   *b0 = alpha;
   *b1 = 0;
   *b2 = -alpha;
  }
  break;

  //calculate BRF-coefficients:
  case RBJ_BRF:
  {
   omega  = 2*PI*freq*sampleRateRec;
   cosine = cos(omega);
   sine   = sin(omega);
   alpha  = sine/(2*q);

   a0  = 1+alpha;
   *a1 = -2*cosine;
   *a2 = 1-alpha;

   *b0 = 1;
   *b1 = -2*cosine;
   *b2 = 1;
  }
  break;

  //calculate APF-coefficients:
  case RBJ_APF:
  {
   omega  = 2*PI*freq*sampleRateRec;
   cosine = cos(omega);
   sine   = sin(omega);
   alpha  = sine/(2*q);

   a0  = 1+alpha;
   *a1 = -2*cosine;
   *a2 = 1-alpha;

   *b0 = 1-alpha;
   *b1 = -2*cosine;
   *b2 = 1+alpha;
  }
  break;

  //calculate peaking-coefficients:
  case RBJ_PEAK:
  {
   A      = pow(10, (0.025*gain) );
   omega  = 2*PI*freq*sampleRateRec;
   cosine = cos(omega);
   sine   = sin(omega);
   alpha  = sine/(2*q);

   a0  = 1+alpha/A;
   *a1 = -2*cosine;
   *a2 = 1-alpha/A;

   *b0 = 1+alpha*A;
   *b1 = -2*cosine;
   *b2 = 1-alpha*A;
  }
  break;

  //calculate low-shelf-coefficients:
  case RBJ_LOW_SHELV:
  {
   A      = pow(10, (0.025*gain) );
   omega  = 2*PI*freq*sampleRateRec;
   sine   = sin(omega);
   cosine = cos(omega);
   beta   = sqrt(A)/q;

   a0  =         (A+1) + (A-1)*cosine + beta*sine;
   *a1 =  -2*  ( (A-1) + (A+1)*cosine              );
   *a2 =         (A+1) + (A-1)*cosine - beta*sine;

   *b0 =   A*  ( (A+1) - (A-1)*cosine + beta*sine  );
   *b1 = 2*A*  ( (A-1) - (A+1)*cosine              );
   *b2 =   A*  ( (A+1) - (A-1)*cosine - beta*sine  );  
  }
  break;

  //calculate high-shelf-coefficients:
  case RBJ_HIGH_SHELV:
  {
   A      = pow(10, (0.025*gain) );
   omega  = 2*PI*freq*sampleRateRec;
   sine   = sin(omega);
   cosine = cos(omega);
   beta   = sqrt(A)/q;

   a0  =          (A+1) - (A-1)*cosine + beta*sine;
   *a1 =    2*  ( (A-1) - (A+1)*cosine              );
   *a2 =          (A+1) - (A-1)*cosine - beta*sine;

   *b0 =    A*  ( (A+1) + (A-1)*cosine + beta*sine  );
   *b1 = -2*A*  ( (A-1) + (A+1)*cosine              );
   *b2 =    A*  ( (A+1) + (A-1)*cosine - beta*sine  );
  }
  break;

  // if mode is not valid, return silence-coefficients by default:
  default:
  {
   *b0 = 0.0;
   *b1 = 0.0;
   *b2 = 0.0;

   a0  = 1.0;
   *a1 = 0.0;
   *a2 = 0.0; 
  }

 } //end switch(mode)

 //scale coefficients:
 a0Rec  = 1.0 / a0;
 *a1    = -(*a1) * a0Rec;
 *a2    = -(*a2) * a0Rec;
 *b0   *= a0Rec;
 *b1   *= a0Rec;
 *b2   *= a0Rec;

}

#endif // BiquadDesigner_h
