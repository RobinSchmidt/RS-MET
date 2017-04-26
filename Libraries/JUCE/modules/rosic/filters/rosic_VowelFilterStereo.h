#ifndef rosic_VowelFilterStereo_h
#define rosic_VowelFilterStereo_h

// rosic-indcludes:
//#include "../basics/rosic_Definitions.h"
#include "../math/rosic_ElementaryFunctionsReal.h"

namespace rosic
{

  /**

  This class implements a cascade of two biquad filters i series, each of which realizes a 
  parameteric equalizer stage to produce a formant for a vowel-like sound.

  */

  class VowelFilterStereo
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    VowelFilterStereo();   
    ///< Constructor.     

    ~VowelFilterStereo();  
    ///< Destructor.

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    void setSampleRate(double newSampleRate);
    ///< Sets up the sample-rate.

    void setVowel(double newVowel);
    /**< Sets the normalized 'position' in the range of vowels. It is assumed to be in between 
    0...1. The range goes like u -> o -> a -> e -> i -> u and 0 corresponds to the 'u' on the left
    and 1 corresponds to the 'u' on the right. The position thus, constitutes a cyclic behaviour. 
    */

    void setAmount(double newAmount);
    /**< Sets the strength of the effect by scaling all the vowels gain-parameters (in the 
    dB-domain). A value of 1 corresponds to the normal vowel gain-values as they were obtained by
    the analysis, 0 means no effect at all and values greater than one lead to overpronounced
    vowel-effects. */

    void setShift(double newShift);
    /**< Shifts the whole formant-spectrum up or down by an amount specified in semitones. */

    INLINE void updateFilterCoefficients(); 
    /**< Calculates filter coefficients from filter parameters. Has to be called each time the
    vowel or the amount changes -> it is NOT called automatically by setVowel or setAmount to 
    avoid multiple calculations of the coeffs when more than one of these parameters is changed at
    the same time (for example by an envelope). */

    //---------------------------------------------------------------------------------------------
    // inquiry (get-, is-, etc. functions):

    double getVowel();      
    ///<  Sets the normalized 'position' in the range of vowels.

    double getAmount();     
    ///< Returns the amount of the filter.

    //---------------------------------------------------------------------------------------------
    // audio processing:

    INLINE void getSampleFrameStereo(double *inL, double *inR, double *outL, double *outR);
    ///< Calculates one output stereo sample-frame at a time. */

    //---------------------------------------------------------------------------------------------
    // others:

    void resetBuffers (); 
    /**< Sets the buffers for the previous input and output samples of all biquad
         stages to zero. */

  //===============================================================================================

  protected:

    // direct form coefficients for the first and second biquad stage:
    doubleA a1_s1, a2_s1, b0_s1, b1_s1, b2_s1;
    doubleA a1_s2, a2_s2, b0_s2, b1_s2, b2_s2;

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

    // filter parameters:
    doubleA vowel;
    doubleA amount;
    doubleA shiftFactor;

    doubleA sampleRateRec;  // reciprocal of the sample-rate
    doubleA sampleRate;

    double gains[6];
    /**< Global gainfactors for each vowel. */

    double vowelData[6][2][3];
    /**< This array will be initialized in the constructor with the frequencies (in Hz), 
    bandwidths (in octaves) and levels (in dB) of the 6 stored formants (u,o,a,e,i,u) where the 
    'u' is redundant. Index 1 is the vowel, index 2 the formant, index 3 is for frequency, 
    bandwidth and gain. */

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions
  // which are supposed to be called at audio-rate (they can't be put into
  // the .cpp file):

  INLINE void VowelFilterStereo::updateFilterCoefficients()
  {
    double p;           // position in between two vowels (normalized to 0...1)
    double f1, b1, g1;  // frequency, bandwidth and gain of the first formant
    double f2, b2, g2;  // frequency, bandwidth and gain of the second formant

    // obtain the design parameters for the two formant filters by interpolating linearly between
    // the measurement values:
    if( vowel < 0.2 )
    {
      p  = 5.0 * vowel;

      f1 = (1.0-p) * vowelData[0][0][0]  +  p * vowelData[1][0][0];
      b1 = (1.0-p) * vowelData[0][0][1]  +  p * vowelData[1][0][1];
      g1 = (1.0-p) * vowelData[0][0][2]  +  p * vowelData[1][0][2];

      f2 = (1.0-p) * vowelData[0][1][0]  +  p * vowelData[1][1][0];
      b2 = (1.0-p) * vowelData[0][1][1]  +  p * vowelData[1][1][1];
      g2 = (1.0-p) * vowelData[0][1][2]  +  p * vowelData[1][1][2];
    }
    else if( vowel < 0.4 )
    {
      p  = 5.0 * (vowel-0.2);

      f1 = (1.0-p) * vowelData[1][0][0]  +  p * vowelData[2][0][0];
      b1 = (1.0-p) * vowelData[1][0][1]  +  p * vowelData[2][0][1];
      g1 = (1.0-p) * vowelData[1][0][2]  +  p * vowelData[2][0][2];

      f2 = (1.0-p) * vowelData[1][1][0]  +  p * vowelData[2][1][0];
      b2 = (1.0-p) * vowelData[1][1][1]  +  p * vowelData[2][1][1];
      g2 = (1.0-p) * vowelData[1][1][2]  +  p * vowelData[2][1][2];
    }
    else if( vowel < 0.6 )
    {
      p  = 5.0 * (vowel-0.4);

      f1 = (1.0-p) * vowelData[2][0][0]  +  p * vowelData[3][0][0];
      b1 = (1.0-p) * vowelData[2][0][1]  +  p * vowelData[3][0][1];
      g1 = (1.0-p) * vowelData[2][0][2]  +  p * vowelData[3][0][2];

      f2 = (1.0-p) * vowelData[2][1][0]  +  p * vowelData[3][1][0];
      b2 = (1.0-p) * vowelData[2][1][1]  +  p * vowelData[3][1][1];
      g2 = (1.0-p) * vowelData[2][1][2]  +  p * vowelData[3][1][2];
    }
    else if( vowel < 0.8 )
    {
      p  = 5.0 * (vowel-0.6);

      f1 = (1.0-p) * vowelData[3][0][0]  +  p * vowelData[4][0][0];
      b1 = (1.0-p) * vowelData[3][0][1]  +  p * vowelData[4][0][1];
      g1 = (1.0-p) * vowelData[3][0][2]  +  p * vowelData[4][0][2];

      f2 = (1.0-p) * vowelData[3][1][0]  +  p * vowelData[4][1][0];
      b2 = (1.0-p) * vowelData[3][1][1]  +  p * vowelData[4][1][1];
      g2 = (1.0-p) * vowelData[3][1][2]  +  p * vowelData[4][1][2];
    }
    else //
    {
      p  = 5.0 * (vowel-0.8);

      f1 = (1.0-p) * vowelData[4][0][0]  +  p * vowelData[5][0][0];
      b1 = (1.0-p) * vowelData[4][0][1]  +  p * vowelData[5][0][1];
      g1 = (1.0-p) * vowelData[4][0][2]  +  p * vowelData[5][0][2];

      f2 = (1.0-p) * vowelData[4][1][0]  +  p * vowelData[5][1][0];
      b2 = (1.0-p) * vowelData[4][1][1]  +  p * vowelData[5][1][1];
      g2 = (1.0-p) * vowelData[4][1][2]  +  p * vowelData[5][1][2];
    }

    // apply the formant shift and gain adjustment:
    g1 *= amount;
    g2 *= amount;
    f1 *= shiftFactor;
    f2 *= shiftFactor;


    // calculate the coefficients for the first formant filter:

    double sqrtg, gamma, w, s, c, gamog, denRec;

    // intermediate variables:
    sqrtg  = sqrt(dB2amp(g1));                  // square root of the linear gain factor
    w      = 2.0*PI*f1*sampleRateRec;           // normalized radian center frequency (omega)
    sinCos(w, &s, &c);                          // sine and cosine of omega
    gamma  = sinh(0.5 * LN2 * b1 * w / s) * s;  // bandwidth parameter
    gamog  = gamma / sqrtg;                     // gamma over square root of gain
    denRec = 1.0 / (1.0 + gamog);               // reciprocal of common denominator

    // the coefficients themselves:
    b0_s1 = (1.0+gamma*sqrtg) * denRec;
    b1_s1 = -2.0*c            * denRec;
    b2_s1 = (1.0-gamma*sqrtg) * denRec;
    a1_s1 = -b1_s1;
    a2_s1 = -(1.0-gamog)      * denRec;
    
    // calculate the coefficients for the second formant filter:

    sqrtg  = sqrt(dB2amp(g2));                  // square root of the linear gain factor
    w      = 2.0*PI*f2*sampleRateRec;           // normalized radian center frequency (omega)
    sinCos(w, &s, &c);                          // sine and cosine of omega
    gamma  = sinh(0.5 * LN2 * b2 * w / s) * s;  // bandwidth parameter
    gamog  = gamma / sqrtg;                     // gamma over square root of gain
    denRec = 1.0 / (1.0 + gamog);               // reciprocal of common denominator

    // the coefficients themselves:
    b0_s2 = (1.0+gamma*sqrtg) * denRec;
    b1_s2 = -2.0*c            * denRec;
    b2_s2 = (1.0-gamma*sqrtg) * denRec;
    a1_s2 = -b1_s2;
    a2_s2 = -(1.0-gamog)      * denRec;

    // normalize the gain:
    /*
    //double normalizer  = 1.0 / sqrt(dB2amp(max(g1, g2)));
    double normalizer  = 0.25;
    b0_s1 *= normalizer;
    b1_s1 *= normalizer;
    b2_s1 *= normalizer;
    b0_s2 *= normalizer;
    b1_s2 *= normalizer;
    b2_s2 *= normalizer;
    */

  }

  INLINE void VowelFilterStereo::getSampleFrameStereo(double *inL, double *inR, 
                                                      double *outL, double *outR)
  {
    // we need some temporary variables for the filtering because we cannot take for granted that
    // *inL is distinct from *outL (the same for right channel)
    double tmpL, tmpR;

    // calculate output of the first stage:
    tmpL = b0_s1*(*inL) + b1_s1*xL_s1_d1 + b2_s1*xL_s1_d2 + a1_s1*yL_s1_d1 + a2_s1*yL_s1_d2;
    tmpR = b0_s1*(*inR) + b1_s1*xR_s1_d1 + b2_s1*xR_s1_d2 + a1_s1*yR_s1_d1 + a2_s1*yR_s1_d2;

    // update buffers of the first stage:
    xL_s1_d2 = xL_s1_d1;
    xL_s1_d1 = *inL;
    yL_s1_d2 = yL_s1_d1;
    yL_s1_d1 = tmpL;

    xR_s1_d2 = xR_s1_d1;
    xR_s1_d1 = *inR;
    yR_s1_d2 = yR_s1_d1;
    yR_s1_d1 = tmpR;
 
    // apply the second stage:
    tmpL = b0_s2*tmpL + b1_s2*xL_s2_d1 + b2_s2*xL_s2_d2 + a1_s2*yL_s2_d1 + a2_s2*yL_s2_d2;
    tmpR = b0_s2*tmpR + b1_s2*xR_s2_d1 + b2_s2*xR_s2_d2 + a1_s2*yR_s2_d1 + a2_s2*yR_s2_d2;

    // update buffers of the second stage:
    xL_s2_d2 = xL_s2_d1;
    xL_s2_d1 = yL_s1_d1;  // this contains still the sample which was input to this stage
    yL_s2_d2 = yL_s2_d1;
    yL_s2_d1 = tmpL;
    xR_s2_d2 = xR_s2_d1;
    xR_s2_d1 = yR_s1_d1;  
    yR_s2_d2 = yR_s2_d1;
    yR_s2_d1 = tmpR;

    // apply some nonlinear distortion:
    //tmpL = clip(tmpL, -1.0, 1.0);
    //tmpR = clip(tmpR, -1.0, 1.0);

    // store the result in the output slots:
    *outL = tmpL;
    *outR = tmpR;
  }

} // end namespace rosic

#endif // rosic_VowelFilterStereo_h
