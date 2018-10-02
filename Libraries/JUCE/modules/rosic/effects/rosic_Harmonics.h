#ifndef rosic_Harmonics_h
#define rosic_Harmonics_h

//// rosic-indcludes:
//#include "../filters/rosic_EllipticSubBandFilterDirectForm.h"
//#include "../filters/rosic_LowpassHighpass.h"

namespace rosic
{

  /**

  This class adds harmonics by means of a parallel connection of polynomial waveshapers, each with 
  its own, appropriately chosen, subband filter in front of it so as to avoid aliasing.

  */

  class Harmonics
  {

  public:

    //---------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    Harmonics();

    /** Destructor */
    ~Harmonics();

    //---------------------------------------------------------------------------------------------
    // parameter settings (set-functions):

    /** Sets the samplerate. */
    void setSampleRate(double newSampleRate);

    /** Sets the input signal gain / drive in dB. */
    void setDrive(double newDrive) { driveFactor = RAPT::rsDbToAmp(newDrive); }

    /** Sets the ratio between the dry and wet signal. */
    void setDryWetRatio(double newRatio) { wet = newRatio; dry = 1.0-wet; }

    /** Sets the cutoff frequency of the highpass filter at the input stage. */
    void setInputHighpassCutoff(double newCutoff)
    { inFilterL.setHighpassCutoff(newCutoff); inFilterR.setHighpassCutoff(newCutoff); }

    /** Sets the cutoff frequency of the lowpass filter at the input stage. */
    void setInputLowpassCutoff(double newCutoff)
    { inFilterL.setLowpassCutoff(newCutoff); inFilterR.setLowpassCutoff(newCutoff); }

    /** Sets the cutoff frequency of the highpass filter at the output stage. */
    void setOutputHighpassCutoff(double newCutoff)
    { outFilterL.setHighpassCutoff(newCutoff); outFilterR.setHighpassCutoff(newCutoff); }

    /** Sets the cutoff frequency of the lowpass filter at the output stage. */
    void setOutputLowpassCutoff(double newCutoff)
    { outFilterL.setLowpassCutoff(newCutoff); outFilterR.setLowpassCutoff(newCutoff); }

    /** Sets a gain factor for one of the harmonics. */
    void setHarmonicGainFactor(int harmonicNumber, double newGain);

    //---------------------------------------------------------------------------------------------
    // audio processing:

    /** Calculates one output stereo sample-frame at a time. */
    INLINE void getSampleFrameStereo(double *inOutL, double *inOutR);

    //---------------------------------------------------------------------------------------------
    // others:

    /** Resets the internal states of the subband filters. */
    void reset();

    //=============================================================================================

  protected:

    static const int numShapers = 12;

    double dry, wet;
    double driveFactor;              // global gain factor at the input stage
    double gainFactors[numShapers];  // gain factors for the individual harmonics
    bool   gainIsZero[numShapers];   // array of flags to indicate zero gain
    bool   chebychevMode;            // flag to indicate that chebychev polynomials should be used

    rsEllipticSubBandFilterDirectForm subbandFiltersL[numShapers], subbandFiltersR[numShapers];
    LowpassHighpass inFilterL, inFilterR, outFilterL, outFilterR;

  };

  //-----------------------------------------------------------------------------------------------
  // from here: definitions of the functions to be inlined, i.e. all functions which are supposed
  // to be called at audio-rate (they can't be put into the .cpp file):

  INLINE void Harmonics::getSampleFrameStereo(double *inOutL, double *inOutR)
  {
    double inL = RAPT::rsClip(inFilterL.getSample(driveFactor * (*inOutL)), -1.0, 1.0);
    double inR = RAPT::rsClip(inFilterR.getSample(driveFactor * (*inOutR)), -1.0, 1.0);  

    double tmpL  = 0.0;
    double tmpR  = 0.0;
    double accuL = 0.0;
    double accuR = 0.0;
    if( chebychevMode == true )
    {
      for(int i=0; i<numShapers; i++)
      {
        if( !gainIsZero[i] )
        {
          tmpL   = subbandFiltersL[i].getSample(inL);
          tmpR   = subbandFiltersR[i].getSample(inR);
          tmpL   = RAPT::rsCheby(tmpL, i+2); 
          tmpR   = RAPT::rsCheby(tmpR, i+2);
          accuL += gainFactors[i]*tmpL;
          accuR += gainFactors[i]*tmpR;
        }
      }
    }
    else
    {
      for(int i=0; i<numShapers; i++)
      {
        if( !gainIsZero[i] )
        {
          tmpL   = subbandFiltersL[i].getSample(inL);
          tmpR   = subbandFiltersR[i].getSample(inR);
          tmpL   = integerPower(tmpL, i+2); 
          tmpR   = integerPower(tmpR, i+2);
          accuL += gainFactors[i]*tmpL;
          accuR += gainFactors[i]*tmpR;
        }
      }
    }

    *inOutL  = dry*inL + outFilterL.getSample(wet*accuL); 
    *inOutR  = dry*inR + outFilterR.getSample(wet*accuR); 
  }

} // end namespace rosic

#endif // #ifndef rosic_Harmonics_h
