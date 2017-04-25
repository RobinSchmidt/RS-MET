#ifndef FuncShaper_h
#define FuncShaper_h

#include "AudioModule.h"
#include "TwoPoleBandpass.h"
#include "Upsampler.h"
#include "TabulatedFunction.h"
#include "AntiAliasFilter.h"

/**

This is a waveshaping  distortion unit which lets the user specify the 
mapping between the input and output as a mathematical formula.

*/

class FuncShaper : public AudioModule 
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          FuncShaper();  ///< Constructor.
 virtual ~FuncShaper();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual bool setFunctionString(char* newFunctionString);
 /**< Sets a new function string in the embedded TabulatedFunction object
      and returns true, if the string is valid */

 virtual void setA(double newA);
 /**< Assigns the a-constant in the expression for the distortion curve to
      a numeric value */

 virtual void setB(double newB);
 /**< Assigns the b-constant in the expression for the distortion curve to
      a numeric value */

 virtual void setC(double newC);
 /**< Assigns the c-constant in the expression for the distortion curve to
      a numeric value */

 virtual void setD(double newD);
 /**< Assigns the d-constant in the expression for the distortion curve to
      a numeric value */

 virtual void setDrive(double newDrive);
 ///< Sets the drive (gain) for the input (expected in dB).

 virtual void setInLpfCutoff(double newInLpfCutoff);
 ///< Sets the cutoff frequency of the input lowpass filter (expected in Hz).

 virtual void setInHpfCutoff(double newInHpfCutoff);
 ///< Sets the cutoff frequency of the input highpass filter (expected in Hz).

 virtual void setDcOffset(double newDcOffset);
 ///< Sets the DC offset to be added to the signal (expected as raw value).

 virtual void setOutLpfCutoff(double newOutLpfCutoff);
 ///< Sets the cutoff frequency of the output lowpass filter (expected in Hz).

 virtual void setOutHpfCutoff(double newOutHpfCutoff);
 ///< Sets the cutoff frequency of the output highpass filter (expected in Hz).

 virtual void setOutVol(double newOutVol);
 ///< Sets the global output volume (expected in dB).

 virtual void setDryWet(double newDryWet);
 /**< Sets the ratio between the original and the effect signal (expected 
      in % wet). */

 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE double getSample(double in);
 ///< Calculates one output sample at a time.

 virtual INLINE void getSampleFrameStereo(double* inL, 
                                          double* inR,
                                          double* outL,
                                          double* outR);
 ///< Calculates one streo output frame at a time.

 //---------------------------------------------------------------------------
 // others:

 //---------------------------------------------------------------------------
 // embedded objects (in the order of the signal flow):
 TwoPoleBandpass   inputFilterL, inputFilterR;
 Upsampler         upsamplerL, upsamplerR;
 TabulatedFunction distortionCurve;
 AntiAliasFilter   antiAliasFilterL, antiAliasFilterR;
 TwoPoleBandpass   outputFilterL, outputFilterR;

 //===========================================================================

protected:

 // parameter variables:
 doubleA a, b, c, d;
 doubleA drive, dcOffset;
 doubleA outVol;
 doubleA dryVol, wetVol;

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE double FuncShaper::getSample(double in)
{
 return distortionCurve.getValueLinear(in);
}

INLINE void FuncShaper::getSampleFrameStereo(double* inL, 
                                             double* inR,
                                             double* outL,
                                             double* outR)
{
 static doubleA tmpL, tmpR;
 static intA i;

 // apply the drive, input filter and dc-offset:
 tmpL = drive * inputFilterL.getSample(*inL) + dcOffset;
 tmpR = drive * inputFilterR.getSample(*inR) + dcOffset;


 /* */

 // feed the samples into the upsampling unit:
 upsamplerL.setInputSample(tmpL);
 upsamplerR.setInputSample(tmpR);

 // calculate 4 frames of the oversampled, distorted, and filtered signal
 // (we do oversampling by a factor of 4):
 for(i=1; i<=4; i++)
 {
  tmpL = upsamplerL.getSample();
  tmpR = upsamplerR.getSample();

  tmpL = distortionCurve.getValueLinear(tmpL);
  tmpR = distortionCurve.getValueLinear(tmpR);

  tmpL = antiAliasFilterL.getSample(tmpL);
  tmpR = antiAliasFilterL.getSample(tmpR);
 }

 //*outL = tmpL;
 //*outR = tmpR;


 /*
 // apply the nonlinear distortion curve:
 tmpL = distortionCurve.getValueLinear(tmpL);
 tmpR = distortionCurve.getValueLinear(tmpR);

 *outL = tmpL;
 *outR = tmpR;
 */ 

 // apply the output filter:
 tmpL = outputFilterL.getSample(tmpL);
 tmpR = outputFilterL.getSample(tmpR);

 // mix dry and wet signal:
 *outL = outVol*wetVol*tmpL + dryVol*(*inL);
 *outR = outVol*wetVol*tmpR + dryVol*(*inR);

 return;
}

#endif // FuncShaper_h