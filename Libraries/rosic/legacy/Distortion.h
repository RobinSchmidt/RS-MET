#ifndef Distortion_h
#define Distortion_h

#include "AudioModule.h"
#include "AntiAliasFilter.h"
#include "float.h"
#include "OnePoleFilter.h"

/** 

This is a memoryless distortion unit which takes an input sample and maps it
to a corresponding output sample. The mapping can be specified for the negative
and the positive half-wave seperately to allow asymmetric distortion 
(symmetric distortion results only in odd harmonics whereas asymmetric
distortion produces also even harmonics). There are three parameters for each
of the both half-wave mappings: threshold specifies the amplitude at which the
distortion begins, saturation is the input amplitude which is mapped to the
clamping level - between threshold and saturation there is a linear mapping
which has in general a different slope than below the threshold. When setting
saturation and clamping to equal values the distortion produces a simple
clipping. There is also an option of rectification of the signal. Because of
the assymetric distortion and rectification capabilities, this unit can also
cause dc-offsets. But these are removed by a first order high-pass filter with
a cutoff frequency of 20 Hz.

*/

class Distortion : public AudioModule
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          Distortion();  ///< Constructor.
 virtual ~Distortion();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

	virtual void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual void setOverSampling(int newOverSampling); 
 ///< Sets an oversampling factor for anti-aliasing.

 virtual void setMode(int newMode);         
 /**< Sets the mapping between "thresh" and "clamp". */

 virtual void setDryWet(double newDryWet);       
 ///< Sets the ratio between dry and wet signal.

 virtual void setDrive(double newDrive);        
 ///< Sets the input gain or "drive" (in dB).

 virtual void setDcOffset(double newDcOffset);
 ///< Adds a dc-offset to the input-signal before distorting.

 virtual void setLowThresh(double newLowThresh);    
 /**< Threshold for the lower half-wave at which the nonlinear part of the 
      curve begins (in dB). */

 virtual void setHighThresh(double newHighThresh);   
 ///< The same for the upper half-wave.

 virtual void setLowSat(double newLowSat);       
 /**< Input level at which the output has reached the clamp level for
      lower half-wave (in dB). */

 virtual void setHighSat(double newHighSat);
 ///< The same for the upper half-wave.

 virtual void setLowClamp(double newLowClamp);
 /**< Level to which the output is forced, when input amplitude reaches
      or exceeds the saturation level. */

 virtual void setHighClamp(double newHighClamp);    
 ///< The same for the upper half-wave.

 virtual void setRectify(double newRectify);      
 /**< Amount of rectification (0.0: no rectification,  
      0.5: half-wave rectification, 1.0: full-wave rectification */

 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE double getSample(double in);
 /**< Calculates a single output sample. */

 virtual INLINE double getSampleAntiAliased(double in);
 /**< Calculates a single output sample while using oversampling by a factor
      adjusted with setOverSampling() and using an AntiAliasFilter before 
      decimating */

 //===========================================================================

protected:

 //embedded objects:
	OnePoleFilter        dcBlocker;
 AntiAliasFilter antiAliasFlt;

	//parameter variables:
 intA    overSampling;
 doubleA overSamplingRec;
 intA    mode;
 doubleA dryVolume, wetVolume;
 doubleA drive;               
 doubleA dcOffset;
 doubleA lowThresh, highThresh;
 doubleA lowClamp, highClamp;   
 doubleA lowSat, highSat;
 doubleA rectify;
 doubleA lowSlope, highSlope;    //slope with which the output is related to
                                 //the input when thresh<=input<sat

 doubleA lowAlpha, lowBeta,      //coeffs in the tanh-function for mode 2
         highAlpha, highBeta;

 bool    lowSlopeInf, highSlopeInf; //when thresh==sat the slope is infinite
                                     //- this case must be treated seperately

 doubleA buffer[64];                //buffer for oversampled signal

 doubleA previousIn;                    //stores previous input sample for 
                                     //interpolation
};

//----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

//calculate a distorted sample:
INLINE double Distortion::getSample(double in)
{
 static doubleA out, temp;

 //scale input sample by drive:
 temp = drive*in;

 //add an dc-offset:
 temp += dcOffset;

 switch(mode)
 {
  case 1:   //linear mapping between "thresh" and "clamp"
  {
   //calculate output for positive inputs:
   if(temp>=0)
   {
    if(highSlopeInf==true)
    {
     if(temp<=highThresh)
      out = temp;
     else //sample is above threshhold
      out = highClamp;
    }
    else  //highSlope is finite
    {
     if(temp<=highThresh)
      out = temp;
     else if(temp<=highSat)
      out = highThresh + highSlope*(temp-highThresh);
     else
      out = highClamp;
    }
   }
   //calculate output for negative inputs:
   else
   {
    temp = -temp;  //invert the sample and do the same procedure as above:

    if(lowSlopeInf==true)
    {
     if(temp<=lowThresh)
      out = temp;
     else //sample is above threshhold
      out = lowClamp;
    }
    else  //lowSlope is finite
    {
     //if(InSamp<=lowThresh)    //<-- falsch! hat aber fetten Sound
     if(temp<=lowThresh)    //<-- richtig
      out = temp;
     else if(temp<=lowSat)
      out = lowThresh + lowSlope*(temp-lowThresh);
     else
      out = lowClamp;
    }

    out = -out; //undo the inversion

    //rectify:
    out -= rectify*out;
   }
  }//end of case 1 (linear mode)
  break;

  case 2: //tanh-mapping between "mThresh" and "mClamp"
  {
   if(temp>=0)    //positive half-wave
   {
    if(temp<=highThresh)
     out = temp;
    else
     out = highThresh + highAlpha*tanh(highBeta*(temp-highThresh));
   }
   else              //negative half-wave
   {
    temp = -temp;
    if(temp<=lowThresh)
     out = temp;
    else
     out = lowThresh + lowAlpha*tanh(lowBeta*(temp-lowThresh));
    //undo the inversion:
    out = -out;
    //rectify:
    out -= rectify*out;
   }
  }
 }//end of switch(mode)

 //mix dry/wet (is done before the dc-blocker, to block dc-components which
 //come from the osc-section (due to waveform-assymetry) too):
 out = dryVolume*in + wetVolume*out;

 //assymetric distortion can produce a dc-offset - this will be removed now:
 //out = dcBlocker.getSample(out); //commented out for test

 return out;
}

//----------------------------------------------------------------------------
//calculate one sample by first linearly interpolating the input signal, then 
//distorting this oversampled signal and applying the anti-aliasing filter:
INLINE double Distortion::getSampleAntiAliased(double in)
{
 static intA  i;
 static doubleA  tmp, 
                out, 
                slope;  //slope for the linear intepolator

 //calculate the interpolation slope:
 slope = (in-previousIn)*overSamplingRec;

 //fill the buffer with the oversampled segment:
 for(i=0; i<overSampling; i++)
  buffer[i] = previousIn + slope*i;

 //distort the buffer and apply anti-alias filter:
 for(i=0; i<overSampling; i++)
 {
  out = getSample(buffer[i]);
  out = antiAliasFlt.getSample(out);
 }

 //store the input sample for next call:
 previousIn = in;

 return out;
}

#endif // Distortion_h
