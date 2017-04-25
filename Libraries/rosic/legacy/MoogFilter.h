#ifndef MoogFilter_h
#define MoogFilter_h

#include "AudioModule.h"
#include "Elliptic8thBandFilter.h"
#include "AntiAliasFilter.h"

/** 

This class makes use of four RC filter stages to model a moog-like ladder
structure consisting of four first order filters in series. The 303-filter
has basically the same structure with the only difference that one of the
first order units has its cutoff at a higher frequency than the others.
The resonance is implemented by inverting the output, multiplying it with
a factor between 0-4 and feeding it back to the input. This results in a
resonant peak, because the signal at the filter output is phase-shifted by
180° and reduced in amplitude by -12 dB at the cutoff frequency. By
multiplication with a factor of 4 and inverting, the output signal at the
cutoff-frequency looks exactly like the input signal at this frequency - so
the interference between the two is constructive. However - in the digital
domain this is not the whole truth: because there must be a unit delay in
the feedback-path the phase shift is not exactly 180° at the cutoff frequency.
This is especially true for high cutoff-frequencies. And thats the reason why
the resonance tends to degrade at high cutoff-frequencies. This effect can be
overcome by effectively running the filter at a higher samplerate such that the
unit delay becomes a shorter period of time. For this reason i have included
an overSampling parameter. Oversampling of 10 works well @ 44100 Hz
samplerate. For feedback factors of 4 or higher the whole filter structure
becomes unstable - that is: it would produce infinite output amplitudes. 
To prevent this, a clipper has been put into the feedback-loop. This clipper
also makes the whole filter a nonlinear system.

*/

class MoogFilter : public AudioModule
{

public:

 //---------------------------------------------------------------------------
 //construction/destruction:

          MoogFilter();  ///< Constructor.        
 virtual ~MoogFilter();  ///< Destructor.

 //---------------------------------------------------------------------------
 //parameter settings:

 virtual void setSampleRate(double newSampleRate);
 ///< Overrides the setSampleRate() method of the AudioModule base class.

 virtual void setOverSamp(int newOverSamp);
 ///< Sets the oversampling factor for the filter. DEPRECATED!!

 virtual void setMode(int newMode);       
 /**< Chooses the mode of the filter. 0:bypass, 1:Moog Lowpass, 
      2:TB-303 Lowpass */

 virtual void setOutputTap(int newOutputTap);  
 /**< Decide after which filter stage the output is taken (range (1-4)). 
      The feedback-loop, however, is always around all four stages. */

	virtual void setDrive(double newDrive);      
 /**< Gain factor for the filter input. The filter input together with the
      feedback signal is clipped at a level of 2.0 before passing through the
      4 filter stages. This is mainly to avoid infinite output amplitudes
      in a self-oscillating state. */

 // parameters that are supposed to be updated at sample-rate 
 // (inlined access-functions):

 virtual INLINE void setCutoff(double newCutoff);
 ///< Sets the cutoff frequency of the 4 filter stages.

 virtual INLINE void setFeedback(double newFeedback);   
 /**< Sets the amount of inverted feedback from the output to the input for
      the resonant peak. */

 virtual INLINE void calcCoeffs();                  
 /**< Causes the filter to re-calculate the coeffiecients - is not called 
      automatically by setCutoff(). ...why not?... */

 //---------------------------------------------------------------------------
 //audio processing:

 virtual INLINE double getSample(double in);
 ///< Calculates one output sample at a time.

 //---------------------------------------------------------------------------
 //others:

 void   resetBuffers();

 //===========================================================================

protected:

 // embedded objects:
 Elliptic8thBandFilter interpolationFilter, 
                       decimationFilter;
 AntiAliasFilter antiAliasFlt;


 //parameter variables:
 double teeBeeFactor;

 static const long maxOverSampling = 64; 
 /**< The maximum oversampling factor. The filter operates on a higher
      sample-rate than the outlying system for closer approximation of a
      delay free feedback-loop and to minmize aliasing artifacts caused by
      the embedded nonlinearity in the feedback loop. */

 doubleA sampleRateRec;         //reciprocal of the sampleRate

 intA mode;  
 doubleA cutoff;
 doubleA feedback;
	doubleA drive;          // drive as amplitude value

 intA outputTap;      // the number of the filter stage after which
                        // the output is taken

 doubleA previousIn;     // previous input sample which is needed for
                        // interpolation (->oversampling); the filter
                        // interpolates with a straight line between the
                        // previous and the current input sample  

 doubleA previousOut;    // previous output sample which is fed back to 
                        // the input

 intA overSampling;   // amount of oversampling (low oversampling rates
                        // cause the resonance to degrade on high cutoff
                        // frequencies, because the phase-shift departs
                        // from 180 degrees at high frequencies)

 doubleA overSamplingRec; // reciprocal

 doubleA buffer[maxOverSampling]; // holds the straight line between the
                                 // previous and the current input sample

 //doubleA b0, b1, a1;              //RC-Filter coefficients

 doubleA b0[4];          // coefficients of the 4 1st order stages
 doubleA b1[4];
 doubleA a1[4];

 //buffering:
 doubleA x_1[4];   // the x[n-1] samples for the 4 filter stages
 doubleA y_1[4];   // the y[n-1] samples for the 4 filter stages

 //the outputs of the 4 filter stages (taps):
 doubleA tapOut[4];

	//protected methods:
 INLINE double calcOverSampledOutSamp(double in);
};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE void MoogFilter::setCutoff(double newCutoff)
{
 //restrict cutoff frequency to the range between 20 and 20000 Hz:
 if( newCutoff <= 20.0 )
  cutoff = 20.0;
 else if( newCutoff >= 20000.0 )
  cutoff = 20000.0;
 else cutoff = newCutoff;
}

INLINE void MoogFilter::setFeedback(double newFeedback)
{
 if( newFeedback >= 0 )
  feedback = newFeedback;
}

INLINE void MoogFilter::calcCoeffs()
{
 static doubleA x;
 static long   i; //loop index

 //intermediate variable for calculation (x is the amount of decay
 //between adjacent samples):
 x = exp( -2 * PI * cutoff*sampleRateRec*overSamplingRec);

 //calculate 1st order low-pass coefficients:
 b0[0] = 1-x;
 b1[0] = 0.0;
 a1[0] = x;

 switch(mode)
 {

  case 1:   //moog lowpass-filter (all rc-filters have the same cutoff frequency)
  {
   for(i=1; i<4; i++)
   {
    b0[i] = b0[0];
    b1[i] = b1[0];
    a1[i] = a1[0];
   }
  }break;

  case 2:   //303 lowpass-filter (one rc-filter has a different cutoff frequency)
  {
   for(i=1; i<3; i++)
   {
    b0[i] = b0[0];
    b1[i] = b1[0];
    a1[i] = a1[0];
   }
   x     = exp( -2 * PI * teeBeeFactor*cutoff*sampleRateRec*overSamplingRec);
   b0[3] = 1-x;
   b1[3] = 0.0;
   a1[3] = x;
  }break;

 }//end of switch
}

INLINE double MoogFilter::getSample(double in)
{
 static doubleA out, intSlope;
 static doubleA tmp;
 static long   i;
/*
 // interpolate with the interpolationFilter, apply filter to each of the 
 // interpolated samples and apply decimationFilter:
 tmp = 8.0*interpolationFilter.getSample(in);
 tmp = calcOverSampledOutSamp(tmp);
 tmp = decimationFilter.getSample(tmp);
 for(i=1; i<=7; i++)
 {
  tmp = interpolationFilter.getSample(0.0);
  tmp = calcOverSampledOutSamp(tmp);
  tmp = decimationFilter.getSample(tmp);
 }
 out = tmp;
*/

 // fill the buffer with the linear interpolated segment between the previous and
 // the current input sample:
 // calculate the interpolation slope:
 intSlope = (in-previousIn)*overSamplingRec;

 // fill the buffer with the oversampled segment:
 for(i=0; i<overSampling; i++)
  buffer[i] = previousIn + intSlope*i;

 // filter the buffer:
 for(i=0; i<overSampling; i++)
  out = calcOverSampledOutSamp(buffer[i]);


 // store the current input sample for the next call:
 previousIn  = in;
 //prevOutSamp = outSamp; //is already done in calcOverSampledOutSamp()

 return out;
}

INLINE double MoogFilter::calcOverSampledOutSamp(double in)
{
 static doubleA tmp1, tmp2; 
 static intA k;        //loop counter for the 4 stages
 static doubleA antiAliasedOut;

 //filter with distortion:
	//apply negative feedback and drive:
	tmp1 = (drive*in - feedback*previousOut);

	//clip the sample:
	CLIP(tmp1, 2.0);  //2.0 is the clipping threshold

 //test:
 //tmp1 = 2.0*tanh(0.5*tmp1);

 // apply the 4 RC-Filter stages:
 tapOut[0] = b0[0]*tmp1 + b1[0]*x_1[0] + a1[0]*y_1[0];
 y_1[0]    = tapOut[0];
 x_1[0]    = tmp1;

 tapOut[1] = b0[1]*tapOut[0] + b1[1]*x_1[1] + a1[1]*y_1[1];
 y_1[1]    = tapOut[1];
 x_1[1]    = tapOut[0];

 tapOut[2] = b0[2]*tapOut[1] + b1[2]*x_1[2] + a1[2]*y_1[2];
 y_1[2]    = tapOut[2];
 x_1[2]    = tapOut[1];

 tapOut[3] = b0[3]*tapOut[2] + b1[3]*x_1[3] + a1[3]*y_1[3];
 y_1[3]    = tapOut[3];
 x_1[3]    = tapOut[2];

 previousOut = tapOut[3];  //used for feedback in the next call

 /*
 //a more general implementation, using a loop (can easily be extended
 //to an arbitrary number of taps):
 tapOut[0] = b0[0]*tmp1 + b1[0]*x_1[0] + a1[0]*y_1[0];
 y_1[0]    = tapOut[0];
 x_1[0]    = tmp1;
 for(k=1; k<4; k++)
 {
  tapOut[k] = b0[k]*tapOut[k-1] + b1[k]*x_1[k] + a1[k]*y_1[k];

  //buffering:
  y_1[k] = tapOut[k];
  x_1[k] = tapOut[k-1];
 }
 prevOutSamp = tapOut[3];  //used for feedback in the next call
 */

 //return tapOut[outputTap];
 return antiAliasFlt.getSample(tapOut[outputTap]);
 //return decimationFilter.getSample(tapOut[outputTap]);
}

#endif // MoogFilter_h
