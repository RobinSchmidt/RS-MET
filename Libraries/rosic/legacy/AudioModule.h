#ifndef AudioModule_h
#define AudioModule_h

#include <stdlib.h>
#include "MoreMath.h"
using namespace MoreMath;

/**

This is the base class for audio modules. It defines only the common interface
and some basic data-types (via include). It is not meant to be instantiated
all by itself. However, it's not made as an abstract class (i.e. class, 
containing purely virtual methods) because some subclasses may want to
override the getSample(void) function and don't care about the
getSample(double in) function while for other subclasses it is the other way
around. That's why for both methods a standard ("do nothing") implementation
is provided here in this class.

*/

class AudioModule
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          AudioModule();  ///< Constructor.
 virtual ~AudioModule();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 virtual void setSampleRate(double newSampleRate);
 ///< Set a new sample-rate.

 //---------------------------------------------------------------------------
 // audio processing:

 virtual INLINE double getSample()          { return 0.0; }
 /**< Calculates one sample at a time. This function should be overriden by
      generator-type audio modules (modules that create sound without an 
      input signal). */

 virtual INLINE double getSample(double in)  { return in; }
 /**< Calculates one sample at a time. This function should be overriden by
      effect-type audio modules (modules that modify the input signal). */

 virtual INLINE void getSampleFrameStereo(double* inL, 
                                          double* inR,
                                          double* outL,
                                          double* outR)                                    
 { return; }
 /**< Calculates one stereo sample frame at a time. */

 virtual INLINE void getSampleFrame(double* in, 
                                    double* out,
                                    int     numInputChannels = 2,
                                    int     numOutputChannels = 2)
 { return; }
 /**< Calculates one sample frame at a time.
      This should be used for multichannel modules */

 //===========================================================================

protected:

 doubleA sampleRate; ///< The sample-rate, aligned at 64 bit boundaries.

 INLINE double getAntiDenormAtDc();
 INLINE double getAntiDenormAtNyquist();
 INLINE double getAntiDenormAtDcAndNyquist();
 INLINE double getAntiDenormSaw(int period, double amplitude);
 INLINE double getAntiDenormNoise();
};

INLINE double AudioModule::getAntiDenormAtDc()
{
 return 1.0e-035;
}
INLINE double AudioModule::getAntiDenormAtNyquist()
{
 static int callCounter = 0;

 // let this varibale alternate between 0 and 1 between succesive 
 // function-calls:
 callCounter += 1;
 callCounter %= 2;

 return (callCounter - 0.5) * 1.0e-035;
}
INLINE double AudioModule::getAntiDenormAtDcAndNyquist()
{
 static int callCounter = 0;

 // let this varibale alternate between 0 and 1 between succesive 
 // function-calls:
 callCounter += 1;
 callCounter %= 2;

 return callCounter * 1.0e-035;
}
INLINE double AudioModule::getAntiDenormSaw(int period, double amplitude)
{
 static intA callCounter = 0;

 // let this varibale alternate between 0 and 1 between succesive 
 // function-calls:
 callCounter += 1;
 callCounter %= period;

 return callCounter * amplitude;
}
INLINE double AudioModule::getAntiDenormNoise()
{
 static doubleA randomNumber;

 randomNumber = (rand()/RAND_MAX)-0.5;

 return randomNumber * 1.0e-035;
}

#endif // AudioModule_h
