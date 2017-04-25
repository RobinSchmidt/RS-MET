#ifndef rosic_DemoVersionNoiseEmitter_h
#define rosic_DemoVersionNoiseEmitter_h

// standard-library indcludes:
#include <stdlib.h>

// rosic-indcludes:
#include "../modulators/rosic_AmpEnvRc.h"

namespace rosic
{

 /**

 This class adds a burst of noise every once in a while to a signal, to be 
 used to restrict the use of demo-versions.

 */

 class DemoVersionNoiseEmitter
 {

 public:

  //---------------------------------------------------------------------------
  // construction/destruction:

  DemoVersionNoiseEmitter();   ///< Constructor.
  ~DemoVersionNoiseEmitter();  ///< Destructor.

  //---------------------------------------------------------------------------
  // parameter settings:

  void setSampleRate(double newSampleRate);
  ///< sets the sample-rate.

  void setBurstInterval(double newBurstInterval);
  /**< Sets the time interval between the noise-bursts in seconds. **/

  //---------------------------------------------------------------------------
  // audio processing:

  INLINE void addNoiseToSampleFrame(double *inL, double *inR);
  ///< adds the noise one output sample-frame at a time.


  //===========================================================================

 protected:

  // embedded audio-modules:
  AmpEnvRc       ampEnv;

  double sampleRate;
  double burstIntervalInSeconds;
  double burstScaler; // scales the volume of the burst inversely proportional 
                      // to the sample-rate to keep the noise-energy in the 
                      // audible band independent from the sample-rate

  double scale, offset;

  unsigned int burstIntervalInSamples;
  unsigned int sampleCounter;

 };

 //-----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):

 INLINE void DemoVersionNoiseEmitter::addNoiseToSampleFrame(double *inL, 
                                                            double *inR)
 {
  // generate the enveloped noise:
  double tmp;
  tmp   = burstScaler * scale * ((double) rand()) + offset;
  tmp  *= ampEnv.getSample();

  // add it to the input-frame:
  *inL += tmp;
  *inR += tmp;

  // keep track of the elapsed time sice the previous burst:
  sampleCounter++;
  if( sampleCounter >= burstIntervalInSamples ) 
  {
   sampleCounter = 0;
   ampEnv.trigger();
  }
 }

} // end namespace rosic

#endif // rosic_DemoVersionNoiseEmitter_h
