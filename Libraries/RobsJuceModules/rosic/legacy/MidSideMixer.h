#ifndef MidSideMixer_h
#define MidSideMixer_h

#include "MoreMath.h"
using namespace MoreMath;

/**

This is a simple Mid/Side mixer which can be used to widen or narrowing the
stereo-image. L and R are converted to M and S channels, individual gain is
applied to those M and S signals and then they are re-converted to L and R.
A global gain is applied afterwards. So it has 3 Parameters: midGain, sideGain
and globalGain.

*/

class MidSideMixer
{

public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 MidSideMixer();  ///< Constructor.
 ~MidSideMixer();  ///< Destructor.

 //---------------------------------------------------------------------------
 // parameter settings:

 void setMidGain(double newMidGain);
 ///< Adjusts the gain for the mid-signal - value is expected in dB. 

 void setSideGain(double newSideGain);
 ///< Adjusts the gain for the side-signal - value is expected in dB.  

 void setGlobalGain(double newGlobalGain);
 ///< Adjusts the global gain for the signal - value is expected in dB.  

 //---------------------------------------------------------------------------
 // audio processing:

 INLINE void getSampleFrameStereo(double* inL,  
                                  double* inR, 
                                  double* outL, 
                                  double* outR);
 /**< Calculates a stereo-ouput frame. */

 //===========================================================================

protected:

 // gain values (as raw factors):
 doubleA midGain;
 doubleA sideGain;
 doubleA globalGain;

};

//-----------------------------------------------------------------------------
//from here: definitions of the functions to be inlined, i.e. all functions
//which are supposed to be called at audio-rate (they can't be put into
//the .cpp file):

INLINE void MidSideMixer::getSampleFrameStereo(double* inL,  
                                               double* inR, 
                                               double* outL, 
                                               double* outR)
{
 static doubleA mid, side;

 mid   = midGain    * SQRT2_INV * ( *inL + *inR );
 side  = sideGain   * SQRT2_INV * ( *inL - *inR ); 
 *outL = globalGain * SQRT2_INV * (  mid + side );
 *outR = globalGain * SQRT2_INV * (  mid - side );
}

#endif // MidSideMixer_h
