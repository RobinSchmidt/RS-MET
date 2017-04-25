#ifndef rosic_WhiteToPinkFilter_h
#define rosic_WhiteToPinkFilter_h

// rosic-indcludes:
#include "../basics/GlobalDefinitions.h"

namespace rosic
{

 /**

 This is a filter-class which approximates a 3 dB/oct lowpass characteristic. As
 such, it turns white noise into pink noise.

 */

 class WhiteToPinkFilter
 {

 public:

  //---------------------------------------------------------------------------
  // construction/destruction:

	 WhiteToPinkFilter();   ///< Constructor.
	 ~WhiteToPinkFilter();  ///< Destructor.

  //---------------------------------------------------------------------------
  // audio processing:

  INLINE double getSample(double in);

  //---------------------------------------------------------------------------
  // others:
  void reset();

  //===========================================================================

 protected:

  doubleA b[7]; // temporary variables

};

 //-----------------------------------------------------------------------------
 // from here: definitions of the functions to be inlined, i.e. all functions
 // which are supposed to be called at audio-rate (they can't be put into
 // the .cpp file):

 INLINE double WhiteToPinkFilter::getSample(double in)
 {
  static doubleA white, pink;

  white = in;

  b[0] = 0.99886 * b[0] + white * 0.0555179;
  b[1] = 0.99332 * b[1] + white * 0.0750759;
  b[2] = 0.96900 * b[2] + white * 0.1538520;
  b[3] = 0.86650 * b[3] + white * 0.3104856;
  b[4] = 0.55000 * b[4] + white * 0.5329522;
  b[5] = -0.7616 * b[5] - white * 0.0168980;
  pink = (b[0] + b[1]) + (b[2] + b[3]) + (b[4] + b[5]) + (b[6] + white*0.5362);
  b[6] = white * 0.115926;

  return 0.25*pink;
 }

} // end namespace rosic

#endif // rosic_WhiteToPinkFilter_h
