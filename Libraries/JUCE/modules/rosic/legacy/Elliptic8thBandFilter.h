#ifndef Elliptic8thBandFilter_h
#define Elliptic8thBandFilter_h

#include "BiquadCascade.h"

/**

This class is derived form the BiquadCascade-class. It implements a 
biquad-cascade consisting of 6 biquad-sections with fixed coeeficients. The 
coefficients were obtained with MatLabs fdatool with the following
design-specifications:

Filter-Type:     Elliptic
Order:           12
Frequency-Units: Normalized to 0...1
wpass:           0.9*0.25
Magnitude-Units: dB
Apass:           0.1
Astop:           96

*/

class Elliptic8thBandFilter : public BiquadCascade
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          Elliptic8thBandFilter();  ///< Constructor.
 virtual ~Elliptic8thBandFilter();  ///< Destructor.

};

#endif // Elliptic8thBandFilter_h
