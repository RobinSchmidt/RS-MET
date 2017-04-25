#ifndef EllipticQuarterbandFilter_h
#define EllipticQuarterbandFilter_h

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

class EllipticQuarterbandFilter : public BiquadCascade
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

          EllipticQuarterbandFilter();  ///< Constructor.
 virtual ~EllipticQuarterbandFilter();  ///< Destructor.

};

#endif // EllipticQuarterbandFilter_h
