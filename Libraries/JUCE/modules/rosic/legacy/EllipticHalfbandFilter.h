#ifndef EllipticHalfbandFilter_h
#define EllipticHalfbandFilter_h

#include "BiquadCascade.h"

/**

This class is derived form the BiquadCascade-class. It implements a 
biquad-cascade consisting of 6 biquad-sections with fixed coeeficients. The 
coefficients were obtained with MatLabs fdatool with the following 
design-specifications:

Filter-Type:     Elliptic
Order:           Minimum Order
Frequency-Units: Normalized to 0...1
wpass:           0.45
wstop:           0.5
Magnitude-Units: dB
Apass:           0.1
Astop:           96

it thus has a passband-frequency of 0.45 * Nyquist-Frequency, a 
stoppand-frequency of 0.5 * Nyquist-Frequency, a magnitude response with
ripples of 0.1 dB in the passband and a stopband attenuation of (at least)
96 dB

*/

class EllipticHalfbandFilter : public BiquadCascade
{
public:

 //---------------------------------------------------------------------------
 // construction/destruction:

 EllipticHalfbandFilter();  ///< Constructor.     
 ~EllipticHalfbandFilter();  ///< Destructor.

};

#endif // EllipticHalfbandFilter_h
