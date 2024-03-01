#pragma once


/** A brickwall lowpass filter class that is made from a chain of 3 filter parts: (1) a lowpass, 
(2) a notch/bandstop, (3) an allpass. The lowpass is responsible for the general lowpass nature of 
the filter. The notch is responsible for reducing the ringing at cutoff frequency by notching some
band around that frequecy out. The allpass is responsible for moving a part of the ringing over to 
the left side of the edge, if you think in terms of the step response or a square wave input. The 
settings of these partial filters have been hand tuned to strike an optimal balance between the
desirable steepness of the filter in the frequency domain and undesirable ringing of the filter in
the time domain. ...TBC...  */

template<class TSig, class TPar>
class rsBrickwallFilter
{

public:

protected:

};