
void rsSamplerFilter::setup(rsSamplerFilter::Type type, float cutoff, float resonance)
{


  int dummy = 0;
}

void rsSamplerFilter::processFrame(float& L, float& R)
{


  int dummy = 0;
}


/*

rsSamplerFilter:
-Maybe make also a struct for very basic 1-pole/1-zero filter. It can be realized by the biquad 
 (and by the other structures, too), but maybe it's more efficient to do it like that. I expect 
 that a simple 1st order lowpass is a quite common thing to use.
-why actually split the variables into coeffs and state - maybe lump them together completely
 maybe that saves some amount of padding as well (declare the "big" TSign variables first). 
 That also eases the implementation. The filters can have their own getSample methods
 ...actually, the we are around the full circle and could just as well use the existing filters
  in RAPT and make a union of those - yeah - maybe that's a good idea we'll see....

*/