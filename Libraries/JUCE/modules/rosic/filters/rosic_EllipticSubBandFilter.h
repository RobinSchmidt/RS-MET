#ifndef rosic_EllipticSubBandFilter_h
#define rosic_EllipticSubBandFilter_h

namespace rosic
{

/** This class implements a 12th order elliptic lowpass filter which is intended to be used for 
situations in which a very steep filter is desirable such as anti-aliasing and anti-imaging in 
resampling applications. */

class EllipticSubBandFilter : public rsBiquadCascade
{

public:

  /** Constructor. Calculates the poles and zeros for the analog unit cutoff prototype filter and 
  stores them in the member arrays prototypePoles and prototypeZeros. */
  EllipticSubBandFilter();

  /** Sets the subdivision factor, for example 2 for a halfband filter (which passes everything 
  below half the Nyquist frequency and stops everything above) or 4 for a quarterband filter. */
  void setSubDivision(double newSubDivision);

protected:

  double  subDivision;
  Complex prototypePoles[6];
  Complex prototypeZeros[6];

};

}

#endif
