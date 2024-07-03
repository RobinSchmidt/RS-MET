#ifndef RAPT_ELLIPTICSUBBANDFILTER_H_INCLUDED
#define RAPT_ELLIPTICSUBBANDFILTER_H_INCLUDED

/** This class implements a 12th order elliptic lowpass filter which is intended to be used for 
situations in which a very steep filter is desirable such as anti-aliasing and anti-imaging in 
resampling applications. */

template<class TSig, class TPar>
class rsEllipticSubBandFilter : public rsBiquadCascade<TSig, TPar>
{

  typedef std::complex<TPar> Complex; // preliminary

public:

  /** Constructor. Calculates the poles and zeros for the analog unit cutoff prototype filter and 
  stores them in the member arrays prototypePoles and prototypeZeros. */
  rsEllipticSubBandFilter();

  /** Sets the subdivision factor, for example 2 for a halfband filter (which passes everything 
  below half the Nyquist frequency and stops everything above) or 4 for a quarterband filter. */
  void setSubDivision(TPar newSubDivision);

protected:

  TPar subDivision;

  Complex prototypePoles[6];
  Complex prototypeZeros[6];
  // What are these arrays used for? Document that.

};


#endif
