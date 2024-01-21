#ifndef RAPT_PROTOTYPEDESIGNER_H_INCLUDED
#define RAPT_PROTOTYPEDESIGNER_H_INCLUDED

/** This class determines the locations of poles and zeros in the s-plane for a continuous time, 
unit cutoff lowpass or low-shelving prototype filter. It supports Butterworth, Chebychev, inverse 
Chebychev, elliptic, Bessel and Papoulis designs. The low-shelving design is a generalization of a 
unit cutoff lowpass filter and its magnitude response can be seen as the lowpass-response raised to 
a pedestal. The height of this pedestal is called the reference gain (which is zero in the 
lowpass-case). The gain the 'passband' (represented by the member variable G) is either a 
boost (when G > 1) or an attenuation/cut (when G < 1) - for standard lowpass designs, this gain 
would be unity.

References:
  (1) Sophocles J. Orfanidis: Lecture Notes on Elliptic Filter Design
  (2) Sophocles J. Orfanidis: High-Order Elliptical Equalizer Design
  (3) Larry D. Paarmann: Design and Analysis of Analog Filters

\todo
 -test the Gaussian design (compare poles with reference design, like Paarmann's GAUSSDE.m)
 -check the gain calculation, include gain factor
 -solve the degree equation for other values than the filter order
 -check getFilterResponseAt() - seems to be buggy - obsolete?
 -get rid of needsSpecialHighShelvTransform() - rewrite the prototype design code in such a way 
  that all approximations can be treated uniformly later on
 -lot's of stinky code duplication -> eliminate it */

template<class T>
class rsPrototypeDesigner
{

  typedef std::complex<T> Complex; 
  // preliminary - maybe we should switch back to using rsComplex because std::complex only 
  // supports float, double and long double for the type T but we may need it for SIMD types or
  // multi-precision types, too. rsComplex does not have that restriction.

public:

  /** This is an enumeration of the available approximation methods. */
  enum approximationMethods
  {
    BUTTERWORTH = 0,   /**< maximally flat at DC: all derivatives zero, up to filter order */
    CHEBYCHEV,         /**< equiripple in passband, monotonic in stopband */
    INVERSE_CHEBYCHEV, /**< equiripple in stopband, monotonic in passband */
    ELLIPTIC,          /**< equiripple in passband and stopband, maximally steep transition */
    BESSEL,            /**< approximates linear phase */
    PAPOULIS,          /**< maximizes steepness at cutoff (selectivity) under constraint of 
                            monotonicity */
    HALPERN,           /**< minimizes ratio of bandwidths at specified magnitudes (shaping factor)
                            under constraint of monotonicity. similar to Papoulis but less steep 
                            at cutoff and steeper in stopband */
    GAUSSIAN,          /**< smallest ringingTime*bandWidth product, good time response (very little
                            overshoot), isotropic in 2D
                            
                            step response rises faster than Bessel. but wait - do we 
                            normalize the cutoff the same way as in bessel? ..with respect to 
                            the asymptote? -> nope - it's not! it's not comparable  */
    NUM_APPROXIMATION_METHODS
  };
  // re-order: COINCINDENT_POLE, GAUSS, BESSEL, BUTTERWORTH, PAPOULIS <-?-> HALPERN, CHEBY1 <-?-> 
  // CHEBY2, ELLIPTIC ->sorted by desirability of time response vs. frequency response (roughly)
  // or maybe sort by ringing time? or maybe list allpole types first, so we may implement an easy
  // isAllpole() check?
  // however, currently the ordering of the choices in EngineersFilter (in ToolChain) depends on
  // the order being the way it is...

  // ideas: try other polynomials, for example non-reversed Bessel, Laguerre, etc. - if roots occur 
  // in the right half-plane, reflect them, maybe try Power-Series expansions of various 
  // interesting functions as it is done with the Gaussian filter. maybe use transfer functions 
  // based on Pade approximations. Maybe try 1 / (1 + cosh(w)) as target function for the 
  // magnitude-squared response

  // in (3), page 16, it says: "Scanlan introduced filters with poles that fall on an ellipse with 
  // equal frequency spacing and noted the tradeoff between magnitude response characteristics and 
  // time-domain response characteristics as the eccentricity of the ellipse is varied" -> try it
  // ...maybe call them "ellipsoid" filters...because "elliptic" already means something else 
  // -> includes butterworth as special case, when the ellipse becomes a circle - maybe the user 
  // parameter should be the log of the ratio of major and minor axis? how about other geometric
  // definitions of pole/zero locations?

  /** This enumerates the two possible prototype filter characterisitics. */
  enum prototypeModes
  {
    LOWPASS_PROTOTYPE = 1,
    LOWSHELV_PROTOTYPE
  };

  //-----------------------------------------------------------------------------------------------
  /** \name Construction/Destruction */

  /** Constructor. */
  rsPrototypeDesigner();

  /** Destructor. */
  ~rsPrototypeDesigner();

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets ups the order of the filter - that is the number of first order sections. Each section 
  increases the slope by 6 dB/oct (for lowpass-designs). If the order is odd, a real pole-/zero 
  pair will be present, otherwise all poles and zeros will be complex
  conjugate. */
  void setOrder(int newOrder);

  /** Chooses one of the approximation methods as enumerated above. */
  void setApproximationMethod(int newApproximationMethod);

  /** Chooses a lowpass or low-shelving prototype as enumerated above. */
  void setPrototypeMode(int newPrototypeMode);

  /** Sets the ripple in the passband for lowpass designs in decibels. */
  void setPassbandRipple(T newPassbandRipple);

  /** Sets the rejection in the stopband for lowpass designs in decibels. */
  void setStopbandRejection(T newStopbandRejection);

  /** Sets the gain for shelving filters in the passband (in dB). This will be positive for a 
  low-boost and negative for low-cut responses. */
  void setGain(T newGain);

  /** Sets up the reference gain for shelving filters (in dB). This will be usually unity. */
  void setReferenceGain(T newReferenceGain);

  /** Selects the fraction of the maximum gain, which the magnitude response assumes at the 
  passband-frequency. This parameter is relevant only for Chebychev and elliptic shelving filters.
  It should be chosen close to unity (for example 0.95) in order to prevent excessive ripple in the
  passband. */
  void setPassbandGainRatio(T newPassbandGainRatio);

  /** Selects the fraction of the maximum gain, which the magnitude response assumes at the 
  stopband-frequency. This parameter is only relevant for inverse Chebychev and elliptic shelving 
  filters. It should be chosen close to zero (for example 0.05) in order to prevent excessive 
  ripple in the stopband. */
  void setStopbandGainRatio(T newBandwidthGainRatio);

  /** Assigns the poles and zeros such that the resulting filter will just pass the signal 
  through. */
  void makeBypass();

  //-----------------------------------------------------------------------------------------------
  /** \name Static Functions */

  /** Given desired the order "N" of the prototype filter, this function returns the number of 
  required 2nd order sections in "L" and 1st order sections in "r" (either 0 or 1). */
  static void getNumBiquadsAndFirstOrderStages(int N, int &L, int &r);

  /** Solves the elliptic degree equation (k from N, k1). */
  static T ellipdeg(int N, T k_1);
  // move to protected

  /** Solves the elliptic degree equation (k_1 from N, k1). */
  static T ellipdeg1(int N, T k);
  // move to protected

  /** Solves the elliptic degree equation (k_1 from N, k) using nomes. */
  static T ellipdeg2(T N, T k);
  // move to protected

  /** Calculates the order required for a Butterworth filter to fullfill given design 
  specifications. The actual order should be chosen to be the next integer. */
  static T getRequiredButterworthOrder(T passbandFrequency, T passbandRipple, T stopbandFrequency, 
    T stopbandRipple);

  /** Calculates the order required for a Chebychev filter (inverse or not) to fullfill given 
  design specifications. The actual order should be chosen to be the next integer. */
  static T getRequiredChebychevOrder(T passbandFrequency, T passbandRipple, T stopbandFrequency, 
    T stopbandRipple);

  /** Calculates the order required for an elliptic filter to fullfill given design 
  specifications. The actual order should be chosen to be the next integer. */
  static T getRequiredEllipticOrder(T passbandFrequency, T passbandRipple, T stopbandFrequency, 
    T stopbandRipple);

  /** Computes the total area under magnitude-squared response of a Butterworth filter of given 
  order. You may also specify the number of passes, if the filter is to be applied multiple times
  (for example 2x in bidirectional filtering or in Linkwitz/Riley filters, but you may also run the
  filter several times to approximate a Gaussian filter. In such cases, it may make sense to scale 
  the cutoff frequency inversely with the factor returned by this function in order to have the 
  filter output energies comparable. */
  static T butterworthEnergy(int order, int numPasses = 1);

  /** Given the arrays of polynomial coefficients "b" and "a" of a transfer function 
  H(s) = N(s)/D(s), this function returns the polynomial coefficients of the corresponding 
  magnitude-squared function H(s)*H(-s) = (N(s)*N(-s)) / (D(s)*D(-s)). "N" is the order of the 
  filter, so "b", "a" should be of length N+1 and "b2", "a2" of length 2*N+1. */
  static void magSquaredNumAndDen(T* b, T* a, T* b2, T* a2, int N);

  /** Given the two arrays "b2" and "a2" of polynomial coefficients (both of length 2*N+1) for 
  numerator and denominator of an s-domain Nth order lowpass prototype magnitude-squared function 
  and gain a constant k, this function computes the numerator coefficients of the magnitude-squared 
  function of the corresponding low-shelving filter with reference-gain G0 and low-frequency gain G 
  and writes them into bS. */
  static void shelvingMagSqrNumFromLowpassMagSqr(T* b2, T* a2, T k, int N, T G0, T G, T* bShelf);

  /** Given the two arrays "b" and "a" of polynomial coefficients (both of length N+1) for 
  numerator and denominator of an s-domain lowpass prototype transfer function and gain a 
  constant k, this function computes the numerator coefficients of the magnitude-squared function 
  of the corresponding low-shelving filter with reference-gain G0 and low-frequency gain G and 
  writes them into bS. Because the magnitude-squared function has twice the order of the transfer 
  function itself, bS will be of length 2*N+1. The left halfplane roots of bS will be the zeros of 
  the shelving filter. The denominator coefficients (and hence, the poles of the filter) are 
  unaffected by the lowpass-to-lowshelf transform, so you may re-use your "a" array in the 
  low-shelving filter. */
  static void shelvingMagSqrNumeratorFromLowpassTransfer(T* b, T* a, T k, int N, T G0, T G, T* bS);

  /** Scales zeros, poles and gain factor, such that the magnitude response at unit frequency 
  equals "g". */
  static void scaleToMatchGainAtUnity(Complex* z, Complex* p, T* k, Complex* zNew, 
    Complex* pNew, T* kNew, int N, T g);

  /** Returns zeros, poles and gain in "zNew", "pNew", "kNew" of a filter that is inverse to the 
  filter with zeros, poles and gain of "z", "p", "k". \todo maybe move to PoleZeroMapper. zNew, 
  pNew may point to the same arrays as z, p. */
  static void getInverseFilter(Complex* z, Complex* p, T* k, Complex* zNew, Complex* pNew, 
    T* kNew, int N);

  /** Given an array of N+1 cofficients for a polynomial of order N, this function returns the 
  left halfplane roots in r and returns the number of such roots in the return-value. The rest of
  the array "r" is left as is - in most cases, you should assume that it contains garbage.
  \todo maybe move to PolynomialAlgorithms */
  static int getLeftHalfPlaneRoots(T* a, Complex* r, int N);


  /** Generates the coefficients of the N-th order reverse Bessel polynomial that occurs in the 
  denominator of the transfer function of Bessel filters. Note that in contrast to Papoulis, 
  Halpern and Gauss filters, thes are the coefficients of the transfer function itself, not the
  magnitude-squared response coeffs. */
  static void besselDenominator(T* a, int N); 
  // allocates

  /** Computes zeros, poles and gain factor for an analog low-shelving Bessel prototype filter.
  @see getBesselLowpassZerosPolesAndGain */
  static void besselZPK(Complex* z, Complex* p, T* k, int N, T G, T G0);

  /** Generates coefficients of a polynomial of order 2*N for the squared polynomial that occurs
  in the denominator of N-th order Papoulis filters. It's the L^2(w) polynomial in Eq. 8.14 in 
  Paarmann: Design and Analysis of Analog Filters.  */
  static void papoulisPolynomial(T *a, int N);
  // allocates

  /** Constructs the denominator polynomial of the magnitude-squared function for Papoulis filters 
  where "N" is the filter order and "a" is of length 2*N+1. */
  static void papoulisDenominator(T* a, int N);
  // allocates

  /** Generates coefficients of a polynomial of order 2*N for the squared polynomial that occurs
  in the denominator of N-th order Halpern filters. It's the T^2(w) polynomial in Eq. 8.18 in 
  Paarmann: Design and Analysis of Analog Filters.   */
  static void halpernPolynomial(T *a, int N);

  /** Analogous to papoulisDenominator */
  static void halpernDenominator(T *a, int N);

  /** Generates coefficients of a polynomial of order 2*N for the squared polynomial that occurs
  in the denominator of Gaussian filters. It's the polynomial in the denominator of Eq. 8.7 in 
  Paarmann: Design and Analysis of Analog Filters. */
  static void gaussianPolynomial(T *a, int N, T wc);

  static void gaussianDenominator(T *a, int N);


  static void papoulisZPK(Complex* z, Complex* p, T* k, int N, T G, T G0);
  static void halpernZPK( Complex* z, Complex* p, T* k, int N, T G, T G0);
  static void gaussianZPK(Complex* z, Complex* p, T* k, int N, T G, T G0);



  /** Computes poles and zeros for one of the lowpass filter types that is determined by it 
  magnitude-squared denominator polynomial, which can be any of the functions
  papoulisDenominator, halpernDenominator, gaussianDenominator. The matchButterworth parameter
  decides whether or not the Butterworth magnitude response should be matched asymptotically. */
  static void zpkFromMagSquaredCoeffsLP(Complex* z, Complex* p, T* k, int N,
    void (*denominatorCoeffsFunction)(T* a, int N), bool matchButterworth);
  // todo: get rid of the z parameter array - fill it with infinity in the calling function
  // zpkFromMagSquaredCoeffsLS

  /** Similar to zpkFromMagSquaredCoeffsLP, but for low-shelving prototypes. The matchButterworth
  parameter is relevant only for the special case G0 == 0 (a lowpass filter) in which case the 
  design will be delegated to zpkFromMagSquaredCoeffsLP. Otherwise, the poles and zeros will be
  scaled according to having the bandwidth gain halfway between G and G0 (verify this) */
  static void zpkFromMagSquaredCoeffsLS(Complex* z, Complex* p, T* k, int N, T G, T G0,
    void (*denominatorCoeffsFunction)(T* a, int N), bool matchButterworth);

  /** Computes poles and zeros for one of the lowpass filter types that is determined by it 
  transfer-function denominator polynomial, which can currently only be the function 
  besselDenominator (maybe more to come later). The matchButterworth parameter decides whether or 
  not the Butterworth magnitude response should be matched asymptotically (unmatched case needs 
  testing). */
  static void zpkFromTransferCoeffsLP(Complex* z, Complex* p, T* k, int N,
    void (*denominatorCoeffsFunction)(T* a, int N), bool matchButterworth);

  /** Given a "denominatorCoeffsFunction" that generates polynomial coefficients for a lowpass prototype 
  transfer function, this function creates the zeros, poles and gain for the corresponding 
  low-shelving prototype with given shelving-gain G and reference-gain G0. */
  static void zpkFromTransferCoeffsLS(Complex* z, Complex* p, T* k, int N, T G, T G0,
    void (*denominatorCoeffsFunction)(T* a, int N), bool matchButterworth);
    // maybe have an optional numeratorCoeffsFunction (defaulting to a nullptr in which case the 
    // numerator is taken to be 1

  //-------------------------------------------------------

  /** Computes zeros, poles and gain factor for an analog elliptic prototype filter of order "N" 
  with passband gain variation (ripple) "Gp" and maximum stopband amplitude "Gs", and stores them 
  in "z", "p" and "k", respectively. */
  static void getEllipticLowpassZerosPolesAndGain(Complex* z, Complex* p, T* k, int N, 
    T Gp, T Gs);

  //static void getLowpassZerosPolesAndGain(Complex* z, Complex* p, T* k, int N, 
  //  int approximationMethod);

  //static void getBesselLowshelfZerosPolesAndGain(Complex* z, Complex* p, T* k, int N, T G, T G0); 

  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Re-calculates the poles and zeros (if necesarry) an writes them into the respective arrays. 
  For complex conjugate pairs, it will only write one representant of the pair into the array. If 
  N denotes the order as passed to setOrder(), let L be the number of second order sections and r 
  be either one or zero indicating presence or absence of a first order section. Then, the array 
  will be filled with one representant for each pole pair into poles[0...L-1] and the real pole 
  will be written into poles[L+r-1], if present. The same goes for the zeros. */
  void getPolesAndZeros(Complex* poles, Complex* zeros);
    // \todo: rename into getNonRedundantPolesAndZeros - provide also a function getPolesAndZeros 
    // that writes all poles and zeros "as is" into the arrays.

  /** Calculates and returns the complex frequency response at some value of the Laplace transform
  variable s - seems to be buggy. */
  Complex getFilterResponseAt(Complex s);

  /** Returns the normalized (divided by the DC gain - ...really?!) magnitude at the given radian
  frequency w. */
  T getMagnitudeAt(T w);


  /** Finds the radian frequency w inside the interval wLow...wHigh at which the specified 
  magnitude (as raw amplitude value, not in dB) occurs. For the function to work properly, the 
  magnitude response inside the given interval should be monotonic and the specified magnitude 
  should occur somewhere inside the given interval. */
  T findFrequencyWithMagnitude(T magnitude, T wLow, T wHigh);

  /** Returns the number of finite poles. */
  int getNumFinitePoles();

  /** Returns the number of finite zeros. */
  int getNumFiniteZeros();

  /** Returns the number of non-redundant finite zeros (each pair of complex conjugate zeros counts
  only once). */
  int getNumNonRedundantFiniteZeros();

  /** Returns the number of non-redundant finite poles (each pair of complex conjugate poles counts 
  only once). */
  int getNumNonRedundantFinitePoles();

  /** Returns the passbandRipple in dB. */
  T getPassbandRipple() const { return Ap; }

  /** Returns the ratio between the peak gain (in dB) and the ripples inside the boosted/cutted 
  band for shelving modes. */
  T getPassbandGainRatio() const { return Rp; }

  /** Returns the approximation method to be used @see enum approximationMethods. */
  int getApproximationMethod() const { return approximationMethod; }

  /** Returns the mode of the prototype @see enum modes. */
  int getPrototypeMode() const { return prototypeMode; }

  /** Returns the filter order. */
  int getOrder() const { return N; }

  /** Returns true if the currently selected approximation method supports a ripple parameter. */
  bool hasCurrentMethodRippleParameter();
  //{ return (approximationMethod == ELLIPTIC) || (approximationMethod == CHEBYCHEV); }

  /** Returns true if the currently selected approximation method supports a rejection 
  parameter. */
  bool hasCurrentMethodRejectionParameter();
  //{ return (approximationMethod == ELLIPTIC) || (approximationMethod == INVERSE_CHEBYCHEV); }

  /** Returns true when the currently selected approximation method requires a special 
  lowshhelf-to-highshelf s-plane frequency transformation. This is a rather ugly kludge - we should 
  rewrite the prototype design code such that all approximation methods can be treated 
  uniformly. */
  bool needsSpecialHighShelvTransform();


  //-----------------------------------------------------------------------------------------------
  /** \name Introspection */

  /** Checks, if the order variables (numBiquads, numPoles, numZeros) are within their sane 
  range. Used mainly for assertions. */
  bool sanityCheckOrderVariables();

  // ToDo: make similar functions for sanityCheckRippleVariables etc.


  //===============================================================================================

protected:

  /** Calculates the poles and zeros according to the selected approximation method and the global 
  gain factor. */
  void updatePolesAndZeros();

  /** Calculates the poles and zeros for a Butterworth lowpass filter. */
  void makeButterworthLowpass();

  /** Calculates the poles and zeros for a Butterworth low-shelving filter. */
  void makeButterworthLowShelv();

  /** Calculates the poles and zeros for a Chebychev lowpass filter. */
  void makeChebychevLowpass();

  /** Calculates the poles and zeros for a Chebychev low-shelving filter. */
  void makeChebychevLowShelv();

  /** Calculates the poles and zeros for an inverse Chebychev lowpass filter. */
  void makeInverseChebychevLowpass();

  /** Calculates the poles and zeros for an inverse Chebychev low-shelving filter. */
  void makeInverseChebychevLowShelv();

  /** Calculates the poles and zeros for an elliptic lowpass filter. */
  void makeEllipticLowpass();

  /** Calculates the poles and zeros for an elliptic low-shelving filter. */
  void makeEllipticLowShelv();

  /** Given gain G, reference gain G0 and a function that computes the zeros, poles and gain
  (such as papoulisZPK, besselZPK, etc.), this function uses the zpk-function to assign our zeros, 
  poles and gain member variables. */
  void makeLowShelfFromZPK(void (*zpkFunc)(Complex* z, Complex* p, T* k, int N, T G, T G0), 
    T G, T G0);

  /** Given the arrays zTmp and pTmp of poles and zeros, this function picks the non-redundant ones
  and copies them into our z and p members, respectively. The arrays zTmp, pTmp are assumed to 
  contain N zeros and poles, where N is the order of the prototype filter. Complex poles/zeros 
  should occur in complex conjugate pairs - this function will select one representant for each 
  such pair.
  ATTENTION: If there is a real zero/pole present, they are assumed to be in zTmp[0] and pTmp[0] 
  respectively - the caller must ensure this.  */
  //void pickNonRedundantPolesAndZeros(Complex *zTmp, Complex *pTmp);
  // obsolete

  // user parameters:
  int N;                   // prototype filter order: N = 2*L + r 
  int approximationMethod; // selected approximation method 
  int prototypeMode;       // selected mode (lowpass or low-shelv)
  T Ap;                    // passband ripple/attenuation for lowpass designs in dB
  T As;                    // stopband rejection for lowpass designs in dB
  T A;                     // boost/cut gain for low-shelv designs in dB
  T A0;                    // reference gain for low-shelv designs in dB
  T Rp;                    // ripple in the boosted/cutted band as fraction of dB-peak-gain for 
                           // shelvers (elliptic and chebychev)
  T Rs;                    // ripple outside the boosted/cutted band as fraction of dB-peak-gain 
                           // for shelvers (elliptic and inverse chebychev)
  // Maybe rename Ap to passRippleDb, As to stopRippleDb, A to shelfGainDb, A0 to refGainDb, Rp to
  // shelfInnerRipplePercent, Rs to shelfOuterRipplePercent. But these variables are used a lot
  // in complex formulas, so maybe it's better to use shorter names which also aggree with the 
  // math notation used in the papers, so maybe it's fine as it is. Not sure.

  // internal variables:
  int L;                   // number of second order sections (biquads)
  int r;                   // number of first order sections (either zero or one)
  int numFinitePoles;      // number of poles (excluding those at infinity)
  int numFiniteZeros;      // number of zeros (excluding those at infinity). 
  // Maybe rename L to numBiquads and r to numOnePoles/numSinglePoles/numFirstOrderSections/
  // numBilins/numBilinears

  static const int maxBiquads = 12;               // maximum number of biquad sections
  static const int maxOrder   = 2 * maxBiquads;   // maximum filter order
  static const int maxCoeffs  = 2 * maxOrder + 1; // maximum number of polynomial coeffs, * 2 
                                                  // because we need mag-squared polynomials
                                                  // todo: maybe use member for temporary coeff arrays
                                                  // -> functions will need less stack memory
  // The maxBiquads value was formerly 10 (until 2024/01/21) but then I discovered that this is not 
  // enough and in certain cases led to writing into memory locations beyond the allowed ones. 

  // arrays for nonredundant poles and zeros:
  Complex z[maxBiquads];   // zeros
  Complex p[maxBiquads];   // poles
  T k = 1;                 // overall gain factor - not yet used
  // ToDo: document why the array sizes are not 2*maxBiquads. Maybe it's because we only store the
  // non-redundant poles and zeros and not their complex conjugates?

  bool stateIsDirty;   // this flag indicates, whether the poles, zeros and gain need to be 
                       // re-calculated or are still valid from a previous calculation 
                       // (todo: use atomic_bool)


  template<class TSig, class TPar> friend class rsEngineersFilter;
  // We want rsEngineersFilter to be able to access our protected maxBiquads member.

};

#endif
