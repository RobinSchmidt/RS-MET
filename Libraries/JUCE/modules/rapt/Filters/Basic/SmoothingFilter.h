#ifndef RAPT_SMOOTHINGFILTER_H_INCLUDED
#define RAPT_SMOOTHINGFILTER_H_INCLUDED


/** A filter for smoothing signals in the time domain. It consists of a series connection of first
order lowpass filters (aka leaky integrators). The user can set the order, i.e. the number of 
lowpasses. 

(todo: The filter will scale the time constants according to the order so as to maintain
comparable transition times, so the transition time will not depend on the order. The transition
doesn't get more and more sluggish, when increasing the order (which would happen without the 
automatic downscaling of the time constant). 
...i should use a table that gets allocated and filled when the first object is created


*/

template<class TSig, class TPar> // signal, parameter types
class rsSmoothingFilter
{

public:

  enum shapes
  {
    SIGMOID = 0,
    FAST_ATTACK,

    NUM_SHAPES
  };

  /** Constructor. */
  rsSmoothingFilter();

  /** Sets the time constant (in seconds) which is the time it takes to reach around 
  1 - 1/e = 63% of the target value when starting from zero. You must also pass the samplerate
  at which the smoother should operate here. */
  //void setTimeConstantAndSampleRate(TPar timeConstant, TPar sampleRate);
  // maybe we should use the half-time, i.e. the time, it takes to reach 0.5
  // setTimeToReachHalf, setNumSamplesToReachHalf...or something, maybe setHalfLifeTime
  // https://en.wikipedia.org/wiki/Exponential_decay#Half-life


  /** Sets up the number of samples that it takes for the unit-step response to reach a value
  of 1/2. */
  void setNumSamplesToReachHalf(TPar numSamples);


  /** Sets the order of the filter, i.e. the number of first order lowpass stages. */
  void setOrder(int newOrder);

  /** Chooses the shape of the step-response. See shapes enum. When all first order lowpass stages
  have the same cutoff frequency, higher order filters will approach a more or less symmetrical 
  sigmoid shape. This is sometimes undesirable because it implies a kind of delayed reaction. With 
  FAST_ATTACK, we scale the time-constants of the successive filters which has the effect that the 
  transition is faster initially. */
  void setShape(int newShape);
    // maybe call the parameter asymmetry and get rid of the setShape function - just check
    // if asymmetry == 0 - if so, use the cheaper computations in updateCoeffs

  /** If some shape other than SIGMOID is chosen, this parameter sets the amount of the shape. If 
  zero, it reduces to the sigmoid shape (i.e. all stages have same time constant). If 1, the time
  constants of the stages follow an 1/n rule. The general rule is 1/n^p where p is the parameter
  set here. */
  void setShapeParameter(TPar newParam);

  /** Returns a smoothed output sample. */
  inline TSig getSample(TSig in)
  {
    //return y1[0] = in + coeff*(y1[0]-in); // this would be the 1st order leaky integrator

    y1[0] = in + coeffs[0]*(y1[0]-in);
    for(int i = 1; i < order; i++)
      y1[i] = y1[i-1] + coeffs[i]*(y1[i] - y1[i-1]);
    return y1[order-1];
  }

  /** Resets the internal filter state to 0. */
  void reset();

protected:

  /** Updates our filter coefficients according to the setting of decay and order. */
  void updateCoeffs();

  std::vector<TSig> y1;     // y[n-1] of the lowpass stages
  std::vector<TPar> coeffs; // lowpass filter coefficients
  TPar decay = TPar(0.1);   // normalized decay == timeConstant * sampleRate ...rename to numSamplesToReachHalf
  int  order = 1;           // number of lowpass stages, now redundant with y1.size()

  int shape = 0;       // remove
  TPar shapeParam = 0; // rename to asymmetry

  // maybe we should instead of "decay" maintain "sampleRate" and "timeConstant" variables and 
  // provide functions to set them separately. That's more convenient for the user. It increases 
  // our data a bit, but the convenience may be worth it.


  // variables for the table of the time-constant scalers:
  static const int  maxOrder = 8;
  static const int  numAsyms = 2;
  static const TPar maxAsym;
  static bool tableIsFilled;
  static rsMatrix<TPar> tauScalers; // 2D table
  void createTauScalerTable();


};

#endif