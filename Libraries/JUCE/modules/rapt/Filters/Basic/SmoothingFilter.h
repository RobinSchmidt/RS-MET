#ifndef RAPT_SMOOTHINGFILTER_H_INCLUDED
#define RAPT_SMOOTHINGFILTER_H_INCLUDED


/** A filter for smoothing signals in the time domain. It consists of a series connection of first
order lowpass filters (aka leaky integrators). The user can set the order, i.e. the number of 
lowpasses. 

(todo: The filter will scale the time constants according to the order so as to maintain
comparable transition times, so the transition time will not depend on the order. The transition
doesn't get more and more sluggish, when increasing the order (which would happen without the 
automatic downscaling of the time constant). */

// todo: templatize the class
template<class TSig, class TPar> // signal, parameter types
class rsSmoothingFilter
{

public:

  /** Constructor. */
  rsSmoothingFilter()
  {
    reset();
  }

  /** Sets the time constant (in seconds) which is the time it takes to reach around 
  1 - 1/e = 63% of the target value when starting from zero. You must also pass the samplerate
  at which the smoother should operate here. */
  void setTimeConstantAndSampleRate(TPar timeConstant, TPar sampleRate)
  {
    decay = sampleRate * timeConstant;
    updateCoeff();
  }
  // maybe we should use the half-time, i.e. the time, it takes to reach 0.5
  // setTimeToReachHalf...or something, maybe setHalfLifeTime
  // https://en.wikipedia.org/wiki/Exponential_decay#Half-life

  /** Sets the order of the filter, i.e. the number of first order lowpass stages. */
  void setOrder(int newOrder)
  {
    order = rsMin(newOrder, maxOrder);
    updateCoeff();
  }

  /** Returns a smoothed output sample. */
  inline TSig getSample(TSig in)
  {
    //return y1[0] = in + coeff*(y1[0]-in); // preliminary - 1st order only

    y1[0] = in + coeff*(y1[0]-in);
    for(int i = 1; i < order; i++)
      y1[i] = y1[i-1] + coeff*(y1[i] - y1[i-1]);
    return y1[order-1];
  }

  /** Resets the internal filter state to 0. */
  inline void reset()
  {
    for(int i = 0; i < maxOrder; i++)
      y1[i] = 0;
  }

protected:

  /** Updates our filter coefficient according to the setting of decay and order. */
  void updateCoeff()
  {
    TPar scaledDecay = decay/order; // preliminary - scale according to order
    coeff = exp(-1 / scaledDecay);
  }


  static const int maxOrder = 10;

  // member variables:
  TSig y1[maxOrder]; // y[n-1] of the lowpass stages - maybe use a std::vector - gets rid of maxOrder
  TPar coeff = 0;    // lowpass filter coefficient
  TPar decay = 0;    // normalized decay == timeConstant * sampleRate
  int  order = 1;    // number of lowpass stages

  // maybe we should instead of "decay" maintain "sampleRate" and "timeConstant" variables and 
  // provide functions to set them separately. That's more convenient for the user. It increases 
  // our data a bit, but the convenience may be worth it.

};

#endif