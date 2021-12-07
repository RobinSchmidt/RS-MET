#ifndef RAPT_DAMPEDSINE_H
#define RAPT_DAMPEDSINE_H

template<class T>
class rsDampedSine
{

public:

protected:

  T a, d, w, p; 
  /**< Parameters for amplitude, decay, frequency and phase for the function:
    f(t) = a * exp(-d*t) * sin(w*t + p)
  The amplitude a is a raw multiplier. The decay d is 1/tau where tau is the time constant, i.e. 
  the time it takes to decay to 1/e. The frequency variable w stands for "omega" and is actually 
  the radian frequency, i.e. w = 2*pi*f. The start phase p is in radians. */

};

//=================================================================================================

template<class T>
class rsDampedSineSum
{

public:

  // implement operators +,*. Addition should just concatenate the arrays. Multiplication will 
  // produce for each pair of factors a pair of sines and the sum and difference frequencies with 
  // amplitudes given by the products, decay factors given by the sum, etc....

protected:

  std::vector<rsDampedSine<T>> sines; /**< The damped sinusoidal components of this function. */

};





#endif