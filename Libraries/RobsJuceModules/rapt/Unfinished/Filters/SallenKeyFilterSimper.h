#ifndef RAPT_SALLENKEYFILTERSIMPER_H
#define RAPT_SALLENKEYFILTERSIMPER_H



template<class T>                   // ToDo: have TSig, TPar template parameters
class rsSallenKeyFilterSimper
{

public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up the filter coefficients ...TBC... */
  inline void setup(T omega, T reso) { setup1(omega, reso); }


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Computes one sample at a time. */
  inline T getSample(T in) { return getSample1(in); }

  /** Resets the internal state. */
  void reset() { ic1eq = ic2eq = 0; }


protected:

  // Setup and process functions for the two different algorithms:
  void setup1(T omega, T reso);
  T getSample1(T in);

  void setup2(T omega, T reso);
  T getSample2(T in);

  // The public setup() and getSample() functions must delegate to setup1()/getSample1() or 
  // setup2()/getSample2() respectively. These are the two variants of the algorithm. The 2nd 
  // algo differs from the first in the computation of a4 (in setup) and of v2 (in getSample). 
  // Algo 2 has less total operations but is also less parallelizable because v1 is used in the 
  // compuation of v2. That reduces the number of operations at the price of introducing a serial 
  // dependency. So, which one should be used should probably be decided by a benchmark. But note 
  // also that for algo 2, there's an additional division in setup. So, when the filter is being
  // modulated at sample rate, the more costly coeff calculation of algo 2 should also be taken 
  // into account. Maybe one or the other algorithm is also more amenable to introducing 
  // nonlinearities? I don't know about that.


  // State:
  T ic1eq = 0;  // Maybe rename to i1
  T ic2eq = 0;

  // Coeffs:
  T a0 = 0, a1 = 0, a2 = 0, a3 = 0, a4 = 0, a5 = 0;
  T k  = 0;

};


template<class T>
void rsSallenKeyFilterSimper<T>::setup1(T omega, T reso)
{
  T g  = tan(0.5*omega);
  T g1 = 1+g;

  k  = 2*reso;
  a0 = 1 / (g1*g1 - g*k);
  a1 = k*a0;
  a2 = g1*a0;
  a3 = g*a2;
  a4 = g*a0;
  a5 = g*a4;
}

template<class T>
T rsSallenKeyFilterSimper<T>::getSample1(T v0)
{
  // Compute node voltages:
  T v1 = a1*ic2eq + a2*ic1eq + a3*v0;
  T v2 = a2*ic2eq + a4*ic1eq + a5*v0;

  // Update state (compute capacitor currents, I guess?):
  ic1eq = 2*(v1 - k*v2) - ic1eq;
  ic2eq = 2*(v2       ) - ic2eq;

  // Return v2 as the lowpass output:
  return v2;

  // See page 3 ("Final Algorithm") here:
  // https://cytomic.com/files/dsp/SkfLinearTrapOptimised2.pdf
}

template<class T>
void rsSallenKeyFilterSimper<T>::setup2(T omega, T reso)
{
  T g  = tan(0.5*omega);
  T g1 = 1+g;

  k  = 2*reso;
  a0 = 1 / (g1*g1 - g*k);
  a1 = k*a0;
  a2 = g1*a0;
  a3 = g*a2;
  a4 = 1 / g1;                             // That's the only difference to setup1()
  a5 = g*a4;
}

template<class T>
T rsSallenKeyFilterSimper<T>::getSample2(T v0)
{
  T v1 = a1*ic2eq + a2*ic1eq + a3*v0;
  T v2 = a4*ic2eq + a5*v1;                // That's the only difference to getSample1()
  ic1eq = 2*(v1 - k*v2) - ic1eq;
  ic2eq = 2*(v2       ) - ic2eq;
  return v2;
}




#endif