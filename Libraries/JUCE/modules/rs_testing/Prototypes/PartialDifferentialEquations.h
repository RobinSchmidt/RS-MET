#pragma once


/** A class to generate sounds based on the heat equation. We consider a one-dimensional medium 
(like a thin rod) that has some initial heat distribution which evolves over time according to the
heat equation. The sound generation is just based on reading out this distribution and interpreting 
it as a cycle of a waveform. At each time-step of the PDE solver, we will get another (simpler) 
shape of the cycle than we had the previous cycle - so time-steps of the PDE solver occur after 
each cycle. The time evolution of the heat equation is such that the distribution becomes simpler 
and simpler over time and eventually settles to a uniform distribution at the mean value of the
intial distribution (which we may conveniently choose to be zero). */

template<class T>
class rsHeatEquation1D
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the maximum length of the cycle in samples. The cycle length is equal to the number of 
  rod-elements in our heat-equation simulation. */
  void setMaxCycleLength(int newLength);

  /** Sets up an (initial) heat distribution for the rod. */
  void setHeatDistribution(T* newDistribution, int rodLength, 
    T targetMean = T(0), T targetVariance = T(1));


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  T getSample()
  {
    return T(0);  // preliminary
  }

protected:

  std::vector<T> rod;

};