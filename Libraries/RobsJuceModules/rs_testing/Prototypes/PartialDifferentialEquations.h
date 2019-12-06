#pragma once

//=================================================================================================

/** A class to generate sounds based on the heat equation. We consider a one-dimensional medium 
(like a thin rod) that has some initial heat distribution which evolves over time according to the
heat equation. The sound generation is just based on reading out this distribution and interpreting 
it as a cycle of a waveform. At each time-step of the PDE solver, we will get another (simpler) 
shape of the cycle than we had the previous cycle - so time-steps of the PDE solver occur after 
each cycle. The time evolution of the heat equation is such that the distribution becomes simpler 
and simpler over time and eventually settles to a uniform distribution at the mean value of the
intial distribution (which we may conveniently choose to be zero). 

the code follows roughly 3blue1brown's explanations here:
Solving the heat equation | DE3  https://www.youtube.com/watch?v=ToIXSwZ1pJU  */

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
  void setHeatDistribution(T* newDistribution, int rodLength);

  void setRandomHeatDistribution(int seed, int rodLength);

  void setTwoValueDistribution(T highFraction, int rodLength);




  /** Normalizes the current heat distribution in the rod to have the given mean and variance. 
  Note that the default mean-value of zero implies that we may have unphysical negative heat 
  values - but that's okay because we want to generate audio signals and don't care about the 
  physical interpretation. You may may want to call this right after calling setHeatDistribution
  in order to remove any DC in your initial distribution and/or in order to normalize the signal
  loudness. */
  void normalizeHeatDistribution(T targetMean = T(0), T targetVariance = T(1));
  // not yet complete

  /** Sets the coefficient by which the heat diffusses. The higher the value, the faster the 
  initial heat distributions evens out over the length of the rod - hence, the faster the 
  complexity and amplitude of the waveform decays. */
  void setDiffusionCoefficient(T newCoeff) { diffusionCoeff = newCoeff; }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Produces one sample at a time and also updates the rod-element corresponding to the 
  currently produced sample. By spreading the rod-update over time in this way, we avoid spikes in 
  CPU load. */
  T getSample()
  {
    //return T(0);  // preliminary

    T out = rodIn[sampleCount];

    // update of a single rod-element:
    T neighborhood;

    
    // cyclic end-handling (better for audio - avoids buzz):
    if(sampleCount == 0)
      neighborhood = T(0.5) * (rodIn[rodLength-1] + rodIn[1]);
    else if(sampleCount == rodLength-1) 
      neighborhood = T(0.5) * (rodIn[rodLength-2] + rodIn[0]);
    else
      neighborhood = T(0.5) * (rodIn[sampleCount-1] + rodIn[sampleCount+1]);


    /*
    if(sampleCount == 0)
      neighborhood = rodIn[1];
    else if(sampleCount == rodLength-1)      // maybe have rodLength-1 as member
      neighborhood = rodIn[rodLength-2];
    else
      neighborhood = T(0.5) * ( rodIn[sampleCount-1] + rodIn[sampleCount+1] );
    */

    



    T delta = neighborhood - out;
    rodOut[sampleCount] = rodIn[sampleCount] + diffusionCoeff * delta;
    // Maybe factor out and/or have other ways of handling the ends (for example cyclically, 
    // clamped, etc.). Currently, we handle the ends by just leaving out the mean-computation and 
    // just adjusting the end-element in the direction to its single neighbor element


    // update counter and return output:
    sampleCount++;
    if(sampleCount == rodLength) {
      sampleCount = 0;
      rsSwap(rodIn, rodOut);
    }
    return out;
  }


  T getSample2()
  {
    T out = rodIn[sampleCount]; // rodIn is currently used for audio output but as input to the update formula
    sampleCount++;

    if(sampleCount == rodLength) {

      // update the whole rod at once:
      T delta   = rodIn[1] - rodIn[0];
      rodOut[0] = rodIn[0] + diffusionCoeff  * delta;
      for(int i = 1; i < rodLength-1; i++)
      {
        T neighborhood = T(0.5) * ( rodIn[i-1] + rodIn[i+1] );
        delta = neighborhood - rodIn[i];
        rodOut[i] = rodIn[i] + diffusionCoeff * delta;
      }
      delta = rodIn[rodLength-2] - rodIn[rodLength-1];
      rodOut[rodLength-1] = rodIn[rodLength-1] + diffusionCoeff * delta;

      // reset counter and swap input and output rods:
      sampleCount = 0;
      rsSwap(rodIn, rodOut);
    }

    if(sampleCount == 0)
      rsPlotArray(rodIn, rodLength); 
      //rsPlotVectors(rodArray1, rodArray2);


    return out;
  }
  // updtae the rod every cycle - should give the same result (needed for test)

  void reset() {
    sampleCount = 0;
    rodIn  = &rodArray1[0];
    rodOut = &rodArray2[0];
  }

protected:

  std::vector<T> rodArray1, rodArray2;
  // array of heat-values in the rod - we need two arrays to alternate between in the evolution
  // of the heat equation

  int rodLength;       // number of (currently used) rod elements

  int sampleCount = 0;

  T diffusionCoeff = 0;


  T *rodIn, *rodOut;
  // Pointers that alternately point to the beginning of rodArray1[0] and rodArray2[0] - used for
  // computing the new heat-distribution from the old

};

//=================================================================================================

/** Implements a numercial solution of the 1D wave equation. This may serve as a model for a string
or tube.

References:
  (1) Numerical Sound Synthesis (Stefan Bilbao)  */

template<class T>
class rsWaveEquation1D
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  /** Sets the number of grid points, i.e. the number of spatial samples along the unit interval
  [0,1]. If this number is N, the spatial sampling interval h will be 1/(N+1). */
  void setNumGridPoints(int newNumGridPoints)
  {
    rsAssert(newNumGridPoints >= 3, "needs at leat 3 grid points (and even that is degenerate)");
    u.resize(newNumGridPoints);
    u1.resize(newNumGridPoints);
    tmp.resize(newNumGridPoints);
  }

  void setWaveSpeed(T newSpeed)
  {
    waveSpeed = newSpeed;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the number of spatial grid points. */
  int getNumGridPoints() const { return (int) u.size(); }

  /** Writes the current state of the string into "state" which is supposed to be "length" long. 
  This "length" is actually supposed to be equal to what was set via setNumGridPoints. */
  void getState(T* state, int length) const
  {
    rsAssert(length == getNumGridPoints(), "array length should match number of grid points");
    RAPT::rsArray::copy(&u[0], state, length);
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Sets the initial positions and velocities. "length" is supposed to be equal to u.size - but 
  client code should pass it in for verification reasons (that's bad API design -> come up with 
  something better) */
  void initPositionsAndVelocities(T* newPositions, T* newVelocities, int length, T timeStep);

  /** Updates the stae of the PDE solver which consists of the current values for the positions and 
  the values one time-step before. */
  void updateState(T timeStep);
  // todo: maybe allow different interpretations of the state such as positions and velocities - 
  // that would be more physical
  // maybe rename to advanceSolution


protected:

  /** Called from updateState. Computes updated values for interior points and stored them in 
  tmp. */
  void computeInteriorPoints(T timeStep);
  // give it a name that says which scheme is being used


  void computeInteriorPointsSimple();
  // uses Eq. 6.54 which is the scheme 6.34 for the special case lambda = 1. in this special case,
  // the recursion formula becomes much simpler



  /** Called from updateState. Computes updated values for boundary points and stored them in 
  tmp. */
  void computeBoundaryPoints(T timeStep);

  /** Moves "u" into "u1" and then "tmp" into "u". */
  void updateStateArrays();


  std::vector<T> u, u1; // current state, state one sample ago
  std::vector<T> tmp;   // temporary buffer used in state updates

  T waveSpeed = T(1);   // which unit? spatial-steps per time-step? -> figure out!

  //int numGridPoints;
  //T gridSpacing;
  //T timeStep;

  // i'm not yet sure, what the most convenient parametrization is


};