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

  /** Sets the coefficient by which the heat diffusss. The higher the value, the faster the 
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

/** Implements a numerical solution of the 1D wave equation. This may serve as a model for a string
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
    rsAssert(newNumGridPoints >= 3, "needs at least 3 grid points (and even that is degenerate)");
    u.resize(newNumGridPoints);
    u1.resize(newNumGridPoints);
    tmp.resize(newNumGridPoints);
  }

  /** Sets the speed by which waves travel along the string/tube. */
  void setWaveSpeed(T newSpeed)
  {
    waveSpeed = newSpeed;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the number of spatial grid points. */
  int getNumGridPoints() const { return (int) u.size(); }

  int getMaxGridIndex() const { return getNumGridPoints()-1; }
  // "N" in (1), see section 5.2.8

  /** Returns the spatial sampling interval.  */
  T getGridSpacing() const { return T(1) / getMaxGridIndex(); }
  // "h" in (1)

  /** Returns the shortest wavelength that can be represented with the current grid settings. */
  T getMinWaveLength() const { return PI / getGridSpacing(); }
  // wavelength is dentoed by "beta" in (1)

  /** Returns the Courant number for the given time-step and our current seetings of number of 
  grid points and wave speed. This number is important for stability analysis of numeric solution 
  schemes for PDEs. If the Courant number is C, then the stability condition is C <= 1. For the 
  special case of C == 1, the finite difference scheme actually produces an exact solution. This is
  only the case for the 1D wave equation and the basis for digital waveguide models. For Courant 
  numbers less than 1, the numerical scheme will produce "numerical dispersion" - waves will 
  propagate with dispersion even though the underlying continuous PDE leads to non-dispersive 
  waves. */
  T getCourantNumber(T timeStep) const;
  // https://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition

  /** Returns the (normalized radian?) frequency for a given wave-number and time-step. */
  T getOmegaForWaveNumber(T waveNumber, T timeStep) const;
  // implements (1), 6.43 - needs test!

  //void rootsOfCharacteristicEquation(T timeStep, T* r1, T* r2);
  // implements (1), Eq. 6.38

  // getKineticEnergy, getPotentialEnergy, getEnergy - mostly for testing energy 
  // conservation/dissipation properties -> plot energies as function of discrete time

  // we need getters for: lambda: Courant number, beta: wavelength, gamma: (normalized) wave-speed
  // h: grid-spacing, k: time-step
  // ...maybe have members for all these variables



  /** Writes the current state of the string into "state" which is supposed to be "length" long. 
  This "length" is actually supposed to be equal to what was set via setNumGridPoints. */
  void getState(T* state, int length) const
  {
    rsAssert(length == getNumGridPoints(), "array length should match number of grid points");
    RAPT::rsArrayTools::copy(&u[0], state, length);
  }

  // it's a bit annoying that we have to pass the time-step to so many fucntions - maybe we should
  // have it as member variable. but this would disallow non-uniform sampling of the solution in 
  // time - which might be convenient thing to have

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Sets the initial positions and velocities. "length" is supposed to be equal to u.size - but 
  client code should pass it in for verification reasons (that's bad API design -> come up with 
  something better) */
  void setInitialConditions(T* newPositions, T* newVelocities, int length, T timeStep);
  // maybe rename to setInitialConditions

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
                        // proportional to frequency of oscillation of single grid point

  //int numGridPoints;
  //T gridSpacing;
  //T timeStep;

  // i'm not yet sure, what the most convenient parametrization is

};

//=================================================================================================

/** Implements a numerical solution of the 2D wave-equation in cartesian coordinates for a 
rectangular membrane.

References:
(1) Numerical Sound Synthesis (Stefan Bilbao) 

*/

template<class T>
class rsRectangularMembrane
{


public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setGridDimensions(int numPointsX, int numPointsY)
  {
    rsAssert(numPointsX == numPointsY, "separate spacings for x- and y not yet supported");
    u.setShape(numPointsX, numPointsY);
    u1.setShape(numPointsX, numPointsY);
    tmp.setShape(numPointsX, numPointsY);
  }
  // maybe rename to setGridResolution

  void setWaveSpeed(T newSpeed) { waveSpeed = newSpeed; }

  void setTimeStep(T newStep) { timeStep = newStep; }



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */


  /** Returns the Courant number for our current settings of the time-step, grid dimensions and 
  wave speed. As opposed to the 1D case, the stability limit is C <= 1/sqrt(2) and there is no 
  special setting for which the numerical scheme produces an exact solution. However, 
  for C == 1/sqrt(2), numerical dispersion in minimized (verify that). */
  T getCourantNumber() const;


  void getState(rsMatrix<T>& state) const { state = u; }
  // check, if this allocates memory (it shouldn't)
  // or maybe return a const reference to our internal state varaiable


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  void setInitialConditions(rsMatrix<T>& displacements, rsMatrix<T>& velocities)
  {
    u  = displacements;
    u1 = u + timeStep * velocities; // why + not -?
  }
  // maybe rename to setInitialValues

  void updateState();

protected:

  void computeInteriorPoints();  // maybe rename to computeInteriorValues
  void computeBoundaryPoints();  // maybe rename to computeBoundaryValues
  void updateStateMatrices();

  void computeInteriorPointsSimple();
  // implements the simplified scheme (1), Eq. 11.12 for special case lambda = 1/sqrt(2)
  // maybe rename to computeInteriorPointsForSpecialCase

  rsMatrix<T> u, u1, tmp; // just like in the 1D case

  T waveSpeed = T(1);
  T timeStep  = T(1);

  // hx = hy = h is natural for isotropic problems ((1), pg. 292)
  // todo: generalize to use separate spatial spacing variables hx, hy

};

// todo: make class rsCircularMembrane


//=================================================================================================

/** Implements numerical solution of the 3D wave-equation in cartesian coordinates for a 
rectangular room. 

References:
(1) Numerical Sound Synthesis (Stefan Bilbao) 

*/

template<class T>  // maybe use TSig, TPar - TSig is used for the pressure
class rsRectangularRoom
{


public:

  //rsRectangularRoom(int numSamplesX, int numSamplesY, int numSamplesZ);


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setGridDimensions(int numSamplesX, int numSamplesY, int numSamplesZ);
  // alternative names setGrid.. Resolution (bad because resolution may be seen as 
  // numSamplesX/lengthX rather than numSamplesX itself

  void setRoomDimensions(T sizeX, T sizeY, T sizeZ)
  {
    Lx = sizeX;
    Ly = sizeY;
    Lz = sizeZ;
    // updateCoeffs()
  }

  void setTimeStep(T newStep) { timeStep = newStep; }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns a const-reference to the current pressure distribution in the room as a 3D array. */
  const rsMultiArray<T>& getState() const { return u; }


  T getPotentialEnergy() const
  {
    return T(0.5) * RAPT::rsArrayTools::sumOfSquares(u.getDataPointerConst(), u.getSize());
    // is this formula correct? we should probably divide by hx*hy*hz in order to approximate the
    // continuous energy integral by a Riemann sum - or maybe use a trapezoidal approximation, 
    // see (1), Eq. 5.20, 5.23, pg.295 bottom

    // what about other inner products like the one used in (1), section 6.2.5 - they use an inner 
    // product of a centered 1st and 2nd time-difference ...what does the "taking the inner product
    // of scheme (6.34) with delta_t. u" mean? is this the way, we arrive at suitable inner 
    // products? taking an inner product of the scheme-equation with some grid-function-difference?

    // implement 5.20 for the 1D case and experiment with it (plot potential/kinetic/total energies
    // over time) - try to generalize from there to 2D and 3D, 
    // try to define other inner products that are closer to physical intuition and plot those 
    // over time, too

  }

  T getKineticEnergy() const
  {
    return T(0.5) * RAPT::rsArrayTools::sumOfSquares(u_t.getDataPointerConst(), u_t.getSize());
    // is this formula correct? ...should also include hx,hy,hz ...maybe also the timeStep?
  }


  T getSecondDerivativeEnergy() const
  {
    return T(0.5) * RAPT::rsArrayTools::sumOfSquares(u_tt.getDataPointerConst(), u_tt.getSize());
    // ad hoc - i don't know, if this as any physical interpretation

    // is this formula correct?
  }

  /** Returns the instantaneous pressure at grid sample location i,j,k. Useful for recording
  output signals. */
  T getPressureAt(int i, int j, int k) const { return u(i,j,k); }

  T getPressureDerivativeAt(int i, int j, int k) const { return u_t(i,j,k); }

  T getPressureSecondDerivativeAt(int i, int j, int k) const { return u_tt(i,j,k); }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  /** Increases the instantaneous pressure at grid sample location i,j,k by the given amount. 
  Useful injecting input signals. */
  void injectPressureAt(int i, int j, int k, T amount) { u(i,j,k) += amount; }
  // maybe assert i,j,k are in valid range...maybe make a function that takes a floating point 
  // position and spreads/de-interpolates the pressure

  /** Updates the state, i.e. the distribution of pressure in the room. */
  void updateState();

  /** Sets the state and all temporary arrays to zero. */
  void reset();

protected:


  void computeLaplacian3D(const rsMultiArray<T>& gridFunction, rsMultiArray<T>& gridLaplacian);
  // maybe factor out - a 3D Laplacian may be useful in many other applications - it may have to
  // take hx,hy,hz variables as parameters ...or maybe 1/hx^2, 1/hy^2, 1/hz^2


  // void updatePressuerAccelerations();
  // void updatePressureVelocities();  // not to be confused with flow of the air (it's a scalar!)
  // void updatePressures;

  int Nx, Ny, Nz;  // redundant but convenient - maybe get rid later
  T   Lx, Ly, Lz;  // room lengths in the coordinate directions

  T timeStep = 1;  // temporal sampling interval

  rsMultiArray<T> u, u_t, u_tt; // pressure, (temporal) pressure-change and change-of-change

  // maybe switch to "scheme"-variables u, u1 later - u1 is the one-sample-delayed pressure - but 
  // for physical intuition, it's easier to represent the temporal derivatives explicitly
  // or maybe use more generic variables like u,t1,t2 where t1,t2 are temporary variables that may
  // be interpreted either as velocity and acceleration or delayed-by-1, delayed-by-2 pressures 
  // the interpretation is up to the solver scheme

};

// how about using spatial pressure gradients and sound-velocities, i.e. vector-fields instead of
// a (scalar) pressure field?