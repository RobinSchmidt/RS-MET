#ifndef RAPT_ALLPASSES_H_INCLUDED
#define RAPT_ALLPASSES_H_INCLUDED


// This file is currently under construction. It should be filled with some classes that are 
// currently located in ReverbStuff.h in the Prototypes folder. This includes allpass chains and 
// nested allpass structures that are used as building blocks for reverbs. Maybe we should also 
// drag the code from rosic::FlatZapper in and call it rsAllpassDisperser or something like that. 
// Maybe we should also implement Thiran allpass filters here. They are used for allpass 
// interpolation.


// This file contains various allpass filter structures such as the Schroeder allpass, a series
// connection of them, a nested allpass structure, etc. These are usefula as building blocks for
// reverb algorithms.

//=================================================================================================

/** An allpass delay that realizes the transfer function and difference equation:

          c +     z^(-M)
  H(z) = ----------------,    y[n] = c * x[n] + x[n-M] - c * y[n-M]
          1 + c * z^(-M)

so it's like a first order allpass filter with coefficient c in which the unit delay was replaced
by a delay line of length M. This is also known as a Schroeder allpass section. Such allpass delays 
can be used as building blocks for reverbs, for example. The implementation of the difference 
equation is not done directly as written down which corresponds to a 1st order direct form 1 
structure with the unit delay replaced by an M sample delay. Instead, we use the equivalent 
difference equation:

  v[n] = x[n] - c * v[n-M]
  y[n] = c * v[n] + v[n-M]

which needs only one delayline and corresponds to a (delay canonical) direct form 2 implementation 
structure.

See:

  https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html
  https://valhalladsp.com/2011/01/21/reverbs-diffusion-allpass-delays-and-metallic-artifacts/  

*/

template<class TSig, class TPar>
class rsAllpassDelay
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsAllpassDelay() {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */


  void setMaxDelayInSamples(int newMaxDelay) { delayLine.setMaximumDelayInSamples(newMaxDelay); }

  void setDelayInSamples(int newDelay) { delayLine.setDelayInSamples(newDelay); }

  void setAllpassCoeff(TPar newCoeff) { allpassCoeff = newCoeff; }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline TSig getSample(TSig x)
  {
    const TPar c = allpassCoeff;         // For convenience.
    TSig vM = delayLine.readOutput();    // Read vM = v[n-M] from the delayline.
    TSig v  = x - c * vM;                // Compute v[n] = x[n] - c * v[n-M].
    delayLine.writeInputAndUpdate(v);    // Write v[n] into the delayline.
    return c * v + vM;                   // Return y[n] = c * v[n] + v[n-M].

    // See: https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html
  }


  void reset() { delayLine.reset(); }



protected:

  TPar allpassCoeff = TPar(0);           // Rename to coeff
  RAPT::rsBasicDelayLine<TSig> delayLine;

};


//=================================================================================================

/** This implements a chain (i.e. series connection) of allpass delays. To achieve this effect, you
could just use a std::vector of rsAllpassDelay which are applied one after the other. This class 
here is a convenience class that does this for you. 

See:

  https://ccrma.stanford.edu/~jos/pasp/Schroeder_Allpass_Sections.html
  https://www.dsprelated.com/freebooks/pasp/Schroeder_Allpass_Sections.html

*/

template<class TSig, class TPar>
class rsAllpassDelayChain
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxNumStages(int newMaxNumStages);

  void setNumStages(int newNumStages);

  void setMaxDelayInSamples(int stageIndex, int newMaxDelay);

  void setDelayInSamples(int stageIndex, int newDelay);

  void setAllpassCoeff(int stageIndex, TPar newCoeff);


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  int getMaxNumStages() const { return (int) allpassDelays.size(); }



  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline TSig getSample(TSig in);

  void reset();


protected:

  std::vector<rsAllpassDelay<TSig, TPar>> allpassDelays;
  int numStages = 0;

};


// Implementation:

template<class TSig, class TPar>
void rsAllpassDelayChain<TSig, TPar>::setMaxNumStages(int newMaxNumStages)
{
  allpassDelays.resize(newMaxNumStages);
}

template<class TSig, class TPar>
void rsAllpassDelayChain<TSig, TPar>::setNumStages(int newNumStages)
{
  RAPT::rsAssert(newNumStages <= getMaxNumStages());
  numStages = newNumStages;
}

template<class TSig, class TPar>
void rsAllpassDelayChain<TSig, TPar>::setMaxDelayInSamples(int stageIndex, int newMaxDelay)
{
  RAPT::rsAssert(stageIndex < getMaxNumStages());
  allpassDelays[stageIndex].setMaxDelayInSamples(newMaxDelay);
}

template<class TSig, class TPar>
void rsAllpassDelayChain<TSig, TPar>::setDelayInSamples(int stageIndex, int newDelay)
{
  RAPT::rsAssert(stageIndex < getMaxNumStages());
  allpassDelays[stageIndex].setDelayInSamples(newDelay);
}

template<class TSig, class TPar>
void rsAllpassDelayChain<TSig, TPar>::setAllpassCoeff(int stageIndex, TPar newCoeff)
{
  RAPT::rsAssert(stageIndex < getMaxNumStages());
  allpassDelays[stageIndex].setAllpassCoeff(newCoeff);
}

template<class TSig, class TPar>
TSig rsAllpassDelayChain<TSig, TPar>::getSample(TSig in)
{
  TSig tmp = in;
  for(int i = 0; i < numStages; i++)
    tmp = allpassDelays[i].getSample(tmp);
  return tmp;
}

template<class TSig, class TPar>
void rsAllpassDelayChain<TSig, TPar>::reset()
{
  for(int i = 0; i < getMaxNumStages(); i++)
    allpassDelays[i].reset();
}


//=================================================================================================

/** Implements a filter structure of nested allpass delays using a lattice ...TBC... */

template<class TSig, class TPar>
class rsAllpassDelayNested
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxNumStages(int newMaxNumStages)
  {
    delayLines.resize(newMaxNumStages);
    coeffs.resize(newMaxNumStages);
    tmp.resize(2*newMaxNumStages + 1);
  }

  void setNumStages(int newNumStages)
  {
    RAPT::rsAssert(newNumStages <= getMaxNumStages());
    numStages = newNumStages;
  }

  void setMaxDelayInSamples(int stageIndex, int newMaxDelay)
  {
    RAPT::rsAssert(stageIndex < getMaxNumStages());
    delayLines[stageIndex].setMaximumDelayInSamples(newMaxDelay);
  }

  void setDelayInSamples(int stageIndex, int newDelay)
  {
    RAPT::rsAssert(stageIndex < getMaxNumStages());
    delayLines[stageIndex].setDelayInSamples(newDelay);
  }

  void setAllpassCoeff(int stageIndex, TPar newCoeff)
  {
    RAPT::rsAssert(stageIndex < getMaxNumStages());
    coeffs[stageIndex] = newCoeff;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  int getMaxNumStages() const { return (int) coeffs.size(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  // Maybe move the implementations out of the class like in rsAllpassDelayChain. They have grown 
  // quite big. Maybe do it for the setters, too.




  /** Needs more tests */
  inline TSig getSample(TSig x)
  {
    switch(numStages)
    {
    case 0:  return x;
    case 2:  return getSample2Stages(x);
    case 3:  return getSample3Stages(x);
    default: return getSampleNStages(x);  // Works in general but may be suboptimal for small N.
    }
    // We should also have a special function for 1 stage and maybe one for 4. Maybe the functions 
    // we dispatch to should be protected. For the unit test, we can then make a subclass that 
    // allows acces to the via delegating public functions
  }



  inline TSig getSampleNStages(TSig x)
  {
    // Shorthands for convenience:
    int N = numStages;
    TSig* t = &tmp[0];

    // Compute the signals in the upper row of the lattice:
    t[0] = x;
    for(int i = 0; i < N; i++)
      t[i+1] = t[i] - coeffs[i] * delayLines[i].readOutput();

    // Compute the signals in the lower row of the lattice:
    for(int i = 0; i < N; i++)
      t[N+i+1] = delayLines[N-i-1].readOutput() + coeffs[N-i-1] * t[N-i];

    // Update the content of the delaylines:
    for(int i = 0; i < N; i++)
      delayLines[i].writeInputAndUpdate(t[2*N-i-1]);

    // The final output is in the 2N-th slot of the temp-buffer:
    return t[2*N];
  }



  /** An unrolled (and therefore potentially optimized) getSample function that can be used 
  alternatively to the general getSample() when there are two allpass stages. It was initially 
  intended to figure out the general algorithm for getSample but I think, it may be worth to keep
  for documentation and optimization reasons. */
  inline TSig getSample2Stages(TSig x)
  {
    // We directly implement the lattice form shown here (in Fig 2.32b "Second-order allpass 
    // filter: (a) Nested direct-form II. (b) Consecutive two-multiply lattice sections"):
    //
    //   https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html
    //   https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
    //
    // but with the unit delays replaced by our delaylines, i.e. the left z^(-1) of the outer 
    // filter becomes z^(-M1) and the right z^(-1) of the inner filter becomes z^(-M2) where
    // M1, M2 are the lengths of our delaylines and the k1, k2 there mapa to our allpass 
    // coefficients. To translate the block diagram into formulas, I assigned names like 
    // t0, t1, t2, ... (t for temporary) to the signals after every adder starting at the top-left 
    // and going around the U-shaped loop (or horseshoe or whatever). Doing this, we get the 
    // difference equations:
    //
    //   Init:
    //   t0[n] = x[n]
    //
    //   Upper row of lattice:
    //   t1[n] = t0[n] - k1 * t3[n-M1]
    //   t2[n] = t1[n] - k2 * t2[n-M2]
    //
    //   Lower row of lattice:
    //   t3[n] = t2[n-M2] + k2 * t2[n]
    //   t4[n] = t3[n-M1] + k1 * t1[n]
    //
    //   Output:
    //   y[n] = t4[n]
    //
    // The equations have been written down in a way that anticipates a generalization to an 
    // arbitrary number of stages where the computations of the upper and lower part can be done in
    // loops. As we see, the first delayline contains t3 and the second contains t2.

    RAPT::rsAssert(numStages == 2, "Function supposes a 2 stage configuration");
    // The function is meant to be called as an unrolled/optimized alternative to the general 
    // getSample() function which works for any number of stages. But, of course, it is a valid 
    // alternative only when the user actually has selected a two stage configuration.

    // Init:
    TSig t0 = x;                                               // t0[n] = x[n]

    // Upper row of lattice:
    TSig t1 = t0 - coeffs[0] * delayLines[0].readOutput();     // t1[n] = t0[n] - k1 * t3[n-M1]
    TSig t2 = t1 - coeffs[1] * delayLines[1].readOutput();     // t2[n] = t1[n] - k2 * t2[n-M2]

    // Lower row of lattice:
    TSig t3 = delayLines[1].readOutput() + coeffs[1] * t2;     // t3[n] = t2[n-M2] + k2 * t2[n]
    TSig t4 = delayLines[0].readOutput() + coeffs[0] * t1;     // t4[n] = t3[n-M1] + k1 * t1[n]

    // Delayline updates:
    delayLines[0].writeInputAndUpdate(t3);                     // t3 goes into 1st delayline
    delayLines[1].writeInputAndUpdate(t2);                     // t2 goes into 2nd delayline

    // Output:
    return t4;                                                 // y[n] = t4[n]
  }


  inline TSig getSample3Stages(TSig x)
  {
    // This uses the same strategy as getSample2Stages. I just extended the block diagram of the 
    // 2-stage lattice to a 3rd stage and did the same thing - assigning names t0,t1,t2,... to the 
    // variables after the adders (except t0 which is the input x itself) and then reading off the
    // difference equations from the diagram. The 3-stage case already shows the general pattern 
    // that is implemented in getSample() using the loops.

    RAPT::rsAssert(numStages == 3, "Function supposes a 3 stage configuration");

    // Init:
    TSig t0 = x;

    // Upper row of lattice:
    TSig t1 = t0 - coeffs[0] * delayLines[0].readOutput();
    TSig t2 = t1 - coeffs[1] * delayLines[1].readOutput();
    TSig t3 = t2 - coeffs[2] * delayLines[2].readOutput();

    // Lower row of lattice:
    TSig t4 = delayLines[2].readOutput() + coeffs[2] * t3;
    TSig t5 = delayLines[1].readOutput() + coeffs[1] * t2;
    TSig t6 = delayLines[0].readOutput() + coeffs[0] * t1;

    // Delayline updates:
    delayLines[0].writeInputAndUpdate(t5);
    delayLines[1].writeInputAndUpdate(t4);
    delayLines[2].writeInputAndUpdate(t3);

    // Output:
    return t6;
  }
  // ToDo:
  // -Write a getSample4Stages (and a unit test for it). Write performance test and check, if it's
  //  better to use our tmp array or stack-allocated variables for the temporary signals.
  //  Check, if it's possible to get away with less temporary variables by overwriting them when
  //  they are not needed anymore.





  void reset()
  {
    for(size_t i = 0; i < delayLines.size(); i++)
      delayLines[i].reset();
    rsSetZero(tmp);            // Not needed for the DSP to work correctly but is cleaner
  }


protected:

  std::vector<RAPT::rsBasicDelayLine<TSig>> delayLines;
  std::vector<TPar> coeffs;
  std::vector<TSig> tmp;
  int numStages = 0;

};











#endif