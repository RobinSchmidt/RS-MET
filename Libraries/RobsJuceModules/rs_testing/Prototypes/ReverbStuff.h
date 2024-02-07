#ifndef RS_REVERBSTUFF_H
#define RS_REVERBSTUFF_H

// ToDo: 
// -Maybe move those classes that eventually go into the RAPT library to the top and those that are
//  really only prototypes (naive implementations that serve as reference for unit tests) to the 
//  bottom


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

  TPar allpassCoeff = TPar(0);
  RAPT::rsBasicDelayLine<TSig> delayLine;

};

// This class actually seems to be production ready when the documentation is completed


//=================================================================================================

/** Under construction....TBC...

obsolete:
This is supposed to be a naive prototype implementation. We just use a std::vector of 
rsAllpassDelay. This is suboptimal because in such a chain, the output delayline of stage i can at 
the same time serve as input delayline for stage i+1 such that we can get rid of almost half of the
delaylines in an optimized implementation. But we don't do that here. Well - when the 
rsAllpassDelay is implemented so as to use only one delayline (e.g. switch from DF1 to DF2 or TDF2, 
see https://www.dsprelated.com/freebooks/filters/Four_Direct_Forms.html), then this optimization 
here won't be needed anymore.


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

/** Implements a filter structure of nested allpass delays using a lattice ...TBC...


This needs clean up and unit tests */

template<class TSig, class TPar>
class rsAllpassDelayNested
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxNumStages(int newMaxNumStages)
  {
    delayLines.resize(newMaxNumStages);
    allpassCoeffs.resize(newMaxNumStages);
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
    allpassCoeffs[stageIndex] = newCoeff;
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  int getMaxNumStages() const { return (int) allpassCoeffs.size(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  // Maybe move the implementations out of the class like in rsAllpassDelayChain. They have grown 
  // quite big. Maybe do it for the setters, too.

  /** Needs more tests */
  inline TSig getSample(TSig x)
  {
    // Shorthands for convenience:
    int N = numStages;
    TSig* t = &tmp[0];

    // Compute the signals in the upper row of the lattice:
    t[0] = x;
    for(int i = 0; i < N; i++)
      t[i+1] = t[i] - allpassCoeffs[i] * delayLines[i].readOutput();

    // Compute the signals in the lower row of the lattice:
    for(int i = 0; i < N; i++)
      t[N+i+1] = delayLines[N-i-1].readOutput() + allpassCoeffs[N-i-1] * t[N-i];

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
    //   https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html
    //   https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
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
    TSig t0 = x;                                                   // t0[n] = x[n]

    // Upper row of lattice:
    TSig t1 = t0 - allpassCoeffs[0] * delayLines[0].readOutput();  // t1[n] = t0[n] - k1 * t3[n-M1]
    TSig t2 = t1 - allpassCoeffs[1] * delayLines[1].readOutput();  // t2[n] = t1[n] - k2 * t2[n-M2]

    // Lower row of lattice:
    TSig t3 = delayLines[1].readOutput() + allpassCoeffs[1] * t2;  // t3[n] = t2[n-M2] + k2 * t2[n]
    TSig t4 = delayLines[0].readOutput() + allpassCoeffs[0] * t1;  // t4[n] = t3[n-M1] + k1 * t1[n]

    // Delayline updates:
    delayLines[0].writeInputAndUpdate(t3);                         // t3 goes into 1st delayline
    delayLines[1].writeInputAndUpdate(t2);                         // t2 goes into 2nd delayline

    // Output:
    return t4;                                                     // y[n] = t4[n]
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
    TSig t1 = t0 - allpassCoeffs[0] * delayLines[0].readOutput();
    TSig t2 = t1 - allpassCoeffs[1] * delayLines[1].readOutput();
    TSig t3 = t2 - allpassCoeffs[2] * delayLines[2].readOutput();

    // Lower row of lattice:
    TSig t4 = delayLines[2].readOutput() + allpassCoeffs[2] * t3;
    TSig t5 = delayLines[1].readOutput() + allpassCoeffs[1] * t2;
    TSig t6 = delayLines[0].readOutput() + allpassCoeffs[0] * t1;

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
  //  Check, if it's possible to get a way with less temporary variables by overwriting them when
  //  they are not needed anymore.





  void reset()
  {
    for(size_t i = 0; i < delayLines.size(); i++)
      delayLines[i].reset();
    rsSetZero(tmp);            // Not needed for the DSP to work correctly but is cleaner
  }


protected:

  std::vector<RAPT::rsBasicDelayLine<TSig>> delayLines;  // Maybe rename to delays
  std::vector<TPar> allpassCoeffs;                       // Maybe rename to coeffs
  std::vector<TSig> tmp;
  int numStages = 0;

};



//#################################################################################################
// From here, we have implementations that are really only for prototyping and as reference for 
// unit testing because they are very suboptimal and/or awkwardly/naively implemented. The 
// implementations here show more clearly, what is going on though, so their value is mostly 
// educational.



//=================================================================================================

/** An allpass delay that realizes the transfer function and difference equation:

          c +     z^(-M)
  H(z) = ----------------,    y[n] = c * x[n] + x[n-M] - c * y[n-M]
          1 + c * z^(-M)

so it's like a first order allpass filter with coefficient c in which the unit delay was replaced
by a delay line of length M. This is also known as a Schroeder allpass section. Such allpass delays 
can be used as building blocks for reverbs, for example. 

See:
https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html  */


template<class TSig, class TPar>
class rsAllpassDelayNaive         // Maybe rename to rsAllpassDelayDF1
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  rsAllpassDelayNaive() {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */


  void setMaxDelayInSamples(int newMaxDelay);

  void setDelayInSamples(int newDelay);

  void setAllpassCoeff(TPar newCoeff) { allpassCoeff = newCoeff; }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline TSig getSample(TSig in);


  void reset();



protected:

  TPar allpassCoeff = 0.0;

  RAPT::rsBasicDelayLine<TSig> inputDelayLine, outputDelayLine;

};


template<class TSig, class TPar>
void rsAllpassDelayNaive<TSig, TPar>::setMaxDelayInSamples(int newMaxDelay)
{
  inputDelayLine.setMaximumDelayInSamples(newMaxDelay);
  outputDelayLine.setMaximumDelayInSamples(newMaxDelay);
}

template<class TSig, class TPar>
void rsAllpassDelayNaive<TSig, TPar>::setDelayInSamples(int newDelay)
{
  inputDelayLine.setDelayInSamples(newDelay);
  outputDelayLine.setDelayInSamples(newDelay);
}

template<class TSig, class TPar>
TSig rsAllpassDelayNaive<TSig, TPar>::getSample(TSig x)
{
  TSig xM = inputDelayLine.getSample(x);                             // x[n-M]
  TSig yM = outputDelayLine.getSampleSuppressTapIncrements(TSig(0)); // y[n-M]
  TSig y  = allpassCoeff * x + xM - allpassCoeff * yM;               // y[n], our current output
  outputDelayLine.addToInput(y);
  outputDelayLine.incrementTapPointers();
  return y;
  // ToDo: verify that this does the right thing with respect to the order of reading, writing and
  // incrementing the taps of the outputDelayLine. Maybe write a unit test that uses a delay of 
  // M = 1 and compare output to a regular first order allpass filter.
  //
  // We want to realize:
  //
  //          c +     z^(-M)
  //  H(z) = ----------------,    y[n] = c * x[n] + x[n-M] - c * y[n-M]
  //          1 + c * z^(-M)
}

template<class TSig, class TPar>
void rsAllpassDelayNaive<TSig, TPar>::reset()
{
  inputDelayLine.reset();
  outputDelayLine.reset();
}

// ToDo:
// -Build a nested allpass in which the z^(-M) term has been replaced by another allpass filter.





//=================================================================================================

template<class TSig, class TPar>
class rsAllpassDelayNestedL1  // L1 mean 1 level of nesting
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Lifetime */

  //rsAllpassDelayNested() {}


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxDelayInSamples(int nestLevel, int newMaxDelay) 
  { 
    if(nestLevel == 0)
      delayLine.setMaximumDelayInSamples(newMaxDelay);
    else
      nestedAllpass.setMaxDelayInSamples(newMaxDelay);
  }

  void setDelayInSamples(int nestLevel, int newDelay) 
  { 
    if(nestLevel == 0)
      delayLine.setDelayInSamples(newDelay);
    else
      nestedAllpass.setDelayInSamples(newDelay);
  }

  void setAllpassCoeff(int nestLevel, TPar newCoeff) 
  { 
    if(nestLevel == 0)
      allpassCoeff = newCoeff;
    else
      nestedAllpass.setAllpassCoeff(newCoeff);
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline TSig getSample(TSig x)
  {
    const TPar c = allpassCoeff;
    TSig vM = nestedAllpass.getSample(delayLine.readOutput());  // Read vM = innerAllpass(v[n-M])
    TSig v  = x - c * vM;
    delayLine.writeInputAndUpdate(v);
    return c * v + vM;

    // The only difference to the implementation of the non-nested case in rsAllpassDelay is that 
    // here we do:
    //   vM = nestedAllpass.getSample(delayLine.readOutput());
    // instead of:
    //   vM = delayLine.readOutput();
  }

  void reset()
  {
    delayLine.reset();
    nestedAllpass.reset();
  }


protected:

  TPar allpassCoeff = TPar(0);
  RAPT::rsBasicDelayLine<TSig> delayLine;
  rsAllpassDelay<TSig, TPar> nestedAllpass;

};

// ToDo:
// -Implement a twice-nested allpass and thrice nested allpass to establish the pattern for how to
//  to it with direct from filters. Then implement a general sturcture for arbitrary many nested 
//  allpasses using a lattice structure - see here:
//  https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
// -Implement unit tests for the arbitrary nesting implementation that compares it to the direct 
//  implementations of 1,2,3 level nesting



//=================================================================================================


template<class TSig, class TPar>
class rsAllpassDelayNestedL2 // L2 means 2 levels of nesting
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxDelayInSamples(int nestLevel, int newMaxDelay) 
  { 
    if(nestLevel == 0)
      delayLine.setMaximumDelayInSamples(newMaxDelay);
    else
      nestedAllpass.setMaxDelayInSamples(nestLevel-1, newMaxDelay);
  }

  void setDelayInSamples(int nestLevel, int newDelay) 
  { 
    if(nestLevel == 0)
      delayLine.setDelayInSamples(newDelay);
    else
      nestedAllpass.setDelayInSamples(nestLevel-1, newDelay);
  }

  void setAllpassCoeff(int nestLevel, TPar newCoeff) 
  { 
    if(nestLevel == 0)
      allpassCoeff = newCoeff;
    else
      nestedAllpass.setAllpassCoeff(nestLevel-1, newCoeff);
  }

  // The only difference to the 1-level nesting case is that the else-branches of the setters now 
  // call the 2-parameter setters of rsAllpassDelayNested instead of the 1-parameter setters of
  // rsAllpassDelay with the nestLevel parameter reduced by one compared to our function argument.
  // This pattern would continue for higher level nesting.


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline TSig getSample(TSig x)
  {
    const TPar c = allpassCoeff;
    TSig vM = nestedAllpass.getSample(delayLine.readOutput());  
    TSig v  = x - c * vM;
    delayLine.writeInputAndUpdate(v);
    return c * v + vM;

    // The only difference to the implementation of the one-level-nested case is that the 
    // nestedAllpass object is now of a different kind - namely itself a 1-level nested allpass
    // delay rather than a regular allpass delay. So, the code here is actually identical to
    // the one-level nesting code - it just means something different because the nestedAllpass
    // member is a different kind of object here.
  }

  void reset()
  {
    delayLine.reset();
    nestedAllpass.reset();
  }

protected:

  TPar allpassCoeff = TPar(0);
  RAPT::rsBasicDelayLine<TSig> delayLine;

  rsAllpassDelayNestedL1<TSig, TPar> nestedAllpass;
  // The only difference to the 1-level nesting case is that this member is now not the simple
  // rsAllpassDelay but itself the 1-level nested allpass delay

};


//=================================================================================================

template<class TSig, class TPar>
class rsAllpassDelayNestedL3
{

public:


  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxDelayInSamples(int nestLevel, int newMaxDelay) 
  { 
    if(nestLevel == 0)
      delayLine.setMaximumDelayInSamples(newMaxDelay);
    else
      nestedAllpass.setMaxDelayInSamples(nestLevel-1, newMaxDelay);
  }

  void setDelayInSamples(int nestLevel, int newDelay) 
  { 
    if(nestLevel == 0)
      delayLine.setDelayInSamples(newDelay);
    else
      nestedAllpass.setDelayInSamples(nestLevel-1, newDelay);
  }

  void setAllpassCoeff(int nestLevel, TPar newCoeff) 
  { 
    if(nestLevel == 0)
      allpassCoeff = newCoeff;
    else
      nestedAllpass.setAllpassCoeff(nestLevel-1, newCoeff);
  }

  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline TSig getSample(TSig x)
  {
    const TPar c = allpassCoeff;
    TSig vM = nestedAllpass.getSample(delayLine.readOutput());  
    TSig v  = x - c * vM;
    delayLine.writeInputAndUpdate(v);
    return c * v + vM;
  }

  void reset()
  {
    delayLine.reset();
    nestedAllpass.reset();
  }

protected:

  TPar allpassCoeff = TPar(0);
  RAPT::rsBasicDelayLine<TSig> delayLine;

  rsAllpassDelayNestedL2<TSig, TPar> nestedAllpass;
  // The only difference to the 1-level nesting case is that this member is now not the 1-level
  // nested allpass delay but the 2-level nested allpass delay. The code of the setters as well as 
  // the code for getSample is literally just copied and pasted from the 2-level implementation.
  // For implementing even more levels of nesting, we can just copy-and-paste the code for one 
  // level of nesting less an replace the nestedAllpass member with an object of the class with 
  // nesting level one less.
};





#endif