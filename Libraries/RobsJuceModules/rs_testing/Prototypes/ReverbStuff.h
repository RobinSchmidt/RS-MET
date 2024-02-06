#ifndef RS_REVERBSTUFF_H
#define RS_REVERBSTUFF_H


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


//=================================================================================================

/** This needs clean up and unit tests */

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

  /** Needs more tests */
  inline TSig getSample(TSig x)
  {
    // Shorthands for convenience:
    int N = numStages;
    TSig* t = &tmp[0];   // maybe use t instead

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
  alternatively to the general getSample() when there are two allpass stages.  */
  inline TSig getSample2(TSig x)
  {
    // We directly implement the lattice form here:
    // https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html
    // https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
    // but with the unit delays replaced by our delaylines

    RAPT::rsAssert(numStages == 2, "Function supposed a 2-stage configuration");

    TSig t0 = x;

    TSig t1 = t0 - allpassCoeffs[0] * delayLines[0].readOutput();
    TSig t2 = t1 - allpassCoeffs[1] * delayLines[1].readOutput();

    TSig t3 = delayLines[1].readOutput() + allpassCoeffs[1] * t2;
    TSig t4 = delayLines[0].readOutput() + allpassCoeffs[0] * t1;

    delayLines[0].writeInputAndUpdate(t3);
    delayLines[1].writeInputAndUpdate(t2);

    return t4;
  }
  // OK - this seems to work - rename to getSample2Stages. This mwas initially intended for 
  // devloping the algo but maybe it should be kept for optimization purposes


  inline TSig getSample3(TSig x)
  {
    TSig t0 = x;

    TSig t1 = t0 - allpassCoeffs[0] * delayLines[0].readOutput();
    TSig t2 = t1 - allpassCoeffs[1] * delayLines[1].readOutput();
    TSig t3 = t2 - allpassCoeffs[2] * delayLines[2].readOutput();

    TSig t4 = delayLines[2].readOutput() + allpassCoeffs[2] * t3;
    TSig t5 = delayLines[1].readOutput() + allpassCoeffs[1] * t2;
    TSig t6 = delayLines[0].readOutput() + allpassCoeffs[0] * t1;

    delayLines[0].writeInputAndUpdate(t5);
    delayLines[1].writeInputAndUpdate(t4);
    delayLines[2].writeInputAndUpdate(t3);

    return t6;
  }
  // rename to getSample3Stages





  void reset()
  {
    for(size_t i = 0; i < delayLines.size(); i++)
      delayLines[i].reset();
    rsSetZero(tmp);            // Not needed for the DSP to work correctly but is cleaner
  }


protected:

  std::vector<RAPT::rsBasicDelayLine<TSig>> delayLines;
  std::vector<TPar> allpassCoeffs;
  std::vector<TSig> tmp;
  int numStages = 0;

};


//=================================================================================================

/** Under construction....TBC...

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

  //std::vector<rsAllpassDelayNaive<TSig, TPar>> allpassDelays;
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



#endif