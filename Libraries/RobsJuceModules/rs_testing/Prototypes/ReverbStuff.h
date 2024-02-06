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

/** This still doesn't work */

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
    t1.resize(newMaxNumStages);
    t2.resize(newMaxNumStages);
    t3.resize(newMaxNumStages);
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

  /** !!!NOT YET TESTED!!! ...and in very early stages of construction */
  inline TSig getSample(TSig x)
  {
    t3[0] = x;

    for(int i = 0; i < numStages; i++)
    {
      const TPar c = allpassCoeffs[i];           // For convenience.
      t1[i] = delayLines[i].readOutput();        // Read vM = v[n-M] from the delayline
      t2[i] = t3[i] - c * t1[i];                 // Compute v[n] = x[n] - c * v[n-M].
      //t3[i] = c * t2[i] + t1[i];                 // Compute y[n] = c * v[n] + v[n-M].
    }

    for(int i = numStages-1; i >= 0; i--)
    {
      const TPar c = allpassCoeffs[i];           // For convenience.
      //t2[i] = t3[i] - c * t1[i];                 // Compute v[n] = x[n] - c * v[n-M].
      t3[i] = c * t2[i] + t1[i];                 // Compute y[n] = c * v[n] + v[n-M].
      delayLines[i].writeInputAndUpdate(t2[i]);  // Write v[n] into the delayline.
    }

    //return t3[numStages];
    return t3[0];
    // Try to get rid of at least one of the temp arrays, i.e. t3. delayLines[i].readOutput can 
    // also be called in the 2nd loop (in t3[i] = ...) and should return the exact same output 
    // there (I think) so we don't need to store it in t1[i]

    /*
    // original allpass delay:    
    const TPar c = allpassCoeff;         // For convenience.
    TSig vM = delayLine.readOutput();    // Read vM = v[n-M] from the delayline.
    TSig v  = x - c * vM;                // Compute v[n] = x[n] - c * v[n-M].
    delayLine.writeInputAndUpdate(v);    // Write v[n] into the delayline.
    return c * v + vM;                   // Return y[n] = c * v[n] + v[n-M].

    // 1 Level nested:
    const TPar c = allpassCoeff;
    TSig vM = nestedAllpass.getSample(delayLine.readOutput());  // Read vM = innerAllpass(v[n-M])
    TSig v  = x - c * vM;
    delayLine.writeInputAndUpdate(v);
    return c * v + vM;
    */
  }


  /*
  // Throw-awy code for testing - implement 2 stages directly
  inline TSig getSample2(TSig x)
  {
    TSig vM[2], v[2], y[3];
    TPar* c = &allpassCoeffs[0];

    y[ 0] = x;

    vM[0] = delayLines[0].readOutput();    // Read vM = v[n-M] from the delayline.
    v[ 0] = y[0] - c[0] * vM[0];           // Compute v[n] = x[n] - c * v[n-M].
    y[ 1] = c[0] * v[0] + vM[0];           // y1[n] = c * v[n] + v[n-M].


    vM[1] = delayLines[1].readOutput();
    v[ 1] = y[1] - c[1] * vM[1];
    y[ 2] = c[1] * v[1] + vM[1]; 



    delayLine.writeInputAndUpdate(v[0]);   // Write v[n] into the delayline.



    return y[0];
  }
  */



  // Try to implement the lattice form here:
  // https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html
  // https://ccrma.stanford.edu/~jos/pasp/Nested_Allpass_Filters.html
  inline TSig getSample2(TSig x)
  {
    TSig y0 = x;                                                    // y0 = x = input
    TSig y1 = y0 - allpassCoeffs[0] * delayLines[0].readOutput();   // y1 = a
    TSig y2 = y1 - allpassCoeffs[1] * delayLines[1].readOutput();   // y2 = b
    TSig y3 = delayLines[1].readOutput() + allpassCoeffs[1] * y2;   // y3 = c
    TSig y4 = delayLines[0].readOutput() + allpassCoeffs[0] * y1;   // y4 = y = output
    delayLines[0].writeInputAndUpdate(y3);
    delayLines[1].writeInputAndUpdate(y2);
    return y4;

    /*
    TSig a = x - allpassCoeffs[0] * delayLines[0].readOutput();
    TSig b = a - allpassCoeffs[1] * delayLines[1].readOutput();
    TSig c = delayLines[1].readOutput() + allpassCoeffs[1] * b;
    TSig y = delayLines[0].readOutput() + allpassCoeffs[0] * a;
    delayLines[0].writeInputAndUpdate(c);
    delayLines[1].writeInputAndUpdate(b);
    return y;
    */

  }
  // OK - this seems to work


  inline TSig getSample3(TSig x)
  {
    /*
    TSig a = x - allpassCoeffs[0] * delayLines[0].readOutput();
    TSig b = a - allpassCoeffs[1] * delayLines[1].readOutput();
    TSig c = a - allpassCoeffs[1] * delayLines[1].readOutput();
    */




    /*
    TSig c = delayLines[1].readOutput() + allpassCoeffs[1] * b;
    TSig y = delayLines[0].readOutput() + allpassCoeffs[0] * a;
    delayLines[0].writeInputAndUpdate(c);
    delayLines[1].writeInputAndUpdate(b);
    return y;
    */
  }


  void reset()
  {
    for(size_t i = 0; i < delayLines.size(); i++)
    {
      delayLines[i].reset();
      t1[i] = TSig(0);        // May not be needed for the DSP but is cleaner
      t2[i] = TSig(0);        // Dito
      t3[i] = TSig(0);        // Dito
    }
    // Maybe call rsSetToZero(delayLines); rsSetToZero(t1); ...
  }


protected:

  std::vector<RAPT::rsBasicDelayLine<TSig>> delayLines;
  std::vector<TPar> allpassCoeffs;
  //std::vector<TSig> tmp;
  std::vector<TSig> t1, t2, t3;  // temp buffers - maybe rename to vM, v, y
  int numStages = 0;


  //TPar allpassCoeff = TPar(0);
  //RAPT::rsBasicDelayLine<TSig> delayLine;

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