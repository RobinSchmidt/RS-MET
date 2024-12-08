#ifndef RS_REVERBSTUFF_H
#define RS_REVERBSTUFF_H

// ToDo: 
// -Maybe move those classes that eventually go into the RAPT library to the top and those that are
//  really only prototypes (naive implementations that serve as reference for unit tests) to the 
//  bottom


//=================================================================================================

/** This implements a chain (i.e. series connection) of allpass delays. To achieve this effect, you
could just use a std::vector of rsAllpassDelay which are applied one after the other. This class 
here is a convenience class that does this for you. Series connections of allpass filters can be 
used for allpass diffusors, for example, as building blocks of a reverb algorithm.

Hmm...maybe it's actually not such a great idea to provide such a convenience class. If we do this,
we may want to have similar convenience classes for other types of allpass filter chains which 
would look very similar - i.e. a lot code duplication and boilerplate. I actually had this class
already in the RAPT library but backed off again and moved it back into the prototypes for this 
reason. 

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

  void setMaxNumStages(int newMaxNumStages)
  {
    allpassDelays.resize(newMaxNumStages);
  }

  void setNumStages(int newNumStages)
  {
    RAPT::rsAssert(newNumStages <= getMaxNumStages());
    numStages = newNumStages;
  }

  void setMaxDelayInSamples(int stageIndex, int newMaxDelay)
  {
    RAPT::rsAssert(stageIndex < getMaxNumStages());
    allpassDelays[stageIndex].setMaxDelayInSamples(newMaxDelay);
  }

  void setDelayInSamples(int stageIndex, int newDelay)
  {
    RAPT::rsAssert(stageIndex < getMaxNumStages());
    allpassDelays[stageIndex].setDelayInSamples(newDelay);
  }

  void setAllpassCoeff(int stageIndex, TPar newCoeff)
  {
    RAPT::rsAssert(stageIndex < getMaxNumStages());
    allpassDelays[stageIndex].setAllpassCoeff(newCoeff);
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  int getMaxNumStages() const { return (int) allpassDelays.size(); }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline TSig getSample(TSig in)
  {
    TSig tmp = in;
    for(int i = 0; i < numStages; i++)
      tmp = allpassDelays[i].getSample(tmp);
    return tmp;
  }

  void reset()
  {
    for(int i = 0; i < getMaxNumStages(); i++)
      allpassDelays[i].reset();
  }


protected:

  std::vector<rsAllpassDelay<TSig, TPar>> allpassDelays;
  int numStages = 0;

};



//#################################################################################################
//
// From here, we have implementations that are really only for prototyping and as reference for 
// unit testing because they are very suboptimal and/or awkwardly/naively implemented. The 
// implementations here show more clearly, what is going on though, so their value is mostly 
// educational. And they can be used for unit testing purposes for producing reference output 
// signals to test the better implementations against.



//=================================================================================================

/** An allpass delay that realizes the transfer function and difference equation:

          c +     z^(-M)
  H(z) = ----------------,    y[n] = c * x[n] + x[n-M] - c * y[n-M]
          1 + c * z^(-M)

so it's like a first order allpass filter with coefficient c in which the unit delay was replaced
by a delay line of length M. This is also known as a Schroeder allpass section. Such allpass delays 
can be used as building blocks for reverbs, for example. A "non-naive" implementation can be found 
in RAPT::rsAllpassDelay. It uses only half of the delay memory. The point to keep this prototype 
class around is to have an implementation that can be verified to be correct by inspection more
easily.

See:
https://www.dsprelated.com/freebooks/pasp/Allpass_Filters.html  */


template<class TSig, class TPar>
class rsAllpassDelayNaive         // Maybe rename to rsAllpassDelayDF1
{


public:

  void setMaxDelayInSamples(int newMaxDelay)
  {
    inputDelayLine.setMaximumDelayInSamples(newMaxDelay);
    outputDelayLine.setMaximumDelayInSamples(newMaxDelay);
  }

  void setDelayInSamples(int newDelay)
  {
    inputDelayLine.setDelayInSamples(newDelay);
    outputDelayLine.setDelayInSamples(newDelay);
  }

  void setAllpassCoeff(TPar newCoeff) { coeff = newCoeff; }

  inline TSig getSample(TSig x)
  {
    TSig xM = inputDelayLine.getSample(x);                             // x[n-M]
    TSig yM = outputDelayLine.getSampleSuppressTapIncrements(TSig(0)); // y[n-M]
    TSig y  = coeff * x + xM - coeff * yM;                             // y[n], our current output
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


  void reset()
  {
    inputDelayLine.reset();
    outputDelayLine.reset();
  }


protected:

  RAPT::rsBasicDelayLine<TSig> inputDelayLine;
  RAPT::rsBasicDelayLine<TSig> outputDelayLine;
  TPar coeff = 0.0;

};

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

  void setMaxDelayInSamples(int newMaxDelay)
  {
    setMaxDelayInSamples(0, newMaxDelay);
    setMaxDelayInSamples(1, newMaxDelay);
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

  void setMaxDelayInSamples(int newMaxDelay)
  {
    setMaxDelayInSamples(0, newMaxDelay);
    setMaxDelayInSamples(1, newMaxDelay);
    setMaxDelayInSamples(2, newMaxDelay);
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

  void setMaxDelayInSamples(int newMaxDelay)
  {
    setMaxDelayInSamples(0, newMaxDelay);
    setMaxDelayInSamples(1, newMaxDelay);
    setMaxDelayInSamples(2, newMaxDelay);
    setMaxDelayInSamples(3, newMaxDelay);
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

/** This is a naive implementation of rsTwoPoleAllpassDelay. The implementation is "correct by 
inspection" but horribly wasteful with memory. Switching to direct form 2 should reduce the 
required delay memory by a factor of two and using a 2-tap delayline instead of two separate 
delaylines for the z^-M and z^-2M terms should reduce it by another another factor of 1.5. So, 
overall, we use 3 times as much delay memory as a sensible implementation should. That's why this 
is a naive prototype. */

template<class TSig, class TPar>
class rsTwoPoleAllpassDelayNaive 
{


public:

  void setMaxDelayInSamples(int newMaxDelay)
  {
    inputDelayLine1. setMaximumDelayInSamples(  newMaxDelay);
    outputDelayLine1.setMaximumDelayInSamples(  newMaxDelay);
    inputDelayLine2. setMaximumDelayInSamples(2*newMaxDelay);
    outputDelayLine2.setMaximumDelayInSamples(2*newMaxDelay);
  }

  void setDelayInSamples(int newDelay)
  {
    inputDelayLine1.setDelayInSamples(   newDelay);
    outputDelayLine1.setDelayInSamples(  newDelay);
    inputDelayLine2.setDelayInSamples( 2*newDelay);
    outputDelayLine2.setDelayInSamples(2*newDelay);
  }

  void setAllpassCoeffs(TPar newCoeff1, TPar newCoeff2) 
  { 
    coeff1 = newCoeff1;
    coeff2 = newCoeff2;
  }

  inline TSig getSample(TSig x)
  {
    // Retrieve delayed inputs and outputs:
    TSig xM  = inputDelayLine1.getSample(x);                             // x[n-M]
    TSig yM  = outputDelayLine1.getSampleSuppressTapIncrements(TSig(0)); // y[n-M]
    TSig x2M = inputDelayLine2.getSample(x);                             // x[n-2*M]
    TSig y2M = outputDelayLine2.getSampleSuppressTapIncrements(TSig(0)); // y[n-2*M]

    // Compute current output:
    TSig y = coeff2 * x + coeff1*xM + x2M - coeff1 * yM - coeff2 * y2M; // y[n], our current output

    // Update the output delaylines and return result:
    outputDelayLine1.addToInput(y);
    outputDelayLine1.incrementTapPointers();
    outputDelayLine2.addToInput(y);
    outputDelayLine2.incrementTapPointers();
    return y;
    // ToDo: verify that this does the right thing with respect to the order of reading, writing and
    // incrementing the taps of the outputDelayLine. Maybe write a unit test that uses a delay of 
    // M = 1 and compare output to a regular first order allpass filter.
    //
    // We want to realize:
    //
    //          c2  +  c1 * z^(-M)  +       z^(-2M)
    //  H(z) = -------------------------------------
    //          1   +  c1 * z^(-M)  +  c2 * z^(-2M)
    //
    //  y[n] = c2 * x[n] + c1 * x[n-M] + x[n-2M] - c1 * y[n-M] - c2 * y[n-2M]
  }


  void reset()
  {
    inputDelayLine1.reset();
    outputDelayLine1.reset();
    inputDelayLine2.reset();
    outputDelayLine2.reset();
  }


protected:

  RAPT::rsBasicDelayLine<TSig> inputDelayLine1;
  RAPT::rsBasicDelayLine<TSig> inputDelayLine2;
  RAPT::rsBasicDelayLine<TSig> outputDelayLine1;
  RAPT::rsBasicDelayLine<TSig> outputDelayLine2;
  TPar coeff1 = 0.0;
  TPar coeff2 = 0.0;
};


//=================================================================================================

/** With this class, we try to generalize from the two pole to the N-pole case. So, what we want 
to realize is:


          c_N  +  c_{N-1} * z^(-M)  +  c_{N-2} * z^(-2M)  + ... +       z^(-NM)
  H(z) = -----------------------------------------------------------------------
          1    +  c_1     * z^(-M)  +  c_2     * z^(-2M)  + ... + c_N * z^(-NM)


  y[n] = c_N * x[n] + c_{N-1} * x[n-M] + c_{N-2} * x[n-2M] + ... +       x[n-NM] 
                    - c_1     * y[n-M] - c_2     * y[n-2M] - ... - c_N * y[n-NM]

This time, we start from the implementation of 2-pole case in rsTwoPoleAllpassDelay and try to 
mold it into an N-pole case.

...TBC... */

template<class TSig, class TPar>
class rsMultiPoleAllpassDelayProto
{


public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxDelayInSamples(int newMaxDelay)
  {
    maxM = newMaxDelay;
    allocateMemory();
  }

  void setDelayInSamples(int newDelay)
  {
    M = newDelay;
    delayLine.setDelayInSamples(N*M);
  }

  void setAllpassCoeffs(const std::vector<TPar>& newCoeffs)
  { 
    N = (int) newCoeffs.size() - 1;
    allocateMemory();
    for(int i = 1; i <= N; i++)
      c[i] = newCoeffs[i];  // Maybe use rsArrayTools::copy
    c[0] = 1;  // This shopuld always be the case
  }
  // ToDo: change this signature later to work with a raw pointer and a length N. But during 
  // development, it's more convenient this way.


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline TSig getSample(TSig x)
  {
    // Read delayed states from delayline:
    for(int i = 1; i <= N; i++)
      v[i] = delayLine.readOutputAt(i*M); // Read v_iM = v[n-i*M] from delayline
    // Actually, that step could be merged with the one below, I think. We don't really need to 
    // copy these out. Optimize this later! ..or maybe we do because we need the cvalues in the 
    // output computation. But actually the values are all still in the delay-line. But maybe it's
    // more efficient that way because readOutput at might be more expensive than a simple array
    // access

    // Compute current state:
    v[0] = x;
    for(int i = 1; i <= N; i++)
      v[0] -= c[i] * v[i];

    // Compute output:
    TSig y = v[N];
    for(int i = 1; i <= N; i++)
      y += c[i] * v[N-i];

    // Write v[n] into delayline, increment taps and return result:
    delayLine.writeInputAndUpdate(v[0]);
    return y;
  }
  // Needs tests! Compare it with N = 2 to the result of rsTwoPoleAllpassDelay.

  /*
  // 2-pole implementation for reference:
  inline TSig getSample(TSig x)
  {
    TSig v1M = delayLine.readOutputAt(M);   // Read v1M = v[n-1*M] from delayline
    TSig v2M = delayLine.readOutput();      // Read v2M = v[n-2*M] from delayline
    TSig v   = x - c1 * vM - c2 * v2M;      // Compute v[n] = x[n] - c1 * v[n-M] - c2 * v[n-2M]
    delayLine.writeInputAndUpdate(v);       // Write v[n] into delayline and increment taps
    return c2 * v + c1 * vM + v2M;          // Return y[n] = c2 * v[n] + c1 * v[n-M] + v[n-2M]
  }
  */

  void reset()
  {
    delayLine.reset();
  }


protected:

  void allocateMemory()
  {
    delayLine.setMaximumDelayInSamples(N * maxM);
    c.resize(N+1);
    v.resize(N+1);

    // We actually need only N, not N+1. We don't really use v[0] and c[0]. But it's easier to
    // write the code this way because we don't have to worry about a lot of -1s.
  }

  RAPT::rsBasicDelayLine<TSig> delayLine;
  std::vector<TPar> c;
  std::vector<TSig> v;

  int N    = 0;  // Prototype order
  int M    = 0;  // Delay amount
  int maxM = 1;  // Maximum for M
};

//=================================================================================================

/** Optimized version */

template<class TSig, class TPar>
class rsMultiPoleAllpassDelay
{


public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxDelayInSamples(int newMaxDelay)
  {
    maxM = newMaxDelay;
    allocateMemory();
  }

  void setDelayInSamples(int newDelay)
  {
    M = newDelay;
    delayLine.setDelayInSamples(N*M);
  }

  void setAllpassCoeffs(const TPar* newCoeffs, int numCoeffsExlcudingC0)
  {
    N = numCoeffsExlcudingC0;
    allocateMemory();
    for(int i = 0; i < N; i++)
      c[i] = newCoeffs[i];
  }
  // needs test



  // Get rid:
  void setAllpassCoeffs(const std::vector<TPar>& newCoeffs)
  { 
    N = (int) newCoeffs.size() - 1;
    allocateMemory();
    for(int i = 1; i <= N; i++)
      c[i] = newCoeffs[i];  // Maybe use rsArrayTools::copy
    c[0] = 1;  
  }
  // ToDo: change this signature later to work with a raw pointer and a length N. But during 
  // development, it's more convenient this way.


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */


  inline TSig getSample(TSig x)
  {
    // Compute current state:
    TSig vNew = x;
    for(int i = 1; i <= N; i++)
      vNew -= c[i-1] * delayLine.readOutputAt(i*M);

    // Compute output:
    TSig y = delayLine.readOutputAt(N*M);
    for(int i = 1; i < N; i++)
      y += c[i-1] * delayLine.readOutputAt((N-i)*M);
    y += c[N-1] * vNew;

    // Write vNew into delayline, increment taps and return result:
    delayLine.writeInputAndUpdate(vNew);
    return y;
  }
  // Needs tests! Compare it with N = 2 to the result of rsTwoPoleAllpassDelay.

  void reset()
  {
    delayLine.reset();
  }


protected:

  void allocateMemory()
  {
    delayLine.setMaximumDelayInSamples(N * maxM);
    c.resize(N+1);
    // Later we will want to use N
  }

  RAPT::rsBasicDelayLine<TSig> delayLine;
  std::vector<TPar> c;

  int N    = 0;  // Prototype order
  int M    = 0;  // Delay amount
  int maxM = 1;  // Maximum for M
};


/*
ToDo: generalize this idea to arbitrary order filters with arbitrary delays, i.e. realize:

  y[n] = c_M x[n] + c_{M-1} x[n-d_1] + c_{M-2} x[n-d_2] + ... +     x[n-d_M]
                  - c_1     y[n-d_1] - c_2     y[n-d_2] - ... - c_M y[n-d_M]

This formula needs to be verified. This could perhaps be realized with a multitap delayline.
Can we then also build nested structure from these units?  
*/




#endif