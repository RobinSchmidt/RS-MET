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
    c[0] = 1;  // This should always be the case in the passed newCoeffs array anyway
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
    // We use lengths of N+1 to better match the indices used in the math equations. This is just 
    // a prototype so it the focus is on ease of recognition of the math concepts and formulas.
  }

  RAPT::rsBasicDelayLine<TSig> delayLine;
  std::vector<TPar> c;
  std::vector<TSig> v;

  int N    = 0;  // Prototype order
  int M    = 0;  // Delay amount
  int maxM = 1;  // Maximum for M
};

//=================================================================================================

/** Production version - Document and move to the library...but first add some more unit tests */

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

  /** The a-array of coefficients should look like as follows:

    a[0]  a[1]  a[2] ... a[N-1]
    c_1   c_2   c_3  ... c_N

  That is, the implicit a[0] = 1 coeff shall *not* be included in the array as a convenience 
  dummy. So the length of the passed array must be equal to the order of the prototype allpass 
  which is equal to N in the table above. */
  void setAllpassCoeffs(const TPar* a, int order)
  {
    N = order;
    allocateMemory();
    for(int i = 0; i < N; i++)
      c[i] = a[i];
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline TSig getSample(TSig x)
  {
    // Compute current state v:
    TSig v = x;
    for(int i = 1; i <= N; i++)
      v -= c[i-1] * delayLine.readOutputAt(i*M);

    // Compute output y:
    TSig y = delayLine.readOutputAt(N*M);
    for(int i = 1; i < N; i++)
      y += c[i-1] * delayLine.readOutputAt((N-i)*M);
    y += c[N-1] * v;

    // Write current state v into delayline, increment taps and return result:
    delayLine.writeInputAndUpdate(v);
    return y;
  }

  void reset()
  {
    delayLine.reset();
  }


protected:

  void allocateMemory()
  {
    delayLine.setMaximumDelayInSamples(N * maxM);
    c.resize(N);
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

//=================================================================================================

/** A super simple unit delay class. Might be convenient for implementing certain prototypes. */

template<class T>
class rsUnitDelay
{

public:

  T getSample(T x)
  {
    T y = x1;   // Output computation
    x1  = x;    // State update
    return y;
  }

  void reset() { x1 = 0; }

protected:

  T x1 = 0;

};

//=================================================================================================

/** We encapsulate into a class the code implemented in  feedbackFilterAllpass()  in 
DelayExperiments.cpp  ...TBC...

*/

template<class TSig, class TPar>
class rsDampedAllpassCombNaive
{

  // Maybe rename to rsDampedAllpassComb. 

public:


  void setMaxDelayInSamples(int newMaxDelay);


  void setup(int delay, TSig feedback, int dampOrder, TPar* dampCoeffsB, TPar* dampCoeffsA, 
             bool predelay);

  // We use the convention that we use  M = delay - 1  for the delayline to compensate for the unit
  // delay. This makes more sense from a user's perspective because then, the spike spacing is 
  // exactly given by delay.


  void reset();

  /** This is the normal getSample funtion to be used when you want to produce the allpass output.
  It calls getSampleComb() and then applyCorrectionFilter() on the result of that. You may
  be interested in using the object without the correction filter to produce only the pure comb 
  filter output. That's why I have split it that way so you can also call getSampleComb if that's
  what you want. You can then just ignore the correction filter - or you can apply it yourself but
  maybe after messing with comb output. I don't know, if that's useful though, but you can do it. 
  ...TBC...   */
  TSig getSample(TSig in);

  /** This implements producing samples for the damped delay feedback loop alone, i.e. without the
  correction filter applied. */
  TSig getSampleComb(TSig in);
  
  /** This function is supposed to be called with the output produced by getSampleComb to apply the
  correction filter. */
  TSig applyCorrector(TSig combOutput);
  // rename to applyCorrector

protected:


  // Objects for implementing the A(z) / (1 + k * z^-1 * F(z) * A(z)), i.e. the uncorrected comb
  // filter with filtered unit delay feedback:
  rsBasicDelayLine<TSig>         mainDelay;
  rsDirectFormFilter<TSig, TPar> damper;

  // Objects for the correction filter:
  rsUnitDelay<TSig>              unitDelay;
  rsBasicDelayLine<TSig>         corDelayM1;
  rsBasicDelayLine<TSig>         corDelayM2;
  rsDirectFormFilter<TSig, TPar> corPoles;
  rsDirectFormFilter<TSig, TPar> invDamper;

  // State for the unit delay feedback loop:
  TSig out = TSig(0);

  // Coefficients:
  TSig k;
  TSig r0, r1, rM1, rM2;

  bool preDelay = false;
};

template<class TSig, class TPar>
void rsDampedAllpassCombNaive<TSig, TPar>::setMaxDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay .setMaximumDelayInSamples(maxM);
  corDelayM1.setMaximumDelayInSamples(maxM+1);
  corDelayM2.setMaximumDelayInSamples(maxM+2);
}

template<class TSig, class TPar>
void rsDampedAllpassCombNaive<TSig, TPar>::setup(
  int delay, TSig feedback, int dampOrder, TPar* b, TPar* a, bool predelay)
{
  rsAssert(dampOrder == 1, "We currently only support 1st order damping filters");
  // The signature allows for higher order damping filters in anticipation of supporting those
  // later.

  rsAssert(a[0] == 1);   
  // We may relax this assumption later. If a[0] != 1, we can just scale all coeffs by 1/a[0]

  k = feedback;
  this->preDelay = predelay;

  // Set up damper and related filters:
  TPar t[5] = { 1,0,0,0,0 };
  damper.setCoefficients(    a, b, dampOrder);
  invDamper.setCoefficients( a, b, dampOrder);
  invDamper.invert();
  corPoles.setCoefficients(  a, t, dampOrder);

  // Set up delaylines:
  int M = delay - 1;                        // -1 corrects for unit delay in feedback path
  mainDelay.setDelayInSamples(M);
  corDelayM1.setDelayInSamples(M+1);
  corDelayM2.setDelayInSamples(M+2);

  // Compute correction coefficients:
  r0  = k * b[1];
  r1  = k * b[0];
  rM1 =     a[1];
  rM2 = 1;
}

template<class TSig, class TPar>
void rsDampedAllpassCombNaive<TSig, TPar>::reset()
{
  mainDelay.reset();
  damper.reset();
  unitDelay.reset();
  corDelayM1.reset();
  corDelayM2.reset();
  corPoles.reset();
  invDamper.reset();
  out = TSig(0);
}

template<class TSig, class TPar>
TSig rsDampedAllpassCombNaive<TSig, TPar>::getSample(TSig in)
{
  return applyCorrector(getSampleComb(in));
}

template<class TSig, class TPar>
TSig rsDampedAllpassCombNaive<TSig, TPar>::getSampleComb(TSig in)
{
  if(preDelay)
  {
    out = mainDelay.getSample(in - k * damper.getSample(out));
    return out;
  }
  else
  {
    out = damper.getSample(in - k * mainDelay.getSample(out));
    return invDamper.getSample(out);
    // It may seem strange that we first apply the damper and then the inverse damper. Doesn't this
    // mean, we could just leave out the damper entirely? No! Because "out" is a state variable 
    // that will be used in the next call. And that state needs to have the damper applied.
  }

  // In the dampedAllpassComb1() experiment where I derived all of this, I actually use a negative
  // sign for the feedback signal. There's some comment about why, but I'm a bit shaky on this. But 
  // this might explain why we have 
}

template<class TSig, class TPar>
TSig rsDampedAllpassCombNaive<TSig, TPar>::applyCorrector(TSig in)
{
  // Apply 1-pole:
  TSig t = corPoles.getSample(in);

  // Apply the FIR part:
  TSig y = 0;
  y += r0  * t;
  y += r1  * unitDelay.getSample(t);
  y += rM1 * corDelayM1.getSample(t);
  y += rM2 * corDelayM2.getSample(t);

  return y;
}


// A free function to set up the object with a more convenient parametrization:
template<class TSig, class TPar>
void rsSetupHighDamp(rsDampedAllpassCombNaive<TSig, TPar>& flt,
  int delay, TSig feedback, TPar dampOmega, TPar dampGain, bool predelay)
{
  // Set up one pole filters:
  TPar b0, b1, a1;
  rsFirstOrderFilterBase<TSig, TPar>::coeffsHighShelfBLT(dampOmega, dampGain, &b0, &b1, &a1);
  a1 = -a1; // We want to use the y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1] sign convention here
  TPar ta[2] = { 1,  a1 };
  TPar tb[2] = { b0, b1 };
  flt.setup(delay, feedback, 1, tb, ta, predelay);
}



//=================================================================================================

/** This is the less naive version meant to go into production someday */

template<class TSig, class TPar>
class rsDampedAllpassComb
{


public:

  rsDampedAllpassComb()
  {


  }


  /** Sets the maximum desired roundtrip delay around the comb. This determines the spacing of the
  spikes in the impulse response in a setting without any decay or damping. In such a case, the 
  first spike appears at delay - 1 and from there, the subsequent ones are spaced apart by delay
  itself. The fact that first spike appears at delay - 1 rather than delay itself has to do with 
  the unit delay in the feedback loop. */
  void setMaxDelayInSamples(int newMaxDelay);



  void setup(int delay, TSig feedback, int dampOrder, TPar* dampCoeffsB, TPar* dampCoeffsA, 
    bool predelay);



  /** Resets the state. */
  void reset();

  /** This is the normal getSample funtion to be used when you want to produce the allpass output.
  It calls getSampleComb() and then applyCorrectionFilter() on the result of that. You may
  be interested in using the object without the correction filter to produce only the pure comb 
  filter output. That's why I have split it that way so you can also call getSampleComb if that's
  what you want. You can then just ignore the correction filter - or you can apply it yourself but
  maybe after messing with comb output. I don't know, if that's useful though, but you can do it. 
  ...TBC...   */
  TSig getSample(TSig in);

  /** This implements producing samples for the damped delay feedback loop alone, i.e. without the
  correction filter applied. */
  TSig getSampleComb(TSig in);

  /** This function is supposed to be called with the output produced by getSampleComb to apply the
  correction filter. */
  TSig applyCorrector(TSig combOutput);


protected:

  TSig applyDelay(TSig x)
  {
    return mainDelay.getSample(x);
  }

  TSig applyDamper(TSig x)
  {
    TSig y = b0 * x + b1 * x1d - a1 * y1d;  // ToDo: maybe use DF2 or TDF2 implementation
    x1d = x;
    y1d = y;
    return y;
  }

  TSig applyInverseDamper(TSig x)
  {
    TSig y = (x + a1 * x1di - b1 * y1di) / b0;
    x1di = x;
    y1di = y;
    return y;
  }

  TSig applyCorrectorOnePole(TSig x)
  {
    TSig y = x - a1 * y1c;
    y1c = y;
    return y;
  }

  void updateDelaysAndCorrectorCoeffs()
  {
    // Set up delaylines:
    mainDelay.setDelayInSamples(M);
    corrDelay.setDelayInSamples(M+2);

    // Compute correction coefficients:
    r0  = k*b1;
    r1  = k*b0;
    rM1 = a1;
  }


  // Objects for implementing the A(z) / (1 + k * z^-1 * F(z) * A(z)), i.e. the uncorrected comb
  // filter with filtered unit delay feedback:
  rsBasicDelayLine<TSig> mainDelay;

  // Objects for the correction filter:
  rsUnitDelay<TSig>      unitDelay;
  rsBasicDelayLine<TSig> corrDelay;

  // State for the unit delay feedback loop:
  TSig combOut = TSig(0);

  // States for the two one pole filters:
  TSig x1d  = 0, y1d  = 0;   // x[n-1], y[n-1] for damper
  TSig x1di = 0, y1di = 0;   // x[n-1], y[n-1] for inverse damper
  TSig y1c  = 0;             // y[n-1] for corrector one pole

  // Coefficients:
  TSig k   = 0;              // Needs to be TSig when we want to use it with complex feedback
  TSig r0  = 0; 
  TSig r1  = 0; 
  TPar rM1 = 0;
  TPar b0 = 0, b1 = 0, a1 = 0;
  // Maybe use arrays b[0], b[1], r[0], r[1], etc. Later they may get longer because we want to
  // support more complex feedback filters
  int  M = 0;

  bool preDelay = false;

};

template<class TSig, class TPar>
void rsDampedAllpassComb<TSig, TPar>::setMaxDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay .setMaximumDelayInSamples(maxM);
  corrDelay.setMaximumDelayInSamples(maxM+2);
}

template<class TSig, class TPar>
void rsDampedAllpassComb<TSig, TPar>::setup(int delay, TSig feedback, int dampOrder,
  TPar* dampCoeffsB, TPar* dampCoeffsA, bool predelay)
{
  M = delay - 1;                            // -1 corrects for unit delay in feedback path
  k = feedback;
  this->preDelay = predelay;

  rsAssert(dampOrder == 1);
  // We do not yet support higher order damping filters. This feature is under constructions

  rsAssert(dampCoeffsA[0] == 1);
  a1 = dampCoeffsA[1];
  b0 = dampCoeffsB[0];
  b1 = dampCoeffsB[1];

  updateDelaysAndCorrectorCoeffs();
}

template<class TSig, class TPar>
void rsDampedAllpassComb<TSig, TPar>::reset()
{
  mainDelay.reset();
  unitDelay.reset();
  corrDelay.reset();
  combOut = TSig(0);
  x1d     = TSig(0);
  y1d     = TSig(0);
  x1di    = TSig(0);
  y1di    = TSig(0);
  y1c     = TSig(0);
}

template<class TSig, class TPar>
TSig rsDampedAllpassComb<TSig, TPar>::getSample(TSig in)
{
  return applyCorrector(getSampleComb(in));
}

template<class TSig, class TPar>
TSig rsDampedAllpassComb<TSig, TPar>::getSampleComb(TSig in)
{
  if(preDelay)
  {
    combOut = applyDelay(in - k * applyDamper(combOut));  // Predelay of M samples
    return combOut;
  }
  else
  {
    combOut = applyDamper(in - k * applyDelay(combOut));  // No predelay
    return applyInverseDamper(combOut);
  }

  // Notes:
  //
  // - This computation has an implicit unit delay applied to the apperance of combOut on the right
  //   hand side. On the right hand side, it's the previous comb output. On the left hand side, 
  //   it's the current comb output.
  //
  // - Applying the mainDelay as inner filter and the damper as outer filter gives us a filter 
  //   without any predelay/latency. Most of the time, this is more desirable, but maybe it could
  //   be useful for something to have predelay built in after all. 
}

template<class TSig, class TPar>
TSig rsDampedAllpassComb<TSig, TPar>::applyCorrector(TSig in)
{
  // Apply 1-pole:
  TSig t = applyCorrectorOnePole(in);

  // Apply the FIR part:
  TSig y = 0;
  y += r0  * t;
  y += r1  * unitDelay.getSample(t);
  y += rM1 * corrDelay.readOutputAt(M+1);
  y +=       corrDelay.readOutputAt(M+2);
  corrDelay.writeInputAndUpdate(t);
  return y;
}

// More Optimization ideas:
// 
// - Maybe implement the unit delay inline. Although I don't think that this causes overhead.
//
// - Maybe the two calls to corrDelay.readOutputAt(M+1); corrDelay.readOutputAt(M+2); can be
//   replaced by a clever arrangement of calling getSample() and readOutput(). The readOutpuAt
//   function is slightly more expensive because it computes the offset taking care of 
//   wrapraounds etc whereas the others just use stored member variables. But I have not yet 
//   figured out how to do it or if it's even possible. Maybe we don't even need the M member
//   anymore then.
//
// - I checked the contents of mainDelay and corrDelay to see if we can use a shared delayline but
//   that doesn't seem to be possible. I've also switched the order of applying FIR part and pole
//   in applyCorrector to see if then the content can be shared. Nope.
//
//
//
// - To do this optimization, we need the following steps:
//
//     (1) Don't use getSample() on mainDelay - do the /write/read/increment manually
//     (2) Use M+2 as length for mainDelay instead of M. That should now be safe to do
//     (3) In applyCorrector(), do the reads from mainDelay rather than corrDelay



// A free function to set up the object with a more convenient parametrization:
template<class TSig, class TPar>
void rsSetupHighDamp(rsDampedAllpassComb<TSig, TPar>& flt,
  int delay, TSig feedback, TPar dampOmega, TPar dampGain, bool predelay)
{
  // Set up one pole filters:
  TPar b0, b1, a1;
  rsFirstOrderFilterBase<TSig, TPar>::coeffsHighShelfBLT(dampOmega, dampGain, &b0, &b1, &a1);
  a1 = -a1; // We want to use the y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1] sign convention here
  TPar ta[2] = { 1,  a1 };
  TPar tb[2] = { b0, b1 };
  flt.setup(delay, feedback, 1, tb, ta, predelay);
}






//=================================================================================================

/** A nonlinear extension of rsDampedAllpassComb */

template<class TSig, class TPar>
class rsDampedAllpassCombNonLin : public rsDampedAllpassComb<TSig, TPar>
{

public:

  TSig getSampleComb(TSig in) // override - but at compile time
  {
    if(preDelay)
      combOut = mainDelay.getSample(fbFunc(in + k * applyDamper(combOut)));
    else
      combOut = applyDamper(fbFunc(in + k * mainDelay.getSample(combOut)));
    return combOut;
  }

protected:

  //std::function<TSig(TSig x)> fbFunc = &tanh;
  //std::function<TSig(TSig x)> fbFunc = tanh;

  std::function<TSig(TSig x)> fbFunc = [](TSig x)
  { 
    return rsTanh(x);
  };


};





#endif