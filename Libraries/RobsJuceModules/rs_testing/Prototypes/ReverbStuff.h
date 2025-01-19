#ifndef RS_REVERBSTUFF_H
#define RS_REVERBSTUFF_H

// ToDo: 
// -Maybe move those classes that eventually go into the RAPT library to the top and those that are
//  really only prototypes (naive implementations that serve as reference for unit tests) to the 
//  bottom





//=================================================================================================

// Maybe remove TPar and replace it by double...but maybe someone wants to use float?


/** This class is a lightweight wrapper around the rsDelay class that gives it an API that is 
compatible with our more advanced interpolating delayline implementations such as rsDelayLinear
and rsDelayAllpass. The idea is that one might start to write a DSP algorithm using the class 
rsDelayRounding (i.e. start simple) and later replace it by some API-compatible other (better) 
delay class as refinement of the algo. */

template<class TSig, class TPar>
class rsDelayRounding
{

public:



  void setMaxDelayInSamples(TPar newMaxDelay)
  { dl.setMaxDelayInSamples(rsRoundToInt(newMaxDelay)); }


  void setDelayInSamples(TPar newDelay)
  { dl.setDelayInSamples(rsRoundToInt(newDelay)); }


  /** Returns the value of the transfer function H(z) at the given value of z. If M is the delay in
  samples, then H(z) = z^-M. */
  rsComplex<TPar> getTransferFunctionAt(rsComplex<TPar> z) const
  {
    int M = dl.getDelayInSamples();        // M is our delay
    return rsPow(z, rsComplex<TPar>(-M));  // H(z) = z^-M
  }

  void getTransferFunction(rsSparseDigitalTransferFunction<TPar>* tf) const
  {
    tf->num.setNumTerms(1); tf->num.setTerm(0, TPar(1), dl.getDelayInSamples());
    tf->den.setNumTerms(1); tf->den.setTerm(0, TPar(1), 0);
  }
  // Needs tests.




  TSig getSample(TSig x) { return dl.getSample(x); }

  void reset() { dl.reset(); }



protected:

  rsDelay<TSig> dl;    // Underlying integer delayline


};

// Hmmm...I think, it might be better to just have a member of type rsDelay<TSig> rather than 
// deriving from it. We do not want expose the API of rsDelay to client code. We want the API to be
// consistent with the ones of rsDelayLinear and rsDelayAllpass.



/** A delayline class that uses linear interpolation to achieve fractional (i.e. non integer) 
delays. */

template<class TSig, class TPar>
class rsDelayLinear
{


public:



  void setMaxDelayInSamples(TPar newMaxDelay)
  {
    dl.setMaxDelayInSamples(rsCeilInt(newMaxDelay) + 1);

    // I think, we really need the +1 because even when the user request an integer delay of 
    // M = maxDelay, we still will access the delayline with a delay of M+1, albeit with a 
    // coeff of zero. Hmm...well...if we just use rsCeilInt without the +1, then trying to access 
    // the delay M+1 may wrap around to another location and read a wrong sample from there. But 
    // since the coeff is zero anyway, it doesn't really matter that we access the wrong sample. 
    // So, it may actually be ok to just use rsCeilInt(newMaxDelay) without the +1. But better safe
    // than sorry. ...maybe do tests. Try it with maxDelay = 7 or 15. The actual maxDelay will always
    // be 2^k-1 for some k to enable the wrapping via bitmasking. If it turns out to be safe to use
    // it without the +1, then maybe do it. I really like being tight with memory allocations, i.e.
    // really allocate exactly as much as you need and nothing extra.
  }

  void setDelayInSamples(TPar newDelay)
  {
    TPar i = rsFloor(newDelay);     // Integer part
    TPar f = newDelay - i;          // Fractional part
    dl.setDelayInSamples(int(i));
    b0 = TPar(1) - f;
    b1 = f;
  }



  TSig getSample(TSig x)
  {
    dl.writeInputNoUpdate(x);
    TSig x0 = dl.readOutput();
    TSig x1 = dl.readOutputWithAdditionalDelay(1);
    dl.incrementTapPointers();

    return b0*x0 + b1*x1;
  }

  void reset() { dl.reset(); }


protected:

  rsDelay<TSig> dl;    // Underlying integer delayline

  TPar b0 = TPar(1);   // Coeff for x[n-M]
  TPar b1 = TPar(0);   // Coeff for x[n-M-1]

  // Maybe store juts one coeff (b1, the fractional part of the delay) and use the one-multiply 
  // form in the implementation...well...maybe try both and benchmark and choose the faster

};
// Needs tests


template<class TSig, class TPar>
class rsDelayAllpass
{


public:

  void setMaxDelayInSamples(TPar newMaxDelay)
  {
    dl.setMaxDelayInSamples(rsCeilInt(newMaxDelay) + 1);
  }

  void setDelayInSamples(TPar newDelay)
  {
    TPar i = rsFloor(newDelay);
    TPar f = newDelay - i;
    dl.setDelayInSamples(int(i));
    c = (1-f) / (1+f);
  }


  TSig getSample(TSig x)
  {
    dl.writeInputNoUpdate(x);
    TSig x0 = dl.readOutput();
    TSig x1 = dl.readOutputWithAdditionalDelay(1);
    dl.incrementTapPointers();

    y1 = c*x0 + x1 - c*y1;   // Verify, maybe optimize to  x1 + c*(x0-y1)
    return y1;

    // Maybe we need to scale the feedback by some number like 0.999 to avoid a parasitic
    // oscillations at the Nyquist freq for certain settings. See the old implemementations. The 
    // oscillation occurs when c is close to 1. This happens the the fractional part f is zero.
  }

  void reset() 
  { 
    dl.reset();
    y1 = 0; 
  }


protected:

  rsDelay<TSig> dl;    // Underlying integer delayline

  TSig y1 = TSig(0);   // Previous output
  TPar c  = TPar(0);   // Allpass filter coefficient

};
// Needs tests





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
    inputDelayLine.setMaxDelayInSamples(newMaxDelay);
    outputDelayLine.setMaxDelayInSamples(newMaxDelay);
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

  RAPT::rsDelay<TSig> inputDelayLine;
  RAPT::rsDelay<TSig> outputDelayLine;
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
      delayLine.setMaxDelayInSamples(newMaxDelay);
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
  RAPT::rsDelay<TSig> delayLine;
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
      delayLine.setMaxDelayInSamples(newMaxDelay);
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
  RAPT::rsDelay<TSig> delayLine;

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
      delayLine.setMaxDelayInSamples(newMaxDelay);
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
  RAPT::rsDelay<TSig> delayLine;

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
    inputDelayLine1. setMaxDelayInSamples(  newMaxDelay);
    outputDelayLine1.setMaxDelayInSamples(  newMaxDelay);
    inputDelayLine2. setMaxDelayInSamples(2*newMaxDelay);
    outputDelayLine2.setMaxDelayInSamples(2*newMaxDelay);
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

  RAPT::rsDelay<TSig> inputDelayLine1;
  RAPT::rsDelay<TSig> inputDelayLine2;
  RAPT::rsDelay<TSig> outputDelayLine1;
  RAPT::rsDelay<TSig> outputDelayLine2;
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
    delayLine.setMaxDelayInSamples(N * maxM);
    c.resize(N+1);
    v.resize(N+1);
    // We use lengths of N+1 to better match the indices used in the math equations. This is just 
    // a prototype so it the focus is on ease of recognition of the math concepts and formulas.
  }

  RAPT::rsDelay<TSig> delayLine;
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
    delayLine.setMaxDelayInSamples(N * maxM);
    c.resize(N);
  }

  RAPT::rsDelay<TSig> delayLine;
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

/** A little helper function to conveniently design a 1st order high shelving filter. This kind of 
filter is used a lot in the damped allpass combs below. The parameters w,g are "omega" and the 
linear shelver gain. respectively. */
template<class T>
void rsMake1stOrderHighShelf(T w, T g, T* b0, T* b1, T* a1)
{
  rsFirstOrderFilterBase<T, T>::coeffsHighShelfBLT(w, g, b0, b1, a1);
  *a1 = -*a1; 
  // We want to design according to the y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1] sign convention here
  // but silly rsFirstOrderFilterBase uses  y[n] = b0*x[n] + b1*x[n-1] + a1*y[n-1], so we need to 
  // flip the sign of a1.
}

template<class T>
void rsMake1stOrderLowShelf(T w, T g, T* b0, T* b1, T* a1)
{
  rsFirstOrderFilterBase<T, T>::coeffsLowShelfBLT(w, g, b0, b1, a1);
  *a1 = -*a1; 
}

// Maybe add rsMake1stOrderAllpass to add dispersion. But we may also add dispersion by using
// maximum-phase low- and high-shelvers



//=================================================================================================

/** We encapsulate into a class the code implemented in  feedbackFilterAllpass()  in 
DelayExperiments.cpp  ...TBC...

*/

template<class TSig, class TPar>
class rsDampedCombAllpassNaive
{

  // Maybe rename to rsDampedCombAllpass. 

public:

  rsDampedCombAllpassNaive()
  {

  }


  void setMaxDelayInSamples(int newMaxDelay);


  void setup(int delay, TPar feedback, int dampOrder, TPar* dampCoeffsB, TPar* dampCoeffsA, 
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
  rsDelay<TSig>         mainDelay;
  rsDirectFormFilter<TSig, TPar> damper;

  // Objects for the correction filter:
  rsUnitDelay<TSig>              unitDelay;

  rsDelay<TSig>         corDelayM1;
  rsDelay<TSig>         corDelayM2;

  rsDirectFormFilter<TSig, TPar> corPoles;
  rsDirectFormFilter<TSig, TPar> invDamper;

  // State for the unit delay feedback loop:
  TSig out = TSig(0);

  // Coefficients:
  TPar k;
  TPar r0, r1, rM1, rM2;

  bool preDelay = false;
};

template<class TSig, class TPar>
void rsDampedCombAllpassNaive<TSig, TPar>::setMaxDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay .setMaxDelayInSamples(maxM);
  corDelayM1.setMaxDelayInSamples(maxM+1);
  corDelayM2.setMaxDelayInSamples(maxM+2);
}

template<class TSig, class TPar>
void rsDampedCombAllpassNaive<TSig, TPar>::setup(
  int delay, TPar feedback, int dampOrder, TPar* b, TPar* a, bool predelay)
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
void rsDampedCombAllpassNaive<TSig, TPar>::reset()
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
TSig rsDampedCombAllpassNaive<TSig, TPar>::getSample(TSig in)
{
  return applyCorrector(getSampleComb(in));
}

template<class TSig, class TPar>
TSig rsDampedCombAllpassNaive<TSig, TPar>::getSampleComb(TSig in)
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

  // In the dampedCombAllpass1() experiment where I derived all of this, I actually use a negative
  // sign for the feedback signal. There's some comment about why, but I'm a bit shaky on this. But 
  // this might explain why we have 
}

template<class TSig, class TPar>
TSig rsDampedCombAllpassNaive<TSig, TPar>::applyCorrector(TSig in)
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

  // Maybe re-implement the FIR part as:
  //
  //   y += corrZerosFront.getSample(t) + corrZerosBack.getSample(corrDelayM1.getSample(t))
  //
  // where corrZerosFront/Back implement the two parts of the FIR filter. One (front) part acts
  // on the non-delayed t and the other (back) on the delayed t. We may then get rid of the 
  // unitDelay and the corDelayM2 because these delays would then be part of the two corZeros
  // filters. This implementation would make it clearer what happens, I think.
}

// A free function to set up the object with a more convenient parametrization:
template<class TSig, class TPar>
void rsSetupHighDamp(rsDampedCombAllpassNaive<TSig, TPar>& flt,
  int delay, TPar feedback, TPar dampOmega, TPar dampGain, bool predelay)
{
  TPar a[2], b[2]; a[0] = 1;
  rsMake1stOrderHighShelf(dampOmega, dampGain, &b[0], &b[1], &a[1]);
  flt.setup(delay, feedback, 1, b, a, predelay);
}

//=================================================================================================

/** This class implements a delayline based allpass filter based on the following block diagram:

  X(z) ---> + -----> z^-M ----------------> C(z) -----> Y(z)
            ^                    |
            |                    |
           -k                   z^-1
            |                    |
            -------- F(z) <-------

There is a delayline of length M around which we have a feedback loop with a feedback filter F(z) 
and a scalar feedback gain k and a unit delay z^-1 in the loop to make it realizable. This part
of the filter so far, taken by itself, implements a comb filter. The strategy is now to apply an 
appropriate compensation filter C(z) to make the overall structure allpass in nature. The 
derivation of this filter C(z) is outlined in Notes/DSP/DampedCombAllpass.txt. The short version 
of it that the comb part has a transfer function:

                   z^-M
  U(z) = ----------------------------
          1 + k * z^-1 * F(z) * z^-M

and we start by defining a preliminary compensation filter C~(z) as follows:

   C~(z) = 1 + k * z^-1 * F(z) * z^-M

This filter would just cancel the denominator and bring us back to a pure delay z^-M. But as, said,
that was just our preliminary compensation filter. The actual compensation filter C(z) is obtained
from the preliminary one by reflecting its zeros about the unit circle. That amounts to just 
reversing its FIR part. The user can specify the feedback filter in terms of its direct form filter 
coefficients and we currently support feedback filters of orders up to 8. 

The class also supports a second mode of operation in which the places of z^-M and F(z) are 
swapped in the block diagram. It turns out that the compensation filter will then need to include 
an inverted damping filter, i.e. F^-1(z) = 1/F(z), but is otherwise the same. This second mode of 
operation has no initial predelay. That is, the first nonzero sample of the impulse response occurs
at sample index n = 0 whereas in the depiction above, the first nonzero sample could clearly not 
occur before n = M because the signal has to go through the delayline before appearing at the 
output. In fact, it appears exactly at n = M. 

The result is an interesting allpass filter with the potential to introduce a frequency dependent 
decay by choosing the feedback filter appropriately. For example, a high shelving filter that 
attenuates high frequencies would make high frequencies decay faster to get Karplus-Strong like
behavior. 


Stability:

For stability, the feedback parameter k should be restricted to -1 <= k <= +1 and the feedback 
filter F(z) should have a magnitude response that is less or equal to one (i.e. |F(w)| <= 1) for 
all frequencies w. Well, strictly speaking, what you actually want is |k * F(w)| <= 1 for all w, so
you could theoretically have F(w) = 2 if you choose |k| <= 0.5, etc. But I really like to normalize
F(w) such that it peaks at 1 and then adjust the overall feedback gain via k. If you want to use 
the mode without predelay, you need to be careful to pass a filter F(z) that has a stable inverse 
(i.e. is minimum phase) and there are no checks and warnings about this. That may mean to stay away
from BLT-based lowpasses as they tend to have zeros at z = -1 which in the inversion will become 
marginally stable poles. I'd rather recommend to go with impulse invariance based allpole lowpasses
or with shelving or peak/bell filters with negative dB gains. 

*/

template<class TSig, class TPar>
class rsDampedCombAllpass
{


public:


  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Standard constructor. Initializes the settings and resets the state to initial conditions. */
  rsDampedCombAllpass()
  {
    // We may want to use delayline interpolation filters of order 1 (e.g. linear or 1st order 
    // allpass) later which need two terms in numerator and denominator:
    A.setNumTerms(2, 2);
    // Maybe factor this out into a protected function setMaxInterpolatorOrder. The purporse of
    // such a function would be mainly for documentation


    initSettings();
    reset();
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup




  /** Sets the maximum desired roundtrip delay around the comb. This total roundtrip delay includes
  the z^-1 unit delay, so the delayline length is actually shorter by one. */
  void setMaxIntDelayInSamples(int newMaxDelay);



  void setMaxDelayInSamples(TPar newMaxDelay)
  {
    setMaxIntDelayInSamples((int)rsCeil(newMaxDelay));


    //setMaxIntDelayInSamples((int)rsCeil(rsReal(newMaxDelay)));

    // Using rsReal() here is needed for enabling instantiating this class also for TPar being a 
    // complex type. This enables complex valued feedback etc. The delay should still be a real 
    // type though, so we just extract the real part.
    //
    // But I think, this still doesn't work. The code for the experiment with a complex feedback 
    // is currently commented out. I think, we need to introduce a 3rd template parameter. Maybe
    // TFdbk (for the feedback) or TDly (for the delay). But I'm not sure, if the transfer function
    // computation will work because it may try to deal with nested complex numbers when TPar is
    // itself complex. Maybe these function need to be templates themselves, introducing their own
    // template parameter TComplex?)
  }





  /** Sets up the filter with the given total roundtrip delay, the scalar feedback gain and the
  coefficients of the damping filter to be used. The filter is supposed to be given in direct form
  and realizes:

    y[n] = b[0] * x[n] + b[1] * x[n-1] + ... + b[P] * x[n-P]
                       - a[1] * y[n-1] - ... - a[P] * y[n-P]

  where P is the filter order, i.e. the dampOrder parameter and the b,a arrays in the formula map 
  t the dampCoeffsB, dampCoeffsA parameters respectively. We assume that the filter coeff arrays 
  are normalized to a[0] = 1. The boolean predelayMode parameter switches between two modes of 
  operation one of which features a predelay of delay-1 samples. The -1 occurs because the 
  delayline length M is given by delay-1. That is: the delay user parameter means the total 
  roundtrip delay which includes the delayline delay and the implicit delay in the feedback 
  loop. */
  void setup(int delay, TPar feedback, 
    int dampOrder, const TPar* dampCoeffsB, const TPar* dampCoeffsA, bool predelayMode);
  // ToDo: should take a TPar for the delay

  /** Initializes all settings to default values. */
  void initSettings();


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  /** Returns the maximum order for the feedback damping filters that is supported. */
  static int getMaxDampingOrder() { return maxDmpOrd; }
    // ToDo: Try to make constexpr. But that seems to incompatible with static. Why?


  /** Evaluates the filter's z-domain transfer function H(z) value at the given value of z. */
  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const;
  // Under construction. ToDo: check, if this works with a complex type for TSig.

  rsComplex<TPar> getCombTransferFunctionAt(const rsComplex<TPar>& z) const;

  rsComplex<TPar> getDamperTransferFunctionAt(const rsComplex<TPar>& z) const;

  rsComplex<TPar> getCorrectorTransferFunctionAt(const rsComplex<TPar>& z) const;



  rsSparseDigitalTransferFunction<TPar> getTransferFunction() const;
  // Allocates!

  rsSparseDigitalTransferFunction<TPar> getCombTransferFunction() const;
  // Allocates! 

  rsSparseDigitalTransferFunction<TPar> getCorrectorTransferFunction() const;
  // Allocates!

  rsSparseDigitalTransferFunction<TPar> getDamperTransferFunction() const;
  // Allocates!




  // Non-allocating versions of the transfer function getters. Well - they *may* allocate - but 
  // will do so only if the passed output parameters and temporary objects have not enough 
  // capacity pre-allocated. Doing so is the responsibility of the caller. They are much less 
  // convenient to use but it's sometimes necessary when one needs to compute these transfer 
  // function in a realtime thread.



  void getCorrectorTransferFunction(rsSparseDigitalTransferFunction<TPar>* tf) const;


  void getCombTransferFunction(rsSparseDigitalTransferFunction<TPar>* tf) const;


  void getDamperTransferFunction(rsSparseDigitalTransferFunction<TPar>* tf) const;


  void getDelayTransferFunction(rsSparseDigitalTransferFunction<TPar>* tf) const;




  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** This is the normal getSample() funtion to be used when you want to produce the allpass 
  output. It calls getSampleComb() and then applyCorrector() on the result of that. You may
  be interested in using the object without the correction filter to produce only the pure comb 
  filter output or you may be interested in shoving some other operations in between the comb and 
  the correction filter. That's why I have split it up that way. */
  TSig getSample(TSig in) { return applyCorrector(getSampleComb(in)); }

  /** This implements producing samples for the damped delay feedback loop alone, i.e. without the
  correction filter applied. */
  TSig getSampleComb(TSig in);

  /** This function is supposed to be called with the output produced by getSampleComb() to apply 
  the correction filter. */
  TSig applyCorrector(TSig combOutput);

  /** Resets the state. */
  void reset();


protected:

  /** Applies the main delay to the input x and updates the state of the main delayline. */
  TSig applyDelay(TSig x)
  {
    return mainDelay.getSample(x);
  }

  /** Applies the feedback damping filter to the signal x and updates the filter's state. */
  TSig applyDamper(TSig x)
  {
    // Compute outputs:
    TSig y = b[0]*x;
    for(int i = 1; i <= dmpOrd; i++)
      y += b[i]*xd[i-1] - a[i] * yd[i-1];

    // Update state and return result:
    rsArrayTools::shiftPushDiscard(xd, dmpOrd, x);
    rsArrayTools::shiftPushDiscard(yd, dmpOrd, y);
    return y;
  }

  /** Applies the inverse feedback damping filter to the signal x and updates the filter's state.
  this filter is need only the "without predelay" mode of operation. */
  TSig applyInverseDamper(TSig x)
  {
    // Compute output:
    TSig y = x;
    for(int i = 1; i <= dmpOrd; i++)
      y += a[i] * xi[i-1] - b[i] * yi[i-1];
    y /= b[0];                                        // ToDo: maybe precompute 1/b[0]

    // Update state and return result:
    rsArrayTools::shiftPushDiscard(xi, dmpOrd, x);
    rsArrayTools::shiftPushDiscard(yi, dmpOrd, y);
    return y;
  }

  /** Applies the poles of the correction filter which are the same as the poles of the damping 
  filter. */
  TSig applyCorrectorPoles(TSig x)
  {
    // Compute output:
    TSig y = x;
    for(int i = 1; i <= dmpOrd; i++)
      y -= a[i] * yc[i-1];
    
    // Update state and return result:
    rsArrayTools::shiftPushDiscard(yc, dmpOrd, y);
    return y;
  }


  static const int maxDmpOrd = 8;      // Maximum damping order

  // Embedded DSP objects:

  //rsDelay<TSig> mainDelay;           // Main delayline for the comb filter    old

  rsDelayRounding<TSig, TPar> mainDelay;  // Main delayline for the comb filter    new
  // With this new code, our dampedCombAllpassComplex() experiment doesn't compile anymore. 
  // Something in it causes rsRoundToInt to get called with a complex argument. Maybe we are trying
  // to set up a complex delay somewhere by having a delay parameter declared as TSig rather than 
  // TPar in some setup function?
  //
  // Hmm - I think, we may need to introduce a 3rd template parameter. We may need:
  // TSig, TCoef, TDly  ...do we really need this? It makes the API more unwieldy. But if that's 
  // what it takes then we have to need to do it....



  rsDelay<TSig>               corrDelay;  // Delayline for the correction filter

  // State:
  TSig combOut = TSig(0);                 // State for the unit delay feedback loop
  TSig xd[maxDmpOrd], yd[maxDmpOrd];      // State for the damping filter
  TSig xi[maxDmpOrd], yi[maxDmpOrd];      // State for the inverse damping filter
  TSig yc[maxDmpOrd];                     // State for the poles of the correction filter

  // Coefficients:
  TPar k = 0;                             // Feedback gain
  TPar b[maxDmpOrd+1];                    // Damping filter feedforward coeffs
  TPar a[maxDmpOrd+1];                    // Damping filter feedback coeffs

  // Settings:
  int  M        = 0;                      // Delayline length
  int  dmpOrd   = 0;                      // Feedback damping filter order
  bool preDelay = false;                  // Switch between with/without predelay mode of operation

  // Temporary object for delay transfer function A(z):
  mutable rsSparseDigitalTransferFunction<TPar> A;
    // This member is needed to support a non-allocating implementation of 
    // getDelayTransferFunction() ...maybe call it D(z) for delay



  // Notes:
  //
  // - The feedback gain k is of type TSig rather than TPar to allow usage with TSig == complex and
  //   then allowing complex feedback factors. There is some experiment that does this. It's 
  //   interesting. ...but maybe make k of type TPar anyway. We can still do this experiment by
  //   just using TPar = complex as well. I think, it makes more sense this way.
  //   ...has been changed back to TPar...Hmm...hmm... I'm not sure about it. Requiring the user
  //   to instantiate it with Par = complex to get complex feedback interferes with using TPar for
  //   the non-integer delay parameter ...and probably also with getTransferFunctionAt. Maybe have
  //   a 3rd template parameter TFdbk for the feedback?
  //
  // - Maybe be more flexible with the order of the damping filter by letting numerator and 
  //   denominator have different orders. Maybe replace dmpOrd by two variables bOrd, aOrd or 
  //   something like that.
};

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::setMaxIntDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay.setMaxDelayInSamples(maxM);

  corrDelay.setMaxDelayInSamples(maxM+maxDmpOrd+1);
  // I think, we may have to add the maximum order of the interpolator filter used in the delayline
}

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::setup(int delay, TPar feedback, int dampOrder,
  const TPar* dampCoeffsB, const TPar* dampCoeffsA, bool predelayMode)
{
  if(dampOrder > maxDmpOrd) 
  {
    rsError("Such high damping order is not supported.");
    initSettings();
    return;
  }

  M        = delay - 1;             // -1 corrects for unit delay in feedback path
  k        = feedback;
  preDelay = predelayMode;
  dmpOrd   = dampOrder;

  rsAssert(dampCoeffsA[0] == TPar(1));  // May be relaxed later by dividing through all coeffs by a[0]
  rsArrayTools::copy(dampCoeffsA, a, dmpOrd+1);
  rsArrayTools::copy(dampCoeffsB, b, dmpOrd+1);

  mainDelay.setDelayInSamples(M);

  corrDelay.setDelayInSamples(M+dmpOrd+1);
  // I think, this may need more delay memory when we have an interpolating delayline. I think, we
  // may have to add the order of the interpolator.
}

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::initSettings()
{
  mainDelay.setDelayInSamples(0);
  corrDelay.setDelayInSamples(0);

  M        = 0;
  k        = 0;
  dmpOrd   = 0;
  preDelay = false;

  using AT = rsArrayTools;
  AT::clear(b, maxDmpOrd+1);
  AT::clear(a, maxDmpOrd+1);
}

template<class TSig, class TPar>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar>::getTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  return getCombTransferFunctionAt(z) * getCorrectorTransferFunctionAt(z);
}

template<class TSig, class TPar>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar>::getCombTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  using Complex = rsComplex<TPar>;
  Complex one(TPar(1));                          // 1 + 0i

  //Complex zM = rsPow(z, Complex(-M));            // z^-M
  Complex zM = mainDelay.getTransferFunctionAt(z);   // A(z)
  // or maybe use A = getDelayTransferFunctionAt(z)

  Complex z1 = one/z;                            // z^-1
  Complex F  = getDamperTransferFunctionAt(z);   // F(z)
  if(preDelay)
    return zM  / (one + k * z1 * F * zM);        // U(z) = z^-M / (1 + k * z^-1 * F(z) * z^-M)
  else
    return one / (one + k * z1 * F * zM);        // U(z) =   1  / (1 + k * z^-1 * F(z) * z^-M)

  // Notes:
  //
  // - The comb transfer function in the 2nd branch (i.e. without predelay) already includes the
  //   inverse damping filter which cancels the effect of the damping filter. That's because
  //   getSampleComb() already applies the inverse damping filter. That's why we don't see an F in 
  //   the numerator. It cancels with the same F that would appear in the denominator due to the
  //   application of the inverse damper. So, overall, the numerator turns out to be just 1.
  //
  // - Rename zM to A (for allpass) and retrieve it from the delaylien via a call like 
  //   mainDelay.getTransferFunctionAt(z) which has to be implemented. We can then replace the 
  //   integer delayline with a fractional one that implements this function also and computes the
  //   correct transfer function for the selected interpolation method (linear, allpass, etc.).
}

template<class TSig, class TPar>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar>::getDamperTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  using Complex = rsComplex<TPar>;
  Complex num = 0, den = 0;
  for(int i = 0; i <= dmpOrd; i++)
  {
    Complex zi = rsPow(z, Complex(-i));          // z^-i
    num += b[i] * zi;
    den += a[i] * zi;
  }
  return num / den;

  // ToDo:
  //
  // - This should be optimized (don't call rsPow - compute the powers on the fly by multiplying by
  //   z) and factored into a library function to compute the transfer function of direct form 
  //   filters. Maybe it should go into rsFilterAnalyzer.
}

template<class TSig, class TPar>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar>::getCorrectorTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  using Complex = rsComplex<TPar>;
  Complex num = 0, den = 0;
  for(int i = 0; i <= dmpOrd; i++)
  {
    den +=     a[i]        * rsPow(z, Complex(-i));
    num += k * b[dmpOrd-i] * rsPow(z, Complex(-i));
    num +=     a[dmpOrd-i] * rsPow(z, Complex(-(M+1+i)));
  }
  return num / den;

  // ToDo:
  //
  // - Optimize! The current implementation is horribly inefficient. But maybe keep it for the 
  //   naive implementation. But before doing this, implement a unit test for the transfer function
  //   computation
}

template<class TSig, class TPar>
rsSparseDigitalTransferFunction<TPar> rsDampedCombAllpass<TSig, TPar>::getTransferFunction() const
{
  return getCombTransferFunction() * getCorrectorTransferFunction();
}

template<class TSig, class TPar>
rsSparseDigitalTransferFunction<TPar> rsDampedCombAllpass<TSig, TPar>
                                      ::getCombTransferFunction() const
{
  using TF = rsSparseDigitalTransferFunction<TPar>;

  TF one; one.num.appendTerm(TPar(1), 0);
  TF z1;  z1.num.appendTerm( TPar(1), 1);
  TF F = getDamperTransferFunction();            // Feedback filter F(z)
  getDelayTransferFunction(&A);                  // Delay filter A(z)
  if(preDelay)
    return A   / (one + k * z1 * F * A);         // U(z) = A(z) / (1 + k * z^-1 * F(z) * A(z))
  else 
    return one / (one + k * z1 * F * A);         // U(z) =   1  / (1 + k * z^-1 * F(z) * A(z))
}

template<class TSig, class TPar>
rsSparseDigitalTransferFunction<TPar> rsDampedCombAllpass<TSig, TPar>
                                      ::getCorrectorTransferFunction() const
{
  using TF = rsSparseDigitalTransferFunction<TPar>;

  TF C = getCombTransferFunction();
  C.invert();
  C.reflectZeros();

  return C;
}

template<class TSig, class TPar>
rsSparseDigitalTransferFunction<TPar> rsDampedCombAllpass<TSig, TPar>
                                      ::getDamperTransferFunction() const
{
  rsSparseDigitalTransferFunction<TPar> H;
  H.setupFromDenseCoeffs(b, dmpOrd+1, a, dmpOrd+1, TPar(0));
  return H;
}

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::getCombTransferFunction(
  rsSparseDigitalTransferFunction<TPar>* tf) const
{
  using Mon = rsMonomial<TPar>;
  getDelayTransferFunction(&A);        // A = A(z) is transfer function of the delay
  getDamperTransferFunction(tf);       // tf = F, F(z) is transfer function in feedback path
  tf->multiplyBy(Mon(k, 1));           // tf = F * k * z^-1
  tf->multiplyBy(A, TPar(0));          // tf = F * k * z^-1 * A
  tf->addConstant(TPar(1), TPar(0));   // tf = 1 + F * k * z^-1 * A
  tf->invert();                        // tf = 1 / (1 + F * k * z^-1 * A)
  if(preDelay)
    tf->multiplyBy(A, TPar(0));        // tf = A / (1 + F * k * z^-1 * A)

  // ToDo:
  //
  // - Verify that all operations above are non-allocating (assuming that tf has enough capacity)
  //   and document that fact. Our member A also needs to have enough capacity to represent the
  //   delay filter including interpolation.
}

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::getCorrectorTransferFunction(
  rsSparseDigitalTransferFunction<TPar>* tf) const
{
  getCombTransferFunction(tf);
  tf->invert();
  tf->reflectZeros();
}

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::getDamperTransferFunction(
  rsSparseDigitalTransferFunction<TPar>* tf) const
{
  tf->setupFromDenseCoeffs(b, dmpOrd+1, a, dmpOrd+1, TPar(0));
}

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::getDelayTransferFunction(
  rsSparseDigitalTransferFunction<TPar>* tf) const
{
  mainDelay.getTransferFunction(tf);
}

template<class TSig, class TPar>
void rsDampedCombAllpass<TSig, TPar>::reset()
{
  mainDelay.reset();
  corrDelay.reset();
  combOut = TSig(0);

  using AT = rsArrayTools;
  AT::clear(xd, maxDmpOrd);
  AT::clear(yd, maxDmpOrd);
  AT::clear(xi, maxDmpOrd);
  AT::clear(yi, maxDmpOrd);
  AT::clear(yc, maxDmpOrd);
}

template<class TSig, class TPar>
TSig rsDampedCombAllpass<TSig, TPar>::getSampleComb(TSig in)
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
  // - This computation has an implicit unit delay applied to the appearance of combOut on the right
  //   hand side. On the right hand side, it's the previous comb output. On the left hand side, 
  //   it's the current comb output.
  //
  // - Applying the mainDelay as inner filter and the damper as outer filter gives us a filter 
  //   without any predelay/latency. Most of the time, this is more desirable, but maybe it could
  //   be useful for something to have predelay built in after all. 
}

template<class TSig, class TPar>
TSig rsDampedCombAllpass<TSig, TPar>::applyCorrector(TSig in)
{
  // Apply the poles:
  TSig t = applyCorrectorPoles(in);

  // Apply the FIR part:
  TSig y = 0;
  corrDelay.writeInputNoUpdate(t);
  for(int i = 0; i <= dmpOrd; i++)
  {
    y += k * b[dmpOrd-i] * corrDelay.readOutputAt(i);
    y +=     a[dmpOrd-i] * corrDelay.readOutputAt(M+1+i);
  }
  corrDelay.incrementTapPointers();
  return y;

  // I think, what's happening in the loop at the bottom can be re-interpreted as reading out the
  // delayline at 0 and at M+1 and applying FIR filters to these outputs and then adding both of 
  // these FIR outputs to y. The coeffs for the FIR for the output at 0 are given by k times the 
  // reversal of the b-coeffs and the FIR coeffs for the output at M+1 are given by the reversal 
  // of the a-coeffs of the feedback filter. Or something along these lines. ToDo: Figure this out 
  // exactly! But I don't think, the implementation should change to implement it literally that 
  // way. The way we do it currently seems more efficient. No additional (redundant) filter objects
  // are used here which is actually a good thing from an economic point of view. But maybe in the 
  // naive prototype, we should do it with the additional filters.
}


// A free function to set up the object with a more convenient parametrization:
template<class TSig, class TPar>
void rsSetupHighDamp(rsDampedCombAllpass<TSig, TPar>& flt,
  int delay, TPar feedback, TPar dampOmega, TPar dampGain, bool predelay)
{
  TPar a[2], b[2]; a[0] = 1;
  rsMake1stOrderHighShelf(dampOmega, dampGain, &b[0], &b[1], &a[1]);
  flt.setup(delay, feedback, 1, b, a, predelay);
}



template<class TSig, class TPar>
void rsSetupFractional_LinViaFb(rsDampedCombAllpass<TSig, TPar>& flt,
  TPar delay, TPar feedback, bool predelay)
{
  int  delayInt  = (int) rsFloor(delay);
  TPar delayFrac = delay - (TPar) delayInt;
  TPar a[2], b[2]; 
  a[0] = 1;

  if(delayFrac == TPar(0))
  {
    b[0] = 1;
    flt.setup(delayInt, feedback, 0, b, a, predelay);
  }
  else
  {
    TPar d = delayFrac;

    //// Allpass interpolation:
    //// See: https://ccrma.stanford.edu/~jos/pasp/First_Order_Allpass_Interpolation.html
    //TPar c = (1-d) / (1+d);
    //b[0] = c;
    //b[1] = 1;
    //a[0] = 1;
    //a[1] = c;
    //// Doesn't work! Gives unstable comb filters! Maybe the meaning of y[n-1] is different in
    //// context here?

    // Linear interpolation:
    b[0] = 1 - d;
    b[1] = d;
    a[0] = 1;
    a[1] = 0;
    // This seems to work.

    flt.setup(delayInt, feedback, 1, b, a, predelay);
  }

}
// Rename to rsSetupFractional_LinViaFb where LinViaFb stands for "linear interpolation via the 
// feedback filter"

// Under construction. Should set up the flt such that it achieves a given overall decay time in 
// samples (in the sense of RT60, i.e. reverb time to decay to -60 dB) and having scaled deacy 
// times for low and high frequencies. The scale factors are given as raw factors for the RT60 and
// crossover frequencies are given as omega.
template<class TSig, class TPar>
void rsSetupDecayTimes_LinViaFb(rsDampedCombAllpass<TSig, TPar>& flt, TPar delay, TPar decayTimeInSamples,
  TPar lowOmega, TPar lowTimeScale, TPar highOmega, TPar highTimeScale, bool predelay)
{
  // Compute desired feedback gains for low, mid and high frequencies:
  TPar a60 = TPar(0.001); // = rsDbToAmp(-60.0). Target amplitude to reach after decayTimeInSamples
  TPar kL  = rsDecayTimeToFeedbackGain(decayTimeInSamples * lowTimeScale , TPar(delay), a60);
  TPar kM  = rsDecayTimeToFeedbackGain(decayTimeInSamples                , TPar(delay), a60);
  TPar kH  = rsDecayTimeToFeedbackGain(decayTimeInSamples * highTimeScale, TPar(delay), a60);
  // These formulas could also be expressed as e.g.:
  //
  //   kM = rsPow(10.0, TPar(-3 * delay) / decayTimeInSamples);
  //
  // which is how they are often seen in the FDN literature. 


  // Compute desired gains for the low and high shelver:
  TPar gL = kL / kM;
  TPar gH = kH / kM;

  // Compute coeffs for low- and high shelver:
  TPar aL[2], bL[2]; aL[0] = 1; rsMake1stOrderLowShelf( lowOmega,  gL, &bL[0], &bL[1], &aL[1]);
  TPar aH[2], bH[2]; aH[0] = 1; rsMake1stOrderHighShelf(highOmega, gH, &bH[0], &bH[1], &aH[1]);

  // Combine low- and high shelver into biquad:
  TPar a[4], b[4];
  rsArrayTools::convolve(aL, 2, aH, 2, a);
  rsArrayTools::convolve(bL, 2, bH, 2, b);

  // Possibly also bake an interpolation filter into the feedback filter to achieve fractional 
  // delay times:
  int  delayInt  = (int) rsFloor(delay);
  TPar delayFrac = delay - (TPar) delayInt;
  if(delayFrac == TPar(0))
  {
    // In the integer delay case, we only need the 2nd order feedback filter that we already have:
    flt.setup(delayInt, kM, 2, b, a, predelay);
  }
  else
  {
    // In the fractional delay case, we create a linear interpolation filter and bake it into the 
    // existing 2nd order feedback filter, thereby turning it into a 3rd order filter:

    TPar d = delayFrac;

    // Design linear interpolation filter:
    TPar aI[2], bI[2];
    bI[0] = 1 - d;
    bI[1] = d;
    aI[0] = 1;
    aI[1] = 0;

    // Bake the interpolation filter into the b,a, arrays:
    rsArrayTools::convolve(a, 3, aI, 2, a);
    rsArrayTools::convolve(b, 3, bI, 2, b);
    flt.setup(delayInt, kM, 3, b, a, predelay);
  }



  // ToDo:
  //
  // - Maybe optionally turn the low- and/or high-shelver into a maximum phase version. Maybe 
  //   optionally let the user also add an allpass filter for additional dispersion in the feedback
  //   path.
  //
  // - Compare to implementation of FeedbackDelayNetwork16::updateDampingAndCorrectionFilters. It 
  //   uses class rosic::DampingFilter and it specifies the gains also at the shelver's crossover
  //   frequencies. The idea is that the linear gain at the crossover freq is not defined to be 
  //   just the geometric mean between the actual shelver gain and unity but instead some gain that
  //   let's the decay time at that frequency be the geometric mean between the mid decay time
  //   and the low (or high) frequency decay time.
  //
  // - Maybe we should use rosic::DampingFilter here for the coefficient calculations, too. But I 
  //   think before that, we should refactor the code in such a way that the damping filter itself
  //   handles the computations of the desired gains at the crossover frequencies - which we 
  //   currently do in FeedbackDelayNetwork16::updateDampingAndCorrectionFilters(). I'm not sure, 
  //   if it's really worth the trouble to do it like this, though. It will just slightly(?) change
  //   the response/feeling of the lowFreq/lowScale, highFreq/highScale parameters. It may be a 
  //   more natural response, though. I think, it may help to decouple the respective freq and 
  //   scale parameters. ...but it's quite complicated to implement... But maybe it doesn't have to
  //   be that complicated. Maybe we could use a general 3-point filter-design routine that lets
  //   the user specify 3 omegas and the 3 corresponding gains. Or maybe we could use the general
  //   5-point biquad design method that takes 5 omegas and 5 magnitudes. The omegas would be
  //   DC, loFreq, sqrt(loFreq*hiFreq), hiFreq, fs/2. But what if loFreq==hiFreq? I guess, we would
  //   get a singular system of equations.
  //
  // - Pass flt by pointer
  //
  // - Maybe implement a filter that realizes a fractional delay by introducing another 1st order 
  //   allpass - like in allpass interpolation. Maybe the allpass should be adjusted to take into
  //   account the effect of the damping filter (which itself may also introduce a frequency 
  //   dependent delay). I think, what we want is to have the correct fractional delay at DC or
  //   maybe at the resonance frequency, so we can tune it exactly.
  //
  // - Can be optimized: design the low shelf directly into a,b, the high shelf into some tmp 
  //   arrays, then bake them into a,b. The allpass can then use the same temp arrays as the 
  //   hi shelf and then also bake them into a,b
}
// Rename to rsSetupDecayTimes_LinViaFb where LinViaFb stands for "linear interpolation via the 
// feedback filter"


// Notes:
// 
// - I checked the contents of mainDelay and corrDelay to see if we can use a shared delayline but
//   that doesn't seem to be possible. I've also switched the order of applying FIR part and pole
//   in applyCorrector to see if then the content can be shared. Nope. But maybe it's possible to
//   share the delayline, if we somehow change the structure of the computations? 
//
//
// ToDo:
//
// - Maybe we can replace the "predelay" parameter (i.e. binary mode switch) with a more general
//   mode switch. I could also imagine to use it in "Schroeder mode" , i.e . with feedforward path
//   around the main delay line. For this, we could repurpose the invDamper for the 2nd damper
//   filter that sits in the feedforward path. But this setup even allows for an implementation
//   with just a single delayline, so maybe it should go into a separate class. But if we can 
//   provide this mode as option here, too then maybe it would be nice to have. I have not yet 
//   worked the math though, so I'm not yet sure if that is even workable. It appears to be 
//   plausible, though.



//=================================================================================================


//=================================================================================================

/** UNDER CONSTRUCTION

A class that creates an allpass filter out of a linear combination of multiple combs. 

...TBC... */


template<class TSig, class TPar>
class rsDampedMultiCombAllpass   // Maybe rename to rsDampedCombBankAllpass
{

public:

  void setFilterOrderLimits(int newMaxDelayInSamples, int newMaxNumCombs)
  {
    // This function is not yet complete and not yet unit tested. We need to pre-allocate enough
    // memory in all the objects here so as to avoid any re-allocations in updateFilters() called
    // by getSample(). Failing to get this right may cause occasional (rare) memory allocations 
    // on the audio thread which may go unnoticed most of the time and cause very rare audio 
    // glitches.


    maxNumCombs = maxNumCombs;
    maxDelay    = newMaxDelayInSamples;
    protoAllpass.setMaxDelayInSamples(maxDelay);
    settings.resize(maxNumCombs);


    int maxBankDelay = maxDelay * maxNumCombs + 4;
    // VERIFY this formula! theoretically and practically (by making a unit test using the maxmimum
    // number of combs each at the maximum possible delay) This is just a first rough guess and 
    // might be wrong!


    combBank.setMaxDelayInSamples(maxBankDelay);
    corrector.setMaxDelayInSamples(maxBankDelay);


    // ToDo:
    //
    // - Preallocate enough memory in U an Ui. For this, we need to first figure out, how much 
    //   could be needed in the worst case.
  }





  void setSampleRate(TPar newSampleRate)        { sampleRate = newSampleRate;     setDirty(); }
  void setFrequency(TPar newFrequency)          { frequency = newFrequency;       setDirty(); }
  void setDecayTimeInSeconds(TPar newDecayTime) { decayTime = newDecayTime;       setDirty(); }
  void setMaxPhaseCombBank(bool useMaxPhase)    { maxPhaseCombBank = useMaxPhase; setDirty(); }

  void setLowCrossoverFreq(TPar newFreq)        { lowCrossFreq = newFreq;        setDirty(); }
  void setLowDecayScale(TPar newTimeScale)      { lowDecayScale = newTimeScale;  setDirty(); }
  void setHighCrossoverFreq(TPar newFreq)       { highCrossFreq = newFreq;       setDirty(); }
  void setHighDecayScale(TPar newTimeScale)     { highDecayScale = newTimeScale; setDirty(); }
  // ToDo compare function names to what we have in the FDN classes and make the consistent

  // Maybe use a setDirty() function instead of dirty = true. Maybe also have a setClean() function
  // or maybe let setDirty have an optional bool parameter






  void setNumCombs(int newNumber) 
  { 
    rsAssert(numCombs <= maxNumCombs);  
    numCombs = rsMin(newNumber, maxNumCombs);
    setDirty();
  }

  void setCombFreqScale(int index, TPar newScale)
  {
    rsAssert(isValidCombIndex(index));
    settings[index].freqScale = newScale;
    setDirty();
  }

  void setCombGain(int index, TPar newGain)
  {
    rsAssert(isValidCombIndex(index));
    settings[index].gain = newGain;
    setDirty();
  }




  
  int getNumCombs() const { return numCombs; }

  bool isValidCombIndex(int index) const { return index <= getNumCombs(); }
  


  rsSparseDigitalTransferFunction<TPar> getTransferFunction() const
  {
    return getCombTransferFunction() * getCorrectorTransferFunction();
  }

  rsSparseDigitalTransferFunction<TPar> getCombTransferFunction() const
  {
    return combBank.getTransferFunction();
    //return U;  should also work
  }
  // allocates - creates copy of the transfer function object.
  // Maybe return a const ref?

  rsSparseDigitalTransferFunction<TPar> getCorrectorTransferFunction() const
  {
    return corrector.getTransferFunction();
  }







  TSig getSample(TSig in)
  {
    if(dirty)
      updateFilters(); // rename to updateCoeffs

    return corrector.getSample(combBank.getSample(in));
  }

  void reset()
  {
    combBank.reset();
    corrector.reset();
  }


protected:

  void setDirty(bool shouldBeDirty = true)
  {
    dirty = shouldBeDirty;   
    // Maybe use dirty.store(shouldBeDirty). I think, it makes no difference but it may add 
    // documentation value. We would document that this is intended to be an atomic operation.
  }

  void updateFilters();
  // Allocates! Not yet realtime ready.




  // Embedded DSP objects:
  rsSparseFilter<TSig, TPar> combBank;
  rsSparseFilter<TSig, TPar> corrector;


  // A damped comb allpass object used prototype to compute the filter coeffs:
  rsDampedCombAllpass<TSig, TPar> protoAllpass; 
  // This is not ideal! It contains itself two delaylines which are not needed here and therefore
  // just waste memory. ToDo: Refactor such that the coefficient calculation can be done without
  // having to use such an object. Maybe the coeff calculation can be done by a static member 
  // function? We'll see...


  /** Struct for the settings that we have per comb */
  struct CombSettings
  {
    TPar freqScale = TPar(1);
    TPar gain      = TPar(1);

    // More Settings to add later:
    //bool onlyOdds  = false;

    //bool maxPhaseLoShelf = false;
    //bool maxPhaseHiShelf = false;

    //bool bypassLoShelf   = false;
    //bool bypassHiShelf   = false;
    // The bypass switches are need because shelves with neutral settings are not really neutral
    // but nontrivial allpasses, I think. ...but figure this out! Or maybe we can come up with a
    // different 1st order shelver design that actually is neutral with neutral settings?
  };

  // Per comb settings:
  std::vector<CombSettings> settings;

  // Global settings:
  TPar sampleRate     = TPar(44100);
  TPar frequency      = TPar(440);
  TPar decayTime      = TPar(0.25);
  TPar lowCrossFreq   = TPar(250);
  TPar lowDecayScale  = TPar(2.0);
  TPar highCrossFreq  = TPar(4000);
  TPar highDecayScale = TPar(0.5);


  // Maximum number of available combs and max delayline length. Must be set up on construction:
  int maxNumCombs = 4;
  int maxDelay    = 16383;  // 16383 = 2^-14 - 1, 2^k - 1 is used anyway for bit-masking
  // ToDo: provide a setter for these. This setter should only be called in suspended state, 
  // though. It may re-allocate

  // Current number of combs:
  int numCombs    = 1;

  // Switch beween min- and max-phase comb bank:
  bool maxPhaseCombBank = false;

  // Flag to indicate that a call to updateFilters() is needed before doing any DSP:
  std::atomic<bool> dirty = true;

  // Transfer function objects used for temporaries in internal computations in updateFilters():
  rsSparseDigitalTransferFunction<TPar> U, Ui;
};


template<class TSig, class TPar>
void rsDampedMultiCombAllpass<TSig, TPar>::updateFilters()
{
  // This is still under construction. It still has allocations and it needs to treat the case 
  // numCombs == 0

  using RatFunc   = rsSparseRationalFunction<TPar>;
  //using TransFunc = rsSparseDigitalTransferFunction<TPar>;

  TPar decaySamples =      decayTime     * sampleRate;
  TPar lowOmega     = 2*PI*lowCrossFreq  / sampleRate;
  TPar highOmega    = 2*PI*highCrossFreq / sampleRate;

  // ToDo: catch special case of numCombs == 0. In this case, we should set up a trivial bypass 
  // filter


  // Accumulate the transfer function of the comb bank:
  U.clear();                                              // Init to U(z) = 0.
  for(int i = 0; i < numCombs; i++)
  {
    const CombSettings& s = settings[i];


    TPar delay = sampleRate / (s.freqScale * frequency);  // Verify!
    // Maybe in case of only odd harmonics, we should scale the freq up by 2? The rationale is that 
    // when using odd harmonics only (by way of the feedback sign), the fundamental frequency 
    // actually goes an octave lower.

    rsSetupDecayTimes_LinViaFb(
      protoAllpass, delay, decaySamples, 
      lowOmega,  lowDecayScale, highOmega, highDecayScale, false);
    // This needs an additional parameter to determine the sign of the feedback, i.e. switch 
    // between all and only odd harmonics


    // Accumulate the i-th comb's tranfer function Ui into our total transfer function U:
    protoAllpass.getCombTransferFunction(&Ui);
    RatFunc::weightedSumDestructive(&U, TPar(1), &Ui, TPar(s.gain), &U, TPar(0));
  }


  // New:
  combBank.setup(U);
  if(maxPhaseCombBank)
    combBank.reflectZeros();
  corrector.setup(U);
  corrector.invert();
  corrector.reflectZeros();

  setDirty(false);
}


//=================================================================================================

/** This is a special trimmed down version of rsDampedCombAllpass that only allows for a first 
order filter in the feedback loop. I think, this is a common case that is worth to have some 
optimized code for. The general version with arbitrary feedback filters needs a much more
complicated implementation. */

template<class TSig, class TPar>
class rsDampedCombAllpass_1p
{

public:

  void setMaxDelayInSamples(int newMaxDelay);
  void setup(int delay, TPar feedback, 
    TPar dampCoeffB0, TPar dampCoeffB1, TPar dampCoeffA1, bool predelay);

  TSig getSample(TSig in);
  TSig getSampleComb(TSig in);
  TSig applyCorrector(TSig combOutput);
  void reset();

protected:

  TSig applyDelay(TSig x)
  {
    return mainDelay.getSample(x);
  }

  TSig applyDamper(TSig x)
  {
    TSig y = b0 * x + b1 * x1d - a1 * y1d; 
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
    mainDelay.setDelayInSamples(M);
    corrDelay.setDelayInSamples(M+2);
    r0  = k*b1;
    r1  = k*b0;
    rM1 = a1;
  }


  rsDelay<TSig> mainDelay;
  rsUnitDelay<TSig>      unitDelay;
  rsDelay<TSig> corrDelay;

  TSig combOut = TSig(0);

  TSig x1d  = 0, y1d  = 0;
  TSig x1di = 0, y1di = 0;
  TSig y1c  = 0;

  TPar k   = 0;  // Use TPar
  TPar r0  = 0; 
  TPar r1  = 0; 
  TPar rM1 = 0;
  TPar b0 = 0, b1 = 0, a1 = 0;
  // Get rid of the r-coeffs! r0 = k*b1, r1 = k*b0, rM1 = a1. Use that directly

  int  M = 0;
  bool preDelay = false;

};

template<class TSig, class TPar>
void rsDampedCombAllpass_1p<TSig, TPar>::setMaxDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay .setMaxDelayInSamples(maxM);
  corrDelay.setMaxDelayInSamples(maxM+2);
}

template<class TSig, class TPar>
void rsDampedCombAllpass_1p<TSig, TPar>::setup(int delay, TPar feedback, 
  TPar dampCoeffB0, TPar dampCoeffB1, TPar dampCoeffA1, bool predelay)
{
  M = delay - 1;
  k = feedback;
  this->preDelay = predelay;

  a1 = dampCoeffA1;
  b0 = dampCoeffB0;
  b1 = dampCoeffB1;

  updateDelaysAndCorrectorCoeffs();
}

template<class TSig, class TPar>
void rsDampedCombAllpass_1p<TSig, TPar>::reset()
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
TSig rsDampedCombAllpass_1p<TSig, TPar>::getSample(TSig in)
{
  return applyCorrector(getSampleComb(in));
}

template<class TSig, class TPar>
TSig rsDampedCombAllpass_1p<TSig, TPar>::getSampleComb(TSig in)
{
  if(preDelay)
  {
    combOut = applyDelay(in - k * applyDamper(combOut));
    return combOut;
  }
  else
  {
    combOut = applyDamper(in - k * applyDelay(combOut));
    return applyInverseDamper(combOut);
  }
}

template<class TSig, class TPar>
TSig rsDampedCombAllpass_1p<TSig, TPar>::applyCorrector(TSig in)
{
  TSig t = applyCorrectorOnePole(in);

  TSig y = 0;
  y += r0  * t;
  y += r1  * unitDelay.getSample(t);
  y += rM1 * corrDelay.readOutputAt(M+1);
  y +=       corrDelay.readOutputAt(M+2);
  corrDelay.writeInputAndUpdate(t);
  return y;
}

template<class TSig, class TPar>
void rsSetupHighDamp(rsDampedCombAllpass_1p<TSig, TPar>& flt,
  int delay, TPar feedback, TPar dampOmega, TPar dampGain, bool predelay)
{
  TPar a[2], b[2]; a[0] = 1;
  rsMake1stOrderHighShelf(dampOmega, dampGain, &b[0], &b[1], &a[1]);
  flt.setup(delay, feedback, b[0], b[1], a[1], predelay);
}


//=================================================================================================

/** Like rsDampedCombAllpass_1p but with two combs in parallel instead of just one. The outputs of
the two combs are scaled by weighting factors and then added together. Then, a compensation filter
is applied to that to make the whole filter allpass..

...TBC... see AllpassStuff.txt in the private repo for more details  */

template<class TSig, class TPar>
class rsDampedAllpassBiComb_1p  // rename to rsDampedBiCombAllpass
{


public:

  rsDampedAllpassBiComb_1p() 
  {
    // Maybe pre-allocate some memory here.
  }


  void setMaxDelayInSamples(int newMaxDelay);

  void setup(
    int delay1, TPar gain1, TPar feedback1, TPar coeffB10, TPar coeffB11, TPar coeffA11,
    int delay2, TPar gain2, TPar feedback2, TPar coeffB20, TPar coeffB21, TPar coeffA21);
  // ToDo: Add a parameter for switching the mode of operation. I'm not yet sure if we should have
  // 2 or 3 modes operation, though. We'll see...



  rsComplex<TPar> getCombTransferFunctionAt(const rsComplex<TPar>& z) const;



  int getCombSumOrder() const { return M1+M2+4; }
  // Verify this! I think, this is the total resulting order of the filter. It can be read off from
  // the line  setA(11, M1+M2+4, b11*b21 * k1*k2);   in  convertCombSumToDirectForm


  /** Converts the weighted sum of the two comb filters into a (sparse) direct form filter. The
  object is passed as pointer - the passed rsSparseFilter object serves as output variable.  */
  void convertCombSumToDirectForm(rsSparseFilter<TSig, TPar>* sparseDirectFormFilter);
  // Not true anymore:
  // The function may trigger a memory allocation in the passed filter object if it doesn't 
  // already have enough memory allocated. You probably wan to avoid calling it on a realtime 
  // thread, or if you do, make very sure that the filter object already has enough delay 
  // capacity.
  // ToDo: document, how much delay memory the filter object needs to have pre-allocated. I think
  // it's M1+M2+4. see setMaxDelayInSamples(). We actually use this function internally with our
  // corrector member



  TSig getSample(TSig in) { return applyCorrector(getSampleCombs(in)); }

  TSig getSampleCombs(TSig in) { return g1 * getSampleComb1(in) + g2 * getSampleComb2(in); }

  TSig applyCorrector(TSig combOutput) { return corrector.getSample(combOutput); }

  void reset();


protected:


  TSig applyDamper1(TSig x)
  {
    TSig y = b10 * x + b11 * x11d - a11 * y11d; 
    x11d = x;
    y11d = y;
    return y;
  }
  TSig applyInverseDamper1(TSig x)
  {
    TSig y = (x + a11 * x11di - b11 * y11di) / b10;
    x11di = x;
    y11di = y;
    return y;
  }
  TSig getSampleComb1(TSig in)
  {
    combOut1 = applyDamper1(in - k1 * mainDelay1.getSample(combOut1));
    if(dampCompensated)
      return applyInverseDamper1(combOut1);
    else
      return combOut1;
  }


  TSig applyDamper2(TSig x)
  {
    TSig y = b20 * x + b21 * x21d - a21 * y21d; 
    x21d = x;
    y21d = y;
    return y;
  }
  TSig applyInverseDamper2(TSig x)
  {
    TSig y = (x + a21 * x21di - b21 * y21di) / b20;
    x21di = x;
    y21di = y;
    return y;
  }
  TSig getSampleComb2(TSig in)
  {
    combOut2 = applyDamper2(in - k2 * mainDelay2.getSample(combOut2));
    if(dampCompensated)
      return applyInverseDamper2(combOut2);
    else
      return combOut2;
  }





  void updateDelaysAndCorrectorCoeffs();



  // The two delayines for the two parallel comb filters:
  rsDelay<TSig> mainDelay1;  // Rename to combDelay1 or delayLine1
  rsDelay<TSig> mainDelay2;

  // Correction filter to turn the whole filter into an allpass:
  rsSparseFilter<TSig, TPar> corrector;

  // The stored comb outputs for use in feedback loop:
  TSig combOut1 = TSig(0);
  TSig combOut2 = TSig(0);

  // Feedback gains:
  TPar k1 = 0;
  TPar k2 = 0;
  // Maybe make them TPar - we don't really want to use it with complex valued feedback...or do we?

  // Gains or weights for the two comb outputs:
  TPar g1 = 1;
  TPar g2 = 1;

  // Coeffs for the two feedback filters in the two combs:
  TPar b10 = 0, b11 = 0, a11 = 0;
  TPar b20 = 0, b21 = 0, a21 = 0;

  // Lengths of the delaylines for the combs:
  int  M1 = 0;
  int  M2 = 0;

  // States of feedback filters:
  TSig x11d  = 0, y11d  = 0;
  TSig x21d  = 0, y21d  = 0;

  // States of the inverse feedback filters:
  TSig x11di  = 0, y11di  = 0;
  TSig x21di  = 0, y21di  = 0;

  // Switch between two modes of operation:
  bool dampCompensated = true;
  //bool dampCompensated = false;
  // ToDo: maybe have 3 modes - include one with predelay where the delays sit in the numerator.

};


template<class TSig, class TPar>
void rsDampedAllpassBiComb_1p<TSig, TPar>::rsDampedAllpassBiComb_1p<TSig, TPar>::reset()
{
  mainDelay1.reset();
  mainDelay2.reset();

  combOut1 = 0;
  combOut2 = 0;

  x11d  = 0;
  y11d  = 0;
  x21d  = 0;
  y21d  = 0;

  x11di = 0;
  y11di = 0;
  x21di = 0;
  y21di = 0;
}

template<class TSig, class TPar>
void rsDampedAllpassBiComb_1p<TSig, TPar>::setMaxDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay1.setMaxDelayInSamples(maxM);
  mainDelay2.setMaxDelayInSamples(maxM);
  corrector.setMaxDelayInSamples(2*maxM+4);
  // See convertCombSumToDirectForm(). The maximum delay that occurs there is: M1+M2+4.
}

template<class TSig, class TPar>
void rsDampedAllpassBiComb_1p<TSig, TPar>::setup(
  int delay1, TPar gain1, TPar feedback1, TPar coeffB10, TPar coeffB11, TPar coeffA11,
  int delay2, TPar gain2, TPar feedback2, TPar coeffB20, TPar coeffB21, TPar coeffA21)
{
  M1 = delay1 - 1;
  M2 = delay2 - 1;

  g1 = gain1;
  g2 = gain2;

  k1 = feedback1;
  k2 = feedback2;

  b10 = coeffB10;
  b11 = coeffB11;
  a11 = coeffA11;

  b20 = coeffB20;
  b21 = coeffB21;
  a21 = coeffA21;

  updateDelaysAndCorrectorCoeffs();
}


template<class TSig, class TPar>
rsComplex<TPar> rsDampedAllpassBiComb_1p<TSig, TPar>::getCombTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  rsError("Not yet implemented");
  return rsComplex<TPar>(0); 

  // ToDo:
  //
  // - Compute the transfer function based on computing the two transfer functions of the two 
  //   combs and forming a weighted sum of them.
}

template<class TSig, class TPar>
void rsDampedAllpassBiComb_1p<TSig, TPar>::convertCombSumToDirectForm(
  rsSparseFilter<TSig, TPar>* sparseFilter)
{
  // Set up feedforward coeffs:
  sparseFilter->setNumNumeratorTerms(9);
  auto setB = [&](int index, int delay, TSig coeff)
  {
    sparseFilter->setNumeratorTerm(index, coeff, delay);
  };
  if(dampCompensated)
  {
    setB(0, 0,    g1 + g2);
    setB(1, 1,    g1*(a11+a21) + g2*(a11+a21));
    setB(2, 2,    g1*a11*a21 + g2*a11*a21);
    setB(3, M1+1, g2*k1*b10);
    setB(4, M1+2, g2*k1*(a21*b10 + b11));
    setB(5, M1+3, g2*k1*a21*b11);
    setB(6, M2+1, g1*k2*b20);
    setB(7, M2+2, g1*k2*(a11*b20 + b21));
    setB(8, M2+3, g1*k2*a11*b21);
  }
  else
  {
    setB(0, 0,    b10*g1 + b20*g2);
    setB(1, 1,    g1*(a21*b10 + b11) + g2*(a11*b20 + b21));
    setB(2, 2,    g1*a21*b11 + g2*a11*b21);
    setB(3, M1+1, g2*k1*b10*b20);
    setB(4, M1+2, g2*k1*(b11*b20 + b10*b21));
    setB(5, M1+3, g2*k1*b11*b21);
    setB(6, M2+1, g1*k2*b10*b20);
    setB(7, M2+2, g1*k2*(b11*b20 + b10*b21));
    setB(8, M2+3, g1*k2*b11*b21);
  }
  // ToDo: Maybe have 3 modes - the third is with the z^-M1, z^-M2 in the numerator. But that 
  // filter is not invertible. But maybe it's invertible up to delay? That would actually be good
  // enough. The modes could be given in an enum: withPredelay, compensated, uncompensated
  // Optimize!


  // Set up feedback coeffs:
  sparseFilter->setNumDenominatorTerms(12);
  auto setA = [&](int index, int delay, TSig coeff)
  {
    sparseFilter->setDenominatorTerm(index, coeff, delay);
  };
  setA( 0, 0,       TPar(1));
  setA( 1, 1,       a11 + a21);
  setA( 2, 2,       a11 * a21);
  setA( 3, M1+1,    b10 * k1);
  setA( 4, M1+2,    (a21*b10 + b11) * k1);
  setA( 5, M1+3,    a21 * b11 * k1);
  setA( 6, M2+1,    b20 * k2);
  setA( 7, M2+2,    (a11*b20 + b21) * k2);
  setA( 8, M2+3,    a11 * b21 * k2);
  setA( 9, M1+M2+2, b10*b20 * k1*k2);
  setA(10, M1+M2+3, (b11*b20 + b10*b21) * k1*k2);
  setA(11, M1+M2+4, b11*b21 * k1*k2);

  // Update the length of the delayline:
  //sparseFilter->ensureEnoughDelayMemory();  // May allocate!
  sparseFilter->updateDelayLineLength();

  // ToDo:
  //
  // - Maybe extract common subexpressions like (a21*b10 + b11), (b11*b20 + b10*b21), g1*k2,
  //   etc.
}


template<class TSig, class TPar>
void rsDampedAllpassBiComb_1p<TSig, TPar>::updateDelaysAndCorrectorCoeffs()
{
  mainDelay1.setDelayInSamples(M1);
  mainDelay2.setDelayInSamples(M2);

  convertCombSumToDirectForm(&corrector);
  corrector.invert();  // This could potentially allocate for the vector swap. Figure this out!
  corrector.reflectZeros();
}


//=================================================================================================

/** A nonlinear extension of rsDampedCombAllpass. At the moment, it's just an experimental stub. */

template<class TSig, class TPar>
class rsDampedCombAllpassNonLin : public rsDampedCombAllpass<TSig, TPar>
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


//=================================================================================================


/** Implements a Schroeder allpass filters with a frequency dependent feedback, i.e. with a filter
in the feedback path. For technical reasons, this feedback filter must be an FIR filter (otherwise
the allpass condition would imply instability, I think). 

...TBC...Elaborate! So far, I think, the filter in the feedforward path must used reversed coeff 
arrays compared to the one in the feedback path. Because this changes the stability for IIR 
filters, we may only be able to use FIR filters. */


template<class TSig, class TPar>
class rsDampedSchroederAllpass
{


public:

  void setMaxDelayInSamples(int newMaxDelay)
  {
    delayLine.setMaxDelayInSamples(newMaxDelay);
  }

  void setMaxDampOrder(int newMaxDampOrder)
  {
    b.reserve(newMaxDampOrder+1);
  }
  // May be called before setup to pre-allocate the memory for the damping filter coeffs to avoid
  // memory allocations in setup()
 
  /** Sets up the delay, feedback gain the damping filter. The damping filter must be an FIR 
  filter and the caller is supposed to pass its coefficients and order. The order is the length of
  the coefficient array plus one. */
  void setup(int delay, TPar feedback, int dampOrder, const TPar* dampCoeffsB)
  {
    delayLine.setDelayInSamples(delay);
    M = delay;
    k = feedback;
    b.resize(dampOrder+1);
    for(int k = 0;  k <= dampOrder; k++)
      b[k] = dampCoeffsB[k];
  }


  /** Evaluates the filter's z-domain transfer function H(z) value at the given value of z. It is 
  given by:

            d^M + sum_{i=0}^P  k * b[i] * d^(i)
    H(z) = ---------------------------------------
             1  + sum_{i=0}^P  k * b[i] * d^(M-i)

  where M is the delay, P is the feedback filter order and d = z^-1 = 1/z. For example, with M = 5
  and P = 2, it would look like:

             d^5 + k*b0     + k*b1*d 
    H(z) =  ---------------------------
              1  + k*b1*d^4 + k*b0*d^5
  */
  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const
  {
    using Complex = rsComplex<TPar>;
    Complex d   = TPar(1) / z;                 // d       = z^-1
    Complex dM  = rsPow(d, Complex(M));        // dM      = z^-M = d^M
    Complex di  = TPar(1);                     // d^i     = z^-i
    Complex dMi = dM;                          // d^(M-i) = z^(-(M-i)) = z^(i-M)
    Complex num = TPar(0);
    Complex den = TPar(0);
    for(size_t i = 0; i < b.size(); i++)
    {
      num += b[i] * di;
      den += b[i] * dMi;
      di  *= d;
      dMi *= z;                                // Equivalent to: dMi /= d; Decrement exponent.
    }
    num = dM      + k*num;
    den = TPar(1) + k*den;
    return num / den;

    // ToDo:
    //
    // - Maybe use rsPowInt(d, M) to compute dM. But the function may not yet work for rsComplex.
    //   Figure that out!
  }



  TSig getSample(TSig in)
  {
    int order = (int) b.size()-1;

    // Apply feedback path:
    TSig tmp = in;
    for(int i = 0; i <= order; i++)
      tmp -= k * b[i] * delayLine.readOutputAt(M-i);
    delayLine.writeInputNoUpdate(tmp);

    // Apply feedforward path:
    tmp = 0;
    for(int i = 0; i <= order; i++)
      tmp += k * b[i] * delayLine.readOutputAt(i);
    tmp += delayLine.readOutputAt(M);

    // Update delayline and return result:
    delayLine.incrementTapPointers();
    return tmp;

    // In order to need only one delayline, we use a direct form 2 implementation. See:
    // https://www.dsprelated.com/freebooks/filters/Direct_Form_II.html

    // ToDo:
    //
    // - Try using size_t for the loop index i. Get rid of the local variable "order". Measure
    //   performance of both variants. Maybe the delayLine needs a member function readOutputAt
    //   that takes a size_t to avoid the conversion.
  }

  void reset()
  {
    delayLine.reset();
  }


protected:

  rsDelay<TSig> delayLine;
  std::vector<TPar> b;
  TPar k;
  int M = 0;

};

/** A naive version of rsDampedSchroederAllpass that implements the transfer function

            k * F(z) + z^-M 
  H(z) = ---------------------
          1 + k * F(z) * z^-M

hmmm but maybe we need to introduce unit delays in front of both k*F*... terms to make it 
implementable? ...noo - I think, that transfer function is wrong anyway. I think, it's more 
something like:

            k * R(z) + z^-M 
  H(z) = -------------------------
          1 + k * F(z) * z^-(M-P)

Where R is the reversal of F and P is the order of F. ...But this also needs verification. We 
assume that F(z) is purely FIR, by the way - because otherwise the whole thing would be unstable,
I think. Maybe try implementing the filter literally like the transfer function above. Use 2 FIR
filters F(z) and R(z) and an input delayline of length M and an output delayline of length M-P. 
Or maybe it's M-P+1 - or M-L where L is the length of the filter (i.e. number of coeffs) */

template<class TSig, class TPar>
class rsDampedSchroederAllpassNaive
{


public:

  void setMaxDelayInSamples(int newMaxDelay)
  {
    inDelay.setMaxDelayInSamples( newMaxDelay);
    outDelay.setMaxDelayInSamples(newMaxDelay);
  }

  void setup(int delay, TPar feedback, int dampOrder, const TPar* dampCoeffsB)
  {
    M = delay;
    k = feedback;
    inDelay.setDelayInSamples( delay);
    outDelay.setDelayInSamples(delay);
    b.resize(dampOrder+1);
    for(int k = 0;  k <= dampOrder; k++)
      b[k] = dampCoeffsB[k];
  }


  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const
  {
    // For a 2-point feedback filter with coeffs b0,b1 and with M = 5 (the delay), the transfer 
    // function should look like this:
    //
    //            k*b0 + k*b1*d + d^M         M=5   k*b0 + k*b1*d + d^5
    //  H(z) = -----------------------------   =   -------------------------
    //          1 + k*b1*d^(M-1) + k*b0*d^M         1 + k*b1*d^4 + k*b0*d^5

    using Complex = rsComplex<TPar>;
    Complex d   = TPar(1) / z;                      // d = z^-1
    Complex num = rsPow(d, Complex(M));
    Complex den = TPar(1);
    for(int i = 0; i < (int)b.size(); i++)
    {
      num += k * b[i] * rsPow(d, Complex(i  ));
      den += k * b[i] * rsPow(d, Complex(M-i));
    }
    return num / den;
  }


  TSig getSample(TSig in)
  {
    int order = (int) b.size()-1;

    // Apply feedforward part:
    inDelay.writeInputNoUpdate(in);
    TSig tmp = inDelay.readOutputAt(M);           // readOutput() should also work (more efficient)
    for(int i = 0; i <= order; i++)
      tmp += k * b[i] * inDelay.readOutputAt(i);

    // Apply feedback path:
    for(int i = 0; i <= order; i++)
      tmp -= k * b[i] * outDelay.readOutputAt(M-i);

    // Update delaylines and return result:
    inDelay.incrementTapPointers();
    outDelay.writeInputAndUpdate(tmp);
    return tmp;
  }

  void reset()
  {
    inDelay.reset();
    outDelay.reset();
  }


protected:

  rsDelay<TSig> inDelay;
  rsDelay<TSig> outDelay;
  std::vector<TPar> b;
  TPar k;
  int M = 0;

};



/** Under Construction.

A variant that implements the transfer function proposed above:

            k * F(z) + z^-M                    k * F(z) + z^-M
  H(z) = ------------------------- = ----------------------------------
          1 + k * R(z) * z^-(M-P)     1 + k * z-^1 * R(z) * z^-(M-P-1)

directly. Mainly to see, if this formula is actually correct. ...TBC...Verify formula  */

template<class TSig, class TPar>
class rsDampedSchroederAllpassNaive2
{


public:

  void setMaxDelayInSamples(int newMaxDelay)
  {
    inDelay.setMaxDelayInSamples( newMaxDelay);
    outDelay.setMaxDelayInSamples(newMaxDelay);
  }

  void setup(int delay, TPar feedback, int dampOrder, const TPar* dampCoeffsB)
  {
    M = delay;
    P = dampOrder;
    k = feedback;

    inDelay.setDelayInSamples( M);
    outDelay.setDelayInSamples(M-P-1);               // -1 compensates for implicit feedback delay
    outFilter.setImpulseResponse(dampCoeffsB, P+1);
    inFilter.setImpulseResponse( dampCoeffsB, P+1);
    outFilter.reverseImpulseResponse();

    // Maybe bake the scaler k into the filter coeffs. I think, instead of reversing the coeffs
    // for the outFilter (i.e. in the feedback path), we should use the given coeffs in the
    // feedback path as is and reverse them in the feedforward path. We need to do this in all the
    // variants then.
  }

  TSig getSample(TSig in)
  {
    TSig tmp;
    tmp = k*inFilter.getSample(in) + inDelay.getSample(in);       // Feedforward path
    out = tmp - k*outFilter.getSample(outDelay.getSample(out));   // Feedback path
    return out;
  }

  void reset()
  {
    inDelay.reset();
    outDelay.reset();
    inFilter.clearInputBuffer();   // Rename to reset
    outFilter.clearInputBuffer(); 
    out = 0;
  }


protected:

  rsDelay<TSig> inDelay;
  rsDelay<TSig> outDelay;

  rosic::ConvolverBruteForce inFilter;
  rosic::ConvolverBruteForce outFilter;

  TPar k;
  int M = 0;   // Delay in samples
  int P = 0;   // Order of feedback- and feedforward filter

  TSig out;

};






#endif