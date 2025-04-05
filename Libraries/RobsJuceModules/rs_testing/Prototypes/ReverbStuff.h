#ifndef RS_REVERBSTUFF_H
#define RS_REVERBSTUFF_H

// ToDo: 
// -Maybe move those classes that eventually go into the RAPT library to the top and those that are
//  really only prototypes (naive implementations that serve as reference for unit tests) to the 
//  bottom


/** Variant of the fast 2x2 Kronecker trafo that uses different coeffs for each stage of the 
transform. The coeffs a,b,c,d are replaced by arrays. They must all be of length log2(N). */
template<class T>
void rsStagedKroneckerTrafo2x2(T* A, int N, T* a, T* b, T* c, T* d)
{
  // UNDER CONSTRUCTION. Not yet tested. 

  rsAssert(rsIsPowerOfTwo(N), "N must be a power of 2");
  int h = 1;
  int L = 0;                 // Transform level
  while(h < N) 
  {
    for(int i = 0; i < N; i += 2*h) 
    {
      for(int j = i; j < i+h; j++) 
      {
        T x = A[j];
        T y = A[j+h];
        A[j]   = a[L]*x + b[L]*y;
        A[j+h] = c[L]*x + d[L]*y;
      }
    }
    h *= 2;
    L++;
  }
}




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
  // Maybe use a template parameter TArg for input and output. See rsDelay

  template<class TTol>
  void getTransferFunction(rsSparseDigitalTransferFunction<TPar, TTol>* tf) const
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


  /** Computes the transfer function H(z) of the delay at the given value of z. */
  template<class TArg>
  TArg getTransferFunctionAt(const TArg& z) const
  {
    int  M  = dl.getDelayInSamples();          // M is our delay
    TArg z1 = TArg(1) / z;                     // z^-1
    TArg zM = rsPow(z, TArg(-M));              // z^(-M)
    return (c*zM + z1*zM) / (TArg(1) + c*z1);  // H(z) = (c*z^(-M) + z^(-M-1)) / (1 + c*z^(-1))
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
    // oscillations at the Nyquist freq for certain settings. See the old implementations. The 
    // oscillation occurs when c is close to 1. This happens the the fractional part f is zero.
    // This is very unfortunate because it means that in the limit of the delay approaching an 
    // integer, the behavior of the allpass interpolated does not approach the one of a simpler 
    // integer delayline. Maybe we should treat f = 0 as special case? But then the question 
    // about an appropriate numeric tolerance arises. Maybe check if higher order Thiran allpass
    // interpolators have the same problem. ...although...maybe it's not actually a problem?
    // It seems that when the delay is an exact integer, the impulse response looks fine. But
    // slightly above an integer (like n + 0.01), we see a long ringing at the Nyquist freq.
    // Slighly below an integer, e.g. n + 0.99, there is not such ringing.

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


/*

Notes

- In general, interpolation is a process that takes in datapoints at discrete points, i.e. pairs
  (x,y) in 1D, and creating a continuous function from these points. It can be sued for various 
  tasks, among them implementing fractional delays and upsampling.

- In the case of fractional delay, interpolation schemes may also be interpreted as filters. The 
  interpolated value at some position n+d will be a linear combination of samples around index n.
  The coefficients of that linear combination depend on the desired fracctional delay d in [0,1)
  and are basically FIR filter coeffs. In this case, we may not care so much about the smoothness
  of the underlying continuous interpolant because we just use one particluar fractional offset 
  anyway (unless the fractional delay is time-varying). In this case, we may be mostly interested
  in the frequency response of our fractional delay filter. Ideally, we want the magnitude response
  to be unity ar all frequencies and the group delay should be equal to our desired fractional 
  delay at all frequencies (Q: why the group delay and not the phase delay?). Lagrange 
  interpolators provide a maximally flat group delay response at DC (see PASP book by Julius 
  Smith), so they might be suitable even though they are not as smooth as e.g. Hermite 
  interpolators of the same order. Thiran allpass interpolators can be seen as an IIR extension to
  the Lagrange interpolators which are FIR (see PASP as well).

- For upsampling, interpolation schemes may be interpreted as so called polyphase filters. For 
  example, to upsample by a factor of 5, we would need 5 different fractional delay filters, namely
  those that implement the delays of: 0.0, 0.2, 0.4, 0.6, 0.8. We would interleave their outputs in
  a round-robin fashion to obtain the upsampled signal. For upsampling, we may care about the 
  smoothness of the underlying continuous interpolant at the datapoints (i.e. where we switch to a
  new segment of the interpolant). This may be especially true if we want to create a nice smooth 
  plot. So, for upsampling for plots, something like Hermite interpolation may be appropriate. 
  Maybe also splines, but they are not suitable for realtime use because each segment depends on
  the whole dataset.

- Maybe it could make sense to design allpass interpolators according to a desired impulse 
  response. For a 1st roder allpass, we would make the ansatz:

     y[n] = c*x[n] + x[n-1] - c*y[n-1]   with   x[0,1,2,3,...] = 1,0,0,0,... and zero for n < 0
     h[0] = c
     h[1] = 1 - c*h[0] = 1 - c^2
     h[2] = c * h[1] = c*(1 - c^2)
     ...
     h[n] = c^(n-1) * (1 - c^2)

  As we have only one coeff, we could prescribe only the first sample of the impulse response. 
  Maybe we should set it to 1-f where f is the desired fractional delay in the interval [0,1). The
  rest would then follow. I don't know if that makes sense. For a 2nd order allpass, we could 
  prescribe two sample h[0], h[1]. Maybe we should set them to (1-f), f. Maybe try it and see what 
  happens. Maybe in the first order case, we could also require: h[0] + h[1] = 1. We want the first
  two samples to sum to the amplitude of the unit impulse. Or maybe the first two samples should 
  have that same eneryg as a unit impulse, i.e. (h[0])^2 + (h[1])^2 = 1? ...some experiments and
  creativity is needed!

*/



//=================================================================================================

/** UNDER CONSTRUCTION...

Implements a universal comb filter. Depending on the coefficients, this filter can be used as
feedforward comb, feedback comb, Schroeder allpass comb, notchpass filter, etc. ...TBC... 

               Feedfwd   Feedback   Blend     Delay    ModDepth     ModShape   
FIR Comb:      g         0          1 or ?
IIR Comb:      0         g=-1..1    1 or c
Allpass:       1         g=-1..1    -g
Delay:         1         0          0
Notchpass:  
Slapback:

The c in the IIR comb can be 1-|g| or sqrt(1-g^2). The former choice ensures peak gain of 1 and 
the latter normalizes the loudness for broadband signals. See (1) pg 70. 

References:

  (1) DAFX (1st Ed., Udo Zoelzer), page 66.
      https://www.dafx.de/DAFX_Book_Page/chapter3.html
      https://www.dafx.de/DAFX_Book_Page_2nd_edition/chapter2.html

*/

template<class TSig, class TPar>
class rsUniversalCombFilter
{

public:

  //-----------------------------------------------------------------------------------------------
  /** \name Setup */

  void setMaxDelayInSamples(int newMaxDelay) { delayLine.setMaxDelayInSamples(newMaxDelay); }

  void setDelayInSamples(int newDelay) { delayLine.setDelayInSamples(newDelay); }

  void setCoeffs(TPar newBlend, TPar newFeedforward, TPar newFeedback)
  {
    bl = newBlend;
    ff = newFeedforward;
    fb = newFeedback;
  }

  void setToAllpass(TPar newAllpassCoeff) 
  { 
    bl =  newCoeff;
    ff =  TPar(1);
    fb = -newCoeff;
  }
  // Needs tests. The filter should be equivalent to rsAllpassDelay with this setting.

  void setToFeedbackComb(TPar newFeedbackCoeff)
  {
    bl = TPar(1);
    ff = TPar(0);
    fb = newFeedbackCoeff;
  }
  // ToDo: Maybe allow the user to set a global gain. This should be assigned to the blend coeff.
  // Make it an optional parameter defaulting to 1.

  void setToFeedforwardComb(TPar newFeedforwardCoeff, TPar newBlendCoeff = TPar(1))
  {
    bl = newBlendCoeff;
    ff = newFeedforwardCoeff;
    fb = TPar(0);
  }
  // ToDo: explain how the ff, bl coeffs relate to a total gain. I think, to achieve a different 
  // overall total gain, we should scale ff and bl by that gain factor

  void setToPureDelay()
  {
    bl = TPar(0);
    ff = TPar(1);
    fb = TPar(0);
  }

  void setToBypass()
  {
    bl = TPar(1);
    ff = TPar(0);
    fb = TPar(0);
  }

  void setToMuted()
  {
    bl = TPar(0);
    ff = TPar(0);
    fb = TPar(0);
  }


  /*
  void setToNotchpass(TPar newPoleCoeff, TPar newZeroCoeff)
  {

  }
  */



  //-----------------------------------------------------------------------------------------------
  /** \name Inquiry */

  /** Returns the amount of delay in samples. */
  int getDelayInSamples() const
  {
    return delayLine.getDelayInSamples();
  }

  /** Returns value of the transfer function H(z) of this filter at the given z. The transfer 
  function is given by:  H(z) = (bl + ff * z^-M) / (1 - fb * z^-M)  where M is the delay and 
  bl, ff, fb are the blend, feedforward and feedback coefficient respectively. */
  template<class TArg>
  TArg getTransferFunctionAt(const TArg& z) const
  {
    int  M  = getDelayInSamples();
    TArg zM = rsPow(z, TArg(-M));              // z^-M
    TArg V  = TArg(1) / (TArg(1) + fb * zM);   // V(z), z-trafo of intermediate signal v[n].
    return bl * V + ff * V * zM;
  }

  /** Assigns the passed tf pointer to our transfer function. We don't use a return value for the
  result to enable pre-allocation of the function object which is important in realtime contexts.*/
  template<class TTol>
  void getTransferFunction(rsSparseDigitalTransferFunction<TPar, TTol>* tf) const
  {
    int M = getDelayInSamples();
    tf->initToZero();
    tf->getNumerator().  _appendTerm(bl, 0);
    tf->getNumerator().  _appendTerm(ff, M);
    tf->getDenominator()._appendTerm(fb, M);

    // We don't use -M because the minus is already baked into the class. It interprets the 
    // function as a rational function in z^-1.
  }


  //-----------------------------------------------------------------------------------------------
  /** \name Processing */

  inline TSig getSample(TSig x)
  {
    TSig vM = delayLine.readOutput();    // Read vM = v[n-M] from the delayline.
    TSig v  = x - fb * vM;               // Intermediate signal: v[n] = x[n] + fb * v[n-M].
    delayLine.writeInputAndUpdate(v);    // Write v[n] into the delayline.
    return bl * v + ff * vM;             // Return y[n] = bl * v[n] + ff * v[n-M].
  }
  // I think, this corresponds to a 1st order filter in DF2:
  //
  //   v[n] = x[n] - a1*v[n-1]
  //   y[n] = b0*v[n] + b1*v[n-1]
  //
  // which in DF1 would be:
  //
  //   y[n] = b0*x[n] + b1*x[n-1] - a1*y[n-1]
  //
  // where the unit delay has been replaced by a delay of M samples. The correspondence between the
  // coeffs is: b0 = bl, b1 = ff, a1 = -fb. Verify this theoretically and numerically! Set up a 
  // rsUniversalCombFilter filter and a rsUniversalCombFilter with delay of 1 and compare the 
  // outputs.

  void reset() { delayLine.reset(); }


protected:

  RAPT::rsDelay<TSig> delayLine;

  TPar bl = TPar(0);   // Blend coeff          maybe rename to b0
  TPar ff = TPar(0);   // Feedforward coeff    maybe rename to b1
  TPar fb = TPar(0);   // Feedback coeff       maybe rename to a1

};



// ToDo:
//
// - Maybe make a class rsNotchpassDelay that combines a normal delayline with a notchpass filter.
//   It may have a setTotalDelay(int) function that sets the delay of both and a function 
//   setDelayRatio(double) that sets the ratio between the lengths. The default should be 0.8 
//   meaning that 80% of the total delay is allocated to the normal delayline and 20% to the 
//   notchpass. That makes the notchpass delay 25% of the normal delay as recommended in Blesser's
//   patent. To implement this, it may make sense to factor out a delayline class that doesn't own 
//   the delay memory and instead just gets a pointer passed in and the memory is managed by the 
//   class that uses the delayline.
//
// - Implement a function getTransferFunction that returns the transfer function as 
//   rsSparseRationalFunction object.
//
// - Maybe implement a function getRingOutTime

// See also:
// https://github.com/isaiahdoyle/universalcombfilter
// https://en.wikipedia.org/wiki/Comb_filter
// https://ccrma.stanford.edu/~jos/pasp/Comb_Filters.html
// https://www.uncini.com/dida/tsa/mod_tsa/Chap_05_special_filters.pdf


//=================================================================================================

/** This implements a chain (i.e. series connection) of allpass delays. To achieve this effect, you
could just use a std::vector of rsAllpassDelay which are applied one after the other. This class 
here is a convenience class that does this for you. Series connections of allpass filters can be 
used for allpass diffusors, for example, as building blocks of a reverb algorithm.

Hmm...maybe it's actually not such a great idea to provide such a convenience class. If we do this,
we may want to have similar convenience classes for other types of allpass filter chains which 
would look very similar - i.e. a lot code duplication and boilerplate. I actually had this class
already in the RAPT library but backed off again and moved it back into the prototypes for this 
reason. On the other hand, chains of Schroeder allpasses are common building block in reverb 
algorithms, so it might be convenient to have that. We'll see...

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


/** Creates a biquad filter made from 2 1st order shelving filters (one low-shelf, one high-shelf)
that can be used as damping filter inside a Karplus-Strong like comb filter algorithm. The b,a, 
arrays are output parameters for the biquad coeffs and  must be (at least) of length 3. The return
value is the overall feedback gain. The parameters are: delay: the length of the delayline, decay: 
the desired decay time (both must be in the same unit, e.g. samples), lo/hiOmega: the corner 
frequencies of the shelvers (as normalized radian frequency omega = 2*pi*f/fs) and lo/hiScale: the 
desired time scaling factors for the decay time for low and high frequenceis. ...TBC...*/
template<class T>
T rsMakeDampBiShelf(T delay, T decay, T loOmega, T loScale, T hiOmega, T hiScale,
  T* b, T* a)
{
  // Compute desired feedback gains for low, mid and high frequencies:
  T a60 = T(0.001);   // = rsDbToAmp(-60.0). Target amplitude to reach after decay (in samples)
  T kL  = rsDecayTimeToFeedbackGain(decay * loScale, delay, a60);
  T kM  = rsDecayTimeToFeedbackGain(decay          , delay, a60);
  T kH  = rsDecayTimeToFeedbackGain(decay * hiScale, delay, a60);
  // These formulas could also be expressed as e.g.:
  //
  //   kM = rsPow(10.0, TPar(-3 * delay) / decay);
  //
  // which is how they are often seen in the FDN literature. 

  // Compute desired gains for the low and high shelver:
  T gL = kL / kM;
  T gH = kH / kM;

  // Compute coeffs for low- and high shelver:
  T aL[2], bL[2]; aL[0] = 1; rsMake1stOrderLowShelf( loOmega, gL, &bL[0], &bL[1], &aL[1]);
  T aH[2], bH[2]; aH[0] = 1; rsMake1stOrderHighShelf(hiOmega, gH, &bH[0], &bH[1], &aH[1]);

  // Combine low- and high shelver into biquad:
  rsArrayTools::convolve(aL, 2, aH, 2, a);
  rsArrayTools::convolve(bL, 2, bH, 2, b);

  // Return the (mid) feedback gain:
  return kM;


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
  // - Can be optimized: design the low shelf directly into a,b, the high shelf into some tmp 
  //   arrays, then bake them into a,b. The allpass can then use the same temp arrays as the 
  //   hi shelf and then also bake them into a,b

}
// Rename to rsMakeFeedbackFilter...maybe


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

/** Under construction - This shopuld eventually go into a pair of files DampedCombFilter.h/cpp in
RAPT/Filters/Musical

This is supposed to factor out some functionality from class rsDampedCombAllpass to facilitate 
re-using it in e.g. rsDampedMultiCombAllpass. We want to decouple the multicomb allpass from the
single comb allpass - mainly because keeping the coupling would require to rewrite some 
functionality in the single comb case to support fractional delays. But that's a complication, I 
don't want to introduce to this class. The fractional delay feature shall be reserved for the
multicomb. We mainly want to factor out all the getTransferFunction stuff that creates the transfer
function objects. ...TBC...  

ToDo: explain intention behind the template parameters TCoef and TDly. TCoef is the type for the
feedback gain and feedback filter coeffs and TDly for the delay. For the former, it can make sense
to have a simd type or maybe even a complex type. For the latter, we can only have scalar real 
number types. It would actually be desirable to have simd vector types for TDly, too - but that's
difficult to implement. The amount of delay in a delayline is not so easily simdified. 

But maybe it can done: With interpolation and damping, we typically do more delayline readouts per 
sample than we do writes. If the read-pointers ("tapIn") are all in-sync but the write pointers 
("tapOut") have different offsets, we could make a meaningful implementation. Try to make one with
rsFloat32x4 for the signal and some rsInt32x4 type (to be written) for the index. The write code 
would extract and use the scalar integer delays stored in the rsInt32x4 vector, i.e. just look like
regular delayline code repeated 4 times. As the write operation is rare, it doesn't matter too much
that the delayline writing can't take advantage of vectorization. The reading code would assume 
that all 4 ints in the vector have the same value and just read out the 4-float vector at that 
index. We could (and probably should) actually even use a scalar int for the read index. If that 
works out, try to templatize it. We may need one template parameter for the signal (e.g. 
rsFloat32x4) and one for the vector index (e.g. rsInt32x4). Maybe we need one for the scalar 
index, too - but maybe we can get away without it.  */

template<class TCoef, class TDly, class TTol>
class rsDampedCombSettings                // Maybe rename to rsDampedCombParams
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  rsDampedCombSettings()
  {
    init();
  }



  //-----------------------------------------------------------------------------------------------
  // \name Setup

   
  enum class InterpolationMode
  {
    sampleHold,     // Sample and hold, truncate/floor read position to int
    nearest,        // Nearest neighbor interpolation, round read position to nearest int
    linear,         // Linear interpolation, connect samples with straight lines
    allpass1        // First order (warped) allpass interpolation (maybe rename to thiran1)
  };
  // Maybe use unsigned char as underlying type for the enum. Maybe offer more interpolation modes
  // like cubic Hermite, cubic Lagrange, 2nd and 3rd order Thiran allpass, etc. See:
  // https://ccrma.stanford.edu/~jos/pasp/Thiran_Allpass_Interpolators.html
  // http://users.spa.aalto.fi/vpv/publications/vesan_vaitos/ch3_pt3_allpass.pdf
  // ...is part of: http://users.spa.aalto.fi/vpv/publications/vesan_vaitos/
  //
  // Maybe add Lagrange and Hermite interpolators. Maybe they can be turned into allpass
  // interpolators by just using the reversed FIR coefficient array for the recursive part? Will 
  // that give meaningful interpolators (i.e. stable and with desirable group delay and/or phase 
  // delay characteristics, etc.)? Maybe implement general functions to compute coeff arrays of
  // N-th order Lagrange, Hermite, Thiran, etc. interpolators
  


  //enum class DampingMode  
  //{
  //  feedbackDamp,        // F(z) is in feedback path, A(z) in feedforward path
  //  forwardDamp,         // F(z) is in feedforward path, A(z) in feedback path
  //  forwardDampComp      // Like forwardDamp but with 1/F(z) compensator in series
  //};
  // Maybe make a mode where both delay and damper are in the feedback path. I think, it should 
  // realize the same transfer function as the forwardDampComp mode but without the need for an
  // inverse damper, so it would be more efficient. Verify this hypothesis theoretically and 
  // numerically. Maybe call that mode feedbackDampDelay
  // Maybe rename to structure, topology, configuration, ...
  //
  // Currently, we only have a boolean preDelay flag. When true, it corresponds to the feedbackDamp
  // case and when false to the forwardDampComp case. That's what we assume in 
  // getCombTransferFunction() 
 

  void init()
  {
    delay         = TDly(0);
    interpolation = InterpolationMode::sampleHold;
    k             = TCoef(0);
    dmpOrd        = 0;
    preDelay      = false;

    using AT = rsArrayTools;
    AT::clear(bD, maxDmpOrd+1); bD[0] = TCoef(1);
    AT::clear(aD, maxDmpOrd+1); aD[0] = TCoef(1);
    AT::clear(bI, maxIntOrd+1); bI[0] = TDly(1);
    AT::clear(aI, maxIntOrd+1); aI[0] = TDly(1);
  }


  void setup(TDly delayInSamples, InterpolationMode interpolationMode,
    TCoef feedback, int dampOrder, const TCoef* dampCoeffsB, const TCoef* dampCoeffsA, 
    bool preDelayMode)
  {
    if(dampOrder > maxDmpOrd) 
    {
      rsError("Such high damping order is not supported.");
      init();
      return;
    }

    delay         = delayInSamples - TDly(1);      // -1 corrects for unit delay in feedback path
    interpolation = interpolationMode;
    k             = feedback;
    preDelay      = preDelayMode;
    dmpOrd        = dampOrder;

    rsAssert(dampCoeffsA[0] == TCoef(1));
    rsArrayTools::copy(dampCoeffsA, aD, dmpOrd+1);
    rsArrayTools::copy(dampCoeffsB, bD, dmpOrd+1);
    //normalizeFilterCoeffs(bD, dmpOrd+1, aD, dmpOrd+1);
      // Maybe call such a function later to relax the dampCoeffsA[0] == 1 assertion. It should 
      // divide all coeffs by aD[0] to normalize the transfer function. 

    updateInterpolatorCoeffs();                    // Assigns the coeffs in the bI, aI arrays
  }
  // Maybe rename to setupFromAlgoParams and have a similar setupFromUserParams function that uses
  // higher level parameters such as decay times at various frequencies, i.e. the currently free
  // function rsSetupDecayTimes_LinViaDly ...Maybe setupViaDecayTimes...the low level function 
  // could be called setupViaCoeffs


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry


  template<class TArg>
  TArg getDamperTransferFunctionAt(const TArg& z) const
  {
    TArg num(0), den(0);
    for(int i = 0; i <= dmpOrd; i++)
    {
      TArg zi = rsPow(z, TArg(-i));          // z^-i
      num += bD[i] * zi;
      den += aD[i] * zi;
    }
    return num / den;

    // ToDo:
    //
    // - This should be optimized (don't call rsPow - compute the powers on the fly by multiplying 
    //   an accumulator by z1 = 1/z, initialized as 1). Maybe it can be factored out into a library
    //   function to compute the transfer function of direct form filters. Maybe it should go into 
    //   rsFilterAnalyzer.
  }
  // We use a template parameter TArg so we don't have to commit to decide between TCoef and TDly 
  // here and we would also allow client code to use std::complex or rsComplex


  /** Multiplies the given transfer function by the transfer function of our delayline. This 
  consists of a factor z^-M for the integer delay of M samples and a factor resulting from the 
  interpolator, for example b0 + b1*z^-1 with b0 = 1-f, b1 = f for the linear interpolator with f 
  being the fractional part of the delay, i.e. delay = M+f. */
  void mulByDelayTransFunc(rsSparseDigitalTransferFunction<TCoef, TTol>* tf) const
  {
    tf->multiplyByDenseCoeffs(bI, intNumOrd+1, aI, intDenOrd+1, TCoef(0));  // Interpolator factor
    tf->addPreDelay((int)delay);                                            // Integer delay factor

    // VERIFY if this is correct! Check and document also, if the order of the calls matters. I 
    // think, it shouldn't. Test it with all the available interpolators.
  }
  // Maybe make protected - it's currently used only internally. But maybe it could be useful for 
  // extenal use, too? Implement also mulByDampTransFunc


  void getCombTransferFunction(rsSparseDigitalTransferFunction<TCoef, TTol>* tf) const
  {
    using Mon = rsMonomial<TCoef>;
    getDamperTransferFunction(tf);        // tf = F, F(z) is transfer function in feedback path
    tf->multiplyBy(Mon(k, 1));            // tf = F * k * z^-1
    mulByDelayTransFunc(tf);              // tf = F * k * z^-1 * A
    //tf->addConstant(TCoef(1), TCoef(0));  // tf = 1 + F * k * z^-1 * A   old - with tol param
    tf->addConstant(TCoef(1));            // tf = 1 + F * k * z^-1 * A
    tf->invert();                         // tf = 1 / (1 + F * k * z^-1 * A)
    if(preDelay)
      mulByDelayTransFunc(tf);            // tf = A / (1 + F * k * z^-1 * A)

    // ToDo:
    //
    // - Verify that all operations above are non-allocating (assuming that tf has enough capacity)
    //   and document that fact.
  }
  // Rename to getTransferFunction


  void getDamperTransferFunction(rsSparseDigitalTransferFunction<TCoef, TTol>* tf) const
  {
    tf->setupFromDenseCoeffs(bD, dmpOrd+1, aD, dmpOrd+1, TCoef(0));

    // ToDo: here and elsewhere in similar calls, do not pass TCoeff(0) for the roundoff error
    // tolerance. Try to come up with a sensible value. Maybe implement a (protected) member 
    // function getExpectedTransFuncRoundoffError() or something like that. It may take into 
    // account the complexity of the filter - or at least provide the infrastructure for doing that
    // later.
  }


  void getDelayTransferFunction(rsSparseDigitalTransferFunction<TCoef, TTol>* tf) const
  {
    tf->setupFromDenseCoeffs(bI, intNumOrd+1, aI, intDenOrd+1, TCoef(0));
    tf->addPreDelay((int)delay);  // VERIFY!
  }
  // Needs more tests with all the different interpolation modes


  /** Returns the maximum possible order of the damping filter F(z). */
  static constexpr int getMaxDampingOrder() { return maxDmpOrd; }

  /** An interpolation scheme that implements a fractional delay can be viewed as a filter. For 
  example, cubic Lagrange interpolation can be viewed as a 3rd order FIR filter and Thiran 
  interpolation uses IIR allpass filters. This function returns the maximum possible order of this 
  interpolator filter. */
  static constexpr int getMaxInterpolationOrder() { return maxIntOrd; }


  static constexpr int getMaxIntPlusDampOrder() { return maxIntOrd + maxDmpOrd; }
  // This is the maximum order increase caused by interpolation and damping taken together. That
  // means, when the maximum nominal delay in samples is M, then the maximum total order of the
  // damped comb filter is M + getMaxIntPlusDampOrder() + 1. The +1 comes from the additional 
  // unit delay in the feedback loop.  ToDo: Verify this and add to documentation!


  /** Returns the maximum possible total order of the full damped comb filter for a given maximum
  delayline length. This order includes delay, interpolation, damping and the implicit unit delay 
  in the feedback loop. */
  int getMaxTotalOrder(int maxDelay) { return maxDelay + getMaxIntPlusDampOrder() + 1; }
  // Verify!






  TDly getDelay() const { return delay; }


  TCoef getFeedbackGain() const { return k; }

  int getDampingOrder() const { return dmpOrd; }

  const TCoef* getDampCoeffsB() const { return bD; }

  const TCoef* getDampCoeffsA() const { return aD; }

  bool isInPreDelayMode() const { return preDelay; }

  int getIntDelay() const { return (int) delay; }

  //int getInterpolationOrder() const { return rsMax(intNumOrd, intDenOrd); }

  //int getTotalOrder() const
  //{  return getIntDelay() + getDampingOrder() + getInterpolationOrder() + 1; }

  // Verify!

  // Maybe we can have more specific functions that return the total numerator and denominator 
  // order separately. I think, this will also depend on the mode/topology. 


  //-----------------------------------------------------------------------------------------------
  // \name Processing. These functions facilitate to implement the actual DSP but leave some 
  // related responsibilities (such as managing the filter states) to the client code. This is a 
  // bit odd and may be refactored later.


  /** Applies the feedback damping filter to the signal "in" and updates the given filter's state 
  x = x[n-1], x[n-2], ... y = y[n-1], y[n-2], ... . ...TBC... */
  template<class TSig>
  inline TSig applyDamper(TSig in, TSig* x, TSig* y)
  {
    // Compute outputs:
    TSig out = bD[0]*in;
    for(int i = 1; i <= dmpOrd; i++)
      out += bD[i] * x[i-1] - aD[i] * y[i-1];

    // Update state and return result:
    rsArrayTools::shiftPushDiscard(x, dmpOrd, in);
    rsArrayTools::shiftPushDiscard(y, dmpOrd, out);
    return out;
  }

  /** Applies the inverse feedback damping filter to the signal x and updates the filter's state.
  this filter is need only the "without predelay" mode of operation. */
  template<class TSig>
  inline TSig applyInverseDamper(TSig in, TSig* x, TSig* y)
  {
    // Compute output:
    TSig out = in;
    for(int i = 1; i <= dmpOrd; i++)
      out += aD[i] * x[i-1] - bD[i] * y[i-1];
    out /= bD[0];                                        // ToDo: maybe precompute 1/b[0]

    // Update state and return result:
    rsArrayTools::shiftPushDiscard(x, dmpOrd, in);
    rsArrayTools::shiftPushDiscard(y, dmpOrd, out);
    return out;
  }

  /** Applies only the poles of the damping filter. */
  template<class TSig>
  inline TSig applyDamperPoles(TSig in, TSig* y)
  {
    // Compute output:
    TSig out = in;
    for(int i = 1; i <= dmpOrd; i++)
      out -= aD[i] * y[i-1];

    // Update state and return result:
    rsArrayTools::shiftPushDiscard(y, dmpOrd, out);
    return out;
  }
  // ToDo: explain where this is needed



  // Maybe factor these out into free functions 
  // rsApplyDirectForm(Inverse)Filter(in, b, x, a, y, order), 
  // rsApplyDirectFormAllpole(in, a, y, order), 
  // Or maybe make them static members of rsDirectFormFilter. Maybe they should be named
  // applyDF1, applyInverseDF1, applyPolesDF1 - the idea being that there could also be 
  // implementations for other (direct) forms.
 






protected:

  /** Updates the arrays of the coefficients for the interpolator filter according to the desired 
  interpolation mode and the fractional part of the delay. */
  void updateInterpolatorCoeffs()
  {
    using IM = InterpolationMode;
    TDly f = delay - rsFloor(delay);        // Fractional part of delay

    switch(interpolation)
    {

    case IM::sampleHold:
    {
      // y[n] = x[n]:
      intNumOrd = 0; bI[0] = TDly(1);
      intDenOrd = 0; aI[0] = TDly(1);
    }
    break;

    case IM::nearest:
    {
      // y[n] = x[n] or x[n-1], depending on f:
      if(f <= TDly(0.5)) { intNumOrd = 0; bI[0] = TDly(1);                  }  // y[n] = x[n]
      else               { intNumOrd = 1; bI[0] = TDly(0); bI[1] = TDly(1); }  // y[n] = x[n-1]
      intDenOrd = 0; aI[0] = TDly(1);
    }
    break;
    // Needs tests

    case IM::linear:
    {
      // y[n] = (1-f)*x[n] + f*x[n-1]:
      intNumOrd = 1; bI[0] = TDly(1)-f; bI[1] = f;
      intDenOrd = 0; aI[0] = TDly(1);
    }
    break;

    case IM::allpass1:
    {
      // y[n] = c*x[n] + x[n-1] - c*y[n-1]  with  c = (1-f) / (1+f):
      TDly c = (TDly(1)-f) / (TDly(1)+f);
      intNumOrd = 1; bI[0] = c;       bI[1] = TDly(1);
      intDenOrd = 1; aI[0] = TDly(1); aI[1] = c;
    }
    break;

    default:
    {
      rsError("Unknown interpolation method.");

      // Use nearest neighbor method in that case:
      intNumOrd = 0; bI[0] = TDly(1);
      intDenOrd = 0; aI[0] = TDly(1);
      // Or maybe we should just output a zero signal by setting bI[0] to zero? aI[0] should 
      // remain 1, though.
    }

    }
    

  }
  // Maybe this should be factored out into a free function rsInterpolatorCoeffs or into a static
  // member function like calcCoeffs() of class rsInterpolator. It should take the fractional delay
  // f, the desired mode and a pointer to the coefficient arrays that will be filled by the 
  // function. We want to call it here like calcCoeffs(f, interpolation, bI, aI). Maybe there 
  // should be a function for FIR interpolators that should not take an a-array. It could be 
  // invoked by the more general function that takes both arrays and fill the a-array with 
  // 1,0,0,0,...


  static const int maxDmpOrd = 8;  // Maximum damping order
  static const int maxIntOrd = 1;  // Maximum interpolation order - ToDo: allow higher orders!

  // Coefficients:
  TCoef k = 0;                     // Feedback gain
  TCoef bD[maxDmpOrd+1];           // Damping filter feedforward coeffs
  TCoef aD[maxDmpOrd+1];           // Damping filter feedback coeffs
  TDly  bI[maxIntOrd+1];           // Interpolation filter feedforward coeffs
  TDly  aI[maxIntOrd+1];           // Interpolation filter feedback coeffs
  // ToDo: Explain why it makes sense to let aI, bI be of type TDly rather than TCoef. ...I'm 
  // actually not quite sure if that is really the right thing to do. Maybe they should be type
  // TCoef. Figure this out - imagine (or better: implement and test) situations where TCoef is a 
  // complex type or a simd type. I actually think, in the simd case, it would make more sense to
  // let the coeffs of the interpolator be of the vector rather than scalar type which amounts to
  // making the interpolator coeffs of type TCoef. Hmm...but maybe then the TDly should also be of
  // the vector type anyway...not sure...


  // Other settings:
  TDly delay  = 0;                 // Delay in samples (may be non integer)
  int  dmpOrd = 0;                 // Feedback damping filter order ToDo: have dmpNumOrd,dmpDenOrd
  int  intNumOrd = 0;              // Interpolator numerator order
  int  intDenOrd = 0;              // Interpolator denominator order


  // Switch between different interpolation modes:
  InterpolationMode interpolation = InterpolationMode::nearest;

  // Switch between with/without predelay mode of operation:
  bool preDelay = false;
    // Maybe rename to something like "mode" or "structure", "configuration", "topology". Maybe 
    // there could be even more modes?


};



// Maybe make that a member of rsDampedCombSettings:
template<class TCoef, class TDly, class TTol>
void rsSetupDecayTimes_LinViaDly(
  rsDampedCombSettings<TCoef, TDly, TTol>& combSettings, 
  TDly delay, TCoef decay, TCoef loOmega, TCoef loScale, TCoef hiOmega, TCoef hiScale, 
  bool preDelay)
{
  // Compute feedback gain and filter coeffs:
  TCoef a[3], b[3];
  TCoef kM = rsMakeDampBiShelf(delay, decay, loOmega, loScale, hiOmega, hiScale, b, a);

  // Set up the rsDampedCombSettings objects:
  using IM = rsDampedCombSettings<TCoef, TDly, TTol>::InterpolationMode;
  combSettings.setup(delay, IM::linear, kM, 2, b, a, preDelay);

  // ToDo:
  //
  // - Pass combSettings by pointer
  //
  // - Maybe the delay should be adjusted to take into account the effect of the damping filter 
  //   (which itself may also introduce a frequency dependent delay). I think, what we want is to
  //   have the correct fractional delay at DC or maybe at the resonance frequency, so we can tune 
  //   it exactly.
}



// Just for comparison/proof-of-concept. We bake the interpolation method not into the delayline 
// but into the feedback damping filter. That can be done only for FIR interpolators, though. In 
// this function, we hardcode the linear interpolation into it. But even in that case, it's 
// questionable if it should be done this way. But we implement it here to verify that it can be 
// done. Well - maybe an IIR interpolator could also be baked into the damping filter when the
// damping filters sits *before* the delay in the topology?
template<class TCoef, class TDly, class TTol>
void rsSetupDecayTimes_LinViaFb(rsDampedCombSettings<TCoef, TDly, TTol>& combSettings, 
  TDly delay, TCoef decay, TCoef loOmega, TCoef loScale, TCoef hiOmega, TCoef hiScale, 
  bool predelay)
{
  // Compute feedback gain and filter coeffs:
  TCoef a[4], b[4];
  TCoef kM = rsMakeDampBiShelf(delay, decay, loOmega, loScale, hiOmega, hiScale, b, a);

  // Possibly also bake an interpolation filter into the feedback filter to achieve fractional 
  // delay times:
  using IM = rsDampedCombSettings<TCoef, TDly, TTol>::InterpolationMode;
  TDly delayInt  = rsFloor(delay);
  TDly delayFrac = delay - delayInt;

  if(delayFrac == TDly(0))
  {
    // In the integer delay case, we only need the 2nd order feedback filter that we already have:
    combSettings.setup(delayInt, IM::nearest, kM, 2, b, a, predelay);
  }
  else
  {
    // In the fractional delay case, we create a linear interpolation filter and bake it into the 
    // existing 2nd order feedback filter, thereby turning it into a 3rd order filter:

    TCoef d = TCoef(delayFrac);

    // Design linear interpolation filter:
    TCoef aI[2], bI[2];
    bI[0] = TCoef(1) - d;
    bI[1] = d;
    aI[0] = TCoef(1);
    aI[1] = 0;

    // Bake the interpolation filter into the b,a feedback damping filter arrays:
    rsArrayTools::convolve(a, 3, aI, 2, a);
    rsArrayTools::convolve(b, 3, bI, 2, b);

    // Set up the delay settings using only the integer part for the delay and trivial sample and 
    // hold interpolation:
    combSettings.setup(delayInt, IM::sampleHold, kM, 3, b, a, predelay);
    //combSettings.setup(delayInt, IM::nearest, kM, 3, b, a, predelay);  // WRONG?!
  }
}
// ToDo: create tests, creating a linearly interpolating damped comb in 2 ways: (1) baking the
// interpolation into the delay, (2) baking the interpolation into the feedback damper. Compare
// the resulting transfer functions.



//=================================================================================================

/** Under Construction. We want to implement a filter based on the following block diagram:


  X(z) ---> + ---> A(z) ----------> Y(z)
            ^                |
            |                |
           -k               z^-1
            |                |
            ------ F(z) <-----

We want to factor out the pure comb filter from rsDampedCombAllpass because I think, it makes sense
to have such a comb filter in it own right. It can be used for Karplus-Strong like plucked string 
synthesis. rsDampedCombAllpass can then either be a subclass of rsDampedCombFilter or have an
object of that class as member. */




template<class TSig, class TPar, class TDly, class TTol>
class rsDampedCombFilter
{


public:

  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  rsDampedCombFilter() { reset(); }




  //-----------------------------------------------------------------------------------------------
  // \name Setup


  /** Initializes all settings to default values. */
  void initSettings();

  /** Sets the maximum desired roundtrip delay around the comb. This total roundtrip delay includes
  the z^-1 unit delay, so the delayline length is actually shorter by one. */
  void setMaxIntDelayInSamples(int newMaxDelay);


  void setMaxDelayInSamples(TDly newMaxDelay)
  { setMaxIntDelayInSamples((int)rsCeil(newMaxDelay)); }


  void setup(TDly delay, TPar feedback, 
    int dampOrder, const TPar* dampCoeffsB, const TPar* dampCoeffsA, bool predelayMode);





  //-----------------------------------------------------------------------------------------------
  // \name Inquiry






  //-----------------------------------------------------------------------------------------------
  // \name Processing

  /** Resets the state. */
  void reset();



protected:

  void updateDelay();

  using Settings = rsDampedCombSettings<TPar, TDly, TTol>;


  // Maximum order of feedback damping filter:
  static const int maxDmpOrd = rsDampedCombSettings<TPar, TDly, TTol>::getMaxDampingOrder(); 

  // Embedded DSP objects:
  rsDelay<TSig> mainDelay;                // Main delayline for the comb filter

  // State:
  TSig combOut = TSig(0);                 // State for the unit delay feedback loop
  TSig xd[maxDmpOrd], yd[maxDmpOrd];      // State for the damping filter
  TSig xi[maxDmpOrd], yi[maxDmpOrd];      // State for the inverse damping filter

  // Settings:
  Settings s;                             // Rename this! ...maybe to settings, params
  int M = 0;                              // Delayline length (redundant but convenient...maybe)


  // Notes:
  //
  // - Maybe the inverse damping filter will become obsoltet at some point. We have it because we 
  //   need it in rsDampedCombAllpass. Maybe we should have it only there - but even there, we may 
  //   not need it anymore at some point. I think, if both filters (damper and delay) sit in the 
  //   feedback path, we may obtain a pre-delay free impulse response that is allpass in nature. We
  //   introduced it in the first place to compensate for the non-allpassness of the transfer
  //   function when the feedforward path isn't allpass. By restructuring the filter, we may get
  //   away without it.
};


template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombFilter<TSig, TPar, TDly, TTol>::initSettings() 
{ 
  s.init();
  updateDelay();
}

template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombFilter<TSig, TPar, TDly, TTol>::setMaxIntDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay.setMaxDelayInSamples(maxM);
  // Verify this!
}


template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombFilter<TSig, TPar, TDly, TTol>::setup(TDly delay, TPar feedback, int dampOrder,
  const TPar* dampCoeffsB, const TPar* dampCoeffsA, bool predelayMode)
{
  s.setup(delay, Settings::InterpolationMode::nearest, 
    feedback, dampOrder, dampCoeffsB, dampCoeffsA, predelayMode);
  // ToDo: Let the user pick the interpolation mode via another parameter.

  updateDelay();
}

template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombFilter<TSig, TPar, TDly, TTol>::reset()
{
  mainDelay.reset();
  combOut = TSig(0);

  using AT = rsArrayTools;
  AT::clear(xd, maxDmpOrd);
  AT::clear(yd, maxDmpOrd);
  AT::clear(xi, maxDmpOrd);
  AT::clear(yi, maxDmpOrd);
}


template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombFilter<TSig, TPar, TDly, TTol>::updateDelay()
{
  M = s.getIntDelay();
  mainDelay.setDelayInSamples(M);
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

template<class TSig, class TPar, class TDly, class TTol>
class rsDampedCombAllpass // ToDo: derive from rsDampedCombFilter ..or use it as member
{


public:


  //-----------------------------------------------------------------------------------------------
  // \name Lifetime

  /** Standard constructor. Initializes the settings and resets the state to initial conditions. */
  rsDampedCombAllpass()
  {
    initSettings();
    reset();
  }


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets the maximum desired roundtrip delay around the comb. This total roundtrip delay includes
  the z^-1 unit delay, so the delayline length is actually shorter by one. */
  void setMaxIntDelayInSamples(int newMaxDelay);
  // Maybe move to protected area. It's confusing to have both available to client code.

  void setMaxDelayInSamples(TDly newMaxDelay)
  { setMaxIntDelayInSamples((int)rsCeil(newMaxDelay)); }
 


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
  static constexpr int getMaxDampingOrder() { return maxDmpOrd; }

  /** Evaluates the filter's z-domain transfer function H(z) value at the given value of z. */
  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const;
  // Under construction. ToDo: check, if this works with a complex type for TSig.

  rsComplex<TPar> getCombTransferFunctionAt(const rsComplex<TPar>& z) const;

  rsComplex<TPar> getDamperTransferFunctionAt(const rsComplex<TPar>& z) const;

  rsComplex<TPar> getCorrectorTransferFunctionAt(const rsComplex<TPar>& z) const;


  rsComplex<TPar> getDelayTransferFunctionAt(const rsComplex<TPar>& z) const
  {
    return rsPow(z, rsComplex<TPar>(-M));  // H(z) = z^-M
  }








  // Non-allocating versions of the transfer function getters. Well - they *may* allocate - but 
  // will do so only if the passed output parameters and temporary objects have not enough 
  // capacity pre-allocated. Doing so is the responsibility of the caller. They are much less 
  // convenient to use but it's sometimes necessary when one needs to compute these transfer 
  // function in a realtime thread.
  template<class TTol>
  void getCorrectorTransferFunction(rsSparseDigitalTransferFunction<TPar, TTol>* tf) const
  {
    getCombTransferFunction(tf);
    tf->invert();
    tf->reflectZeros();
  }

  //template<class TTol>
  void getCombTransferFunction(rsSparseDigitalTransferFunction<TPar, TTol>* tf) const
  {
    s.getCombTransferFunction(tf);
  }

  //void getDamperTransferFunction(   rsSparseDigitalTransferFunction<TPar>* tf) const
  //{
  //  tf->setupFromDenseCoeffs(b, dmpOrd+1, a, dmpOrd+1, TPar(0));
  //}
  //// This may not be needed

  //template<class TTol>
  void getDelayTransferFunction(rsSparseDigitalTransferFunction<TPar, TTol>* tf) const
  {
    tf->getNumerator().  _setNumTerms(1); tf->getNumerator().  _setTerm(0, TPar(1), M);
    tf->getDenominator()._setNumTerms(1); tf->getDenominator()._setTerm(0, TPar(1), 0);
  }



  // Maybe get rid of them - they are merely more convenient versions of the ones above (but they 
  // allocate)...but: the overall getTransferFunction function is missing. I think, to implement 
  // that in a non-allocating way, we would need a temporary transfer function variable:

  //template<class TTol>
  rsSparseDigitalTransferFunction<TPar, TTol> getTransferFunction(TTol tol) const
  {
    return getCombTransferFunction(tol) * getCorrectorTransferFunction(tol);
  }

  //template<class TTol>
  rsSparseDigitalTransferFunction<TPar, TTol> getCombTransferFunction(TTol tol) const
  {
    using TF = rsSparseDigitalTransferFunction<TPar, TTol>;

    TF one; one.getNumerator()._appendTerm(TPar(1), 0); // Use setToOne()       z^-0
    TF z1;  z1.getNumerator()._appendTerm( TPar(1), 1); // Use setToIdentity()  z^-1
    TF F = getDamperTransferFunction(tol);    // Feedback filter F(z)
    TF A; getDelayTransferFunction(&A);       // Delay filter A(z)
    TPar k = s.getFeedbackGain();
    if(s.isInPreDelayMode())
      return A   / (one + k * z1 * F * A);    // U(z) = A(z) / (1 + k * z^-1 * F(z) * A(z))
    else 
      return one / (one + k * z1 * F * A);    // U(z) =   1  / (1 + k * z^-1 * F(z) * A(z))
  }

  //template<class TTol>
  rsSparseDigitalTransferFunction<TPar, TTol> getCorrectorTransferFunction(TTol tol) const
  {
    using TF = rsSparseDigitalTransferFunction<TPar, TTol>;

    TF C = getCombTransferFunction(tol);
    C.invert();
    C.reflectZeros();

    return C;
  }

  //template<class TTol>
  rsSparseDigitalTransferFunction<TPar, TTol> getDamperTransferFunction(TTol tol) const
  {
    rsSparseDigitalTransferFunction<TPar, TTol> H;
    H.setRoundoffTolerance(tol);
    H.setupFromDenseCoeffs(s.getDampCoeffsB(), s.getDampingOrder()+1,
      s.getDampCoeffsA(), s.getDampingOrder()+1, TPar(0));
    return H;

    // Factor out into s.getDamperTransferFunction();
  }
  // Maybe get rid of them. Yes, they make some things more convenient - but they blow up the API
  // to much for that. Maybe, if needed, implement them as free functions somewhere in the
  // rs_testing module



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
  TSig applyDelay(TSig in) { return mainDelay.getSample(in); }

  /** Applies the feedback damping filter to the signal x and updates the filter's state. */
  TSig applyDamper(TSig in) { return s.applyDamper(in, xd, yd); }

  /** Applies the inverse feedback damping filter to the signal x and updates the filter's state.
  this filter is need only the "without predelay" mode of operation. */
  TSig applyInverseDamper(TSig in) { return s.applyInverseDamper(in, xi, yi); }

  /** Applies the poles of the correction filter which are the same as the poles of the damping 
  filter. */
  TSig applyCorrectorPoles(TSig in) { return s.applyDamperPoles(in, yc); }

  /** Updates the lengths of the delaylines according to the settings. */
  void updateDelays();



  // Embedded DSP objects:
  rsDelay<TSig> mainDelay;                // Main delayline for the comb filter
  rsDelay<TSig> corrDelay;                // Delayline for the correction filter

  // Maximum order of feedback damping filter:
  static const int maxDmpOrd = rsDampedCombSettings<TPar, TDly, TTol>::getMaxDampingOrder(); 
  // Maybe replace double with TDly - a 3rd template parameter for this class

  // State:
  TSig combOut = TSig(0);                 // State for the unit delay feedback loop
  TSig xd[maxDmpOrd], yd[maxDmpOrd];      // State for the damping filter
  TSig xi[maxDmpOrd], yi[maxDmpOrd];      // State for the inverse damping filter
  TSig yc[maxDmpOrd];                     // State for the poles of the correction filter

  // Settings:
  rsDampedCombSettings<TPar, TDly, TTol> s;   // Rename this! ...maybe to settings
  int M = 0;                                  // Delayline length (redundant but convenient)

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
  //
  // - The delayline length member M is redundant (it's equal to mainDelay.getDelayInSamples()) but
  //   we keep it here for convenience (it's used in some formulas and filtering algorithms)
  //
  //
  // ToDo:
  //
  // - Maybe factor out a class rsDampedComb that has everything except the stuff related to the 
  //   correction filter. It should have the mainDelay and the "State" and "Settings" stuff. The 
  //   only extra data member in rsCombAllpass would be the corrDelay. Maybe the stuff under State
  //   could also be moved into a class rsDampedCombState. Maybe that should also contain the
  //   mainDelay delayline (the content of this is also part of the overall state). We do not 
  //   really have much use for such a state class, but it would be cleaner to have it from a 
  //   design perspective
};

template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombAllpass<TSig, TPar, TDly, TTol>::setMaxIntDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay.setMaxDelayInSamples(maxM);
  corrDelay.setMaxDelayInSamples(maxM+maxDmpOrd+1);

  // I think, we may have to add the maximum order of the interpolator filter used in the 
  // delayline. Maybe we should use something like 
  //   corrDelay.setMaxDelayInSamples(maxM + s.getMaxDampPlusIntOrder() + 1);
}

template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombAllpass<TSig, TPar, TDly, TTol>::setup(int delay, TPar feedback, int dampOrder,
  const TPar* dampCoeffsB, const TPar* dampCoeffsA, bool predelayMode)
{
  M = delay - 1;             // -1 corrects for unit delay in feedback path
  s.setup(delay, rsDampedCombSettings<TPar, TDly, TTol>::InterpolationMode::nearest, 
    feedback, dampOrder, dampCoeffsB, dampCoeffsA, predelayMode);
  updateDelays();

  // ToDo:
  //
  // - We need to take the delay as TPar to allow for non-integer delays and then pass that value
  //   to mainDelay.setDelayInSamples(). The member variable M may then be obsolete. I think, the 
  //   corrDelay needs to use  (int) (delay-1) + dmpOrd + 1 + interpolationOrder  but I'm not 
  //   totally sure about that, so that needs to be verified and unit tested. Oh! But supporting
  //   fractional delays is actually more complicated and would also require re-implementation of
  //   applyCorrector(). The current implementation is really ony applicable to the case of an 
  //   integer delay. Hmmmm....maybe we should revert this class to support only integer delays and 
  //   factor out the stuff that designs comb allpasses with non-integer delays. Maybe factor out
  //   a class rsDampedCombAllpassSettings containing  k,b,a,M,dmpOrd,predelay  as data members
  //   and all the  getTransferFunction...  functions as member functions. And maybe some 
  //   design/setup functions. This class here should then maintain a settings member of this type.
}

template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombAllpass<TSig, TPar, TDly, TTol>::initSettings()
{
  mainDelay.setDelayInSamples(0);
  corrDelay.setDelayInSamples(0);
  s.init();
  M = 0;
}

template<class TSig, class TPar, class TDly, class TTol>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar, TDly, TTol>::getTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  return getCombTransferFunctionAt(z) * getCorrectorTransferFunctionAt(z);
}

template<class TSig, class TPar, class TDly, class TTol>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar, TDly, TTol>::getCombTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  using Complex = rsComplex<TPar>;
  Complex one(TPar(1));                         // 1 + 0i
  Complex A  = getDelayTransferFunctionAt(z);   // A(z)
  Complex z1 = one/z;                           // z^-1
  Complex F  = getDamperTransferFunctionAt(z);  // F(z)
  TPar k = s.getFeedbackGain();
  if(s.isInPreDelayMode())
    return A   / (one + k * z1 * F * A);        // U(z) = A(z) / (1 + k * z^-1 * F(z) * A(z))
  else
    return one / (one + k * z1 * F * A);        // U(z) =   1  / (1 + k * z^-1 * F(z) * A(z))

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
  //
  //
  // ToDo:
  //
  // - Maybe move this function into rsDampedCombSettings. But in the "without-predelay" mode, the
  //   transfer function represented by this class should probably feature an F in the numerator...
  //   unless we assume already there, that the F filter will be compensated for - which we do 
  //   here. It's a bit messy. Maybe the 2nd branch needs two sub-branches, switched by a boolean
  //   parameter "isCompensated" or something...which is a bit ugly. Or maybe instead of a boolean 
  //   to indicate predelay and another for indicate the feedabck compensation in case of 
  //   no-predelay, have a an enum parameter for the mode. Options: withPreDelay, noPreDelay, 
  //   noPreDelayCompensated ...or damperInFeedback, delayInFeedback, 
  //   delayInFeedbackDampCompensated...or: feedbackDamped, forwardDamped, forwardDampedCompensated
}

template<class TSig, class TPar, class TDly, class TTol>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar, TDly, TTol>::getDamperTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  return s.getDamperTransferFunctionAt(z);
}

template<class TSig, class TPar, class TDly, class TTol>
rsComplex<TPar> rsDampedCombAllpass<TSig, TPar, TDly, TTol>::getCorrectorTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  int dmpOrd = s.getDampingOrder();
  const TPar* b = s.getDampCoeffsB();
  const TPar* a = s.getDampCoeffsA();
  TPar k = s.getFeedbackGain();

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

////template<class TTol>
////template<class TSig, class TPar, class TDly>
//template<class TSig, class TPar, class TDly, class TTol>
//rsSparseDigitalTransferFunction<TPar, TTol> rsDampedCombAllpass<TSig, TPar, TDly>
//                                      ::getTransferFunction(TTol tol) const
//{
//  return getCombTransferFunction(tol) * getCorrectorTransferFunction(tol);
//}

//template<class TTol>
//template<class TSig, class TPar, class TDly>
//rsSparseDigitalTransferFunction<TPar, TTol> rsDampedCombAllpass<TSig, TPar, TDly>
//                                      ::getCombTransferFunction(TTol tol) const
//{
//  using TF = rsSparseDigitalTransferFunction<TPar, TTol>;
//
//  TF one; one.num._appendTerm(TPar(1), 0);
//  TF z1;  z1.num._appendTerm( TPar(1), 1);
//  TF F = getDamperTransferFunction();       // Feedback filter F(z)
//  TF A; getDelayTransferFunction(&A);       // Delay filter A(z)
//  TPar k = s.getFeedbackGain();
//  if(s.isInPreDelayMode())
//    return A   / (one + k * z1 * F * A);    // U(z) = A(z) / (1 + k * z^-1 * F(z) * A(z))
//  else 
//    return one / (one + k * z1 * F * A);    // U(z) =   1  / (1 + k * z^-1 * F(z) * A(z))
//}

//template<class TTol>
//template<class TSig, class TPar, class TDly>
//rsSparseDigitalTransferFunction<TPar, TTol> rsDampedCombAllpass<TSig, TPar, TDly>
//                                      ::getCorrectorTransferFunction(TTol tol) const
//{
//  using TF = rsSparseDigitalTransferFunction<TPar, TTol>;
//
//  TF C = getCombTransferFunction(tol);
//  C.invert();
//  C.reflectZeros();
//
//  return C;
//}

//template<class TTol>
//template<class TSig, class TPar, class TDly>
//rsSparseDigitalTransferFunction<TPar, TTol> rsDampedCombAllpass<TSig, TPar, TDly>
//                                      ::getDamperTransferFunction(TTol tol) const
//{
//  rsSparseDigitalTransferFunction<TPar, TTol> H;
//  H.setRoundoffTolerance(tol);
//  H.setupFromDenseCoeffs(s.getDampCoeffsB(), s.getDampingOrder()+1,
//                         s.getDampCoeffsA(), s.getDampingOrder()+1, TPar(0));
//  return H;
//
//  // Factor out into s.getDamperTransferFunction();
//}

template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombAllpass<TSig, TPar, TDly, TTol>::reset()
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

template<class TSig, class TPar, class TDly, class TTol>
TSig rsDampedCombAllpass<TSig, TPar, TDly, TTol>::getSampleComb(TSig in)
{
  TPar k = s.getFeedbackGain();
  if(s.isInPreDelayMode())
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

template<class TSig, class TPar, class TDly, class TTol>
TSig rsDampedCombAllpass<TSig, TPar, TDly, TTol>::applyCorrector(TSig in)
{
  // Retrieve feedback gain and damping coeffs:
  TPar k        = s.getFeedbackGain();
  int  dmpOrd   = s.getDampingOrder();
  const TPar* b = s.getDampCoeffsB();
  const TPar* a = s.getDampCoeffsA();

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

template<class TSig, class TPar, class TDly, class TTol>
void rsDampedCombAllpass<TSig, TPar, TDly, TTol>::updateDelays()
{
  mainDelay.setDelayInSamples(M);
  corrDelay.setDelayInSamples(M + s.getDampingOrder() + 1); 

  // I think, this may need more delay memory when we have an interpolating delayline. I think, we
  // may have to add the order of the interpolator. Maybe we should use a function like
  // s.getDampPlusIntOrder() which returns the sum of the orders of damper and interpolator. Or we
  // just stick to a non-interpolating delaylije here.
}





// A free function to set up the object with a more convenient parametrization:
template<class TSig, class TPar, class TDly, class TTol>
void rsSetupHighDamp(rsDampedCombAllpass<TSig, TPar, TDly, TTol>& flt,
  int delay, TPar feedback, TPar dampOmega, TPar dampGain, bool predelay)
{
  TPar a[2], b[2]; a[0] = 1;
  rsMake1stOrderHighShelf(dampOmega, dampGain, &b[0], &b[1], &a[1]);
  flt.setup(delay, feedback, 1, b, a, predelay);
}

// A function to set up the object with a fractionla delay by incorporating a linear interpolation 
// filter into the feedback filter. LinViaFb stands for "linear interpolation via the feedback 
// filter" ...ToDo: Explain why this works. Does it actually work, though? ..I think it works for
// FIR interpolation filters but not IIR (like allpass interpolators)
//template<class TSig, class TPar>
template<class TSig, class TPar, class TDly, class TTol>
void rsSetupFractional_LinViaFb(rsDampedCombAllpass<TSig, TPar, TDly, TTol>& flt,
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
//template<class TSig, class TPar>
template<class TSig, class TPar, class TDly, class TTol>
void rsSetupDecayTimes_LinViaFb(rsDampedCombAllpass<TSig, TPar, TDly, TTol>& flt, 
  TPar delay, TPar decay, TPar loOmega, TPar loScale, TPar hiOmega, TPar hiScale, 
  bool predelay)
{
  // Compute feedback gain and filter coeffs:
  TPar a[4], b[4];
  TPar kM = rsMakeDampBiShelf(delay, decay, loOmega, loScale, hiOmega, hiScale, b, a);

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


  // ToDo: Use  rsSetupDecayTimes_LinViaFb(s, ...)  and then maybe call a flt.updateDelays()
  // function if necessary. This function here copies a lot of code from the other rsSetup...
  // function. We want to get rid of that duplication.
}


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


template<class TSig, class TPar, class TTol> // Get rid of TTol - Give member functions that need it their own template param
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
    //protoAllpass.setMaxDelayInSamples(maxDelay);
    settings.resize(maxNumCombs);


    int maxBankDelay = maxDelay * maxNumCombs + 4;
    // VERIFY this formula! Theoretically and practically (by making a unit test using the maxmimum
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

  void setMaxPhaseCombMode(bool shouldBeMaxPhase)
  {
    maxPhaseCombBank = shouldBeMaxPhase;
  }

  void setSerialCombsMode(bool shouldBeSerial)
  {
    serialCombs = shouldBeSerial;
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
  


  rsComplex<TPar> getCombTransferFunctionAt(const rsComplex<TPar>& z) const
  {
    return combBank.getTransferFunctionAt(z);
  }

  rsComplex<TPar> getCorrectorTransferFunctionAt(const rsComplex<TPar>& z) const
  {
    return corrector.getTransferFunctionAt(z);
  }

  rsComplex<TPar> getTransferFunctionAt(const rsComplex<TPar>& z) const
  {
    return getCombTransferFunctionAt(z) * getCorrectorTransferFunctionAt(z);
  }





  // Maybe let them take a tol param
  rsSparseDigitalTransferFunction<TPar, TTol> getTransferFunction() const
  {
    return getCombTransferFunction() * getCorrectorTransferFunction();
  }

  rsSparseDigitalTransferFunction<TPar, TTol> getCombTransferFunction() const
  {
    return combBank.getTransferFunction();
    //return U;  // Should also work, I think.
  }
  // allocates - creates copy of the transfer function object.
  // Maybe return a const ref?

  rsSparseDigitalTransferFunction<TPar, TTol> getCorrectorTransferFunction() const
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
    // Maybe use dirty.store(shouldBeDirty). I think, it makes no difference with regard to the
    // generated assembly code but it may add documentation value. We would document that this is 
    // intended to be an atomic operation.
  }

  void updateFilters(); // Allocates! Not yet realtime ready. ...might be fixed....verify!
                        // Maybe rename to updateCoeffs or calcCoeffs


  // Embedded DSP objects:
  rsSparseFilter<TSig, TPar, TTol> combBank;
  rsSparseFilter<TSig, TPar, TTol> corrector;



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
  std::vector<CombSettings> settings;   // Rename to perCombSettings

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
  // Maybe the name "bank" is not appropriate anymore. It can now also be a "chain". Banks are 
  // parallel connections, chains serial ones.

  bool serialCombs = false;

  // Flag to indicate that a call to updateFilters() is needed before doing any DSP:
  std::atomic<bool> dirty = true;



  // This object will be used compute the filter coeffs of the prototype combs:
  //rsDampedCombSettings<TPar> protoComb;
  rsDampedCombSettings<TPar, double, TTol> protoComb;
  // TODO: Replace double with TDly

  // Transfer function objects used for temporaries in internal computations in our updateFilters()
  // function. They are members rather than locals there to avoid heap allocations in this 
  // function.
  rsSparseDigitalTransferFunction<TPar, TTol> U, Ui;
};


template<class TSig, class TPar, class TTol>
void rsDampedMultiCombAllpass<TSig, TPar, TTol>::updateFilters()
{
  // This is still under construction. It still has allocations and it needs to treat the case 
  // numCombs == 0. ...the allocations might be gone now - but verify this! To treat numCombs == 0,
  // I think, we just need to assign U(z) = 1/1. Maybe we could generally add a direct path with a 
  // given gain to the comb-bank. Then instead of doing U.clear() before the accumulation loop, we
  // could just init it as U.initToConst(directGain). If directGain happens to be 0, we would get
  // back to the currently implemented case. Then we wouldn't need a special treatment and would
  // have an even more flexible filter.

  using RatFunc   = rsSparseRationalFunction<TPar, TTol>;
  //using TransFunc = rsSparseDigitalTransferFunction<TPar>;

  TPar decaySamples =      decayTime     * sampleRate;
  TPar lowOmega     = 2*PI*lowCrossFreq  / sampleRate;
  TPar highOmega    = 2*PI*highCrossFreq / sampleRate;

  // ToDo: catch special case of numCombs == 0. In this case, we should set up a trivial bypass 
  // filter


  // Init transfer function U:
  if(serialCombs == false)
    U.initToZero();            // Init to U(z) = 0 for additive accumulation
  else
    U.initToOne();             // Init to U(z) = 1 for multiplicative accumulation
  
  // Accumulate the transfer function of the comb bank or chain:
  for(int i = 0; i < numCombs; i++)
  {
    const CombSettings& s = settings[i];


    TPar delay = sampleRate / (s.freqScale * frequency);  // Verify!
    // Maybe in case of only odd harmonics, we should scale the freq up by 2? The rationale is that 
    // when using odd harmonics only (by way of the feedback sign), the fundamental frequency 
    // actually goes an octave lower.


    rsSetupDecayTimes_LinViaFb(protoComb, delay, decaySamples, 
      lowOmega, lowDecayScale, highOmega, highDecayScale, false);
    // This needs an additional parameter to determine the sign of the feedback, i.e. switch 
    // between all and only odd harmonics


    // Accumulate the i-th comb's transfer function Ui into our total transfer function U:
    protoComb.getCombTransferFunction(&Ui);
    if(serialCombs == false)
      RatFunc::weightedSumDestructive(&U, TPar(1), &Ui, TPar(s.gain), &U, TPar(0));
    else
      U.multiplyBy(Ui, TPar(0));  // Verify if this works in place!

  }

  // Set up comb bank and corrector:
  combBank.setup(U);
  if(maxPhaseCombBank)
    combBank.reflectZeros();
  corrector.setup(U);
  corrector.invert();
  corrector.reflectZeros();

  // After updating the filters, we are in clean state:
  setDirty(false);


  // ToDo:
  //
  // - Have a member serial (defaulting to false). If true, do a serial connection accumulation 
  //   loop. Init U(z) to 1 instead of 0 and accumulate multiplicatively instead of additively.
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

template<class TSig, class TPar, class TTol>
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
  void convertCombSumToDirectForm(rsSparseFilter<TSig, TPar, TTol>* sparseDirectFormFilter);
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
  rsSparseFilter<TSig, TPar, TTol> corrector;

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


template<class TSig, class TPar, class TTol>
void rsDampedAllpassBiComb_1p<TSig, TPar, TTol>::rsDampedAllpassBiComb_1p<TSig, TPar, TTol>::reset()
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

//template<class TSig, class TPar>
template<class TSig, class TPar, class TTol>
void rsDampedAllpassBiComb_1p<TSig, TPar, TTol>::setMaxDelayInSamples(int newMaxDelay)
{
  int maxM = newMaxDelay - 1;
  mainDelay1.setMaxDelayInSamples(maxM);
  mainDelay2.setMaxDelayInSamples(maxM);
  corrector.setMaxDelayInSamples(2*maxM+4);
  // See convertCombSumToDirectForm(). The maximum delay that occurs there is: M1+M2+4.
}

//template<class TSig, class TPar>
template<class TSig, class TPar, class TTol>
void rsDampedAllpassBiComb_1p<TSig, TPar, TTol>::setup(
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


//template<class TSig, class TPar>
template<class TSig, class TPar, class TTol>
rsComplex<TPar> rsDampedAllpassBiComb_1p<TSig, TPar, TTol>::getCombTransferFunctionAt(
  const rsComplex<TPar>& z) const
{
  rsError("Not yet implemented");
  return rsComplex<TPar>(0); 

  // ToDo:
  //
  // - Compute the transfer function based on computing the two transfer functions of the two 
  //   combs and forming a weighted sum of them.
}

//template<class TSig, class TPar>
template<class TSig, class TPar, class TTol>
void rsDampedAllpassBiComb_1p<TSig, TPar, TTol>::convertCombSumToDirectForm(
  rsSparseFilter<TSig, TPar, TTol>* sparseFilter)
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


//template<class TSig, class TPar>
template<class TSig, class TPar, class TTol>
void rsDampedAllpassBiComb_1p<TSig, TPar, TTol>::updateDelaysAndCorrectorCoeffs()
{
  mainDelay1.setDelayInSamples(M1);
  mainDelay2.setDelayInSamples(M2);

  convertCombSumToDirectForm(&corrector);
  corrector.invert();  // This could potentially allocate for the vector swap. Figure this out!
  corrector.reflectZeros();
}


//=================================================================================================

/** A nonlinear extension of rsDampedCombAllpass. At the moment, it's just an experimental stub. */

template<class TSig, class TPar, class TDly, class TTol>
class rsDampedCombAllpassNonLin : public rsDampedCombAllpass<TSig, TPar, TDly, TTol>
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





//=================================================================================================



/**  Under construction...  */

template<class TSig, class TPar>
class rsProtoFDN
{



public:

  //-----------------------------------------------------------------------------------------------
  // \name Setup

  //** Sets up the lengths of the delaylines. */
  void setDelays(const std::vector<int>& newDelays)
  {
    int N = (int) newDelays.size();
    delays.resize(N);
    for(int i = 0; i < N; i++)
    {
      int d = newDelays[i];
      rsAssert(d >= 1, "Delay must be at leat 1 to avoid delayless feedback loop.");
      d = rsMax(d, 1) - 1;                // -1 to compensate for implicit feedback loop delay
      delays[i].setMaxDelayInSamples(d);  // Reallocates in case of too small capacity
      delays[i].setDelayInSamples(   d);

      // ToDo: Try to solve that problem more elegantly. I'm not sure how, though. Maybe we have to
      // do it like that. ...yeah - I think so. Otherwise, we'd get delay-free feedback loops. In 
      // the extended FDN, we need to do that only for the pre-matrix delaylines. We actually want
      // to be able to set the post-matrix delays to zero because that's how we deactivate them.
      //
      // But: what about the dampFactors? Are they now wrong? I don't think so, though.
    }
  }
  // May allocate. Maybe make it safely non-allocating. Perhaps provide a setMaxDelayInSamples()
  // method. Maybe, for optimization of memory usage, it may also take a vector parameter such that
  // the different delaylines may use different maximal delays. This might be more more economic 
  // when we want to use a wide range of (maximum) delays such that it would be wasteful to give 
  // all delaylines the greatest maximum length. Although - none of that really matters unless 
  // this class is used in a production environment which shouldn't be done anyway. It's purely for
  // research purposes.


  void setFeedbackMatrix(const rsMatrix<TPar>& newFeedbackMatrix)
  {
    rsAssert(newFeedbackMatrix.isSquare());
    feedbackMatrix = newFeedbackMatrix;
    state.resize(feedbackMatrix.getNumRows());
    outs.resize( feedbackMatrix.getNumRows());
  }

  void setInputMatrix(const rsMatrix<TPar>& newInputMatrix)    { inMatrix  = newInputMatrix; }

  void setOutputMatrix(const rsMatrix<TPar>& newOutputMatrix)  { outMatrix = newOutputMatrix; }

  void setDampFactors(const std::vector<TPar>& newDampFactors) { dampFactors = newDampFactors; }

  // ToDo: Maybe provide a setup() function similar to the one in rsStateSpaceFilter. It needs
  // an additional parameter (of type std::vector<int>) for the delays


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getNumDelayChannels() const { return (int) delays.size(); }

  int getNumInputs()        const { return inMatrix.getNumColumns(); }

  int getNumOutputs()       const { return outMatrix.getNumRows(); }

  /** Returns the amount of delay of the i-th delayline. This amount includes the implicit unit 
  delay of the feedback loop. */
  int getDelay(int i) const 
  { 
    rsAssert(i >= 0 && i < (int)delays.size(), "Index out of range in rsProtoFDN::getDelay");
    return delays[i].getDelayInSamples() + 1;
    // +1 for the implicit unit delay in the feedback loop.
  }

  /** Checks, if the lengths and shapes of the various vectors and matrices fit together. */
  bool areSettingsConsistent() const
  {
    bool ok = true;
    int N = getNumDelayChannels();

    ok &= feedbackMatrix.hasShape(N, N);
    ok &= inMatrix.getNumRows()     == N;
    ok &= outMatrix.getNumColumns() == N;
    ok &= (int) delays.size()       == N;
    ok &= (int) dampFactors.size()  == N;
    ok &= (int) state.size()        == N;
    ok &= (int) outs.size()         == N;

    return ok;
  }







  template<class TArg>
  TArg getDelayTransferFunctionAt(const TArg& z, int n)
  {
    TArg H = rsUnityValue(z) / z;                // Init with z^-1 for the implicit unit delay.
    H *= delays[n].getTransferFunctionAt(z);     // Bake in the delay transfer function.
    H *= dampFactors[n];                         // Bake in the damping/decay factor.
    //H *= dampers[n].getTransferFunctionAt(z);  // ToDo: Bake in the damping filter here.
    return H;
  }

  // Maybe rename to getTransferMatrixAt (also in class rsStateSpaceFilter for consistency)
  template<class TArg>
  rsMatrix<TArg> getTransferFunctionAt(const TArg& z) 
  {
    using Mat    = rsMatrix<TArg>;
    using LinAlg = rsLinearAlgebraNew;

    // Compute delay matrix D(z):
    int N = getNumDelayChannels();
    Mat D(N, N, rsZeroValue(z));
    for(int n = 0; n < N; n++)
      D(n, n) = TPar(1) / getDelayTransferFunctionAt(z, n);  // Why reciprocal?

    // Convert feedback-, input- and output-matrices to TArg (which is typically complex):
    Mat A; rsConvert(feedbackMatrix, &A);
    Mat B; rsConvert(inMatrix,       &B);
    Mat C; rsConvert(outMatrix,      &C);

    // Compute and return the transfer function matrix:
    Mat M = LinAlg::inverse(D-A);  // (D-A)^-1
    Mat H = C * M * B;
    return H;


    // According to the DAFX book (1st Ed), page 182, the transfer function of an FDN is given by:
    //
    //   H(z) = c^T * (D - A)^-1 * b + d
    //
    // where
    //
    //   D: delay matrix = D(z^-1) = diag(z^-m1, z^-m2, ..., z^-mN)
    //   A: NxN feedback matrix
    //   b: input vector
    //   c: output vector
    //   d: direct path gain (a scalar in the book - for mono input) - it's zero here
    //
    // The poles and zeros are the solutions of:
    //
    //   det(A - D) = 0
    //   det(A - b*c^T / d - D) = 0 
    //
    // But using the formulas from the DAFX book as is didn't work out. I had to massage them quite
    // a bit by trial and error and compare with the implementation of getTransferFunctionAt() of
    // class rsStateSpaceFilter to make it finally work. Without that working implementation in the
    // state space filter, I wouldn't have had a chance to figure it out. It seems to work now but 
    // we should really do some thorough unit tests. That implies that the formulas for the poles
    // and zeros should also be taken with a grain of salt. Also, in our case here, d is zero. What
    // does that mean for the formula for the zeros which has a division by d? Maybe it means that 
    // we don't have any zeros? Can we tune b, c and d such that the zeros are mirror images of the
    // poles reflected at the unit circle such that we get an overall allpass filter? That might be
    // a nice feature.
    //
    // Try to get rid of computing the inverse by replacing it by a call to a solver of a linear
    // system. Not sure, if that works though. The inverse´matrix is sandwiched between two vectors
    // so we can't just rewrite it as a linear system. In general, whenever we have a matrix 
    // equation like X = A^-1 * B for an unknown matrix X, it's better to write it as A * X = B
    // and give it to a linear system solver. But here, the A^-1 is sandwiched (our matrix M has 
    // the role of A^-1 of the general form) between two other matrices, so I don't know, if we can
    // do something similar here. Maybe we could premultiply both sides with C^-1 to get the form
    // C^-1 * H = C^-1 * C * M * B = M * B, compute the LHS using a linear solver and then do a
    // matrix multiply by C to get H? But C need not to be invertible. It's usually not even a 
    // square matrix. Maybe the pseudoinverse could be used? ...but maybe not.
    //
    //
    // ToDo:
    //
    // - Figure out, why we need to take the reciprocal in the computation of D(n, n). I've figured 
    //   that out by trial and error. Maybe it has to do with the fact that in the DAFX book, 
    //   they call the matrix D(z^-1) rather than D(z)?
    //
    // - Maybe rename the D matrix to Z to be compatible with the notation in class 
    //   rsStateSpaceFilter. This class also has a D matrix but there, it has a different menaing: 
    //   It is the "feedaround" matrix (i.e. the direct path from inputs to outputs) there. Maybe
    //   introduce such a feedaround matrix here, too (defaulting to all zeros) and call *that* D. 
    //   Document that deviation in notation from DAFX (and possibly the wider DSP/FDN literature).
  }
  // Allocates! Not for realtime use!
  

  template<class TTol>
  void getDelayTransferFunction(int n, rsSparseDigitalTransferFunction<TPar, TTol>* tf) const
  {
    int M = getDelay(n);
    tf->getNumerator()._setNumTerms(1); 
    tf->getNumerator()._setTerm(0, TPar(dampFactors[n]), M+1);
    tf->getDenominator()._setNumTerms(1); 
    tf->getDenominator()._setTerm(0, TPar(1), 0);
    // Maybe simplify this by introducing a function tf->setToPower(dampFactors[n], M+1)
  }
  // Needs tests


  // Under construction:
  template<class TTol>
  rsMatrix<rsSparseDigitalTransferFunction<TPar, TTol>> getTransferFunction(TTol tol)
  {
    using TF  = rsSparseDigitalTransferFunction<TPar, TTol>;
    using LA  = rsLinearAlgebraNew;
    using Mat = rsMatrix<TF>;

    int numIns   = getNumInputs();
    int numOuts  = getNumOutputs();
    int numChans = getNumDelayChannels();

    // Create the required matrices with transfer functions as elements:
    Mat D(numChans, numChans);
    for(int n = 0; n < numChans; n++)
    {
      // What about the tolerance? Should we set it here? Like D(n,n).setRoundoffTolerance(tol)?

      getDelayTransferFunction(n, &(D(n,n))); 
      D(n,n).invert();
      // After the loop, the denominators of the diagonal elements of D are almost 1 but not 
      // exactly. Maybe an epsilon below or something? But why? Shouldn't the values be exact? 
      // Figure out and document!
    }
    Mat A; rsConvert(feedbackMatrix, &A);
    Mat B; rsConvert(inMatrix,       &B);
    Mat C; rsConvert(outMatrix,      &C);


    // Compute and return the transfer function matrix:
    // Mat dmA = D-A;  // This is fine

    Mat M = LA::inverse(D-A);  // (D-A)^-1
    // Triggers rsAssert. It occurs when searching for the pivot in the inversion algorithm. 
    // Apparently, one of the matrix elements isn't in canonical representation when we expect it
    // to be. Figure out why and fix this! Maybe implement a unit test that tests linear algebra
    // with rational functions (sparse and non-sparse)

    Mat H = C * M * B;
    return H;


    // Preliminary (and wrong!):
    //rsMatrix<TF> H(numOuts, numIns);
    //return H;





    // I think, the matrix inversion step will be complicated. Maybe as a preliminary, we need to
    // try creating matrices of sparse transfer functions and doing linear algebra with them. 
    // That's actually a pretty complicated thing to do. I'm not sure, how pivoting should work in
    // this case, for example. Maybe in the implementation of the Gauss-Jordan algorithm, when we
    // search for the pivot, we should use a template function rsIsBetterPivot(T lhs, T rhs) that
    // defaults to calling: rsGreaterAbs(lhs, rhs) but can be "overriden" for other datatypes by
    // explicit specialization. For matrices of type rsFraction, the greater-abs criterion may also
    // be inappropriate because in this case, floating point rounding errors are no issue. Instead,
    // we may have to worry about integer overflow. There, we may want to choose the "simplest" 
    // fraction in order to make overflow less likely. ...and with rational functions? Maybe the
    // simplest (nonzero) rational function in terms of degrees of numerator and denominator would 
    // make for the best pivot? Do we even need pivoting or can we get away without it? Maybe the 
    // matrix D-A has special properties that make pivoting superfluous (maybe it's diagonally 
    // dominant or something?)
    // https://math.stackexchange.com/questions/2485574/strictly-column-diagonally-dominant-matrices-and-gaussian-elimination-with-parti?rq=1


    // Notes:
    //
    // - We do not need to explicitly initialize the D matrix by all zeros (even though the 
    //   constructor of rsMatrix does no zero-initialization) because the default 
    //   constructor of class rsSparseDigitalTransferFunction class initializes with zero.
  }





  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void processFrame(const std::vector<TSig>& inputs, std::vector<TSig>& outputs)
  {
    int numIns   = (int) inputs.size();
    int numOuts  = (int) outputs.size();
    int numChans = getNumDelayChannels();

    rsAssert(areSettingsConsistent(), "Inconsistent settings in rsProtoFDN::processFrame()");
    rsAssert(numIns  == getNumInputs());
    rsAssert(numOuts == getNumOutputs()); 

    using Vec = std::vector<TSig>;

    // Form the outputs:
    rsSetZero(outputs);
    for(int i = 0; i < numOuts; i++)
      for(int j = 0; j < numChans; j++)
        outputs[i] += outMatrix(i, j)  * outs[j];

    // Form the FDN input by applying the input matrix to the inputs vector:
    Vec x(numChans); 
    rsSetZero(x);     // rsSetZero is superfluous now but maybe later we use a member for x
    for(int i = 0; i < numChans; i++)
      for(int j = 0; j < numIns; j++)
        x[i] += inMatrix(i, j) * inputs[j];

    // Form the inputs to the delaylines:
    Vec u(numChans);
    for(int i = 0; i < numChans; i++)
      u[i] = x[i] + state[i];

    // Apply the delaylines:
    for(int i = 0; i < numChans; i++)
    {
      outs[i] = dampFactors[i] * delays[i].getSample(u[i]);

      // I think, when damping filters are included, they should be applied here, like so:
      //outs[i] = dampFactors[i] * dampers[i].getSample(delays[i].getSample(u[i]));
    }

    // Apply the feedback matrix:
    rsSetZero(state);
    for(int i = 0; i < numChans; i++)
      for(int j = 0; j < numChans; j++)
        state[i] += feedbackMatrix(i, j) * outs[j];

    // Can we reorder the operations to get rid of the "outs" member?
  }


  void reset()
  {
    rsSetZero(state);
    rsSetZero(outs);
    for(size_t i = 0; i < delays.size(); i++)
      delays[i].reset();
  }



protected:

  using Delay = rsDelay<TSig>;        // Type alias for convenience

  // State:
  std::vector<Delay> delays;          // Delaylines
  std::vector<TSig>  state;           // State of the FDN
  std::vector<TSig>  outs;            // Output signals

  // Settings:
  std::vector<TPar>  dampFactors;     // Damping/decay factors
  rsMatrix<TPar>     feedbackMatrix;  // Feedback matrix
  rsMatrix<TPar>     inMatrix;        // Input matrix
  rsMatrix<TPar>     outMatrix;       // Output matrix
  //rsMatrix<TPar>     thruMatrix;       // Direct throughput/feedaround matrix

  // ToDo:
  //
  // - Use an interpolating delayline class. Using allpass interpolation makes the most sense in 
  //   this context, I think. It doesn't destroy the unitarity of the prototype network (assuming 
  //   a unitary feedback matrix and dampFactors of all 1s), if I'm not mistaken.
  //
  // - Maybe we could even use the class rsUniversalComb instead of the simple delaylines. But 
  //   maybe then wen should have a variant of this class that also incoprorates (allpass) 
  //   interpolation. It may make sense to use the universal comb in notchpass mode. 
  //
  // - Try to make an 1-in/1-out FDN with an overall allpass transfer function. Maybe take the 
  //   feedback matrix and input matrix as given and try to tune the output matrix. Or maybe take
  //   only the feedback matrix as given and tune in and output matrices. Or maybe take feedback- 
  //   input and output matrices as given and try to add an appropriate feedaround matrix that 
  //   turns the whole FDN into an allpass. Or take the whole FDN as given and try to create a 
  //   compensation filter that turns the series of FDN -> compensator into an allpass. Maybe as a
  //   preliminary step, investigate how to create state space filters with allpass characteristic
  //   and then generalize the findings to FDNs.
};








//=================================================================================================

/** Under construction...

Prototype implementation of a feedback delay network that is meant purely for experimentation in
research and development. For this purpose, the focus is deliberately not on efficiency but rather 
on generality, flexibility with respect to the general configuration (with respect to where we 
inject inputs, pick up outputs, etc), the number of delaylines (no restriction to powers or two, 
for example), the feedback matrix (no restriction matrices with efficient implementation of the 
matrix vector product), etc. It should enable convenient experimentation with various 
architectures, settings, etc. It's a vehicle to explore the vast search space of possible reverb 
designs. Even if we restrict ourselves to an FDN architecture, the design space is still vast. 
Actually, many of the classic reverb designs can be recast in terms of an FDN structure as well, so
FDNs indeed provide a general framework to implement and investigate these other structures as well
(albeit in a suboptimal way with regard to efficiency). The goal is to identify promising designs 
which may later be implemented in a more efficient way for use in production. */

template<class TSig, class TPar>
class rsExtendedProtoFDN
{

public:


  //-----------------------------------------------------------------------------------------------
  // \name Setup

  /** Sets up the lengths of the delaylines before and after the feedback matrix and the feedback
  matrix itself. The feedback matrix should be unitary. The damping is taken care of elsewhere. */
  //void setDelaysAndFeedback(
  //  const std::vector<int>& newPreMatrixDelays,
  //  const rsMatrix<TPar>& newFeedbackMatrix,
  //  const std::vector<int>& newPostMatrixDelays);
  // We allow to set these 3 things only all at once because there are consistency constraints that 
  // need to be observed (pre- and post matrix delays must have same size, matrix must be of shape 
  // size x size). That would be messy to ensure with separate setters for each of the 3. Maybe
  // make the last parameter optional. ..but then we can't use references, I think. Or can we? If
  // not, just provide a 2nd function that has only 2 parameters. Inside of it, create a dummy 
  // vector for the newPostMatrixDelays of all zeros and then delegate to the 3-parameter version 
  // of the function.


  void setDelays(
    const std::vector<int>& newPreMatrixDelays,
    const std::vector<int>& newPostMatrixDelays);

  void setFeedbackMatrix(const rsMatrix<TPar>& newFeedbackMatrix)
  {
    rsAssert(newFeedbackMatrix.isSquare());
    feedbackMatrix = newFeedbackMatrix;
    state.resize(feedbackMatrix.getNumRows());
  }

  void setInputMatrices(
    const rsMatrix<TPar>& newInputMatrixPre,
    const rsMatrix<TPar>& newInputMatrixPost)
  {
    inMatrixPre  = newInputMatrixPre;
    inMatrixPost = newInputMatrixPost;
  }

  void setOutputMatrices(
    const rsMatrix<TPar>& newOutputMatrixPre,
    const rsMatrix<TPar>& newOutputMatrixPost)
  {
    outMatrixPre  = newOutputMatrixPre;
    outMatrixPost = newOutputMatrixPost;
  }

  void setDampFactors(
    const std::vector<TPar>& newDampFactorsPre,
    const std::vector<TPar>& newDampFactorsPost)
  {
    dampFactorsPre  = newDampFactorsPre;
    dampFactorsPost = newDampFactorsPost;
  }


  //-----------------------------------------------------------------------------------------------
  // \name Inquiry

  int getNumDelayChannels() const
  {
    return (int) delaysPre.size();

    // ToDo: Maybe assert that delaysPost() and state have the same size and the feedback matrix is 
    // size x size
  }

  /** Checks, if the lengths and shapes of the various vectors and matrices fit together. */
  bool areSettingsConsistent() const;


  //-----------------------------------------------------------------------------------------------
  // \name Processing

  void processFrame(const std::vector<TSig>& inputs, std::vector<TSig>& outputs);
  // We use std::vector rather than raw arrays because we don't really need the flexibility to deal
  // with anything other than std::vector in our experiments and std::vector is more convenient and 
  // safer than raw arrays.



protected:

  //void processDelayChannels(std::vector<TSig>& ioData);
  // Assumes that ioData is of length getNumDelayChannels()
  // Maybe use a member array for the ioData..but maybe it's nicer to pass an array in. Maybe we 
  // want to make this function public later. Then we may want to look into the ins and outs. Maybe
  // we should use separate arrays for ins and outs. It may be convenient and we don't care about 
  // optimizing memory usage here.


  // State of the FDN:
  std::vector<TSig> state;                // State of the FDN
  // Maybe use several vectors for the state at different position in the signal processing 
  // pipeline, e.g. before the delays, after the pre-matrix delays, after the pre-matrix dampers, 
  // after the matrix, after the post matrix dampers, etc. and then provide accessors to allow 
  // client code to peek into these intermediate signals after calling processFrame(). That might 
  // be useful for R&D.

  // Processing elements:
  std::vector<rsDelay<TSig>> delaysPre;   // Delaylines pre feedback matrix
  std::vector<rsDelay<TSig>> delaysPost;  // Delaylines post feedback matrix

  rsMatrix<TPar> feedbackMatrix;

  rsMatrix<TPar> inMatrixPre;
  rsMatrix<TPar> inMatrixPost;

  rsMatrix<TPar> outMatrixPre;
  rsMatrix<TPar> outMatrixPost;


  std::vector<TPar> dampFactorsPre;
  std::vector<TPar> dampFactorsPost;



  // ToDo: 
  //
  //
  // - Maybe use interpolating delayline to enable experimentation with fractional delays. I think,
  //   using allpass interpolation is most appropriate for this purpose.
  //
  // - Maybe rename to rsFeedbackDelayExplorer
};


template<class TSig, class TPar>
void rsExtendedProtoFDN<TSig, TPar>::setDelays(
  const std::vector<int>& newPreMatrixDelays,
  const std::vector<int>& newPostMatrixDelays)
{
  int N = (int) newPreMatrixDelays.size();
  rsAssert((int) newPostMatrixDelays.size() == N);

  delaysPre.resize(N);
  for(int i = 0; i < N; i++)
  {
    delaysPre[i].setMaxDelayInSamples(newPreMatrixDelays[i]);
    delaysPre[i].setDelayInSamples(   newPreMatrixDelays[i]);
  }

  delaysPost.resize(N);
  for(int i = 0; i < N; i++)
  {
    delaysPost[i].setMaxDelayInSamples(newPostMatrixDelays[i]);
    delaysPost[i].setDelayInSamples(   newPostMatrixDelays[i]);
  }

  // Notes:
  //
  // - We call setMaxDelayInSamples() before calling setDelayInSamples() to ensure that the 
  //   delaylines have enough memory allocated.
}

template<class TSig, class TPar>
bool rsExtendedProtoFDN<TSig, TPar>::areSettingsConsistent() const
{
  bool ok = true;

  int N = getNumDelayChannels();

  ok &= feedbackMatrix.hasShape(N, N);

  ok &= inMatrixPre.getNumRows()  == N;
  ok &= inMatrixPost.getNumRows() == N;

  ok &= outMatrixPre.getNumColumns()  == N;
  ok &= outMatrixPost.getNumColumns() == N;

  ok &= (int) delaysPre.size()  == N;
  ok &= (int) delaysPost.size() == N;

  ok &= (int) dampFactorsPre.size()  == N;
  ok &= (int) dampFactorsPost.size() == N;

  ok &= (int) state.size() == N;

  return ok;
}

template<class TSig, class TPar>
void rsExtendedProtoFDN<TSig, TPar>::processFrame(
  const std::vector<TSig>& inputs, std::vector<TSig>& outputs)
{
  rsAssert(areSettingsConsistent(), "Inconsistent settings in rsExtendedProtoFDN::processFrame()");

  int numIns   = (int) inputs.size();
  int numOuts  = (int) outputs.size();
  int numChans = getNumDelayChannels();

  using Vec = std::vector<TSig>;

  // Form the FDN input by applying the pre-feedback input matrix to the inputs vector:
  Vec x(numChans); rsSetZero(x);     // rsSetZero is superfluous now but maybe later we use a member for x
  for(int i = 0; i < numChans; i++)
    for(int j = 0; j < numIns; j++)
      x[i] += inMatrixPre(i, j) * inputs[j];

  // Form the inputs to the pre feedback matrix delaylines:
  Vec u(numChans);
  for(int i = 0; i < numChans; i++)
    u[i] = x[i] + state[i];

  // Apply the pre feedback matrix delaylines:
  Vec y(numChans);
  for(int i = 0; i < numChans; i++)
    y[i] = dampFactorsPre[i] * delaysPre[i].getSample(u[i]);

  // Apply the feedback matrix:
  Vec v(numChans); rsSetZero(v);
  for(int i = 0; i < numChans; i++)
    for(int j = 0; j < numChans; j++)
      v[i] += feedbackMatrix(i, j) * y[j];

  // Inject inputs into v via the post feedback matrix input matrix:
  for(int i = 0; i < numChans; i++)
    for(int j = 0; j < numIns; j++)
      v[i] += inMatrixPost(i, j) * inputs[j];

  // Apply the post feedback matrix delaylines and save their outputs in the state:
  Vec z(numChans);
  for(int i = 0; i < numChans; i++)
    z[i] = state[i] = dampFactorsPost[i] * delaysPost[i].getSample(v[i]);
  // Using a local vector z is actually superfluous

  // Form the outputs:
  rsSetZero(outputs);
  for(int i = 0; i < numOuts; i++)
  {
    for(int j = 0; j < numChans; j++)
    {
      outputs[i] += outMatrixPre(i, j)  * y[j];
      outputs[i] += outMatrixPost(i, j) * z[j];
    }
  }
  // I think, this should be done first...but maybe only the "pre" half of it should be done first 
  // and the "post" half last?
}

// ToDo:
//
// - Maybe get rid of the assertions in the setters and then maybe split the 
//   setDelaysAndFeedbackMatrix into two functions.
//
// - Make a class rsProtoFDN, i.e. a non-extended FDN that gets rid of the post-matrix delays.
//
// - Let the class rsExtendedProtoFDN have a method to convert itself to a non-extended, i.e. basic 
//   FDN. The basic FDN will have to be of twice the order, I think.
//
// - Implement a getTransferFunctionAt method for both. First for the basic one, then for the 
//   extended one. A first implementation can be based on converting to basic and then using the
//   basic implementation. Eventually, we want to implement it without conversion, though.
//
// - Find a way to compute the poles of the FDNs in a reasonably efficient way. We want to use the
//   mode density as a quality measure and for this, I think, we need to compute the poles.





#endif