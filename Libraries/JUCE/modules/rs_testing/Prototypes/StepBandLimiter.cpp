
template<class TSig, class TTim>
rsStepBandLimiter<TSig, TTim>::rsStepBandLimiter()
{
  // rectangular window by default (...use something better later):
  windowCoeffs[0] = TTim(1);
  for(int i = 1; i < maxNumWindowCoeffs; i++)
    windowCoeffs[i] = TTim(0);

  // test:
  windowCoeffs[0] = TTim(0.5);
  windowCoeffs[1] = TTim(0.5);

  setLength(5);
}

template<class TSig, class TTim>
void rsStepBandLimiter<TSig, TTim>::updateTables()
{
  //int L = delayLength * samplesPerLobe; // table length

  int L = 2 * sincLength * samplesPerLobe + 1; // table length
  timeTbl.resize(L);
  blitTbl.resize(L);
  blitDrvTbl.resize(L);
  blepTbl.resize(L);
  blampTbl.resize(L);


  int ic      = (L-1)/2;  // center index
  timeTbl[ic] = TTim(0);  // time axis in samples
  //blitTbl[ic] = TTim(1);
  blitTbl[ic] = window(TTim(0));
  int i;
  for(i = 1; i <= ic; i++) {
    TTim t = TTim(i) / TTim(samplesPerLobe);  // time in samples
    TTim s = sin(PI*t) / (PI*t); // todo: apply window later
    // todo: apply a window function - use a windowed sinc - try to find analytic expressions for 
    // the integral of the windowed sinc - if none can be found, use numeric integration
    // wolfram can evaluate the indefinite integral of a cosine-term times the sinc:
    // https://www.wolframalpha.com/input/?i=integral+cos(n*pi*x)+*+sin(pi*x)%2F(pi*x)
    // input: Integrate[Cos[n Pi x] (Sin[Pi x]/(Pi x)), x]
    // output: (SinIntegral[(1 + n) Pi x] + SinIntegral[Pi x - n Pi x])/(2 Pi)
    // so if we use a window that is a sum of cosine terms, we will be able to derive an expression
    // for the desired integral (it will be a sum of such terms, each weighted by the coefficient 
    // of the cosine term)

    s *= window(t);

    timeTbl[ic + i] =  t;
    timeTbl[ic - i] = -t;
    blitTbl[ic + i] =  s;
    blitTbl[ic - i] =  s;
  }

  // ..and we also need to fill the blitDrv table with the derivative of the blit - can be 
  // computed analytically...the integrals also (in terms of the Si function)

  // numerically integrate the blit to obtain the blep, choose the integration constant such that
  // the center sample is exactly 0.5:
  rsNumericIntegral(&timeTbl[0], &blitTbl[0], &blepTbl[0],  L, TTim(0));
  TTim c = blepTbl[ic] - 0.5;  
  rsArray::add(&blepTbl[0], -c, &blepTbl[0], L);

  // integrate blep to get blamp:
  rsNumericIntegral(&timeTbl[0], &blepTbl[0], &blampTbl[0], L, TTim(0));
  // maybe use better numeric integration later or find analytic expressions
  // maybe we could estimate the integral symmetrically - i.e. one forward pass, one reverse pass
  // and take the average? does that make any sense at all? dunno - try it!

  //// we actually don't want the blit/blep/blamp itself but rather the residual, i.e. the difference
  //// between them and a naive impulse/step/ramp:
  //blitTbl[ic] -= 1;
  //for(int i = ic; i < L; i++) {
  //  blepTbl[i]  -= 1;
  //  blampTbl[i] -= timeTbl[i];
  //}

  // todo: maybe store only the right half and obtain the left half by symmetry when computing 
  // values


  //GNUPlotter plt;
  ////plt.addDataArrays(L, &timeTbl[0], &blitTbl[0], &blepTbl[0]);
  //plt.addDataArrays(L, &timeTbl[0], &blitTbl[0], &blepTbl[0], &blampTbl[0]);
  //plt.plot();
  // they are a bit inexact - compare with the plots in the blamp-paper - that's probably due to 
  // imperfect numeric integration
  
}

template<class TSig, class TTim>
void rsStepBandLimiter<TSig, TTim>::allocateBuffers()
{
  // corrector buffer needs sincLength+1 samples, delay buffer needs sicnLength samples

  tempBuffer.resize(2*sincLength); 
  // for sampled correction signal - to be spread between delayline and corrector

  delayline.resize(bufferSize);
  corrector.resize(bufferSize);


  // maybe apply the corrector directly to stored past samples in addImpulse and use the 
  // future-corrector when writing into the delayline

  reset();
}

template<class TSig, class TTim>
void rsStepBandLimiter<TSig, TTim>::reset()
{
  rsSetZero(delayline);
  rsSetZero(corrector);
  bufIndex = 0;
}

/*

Background/Ideas:

make a class rsBlep (or rsBlepper, BlepApplicator).or rsStepBandLimiter
-double getSample(double in) produces an output sample at a time as usual
-void addStep(double time, double amplitude) may be called immediately before getSample to announce 
 that a step has occurred in between the previous call to getSample and the one to follow, "time"
 gives the fractional sample instant at which the step occurred (which determines at which phase we
 have to sample the corrector to be added) and "amplitude" gives the size of the step (which 
 determines the multiplication factor by which the corrector has to be scaled)
-the class should manage the required delay (i.e. delay the output by whatever number of samples is 
 needed)
-the user should be able to set the number of samples for the blep
-implementation:
 -use a delayline on N samples
  -setStep adds the corrector signal into the delayline (i.e. adds the scaled and shifted 
   step-corrector into the whole delayline
  -getSample adds to the incoming input a sample from the delayline at current index i, clears the 
   delay-element at the current index i and increments the index - the corrector signal as index i
   has been consumed by getSample
-maybe extend it to allow for supressing jumps in the derivative as well:
 -addDerivativeStep(double time, double amplitude)
 -see: http://dafx16.vutbr.cz/dafxpapers/18-DAFx-16_paper_33-PN.pdf
    http://research.spa.aalto.fi/publications/papers/dafx16-blamp/
-to fill the belp/blamp/blarabola/blolynomial buffer (the delayline), we need a function to 
 evaluate the an approximation of the residual (which is the difference between a heaviside-step 
 and a bandlimited step) at arbitrary times
pseudocode (likely wrong):
void addStep(double time, double amplitude)
{
  for(int i = 0; i < L; i++)  // L is the length of the delayline
    blepBuffer[wrap(bufIndex + i)] += amplitude  * blepResidual(time + i);
}
double getSample(double in)
{
  double y = delayBuffer[bufIndex] + blepBuffer[bufIndex];
  blepBuffer[bufIndex] = 0;            // clear blep at this position - it has been consumed
  delayBuffer[wrap(bufIndex+L)] = in;  // write input into delayline
  bufIndex = wrap(bufIndex + 1);
}
-maybe the blepResidual function could use a function oject rsBlepApproximator
-approximate the blep by writing a sampled sinc-function into an array (let the user select an 
 integer that says, how many samples should be taken per sinc-zero-crossing: 
 1: sample at zero-crossings, 2: at zero-crossings and maxima, etc.
 -a 2nd array is filled with the sinc-derivative values at the same instants (this derivative 
  should be evaluated analytically)
 -the sinc values together with the derivative values can be used for a hermite approximation
 -this piecewise polynomial approximation can be integrated to obtain the blep and integrated twice
  to obtain the blamp (alternatively, the blep can be calculated analytically)
 -maybe apply window to the sinc before differentiating/integrating ..or obtain the resisuals first
  from the unwidowed versions and the window the residual - try both, use whichever gives better 
  results
 -maybe, may compute the blep and blamp residual simultaneously (by evaluating a polynomial and its
  derivative simultaneously)

here it says that the *residual* signal should be tabulated (not the blep/blamp itself)...but 
somehow that seems strange in the case of the blit - should we just zero out the center value of 
the table? my gut feeling is that this wouldn't make much sense...
http://metafunction.co.uk/all-about-digital-oscillators-part-2-blits-bleps/
..for the blep, on the other hand, that approach seems sensible ...maybe the blit is a special 
oddball here because its continuous time version is not an actual function but rather a 
distribution?


maybe make a version that only applies the correction to future samples (minblep) - maybe 
use a (windowed) impulse-, step- and ramp- response of an elliptic filter to window step- 
and ramp-response, subtract the naive versions, apply the window and add the naive versions
back

KVR threads:

https://www.kvraudio.com/forum/viewtopic.php?f=33&t=523056
mystran says:
If you go with the LUT approach, make sure you reorder the BLEP data in the LUT in such a way that 
for any given transition you only need to fetch the minimum number of cache lines

...we should probably keep the unoptimized prototype implementation as is for a unit tests with an 
optimized production version...because the details are a bit tricky

Tutorial: BLEPs (using PolyBLEPs but extensible)
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=398553

Windowing BLEPs
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=511023

Triangle hard sync with BLEP's opinions
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=492537

Hard-Sync with MinBlep / How to manage edge cases
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=489691

minBLEPS once and for all
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=364256

Did any one get Minblep hardsync right?
https://www.kvraudio.com/forum/viewtopic.php?f=33&t=461630

*/