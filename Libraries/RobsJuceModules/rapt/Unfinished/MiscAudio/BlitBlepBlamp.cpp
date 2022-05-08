
template<class TSig, class TTim>
rsTableBlep<TSig, TTim>::rsTableBlep()
{
  // rectangular window by default (...use something better by default later):
  windowCoeffs[0] = TTim(1);
  for(int i = 1; i < maxNumWindowCoeffs; i++)
    windowCoeffs[i] = TTim(0);

  // test - Hanning window:
  windowCoeffs[0] = TTim(0.5);
  windowCoeffs[1] = TTim(0.5);
}

//=================================================================================================

template<class TSig, class TTim>
rsTableLinBlep<TSig, TTim>::rsTableLinBlep()
{
  setLength(5);
}

template<class TSig, class TTim>
void rsTableLinBlep<TSig, TTim>::updateTables()
{
  std::vector<TTim> timeTbl;  // needed temporarily

  int L  = 2 * halfLength * this->tablePrecision + 1; // nominal table length
  int Lg = L+1;                                       // table length with guard
  timeTbl.resize(Lg);
  this->blitTbl.resize(Lg);
  this->blepTbl.resize(Lg);
  this->blampTbl.resize(Lg);

  int ic      = (L-1)/2;  // center index
  timeTbl[ic] = TTim(0);  // time axis in samples
  this->blitTbl[ic] = this->window(TTim(0));
  int i;
  for(i = 1; i <= ic; i++) {
    TTim t = TTim(i) / TTim(this->tablePrecision);  // time in samples
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

    s *= this->window(t);

    timeTbl[ic + i] =  t;
    timeTbl[ic - i] = -t;
    this->blitTbl[ic + i] =  s;
    this->blitTbl[ic - i] =  s;
  }

  // ..and we also need to fill the blitDrv table with the derivative of the blit - can be
  // computed analytically...the integrals also (in terms of the Si function)

  // numerically integrate the blit to obtain the blep:
  rsNumericIntegral(&timeTbl[0], &this->blitTbl[0], &this->blepTbl[0], L, TTim(0));
  // do the numerical integration only up to L/2+1 and obtain the rest via symmetry (-> less
  // accumulation of numerical integration error)
  // or maybe even store only the right half and obtain the left half by symmetry when computing
  // values

  // scale it such that the last sample is exactly 1 (the first sample is already exactly zero by
  // construction):
  //rsArrayTools::scale(&blepTbl[0], &blepTbl[0], L, TTim(1)/rsLast(blepTbl));*
  rsArrayTools::scale(&this->blepTbl[0], &this->blepTbl[0], L, TTim(1)/this->blepTbl[L-1]);

  // integrate blep to get blamp:
  rsNumericIntegral(&timeTbl[0], &this->blepTbl[0], &this->blampTbl[0], L, TTim(0));
  // maybe use better numeric integration later or find analytic expressions

  //rsPlotVectors(blitTbl, blepTbl, blampTbl);

  // We actually don't want the blep/blamp itself but rather the residual, i.e. the difference
  // between them and a naive step/ramp ..but i think, that works only for blep and higher order
  // integrals - because the continuous delta distribution is not actually a function ...or
  // something like that:
  for(int i = ic; i < L; i++) {
    //blepTbl[i]  -= 1;
    this->blampTbl[i] -= timeTbl[i];
  }
  //blepTbl[ic] = 0;
  // is this correct? supposed to fix spurious peaks when linearly interpolating
  // the table at the discontinuity (of the table)

  // OK - it seems that also for the blep, we cannot just tabulate the residual because the
  // residual contains a discontinuity which causes artifacts when trying to linearly interpolate
  // the residual-table near this discontinuity - so for the moment, we also store the blep
  // directly and obtain the residual in the discrete time domain


  //rsPlotVectors(blitTbl, blepTbl, blampTbl);

  // the guard samples in the tables just repeat the last actual samples:
  this->blitTbl[L]  = this->blitTbl[L-1];
  this->blepTbl[L]  = this->blepTbl[L-1];
  this->blampTbl[L] = this->blampTbl[L-1];

  //rsPlotVectors(blitTbl, blepTbl, blampTbl);

  //GNUPlotter plt;
  //plt.addDataArrays(L, &timeTbl[0], &blitTbl[0], &blepTbl[0]);
  ////plt.addDataArrays(L, &timeTbl[0], &blitTbl[0], &blepTbl[0], &blampTbl[0]);
  //plt.plot();
  ////// they are a bit inexact - compare with the plots in the blamp-paper - that's probably due to
  ////// imperfect numeric integration

  // the timeTbl is actually not needed anymore now - it can be cleared - or maybe it should be a
  // local variable here - we actually only need it, because the numerical integration routine
  // needs an x-axis to be passed - maybe write another (simpler) integration routine that doesn't
  // need an x axis (and instead assumes equal spacing of the y-values by some distance h)
  // ...ah - and we need it also when creating the blamp residual - but maybe the time-values
  // should be computed on the fly (again - they are already computed in the loop creating the
  // sinc)...hmm - maybe a local array is better indeed to avoid the recomputation...


  // actually, it seems from the plots, that for the blamp, a much shorter kernel length would
  // be sufficient than for the blep - which makes sense - so maybe the blamp should actually be
  // a separate class - maybe make 3 classes rsTableBlit/Blep/Blamp and also rsPolyBlit/Blep/Blamp
  // on the other hand, if a blep and blamp should be applied at the same time (for example, when
  // hard-syncing triangle waves), it *does* make sense to have one object that does both
  // simultaneously - especially because of the consolidation of the delay
}

template<class TSig, class TTim>
void rsTableLinBlep<TSig, TTim>::allocateBuffers()
{
  // corrector buffer needs sincLength+1 samples, delay buffer needs sicnLength samples

  tempBuffer.resize(2*halfLength);
  // for sampled correction signal - to be spread between delayline and corrector

  this->delayline.resize(this->bufferSize);
  this->corrector.resize(this->bufferSize);

  // maybe apply the corrector directly to stored past samples in addImpulse and use the
  // future-corrector when writing into the delayline

  reset();
}

template<class TSig, class TTim>
void rsTableLinBlep<TSig, TTim>::reset()
{
  rsSetZero(this->delayline);
  rsSetZero(this->corrector);
  this->bufIndex = 0;
}

//=================================================================================================

template<class TSig, class TTim>
rsTableMinBlep<TSig, TTim>::rsTableMinBlep()
{
  setLength(32);
}

template<class TSig, class TTim>
void rsTableMinBlep<TSig, TTim>::updateTables()
{
  int L = this->blepLength * this->tablePrecision + 1; // why +1? guard for interpolator? ...nope...but...
  //int L = blepLength * tablePrecision;   // ...without, it crashes for blepLength == 0
  int Lg = L+1;                            // this +1 here is for the guard sample
  this->blitTbl.resize(Lg);
  this->blepTbl.resize(Lg);
  this->blampTbl.resize(Lg);

  // create temporary elliptic subband filter object - we use its impulse/step/ramp response for
  // the tables (maybe later allow the user to set up the EngineersFilter settings):
  rsEllipticSubBandFilter<TTim, TTim> flt;
  flt.setSubDivision(this->tablePrecision);
  // later: tablePrecision * cutoffScaler ..or maybe use flt.setCutoffScaler (that function doesn't
  // exist yet - rsEllipticSubBandFilter uses a fixed scaler of 0.9)

  // create time axis and (unilateral) window:
  std::vector<TTim> timeAxis(L), wnd(L);
  TTim ts = TTim(1) / TTim(this->tablePrecision);  // time axis scaler
  for(int n = 0; n < L; n++) {
    timeAxis[n] = ts * TTim(n);
    wnd[n] = this->window(0.5 * timeAxis[n]);      // 0.5 because of unilaterality
  }

  // create blit:
  this->blitTbl[0] = wnd[0] * (flt.getSample(TTim(this->tablePrecision))); // correct to scale by tablePrecision?
  for(int n = 1; n < L; n++)
    this->blitTbl[n] = wnd[n] * flt.getSample(TTim(0));

  // create blep-resiudal:
  flt.reset();
  for(int n = 0; n < L; n++)
    this->blepTbl[n] = wnd[n] * ( flt.getSample(TTim(1)) - TTim(1) );

  // create blamp-residual:
  flt.reset();
  for(int n = 0; n < L; n++)
    this->blampTbl[n] = wnd[n] * (flt.getSample(timeAxis[n]) - timeAxis[n]);

  // higher order residuals could be computed as:
  // blarabola: flt.getSample( t^2/2 ) - t^2/2;
  // blubic:    flt.getSample( t^3/6 ) - t^3/6;
  // etc.: t^k / k!

  // the guard samples in the tables just repeat the last actual samples:
  this->blitTbl[L]  = this->blitTbl[L-1];
  this->blepTbl[L]  = this->blepTbl[L-1];
  this->blampTbl[L] = this->blampTbl[L-1];

  //rsPlotVectors(blitTbl, blepTbl);
  //rsPlotVectors(blitTbl, blepTbl, blampTbl, wnd);
  // blamp looks weird but i think, it's ok - the lowpassed version is sort of shifted
  // (to the right/bottom) with respect to the perfect ramp
}

template<class TSig, class TTim>
void rsTableMinBlep<TSig, TTim>::allocateBuffers()
{
  this->corrector.resize(this->bufferSize);
  reset();
}

template<class TSig, class TTim>
void rsTableMinBlep<TSig, TTim>::reset()
{
  rsSetZero(this->corrector);
  this->bufIndex = 0;
}















//=================================================================================================

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

also interesting: polynomial transtition regions (PTR):
https://www.yofiel.com/software/cycling-74-patches/antialiased-oscillators
http://research.spa.aalto.fi/publications/papers/spl-ptr/
http://home.mit.bme.hu/~bank/publist/smc13.pdf

"New Perspectives on Distortion Synthesis for Virtual Analog Oscillators" Computer Music Journal, 2010.http://eprints.maynoothuniversity.ie/4104/1/VL_New_perspectives.pdf


maybe make a version that only applies the correction to future samples (minblep) - maybe
use a (windowed) impulse-, step- and ramp- response of an elliptic filter to window step-
and ramp-response, subtract the naive versions, apply the window and add the naive versions
back

Resources:

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

Stefan Stenzel - The amazing usefulness of band limited impulse trains, shown for oscillator banks..
https://www.youtube.com/watch?v=lpM4Tawq-XU




*/
