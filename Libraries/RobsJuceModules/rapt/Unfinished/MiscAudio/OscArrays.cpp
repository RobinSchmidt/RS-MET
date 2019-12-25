template<class T>
rsOscArray<T>::rsOscArray(RAPT::rsRatioGenerator<T>* _ratioGenerator)
{
  ratioGenerator = _ratioGenerator;
  setMaxDensity(8);
}

template<class T>
void rsOscArray<T>::updateIncrements()
{
  typedef rsRatioGenerator<T> RG;

  if(numOscs == 1)
  {
    incs[0] = 1;
  }
  else
  {
    ratioGenerator->fillRatioTable(&incs[0], numOscs);
    // actually, if we know that the min is in incs[0] and the max is in incs[numOscs-1], we don't
    // need to search for min/max -> optimize
    // -> factor out a rangeConverversionCoeffs function from transformRange that computes a, b
    // from given x1,x2,y1,y2

    // we get all 1s in case of prime-power and all 0s in case of prime-power-diff - how should be
    // handle this situation? maybe by letting all increments be the same? i think, that would be
    // the natural limit - or wait! no! it's not! what is the limit? figure out! make plots (7
    // graphs for the resulting normalized ratios as functions of the power p)


    // the stuff like taking 1.f and/or taking reciprocal and/or sorting should all be done in
    // rsRatioGenerator ...or maybe we should do it here - first: transform to interval 1..2
    // the take the reciprocal, giving numbers in the interval 0.5...1 but in descending order, if
    // the original array is sorted ascending
    bool invertRatios = false; // make user parameter
    if(invertRatios)
    {
      rsArrayTools::transformRange(&incs[0], &incs[0], numOscs, T(1), T(2));  // ratios in 1...2
      for(int i = 0; i < numOscs; i++)
        incs[i] = T(1) / incs[i];                                        // ratios in 0.5...1
      rsArrayTools::reverse(&incs[0], numOscs);
    }

    // transform the   increments to their actual target range:
    T minRatio = T(1) - T(0.5) * detune;
    T maxRatio = T(1) + T(0.5) * detune;
    rsArrayTools::transformRange(&incs[0], &incs[0], numOscs, minRatio, maxRatio);
  }

  // whith numOscs==1, we get -inf in incs[0] - we may have to treat that as special case

  // multiply in the reference increment:
  for(int i = 0; i < numOscs; i++)
    incs[i] *= refInc;
  // ...hmm - or maybe it would be better to do this in getSample and therefore avoid updating the
  // ratios on freq-changes? oh yes: definitely do that, now that the ratio computation has become
  // so expensive ...actually, we should probably just keep an array of prototype ratios (or just
  // offsets between 0..1) and apply the mapping to the desired inc-range on the fly as well -
  // makes detune modulation more efficient (i.e. less inefficient)
}

template<class T>
void rsOscArray<T>::updateAmplitudes()
{
  T p = T(0.5) * (stereoSpread + T(1));

  //T s = T(1.0) / sqrt(numOscs);              // scaler
  //T s = T(1.0) / numOscs;
  T s = 2.0;


  T evenAmp = s * rsCubicFadeIn(p);
  T oddAmp  = s * rsCubicFadeOut(p);
  // no - i think the cubic fade (approximating the const-power rule) is wrong - we need a linear
  // fade because the left/right signals of a *single* osc is coherent - the const-power rule is
  // only for uncorrelated signals - try it with a single saw - for which rule there's no gap in
  // the middle - hmm - but sonically, the const power rule seems better - at least with
  // headphones -> try with loudspeakers

  //T evenAmp = s * p;
  //T oddAmp  = s * (1-p);


  if(rsIsOdd(numOscs)) {                         // odd density - pan middle osc to center
    int mid = numOscs/2;
    for(int i = 0; i < mid; i += 2) {
      ampsL[i]   = evenAmp;
      ampsR[i]   = oddAmp;
      ampsL[i+1] = oddAmp;
      ampsR[i+1] = evenAmp;
    }
    if(rsIsOdd(mid)) {
      for(int i = 0; i < mid; i += 2) {
        ampsL[mid+i]   = evenAmp;
        ampsR[mid+i]   = oddAmp;
        ampsL[mid+i+1] = oddAmp;
        ampsR[mid+i+1] = evenAmp;
      }
    } else {
      for(int i = 0; i <= mid; i += 2) {
        ampsL[mid+i+1] = evenAmp;
        ampsR[mid+i+1] = oddAmp;
        ampsL[mid+i]   = oddAmp;
        ampsR[mid+i]   = evenAmp;
      }
    }
    ampsL[mid] = ampsR[mid] = s * 0.5;
  } else {                                       // even density - pan all oscs alternatingly
    for(int i = 0; i < numOscs; i += 2) {
      ampsL[i]   = evenAmp;
      ampsR[i]   = oddAmp;
      ampsL[i+1] = oddAmp;
      ampsR[i+1] = evenAmp;
    }
  }



  // there seems to be something wrong with that - when density = 3 and stereo-spread is +-1
  // we have one saw at one side and two on the other - ok - the middle osc is indeed panned center
  // but the outer ones are wrong - the amp arrays are: L: 2,1,2, R: 0,1,0 but we want
  // L: 0,1,2, R: 2,1,0 - we need to start at the middle and alternate the amplitudes outward
  // instead of starting at 0 and alternate upward


  // todo: apply bell curve, maybe we should treat even densities different from odd densities by
  // having one saw in the middle in the odd case? -but hwo would we handle continuous density
  // then?
}

template<class T>
void rsOscArray<T>::reset()
{
  rsNoiseGenerator<T> ng;
  ng.setSeed(startPhaseSeed);
  ng.setRange(T(0), T(1));
  for(int i = 0; i < numOscs; i++) {
    T pos = startPhaseDist * T(i) / T(numOscs);
    pos += startPhaseRand * ng.getSample();
    resetOsc(i, pos);
  }

}

//=================================================================================================

template<class T, class TOsc, class TBlep>
rsBlepOscArray<T, TOsc, TBlep>::rsBlepOscArray(rsRatioGenerator<T>* _ratioGenerator)
  : rsOscArray<T>(_ratioGenerator)
{
  //ratioGenerator = _ratioGenerator;
  //setMaxDensity(8);
  oscs.resize(this->getMaxDensity());
}




/*
Ideas:

massive oscillator bank (MOB) synthesis:
-when looking at the spectrum of a regular supersaw, one notices that the absolute frequency
 spreading gets larger towards higher harmonics - which it obviously must do, since the relative
 freq spread must remain constant
-that leads to the perhaps undesirable effect that the clusters of sines at higher harmonics start
 to overlap - for example, at a fundamenal of 100Hz with a detune of 20%, we have the lowest
 "fundamental" sine at 90 and the highest at 110 - for the 10th harmonic at 1000 Hz, the lowest
 would be 900 and the highest 1100 - so there's overlap with the 9th and 11th harmonic cluster
-it may be desirable to narrow down the freq-spread of the clusters at the higher harmonics, maybe
 in such a way to keep constant absolute cluster width
-the only way i can think of, to do this, it to generate all sines separately
-so, to synthesize a single saw for a fundamental at 100Hz with spectrum up to 20 kHz, we need
 200 sine oscs - if we want 7 saws like in a supersaw, we need 1400 sine oscs - that's a massive
 number of oscillators to use for a single note - hence the name: massive oscillator bank synthesis
-it goes without saying, that the simplemost trigonometric recursion should be used for the sines,
 pretty much ruling out any realtime modulation of parameters
-even then, we should use the newest (widest) possible SIMD vector instructions and use single
 precision float, and/or implement this stuff on the GPU using its natural massive parallelism
 capabilities
-having this in place, we could adjust the frequency spread for each harmonic cluster separately,
 potentially leading to interesting new variants of the "supersaw" with the regular supersaw as
 the special case where the absoulte widths of the harmonic clusters are proportional to the
 harmonic's frequency
-we can even do inharmonic spacing of the partial clusters
-hmm...maybe harmonic cluster (HC) synthesis would also be a fitting name - actually the MOB
 concept is far more general - we do not necessarily need to arrange the frequencies in clusters
 around the harmonics - that's just one possbility - so maybe HC synthesis could be considered
 as subset of MOB synthesis
-actually, this is just additive synthesis with some specific way of macro-controlling the whole
 array of oscs - so maybe it doesn't really deserve a new name - but "MOB synthesis" just sounds
 too cool to not use the term :-D

*/

