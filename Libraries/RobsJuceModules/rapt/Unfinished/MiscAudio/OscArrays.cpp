

template<class T, class TOsc, class TBlep>
rsBlepOscArray<T, TOsc, TBlep>::rsBlepOscArray(rsRatioGenerator<T>* _ratioGenerator)
{
  ratioGenerator = _ratioGenerator;
  setMaxNumOscillators(8);
}


template<class T, class TOsc, class TBlep>
void rsBlepOscArray<T, TOsc, TBlep>::updateIncrements()
{
  typedef rsRatioGenerator<T> RG;

  if(numOscs == 1)
    incs[0] = 1;
  else
    ratioGenerator->fillRatioTable(&incs[0], numOscs);
    // actually, if we know that the min is in incs[0] and the max is in incs[numOscs-1], we don't
    // need to search for min/max -> optimize
    // -> factor out a rangeConverversionCoeffs function from transformRange that computes a, b
    // from given x1,x2,y1,y2

  // the stuff like taking 1.f and/or taking reciprocal and/or sorting should all be done in
  // rsRatioGenerator ...or maybe we should do it here - first: transform to interval 1..2
  // the take the reciprocal, giving numbers in the interval 0.5...1 but in descending order, if 
  // the original array is sorted ascending
  bool invertRatios = false; // make user parameter
  if(invertRatios)
  {
    rsArray::transformRange(&incs[0], &incs[0], numOscs, T(1), T(2));  // ratios in 1...2
    for(int i = 0; i < numOscs; i++) 
      incs[i] = T(1) / incs[i];                                        // ratios in 0.5...1
    rsArray::reverse(&incs[0], numOscs); 
  }

  // transform the   increments to their actual target range:
  T minRatio = T(1) - T(0.5) * detune;
  T maxRatio = T(1) + T(0.5) * detune;
  rsArray::transformRange(&incs[0], &incs[0], numOscs, minRatio, maxRatio);

  // multiply in the reference increment:
  for(int i = 0; i < numOscs; i++)
    incs[i] *= refInc;
  // ...hmm - or maybe it would be better to do this in getSample and therefore avoid updating the
  // ratios on freq-changes? oh yes: definitely do that, now that the ratio computation has become
  // so expensive ...actually, we should probably just keep an array of prototype ratios (or just 
  // offsets between 0..1) and apply the mapping to the desired inc-range on the fly as well - 
  // makes detune modulation more efficient (i.e. less inefficient)
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
 would be 900 and the highst 1100 - so there's overlap with the 9th and 11th harmonic cluster
-it may be desirable to narrow down the freq-spread of the clusters at the higher harmonics
-the only way i can think of, to do this, it to generate all sines separately
-so, to synthesize a single saw for a fundamental at 100Hz with spectrum up to 20 kHz, we need
 200 sine oscs - if we want 7 saws, we need 1400 sine oscs - that's a massive number of 
 oscillators to use for a single note - hence the name: massive oscillator bank synthesis
-it goes without saying, that the simplemost trigonometric recursion should be used, pretty much
 ruling out any realtime modulation of parameters
-even then, we should use the newest (widest) possible SIMD vector instructions and use single 
 precision float, and/or implement this stuff on the GPU using its natural massive parallelism 
 capabilities
-having this in place, we could adjust the frequency spread for each harmonic cluster separately,
 potentially leading to interesting new variants of the "supersaw" with the regular supersaw as
 the special case where the absoulte widths of the harmonic clusters are proprtional to the
 harmonic's frequency
-we can even do inharmonic spacing of the partial clusters
-hmm...maybe harmonic cluster (HC) synthesis would also be a fitting name - actually the MOB 
 concept is far more general - we do not necessarily need to arrange the frequencies in clusters
 around the harmonics - that's just one possbility - so maybe HC synthesis could be considered
 as subset of MOB synthesis








*/

