#ifndef RAPT_AUDIOFUNCTIONS_H_INCLUDED
#define RAPT_AUDIOFUNCTIONS_H_INCLUDED

/** Functions that are specific to audio-processing and musical applications such as 
conversion between amplitudes and decibels, frequencies and midi-pitches, beats and seconds, 
etc. */

// todo: merge with the file in "Unfinished"

/** Converts a raw amplitude value/factor to a value in decibels. */
template<class T>
inline T rsAmpToDb(T amp)
{
  return T(8.6858896380650365530225783783321) * log(amp);
}

/** Converts a raw amplitude value/factor to a value in decibels with a check, if the amplitude
is close to zero (to avoid log-of-zero
and related errors). */
template<class T>
inline T rsAmpToDbWithCheck(T amp, T lowAmplitude)
{
  if( amp >= lowAmplitude )
    return rsAmpToDb(amp);
  else
    return rsAmpToDb(lowAmplitude);
}

/** Returns true, if p1 and p2 are a multiple of 2*pi apart (within some given tolerance). */
template<class T>
bool rsArePhasesConsistent(T p1, T p2, T tol = 1.e-13);

/** Converts a time-stamp given in beats into seconds acording to a tempo measured in beats per
minute (bpm). */
template<class T>
inline T rsBeatsToSeconds(T beat, T bpm)
{
  return (60 / bpm)*beat;
}

/** Converts a value in decibels to a raw amplitude value/factor. */
template<class T>
inline T rsDbToAmp(T dB)
{
  return exp(dB * T(0.11512925464970228420089957273422));
}

/** Given a value x between 0 and 1, this function returns a value of a cubic polynomial that 
can be used as fade-in function. The polynomial satisfies: y(0)=0, y'(0)=pi/2, y(1)=1, y'(1)=0.
Together with rsCubicFadeOut, the curve approximates a constant power (cross)fade. */
template<class T>
inline T rsCubicFadeIn(T x)
{
  //return x*(1+x*(1-x)); // slope = 1 at x = 0
  return x*(x*((PI/2-2)*x+(3-PI))+(PI/2));
}

/** Given a value x between 0 and 1, this function returns a value of a cubic polynomial that 
can be used as fade-out function. The polynomial satisfies: y(0)=1, y'(0)=0, y(1)=0, y'(1)=-pi/2.
Together with rsCubicFadeIn, the curve approximates a constant power (cross)fade. */
template<class T>
inline T rsCubicFadeOut(T x)
{
  //return x*x*(x-2)+1;   // slope = -1 at x = 1
  return x*x*((2-PI/2)*x+(PI/2-3))+1;
}

/** NOT YET USABLE (buggy)!

Given a tentative unwrapped phase value (in 0..inf), and a wrapped target phase value 
(in 0..2pi), this function computes a phase that is in the neighbourhood of the tentative value but 
also a multiple of 2*pi above the targetPhase value, such that it is consistent with the target 
phase value. The actual returned value will be the one that is closest to the original tentative 
phase that also satisfies the consistency criterion of being k*2*pi away from targetPhase. */
//template<class T>
//T rsFindCosistentPhase(T targetPhase, T tentativePhase);
// maybe switch argument order, let the targetPhase be in the range -PI..+PI instead of 0..2pi


/** Assume to have some interval of values, where values that fall outside the interval are 
identified with some value from within the interval - for example, a range of angles from 0 to 
2*pi - any angle x above 2*pi or below 0 would be identified with some number y = x + 2*k*pi for 
some suitably chosen k. What this function does, is to take a tentative unwrapped values (in the 
range -inf...+inf and returns another unwrapped value that is consistent with some target value
from the base-range - consistent in the sense that is a multiple of the range-size away from it.
It is also consistent with the tentative unwrapped value in the sense that is at most 
(rangeMax-rangeMin)/2 away from that. So, the returned value is in the range 
preliminaryUnwrappedValue +- (rangeMax-rangeMin)/2 but also a multiple of (rangeMax-rangeMin) away
from targetWrappedValue. This function is useful for finding unwrapped phase values that should be 
consistent with some measured phase-value from the base-range (of 0..2pi, say) but should also be 
in the neighbourhood of some unwrapped value that may be the result of integrating/accumulating a 
frequency over time. */
template<class T>
T rsConsistentUnwrappedValue(T preliminaryUnwrappedValue, T targetWrappedValue, 
  T rangeMin, T rangeMax);

/** ..as above but supposes rangeMin to be zero. Uses a different algorithm based on fmod */
template<class T>
T rsConsistentUnwrappedValue0(T preliminaryUnwrappedValue, T targetWrappedValue, T rangeSize);
// rename to rsUnwrap

/** Computes the integer factor k such that abs(v+k*p - t) becomes smallest. It uses a loop trying
different values of k, so to keep the loop short, you should pass an initial guess for k that is 
likely close to the actually desired k. Parameters: v: value, t: target, p: period, k: guess. */
template<class T>
int rsUnwrapFactor(T v, T t, T p, int k)
{
  while(rsAbs((v+(k*p))-t) > rsAbs((v+((k+1)*p))-t))  k++;
  while(rsAbs((v+(k*p))-t) > rsAbs((v+((k-1)*p))-t))  k--;
  return k;
}

/** Phase distance between given phase p and given target phase t. That's the absolute value of
p-t or (p+2pi)-t or (p-2pi)-t or generally of (p+2*k*pi)-t for a k chosen such that the absolute
value of that difference becomes a minimum. */
template<class T>
T rsPhaseDistance(T p, T t)
{
  int k = rsUnwrapFactor(p, t, 2*PI, 0);
  return rsAbs(p+2*k*PI - t);
}
// needs test


/** Converts a frequency in Hz into a MIDI-note value. It can be used also for tunings different 
than the default the 440 Hz. */
template<class T>
inline T rsFreqToPitch(T freq, T masterTuneA4 = T(440))
{
  //return 12.0 * rsLog2(freq / masterTuneA4) + 69.0;  // old - produced compiler warnings
  return T(12) * rsLog2(freq / masterTuneA4) + T(69);  // new - may be less precise for T=float
}
// ToDo: templatize rsLog2. It currently uses double for the data type

/** Returns, how far two phase values are apart after wrapping them both into the interval
0..2pi. The returned value will be in 0..2pi. */
template<class T>
T rsPhaseError(T p1, T p2);
// hmm - shouldn't it be in 0...pi? when the error above pi, we could use 2*PI - error. this would 
// mean, we look into the otherdirection and it would make the function symmetrical in its 
// arguments

/** Converts a pitch-offset in semitones value into a frequency multiplication factor. Uses one 
multiplication and one call to exp.*/
template<class T>
inline T rsPitchOffsetToFreqFactor(T pitchOffset)
{
  return exp(T(0.057762265046662109118102676788181) * pitchOffset);
  //return pow(2.0, pitchOffset/12.0); // naive, slower but numerically more precise
}
// rename to rsPitchShiftToFreqRatio

/** Converts a MIDI-note value into a frequency in Hz assuming A4 = 440 Hz. Uses two 
multiplications and one call to exp. */
template<class T>
inline T rsPitchToFreq(T pitch)
{
  return T(8.1757989156437073336828122976033 * exp(0.057762265046662109118102676788181 * pitch));
  //return 440.0*( pow(2.0, (pitch-69.0)/12.0) ); // naive, slower but numerically more precise
}
// todo: make a function rsPitchToOmega(T pitch, T sampleRate)
// maybe compute the coefficients as constexpr here to document the formulas used - something with
// a * pow(b, c) = a * exp(log(b) * c) = ...something -> look it up!
// Maybe it can be further optimized by replacing exp by exp2? Maybe it's cheaper to compute 2^x 
// than e^x? Do some benchmarks! If so, then there are a lot of places where this optimization can
// be applied. See: https://en.cppreference.com/w/cpp/numeric/math/exp2

/** Converts a MIDI-note value into a frequency in Hz for arbitrary master-tunings of A4. Uses 
three multiplications and one call to exp. */
template<class T>
inline T rsPitchToFreq(T pitch, T masterTuneA4)
{
  return masterTuneA4 * 0.018581361171917516667460937040007
    * exp(0.057762265046662109118102676788181 * pitch);
}

/** Converts a time value in seconds into a time value measured in beats. */
template<class T>
inline T rsSecondsToBeats(T timeInSeconds, T bpm)
{
  return timeInSeconds * (bpm/60.0);
}

/** Computes a frequency multiplication factor (for e.g. filter cutoff frequencies, oscillators, 
etc.) based on midi note key and velocity info and corresponding key- and vel-tracking parameters.
It is meant to consolidate the used formula for such tracking functionality into one place. The key
and vel parameters are midi values (i.e. 0..127) and keytrack/veltrack are in percent with the 
following meaning: At the reference key and vel (both at 64), the factor is 1.0, i.e. the frequency
should be unchanged. At 100% keytrack will respond with a 1 octave per octave characteristic, i.e.
1 octave above (i.e. 12 semitones above refKey=64, i.e. at key = 76) will give rise to a factor of 
2. For veltrack, 100% means that a factor 2 will be produced at full velocity (vel=127) and 0.5 
will be produced at the lowest nonzero velocity (vel=1). That makes sense because vel=0 is used for
noteOff anyway. VERIFY!  ...TBC... */
template<class T>
inline T rsMidiKeyAndVelToFreqFactor(int key, int vel, T keytrack, T veltrack, 
  int refKey = 69, int refVel = 64)
{
  return   pow( T(2), (T(0.01/12.0)*keytrack)*(key-T(refKey)) ) 
         * pow( T(2), (T(0.01/63.0)*veltrack)*(vel-T(refVel)) );

  // Notes:
  // -refKey = 69 was chosen as default because it corresponds to A4 with 440 Hz which is used as
  //  reference pitch for musical tuning. Using that as reference implies that when we dial in 
  //  440 Hz for some oscillator or filter frequency then it will indeed still be tuned to 440 Hz
  //  even with nonzero keytrack when A4 is being played. That is convenient for dialing in 
  //  frequencies. Other options would be C4 = 60 because this is commonly used as reference or 
  //  E4 = 64 because it's in the middle. But then we would have to dial in ugly noninteger 
  //  frequencies to get the tuning right.
  // -refVel = 64 was chosen as default because it is right in the middle.
}
// Needs unit tests (seems to work fine in practice, though)
//
// The formula was adapted from  rosic::MultiModeFilter::updateFreqWithKeyAndVel  and factored out
// into into a RAPT library function with the intention to consolidate the formula used for 
// key/vel-tracking of frequency parameters into one central place - here.
// ToDo: 
// -Call it from there and find other places where the formula is currently implemented 
//  directly and it can be replaced by a call to this. But there, the reference key is 64. Maybe
//  let the caller pass it as parameter
// -Try to optimize to make use of exp instead of pow using a^x = e^(x*log(a)) and use also 
//  2^x * 2^y = 2^(x+y). We can replace 2 pow-calls with one exp-call. That's a good optimization.


/** Returns the frequency ratio (with respect to the fundamental) of a stiff string, such as a 
piano string, for the given harmonic number. The stiffness parameter controls the amount of 
inharmonicity. At zero, the ratio is strictly harmonic, i.e. an integer.  
References:
The Physics of Musical Instruments, 2nd Edition, Eq. 2.67a (p.64)
 see also 2.67b (p.65), Eq. 12.5 (p.363) ...maybe that last one is better - it normalizes the 
 fundamental back to 1 (i think)
http://www.simonhendry.co.uk/wp/wp-content/uploads/2012/08/inharmonicity.pdf  Eq.10
http://www.jbsand.dk/div/StivStreng.pdf  
*/
template<class T>
inline T rsStiffStringFreqRatio(T harmonicNumber, T stiffness)
{
  T n = harmonicNumber;
  T B = stiffness;
  return n*sqrt(1+B*n*n);
}
// todo: maybe implement formula 24 from http://www.jbsand.dk/div/StivStreng.pdf - fix the 1st
// partial at 1 and let the user set the frequency of the 2nd partial - all other partial 
// frequencies follow from these choices

/** Converts a time-constant (typically denoted as "tau") of an exponential decay function to
the time-instant at which the decay reaches the given level in decibels. For example
tauToDecayTime(tau, -60.0) returns the time required to fall to -60 dB for the given value of
tau. */
template<class T>
inline T rsTauToDecayTime(T tau, T decayLevel)
{
  static const T k = LN10 / 20.0;
  return -k*tau*decayLevel;
}

/** A waveform that morphs between saw-down (p = -1), triangle (p = 0) and saw-up (p = +1) - but 
note that the limiting values (+-1) should not actually be used because of divisions by zero, so
use a range -0.99...+0.99 instead for example. It expects the phase argument x in 0..2*pi such that
it has the same API as std::sin, etc. */
template<class T>
T rsTriSaw(T x, T p)
{
  x *= T(1/(2*PI));   // We expect the phase argument x in 0..2*pi but the code below needs 0..1.


  //x  = fmod(x, 1);    // We also allow x to be off by a multiple of the period. ..old - wrong
  x = rsWrapAround(x, 1.0); // New


  T r2 = 0.25*(p+1);
  if(x < r2)
    return x / r2;
  else if(x > 1-r2) 
    return (x-1)/r2;
  else
    return 2*(r2-x)/(1-2*r2) + 1;
}
// todo: introduce flat top and/or bottom - maybe it's simplest to just amplify and clip the 
// waveform, maybe make an oscillator class that pre-computes 1/r2 and maybe also allows for a 
// p-range -1..+1 by somehow taking care of div-by-0 itself
// maybe make a simplified version:
// if(x < h)
//   return a0 + a1*t;
// else
//   return b0 + b1*t;   
// -> no divisions at all, h should the be turning point in the cycle, i.e. a number 0..1 that 
// separates upward and downward "half"-waves (they are not actually half a cycle long) and the
// coeffs are computed according to the desired output range - it can be nicely extendend to higher
// order polynomial shapes as well (for example to introduce smoothness constraints)
// -Make a variant of the function that expects the argument in 0..1. Maybe call it rsTriSaw_01.
//  We can then get rid of the x *= (1/(2*PI));
// -Verify if fmod behaves correctly for negative x.



/** Converts a time-stamp given in whole notes into seconds according to a tempo measured
in beats per minute (bpm). */
template<class T>
inline T rsWholeNotesToSeconds(T noteValue, T bpm)
{
  return (240.0/bpm) * noteValue;
}

#endif