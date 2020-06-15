#ifndef RAPT_AUDIOFUNCTIONS_H
#define RAPT_AUDIOFUNCTIONS_H

// todo: maybe wrap into class, get rid of redundant implementations
// merge with other AudioFunctions.h file make another file AudioAnalysisFunctions

  /** Converts a raw amplitude value/factor to a value in decibels. */
template<class T>
RS_INLINE T rsAmp2dB(T amp);

/** Converts a raw amplitude value/factor to a value in decibels with a check, if the amplitude
is close to zero (to avoid log-of-zero
and related errors). */
template<class T>
RS_INLINE T rsAmp2dBWithCheck(T amp, T lowAmplitude = 0.000001);


/** Converts a value in decibels to a raw amplitude value/factor. */
template<class T>
RS_INLINE T rsDB2amp(T x);

/** Computes the centroid of the length N array x. */
template<class T>
T rsCentroid(T *x, int N); // maybe move to rsArrayTools

/** Computes the centroid of energy the length N array x. */
template<class T>
T rsCentroidOfEnergy(T *x, int N); // maybe move to rsArrayTools

/** Given 2 successive sample value of a sinusoidal function with known (normalized radian)
frequency w, such that:
y0 = a * sin(p)
y1 = a * sin(p+w)
this function computes the amplitude a and initial phase p which are passed as pointer
parameters to the function. The computed phase is in the range -PI...PI */
template<class T>
void rsSineAmplitudeAndPhase(T y0, T y1, T w, T *a, T *p);

/** Given 3 successive sample values of a sinusoidal function with arbitrary amplitude a and
initial phase p, such that
y0 = a * sin(p)
y1 = a * sin(p+w)
y2 = a * sin(p+2*w)
this function computes the (normalized radian) frequency w of the sinusoid. The problem is ill
defined, if y1 is (numerically) equal to zero, so the caller is responsible for making sure, that
y1 != 0. In the computations, there's a division by y1 and from the symmetry of the sine, it can
be seen that with y1 = 0, we would also have y0 = -y2 and no matter what the values of y0 and y2
actually are, we could assign any arbitrary value to the frequency w and adjust the amplitude a
in a way, to compensate - so, indeed, there's no unique solution in this case, unless the
amplitude is known in which case w = asin(|y2|/a) = asin(|y0|/a) [verify this]. The last
parameter "small" determines how close y1 may be to 0 until the computations are considered to be
too error-prone. For a sine with changing frequency, the estimate should be considered to be the
instantaneous value at the center sample, i.e. at y1. */
template<class T>
T rsSineFrequency(T y0, T y1, T y2, T smalll = 1.e-8);
// "smalll" because "small" is defined as "char" somewhere leading to weird compiler errors
// maybe rename y0,y1,y2 to yL,yC,yR (for left,center,right)

/** Assuming a sinusoidal input signal x of length N (and unknown frequency), this function
computes the sinusoid's instantaneous normalized radian frequency at sample-instant n0. The value
is more accurate, if the frequency varies less over 3 successive samples centered around n0 and
it's also more accurate, if n0 corresponds to a peak or valley in the signal. If it's not
possible to use 3 values centered around n0 (because n0 = 0 or n0 = N-1 and/or the sample value
of n0 is close to zero), the actual measurement point might be shifted (by at most 2 samples)
which should be inconsequential, if the sinuosoid has a stable frequency. The function is used
internally by rsSineFrequencyAt which is probably the function, you want to use instead. */
template<class T>
T rsSineFrequencyAtCore(T *x, int N, int n0, T smalll = 1.e-8);
  // find a better name

/** Assuming a sinusoidal input signal x of length N (and unknown frequency), this function
computes the sinusoid's instantaneous normalized radian frequency at sample-instant n0. It will
use the function rsSineFrequencyAtCore at the peak/valley before and after the sample n0 (because at
at peaks and valleys rsSineFrequencyAtCore returns the most accurate values) and return a
weighted average of these values. If the sample index n0 is not surrounded by two peaks/valleys
(because it's close the start or end of the signal) two peaks/valleys before and after will be
used and the frequency value will be linearly extrapolated. This implies, for this function to
work, the caller must make sure that there are at least one peak and one valley in the input
signal. */
template<class T>
T rsSineFrequencyAt(T *x, int N, int n0, bool refine = false);
 // todo: give the function another parameter that crosfades between arithmetic and geometric 
 // mean/extrapolation

/** Assuming a sinusoidal input signal x of length N and normalized radian frequency w, this
function computes the sinusoid's instantaneous phase at sample-instant n0. */
template<class T>
T rsSinePhaseAt(T *x, int N, int n0, T w);

/** Assuming a sinusoidal input signal x of length N (and unknown frequency), this function
computes the sinusoid's instantaneous phase at sample-instant n0. */
template<class T>
T rsSinePhaseAt(T *x, int N, int n0);


/** Finds an index in the length N array x close to n0, where the signal x crosses zero. The
searchDirection parameter determines whether we search into the leftward (-1) or rightward (+1)
direction and the upward/downward parameters determine, whether upward and/or downward crossings
should count. if no zero crossing could be found, -1 will be returned. */
template<class T>
int rsFindZeroNear(T *x, int N, int n0, int searchDirection, bool upward, bool downward);
// move to rsZeroCrossingFinder


/**
NOTE: this function works only if the sample index n0 is between two zero-crossings of the signal
x
*/
template<class T>
T rsSinePhaseAtViaZeros(T *x, int N, int n0, int precision);


/** Assuming a sinusoidal input signal x of length N with normalized radian frequency w, this
function computes an offset in samples by which the sinusoid must be shifted, such that the
sinusoid's phase at sample instant n0 equals p0. The shift will be either forward or backward,
whichever is smaller. x: signal, N: signal length, n0: sample-instant for target phase p0,
p0: target phase at n0, return value: number of samples to shift (might be positive for
rightshift or negative for leftshift). */
template<class T>
T rsSineShiftAmount(T *x, int N, int n0, T p0, T w);

/** Similar to rsSineShiftAmount(T *x, int N, int n0, T p0, T w) but it works for
sinusoids with unknown frequency. It measures the sinusoids frequency first using
rsSineFrequencyAt and then uses the function that works for known frequencies. */
template<class T>
T rsSineShiftAmount(T *x, int N, int n0, T p0);


//-----------------------------------------------------------------------------------------------
// implementation:

template<class T>
RS_INLINE T rsAmp2dB(T amp)
{
  return T(8.6858896380650365530225783783321) * log(amp);
}

template<class T>
RS_INLINE T rsAmp2dBWithCheck(T amp, T lowAmplitude)
{
  if(amp >= lowAmplitude)
    return rsAmp2dB(amp);
  else
    return rsAmp2dB(lowAmplitude);
}


template<class T>
RS_INLINE T rsDB2amp(T dB)
{
  return exp(dB * T(0.11512925464970228420089957273422));
}

#endif
