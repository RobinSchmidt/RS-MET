#ifndef RS_AUDIOFUNCTIONS_H
#define RS_AUDIOFUNCTIONS_H

namespace RSLib
{

  /** Converts a raw amplitude value/factor to a value in decibels. */
  RS_INLINE double rsAmp2dB(double amp);

  /** Converts a raw amplitude value/factor to a value in decibels with a check, if the amplitude
  is close to zero (to avoid log-of-zero
  and related errors). */
  RS_INLINE double rsAmp2dBWithCheck(double amp, double lowAmplitude = 0.000001);

  /** Converts a time-stamp given in beats into seconds acording to a tempo measured in beats per
  minute (bpm). */
  RS_INLINE double rsBeatsToSeconds(double beat, double bpm);

  /** Converts a value in decibels to a raw amplitude value/factor. */
  RS_INLINE double rsDB2amp(double x);

  /** Computes the centroid of the length N array x. */
  RSLib_API double rsCentroid(double *x, int N);

  /** Computes the centroid of energy the length N array x. */
  RSLib_API double rsCentroidOfEnergy(double *x, int N);

  /** Given a value x between 0 and 1, this function returns a value of a cubic polynomial that 
  can be used as fade-in function. The polynomial satisfies: y(0)=0, y'(0)=pi/2, y(1)=1, y'(1)=0.
  Together with rsCubicFadeOut, the curve approximates a constant power (cross)fade. */
  RS_INLINE double rsCubicFadeIn(double x);

  /** Given a value x between 0 and 1, this function returns a value of a cubic polynomial that 
  can be used as fade-out function. The polynomial satisfies: y(0)=1, y'(0)=0, y(1)=0, y'(1)=-pi/2.
  Together with rsCubicFadeIn, the curve approximates a constant power (cross)fade. */
  RS_INLINE double rsCubicFadeOut(double x);

  /** Converts a frequency in Hz into a MIDI-note value. It can be used also for tunings different 
  than the default the 440 Hz. */
  RS_INLINE double rsFreqToPitch(double freq, double masterTuneA4 = 440.0);

  /** Converts a picth-offset in semitones value into a frequency multiplication factor. */
  RS_INLINE double rsPitchOffsetToFreqFactor(double pitchOffset);

  /** Converts a MIDI-note value into a frequency in Hz assuming A4 = 440 Hz. */
  RS_INLINE double rsPitchToFreq(double pitch);

  /** Converts a MIDI-note value into a frequency in Hz for arbitrary master-tunings of A4. */
  RS_INLINE double rsPitchToFreq(double pitch, double masterTuneA4);

  /** Converts a time value in seconds into a time value measured in beats. */
  RS_INLINE double rsSecondsToBeats(double timeInSeconds, double bpm);

  /** Given 2 successive sample value of a sinusoidal function with known (normalized radian) 
  frequency w, such that:
  y0 = a * sin(p) 
  y1 = a * sin(p+w)
  this function computes the amplitude a and initial phase p which are passed as pointer 
  parameters to the function. The computed phase is in the range -PI...PI */
  RSLib_API void rsSineAmplitudeAndPhase(double y0, double y1, double w, double *a, double *p);

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
  too error-prone. */
  RSLib_API double rsSineFrequency(double y0, double y1, double y2, double smalll = 1.e-8);

  /** Assuming a sinusoidal input signal x of length N (and unknown frequency), this function 
  computes the sinusoid's instantaneous normalized radian frequency at sample-instant n0. The value 
  is more accurate, if the frequency varies less over 3 successive samples centered around n0 and 
  it's also more accurate, if n0 corresponds to a peak or valley in the signal. If it's not 
  possible to use 3 values centered around n0 (because n0 = 0 or n0 = N-1 and/or the sample value 
  of n0 is close to zero), the actual measurement point might be shifted (by at most 2 samples) 
  which should be inconsequential, if the sinuosoid has a stable frequency. The function is used 
  internally by rsSineFrequencyAt which is probably the function, you want to use instead. */
  RSLib_API double rsSineFrequencyAtCore(double *x, int N, int n0, double smalll = 1.e-8);
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
  RSLib_API double rsSineFrequencyAt(double *x, int N, int n0, bool refine = false);
   // todo: give the function another parameter that crosfades between arithmetic and geometric 
   // mean/extrapolation

  /** Assuming a sinusoidal input signal x of length N and normalized radian frequency w, this 
  function computes the sinusoid's instantaneous phase at sample-instant n0. */
  RSLib_API double rsSinePhaseAt(double *x, int N, int n0, double w);

  /** Assuming a sinusoidal input signal x of length N (and unknown frequency), this function 
  computes the sinusoid's instantaneous phase at sample-instant n0. */
  RSLib_API double rsSinePhaseAt(double *x, int N, int n0);


  /** Finds an index in the length N array x close to n0, where the signal x crosses zero. The 
  searchDirection parameter determines whether we search into the leftward (-1) or rightward (+1)
  direction and the upward/downward parameters determine, whether upward and/or downward crossings
  should count. if no zero crossing could be found, -1 will be returned. */
  RSLib_API int rsFindZeroNear(double *x, int N, int n0, int searchDirection, bool upward, 
    bool downward);

  
  /**
  
  NOTE: this function works only if the sample index n0 is between two zero-crossings of the signal
  x
  */
  RSLib_API double rsSinePhaseAtViaZeros(double *x, int N, int n0, int precision);



  /** Assuming a sinusoidal input signal x of length N with normalized radian frequency w, this
  function computes an offset in samples by which the sinusoid must be shifted, such that the 
  sinusoid's phase at sample instant n0 equals p0. The shift will be either forward or backward, 
  whichever is smaller. x: signal, N: signal length, n0: sample-instant for target phase p0, 
  p0: target phase at n0, return value: number of samples to shift (might be positive for 
  rightshift or negative for leftshift). */
  RSLib_API double rsSineShiftAmount(double *x, int N, int n0, double p0, double w);

  /** Similar to rsSineShiftAmount(double *x, int N, int n0, double p0, double w) but it works for 
  sinusoids with unknown frequency. It measures the sinusoids frequency first using 
  rsSineFrequencyAt and then uses the function that works for known frequencies. */
  RSLib_API double rsSineShiftAmount(double *x, int N, int n0, double p0);



  /** Converts a time-constant (typically denoted as "tau") of an exponential decay function to
  the time-instant at which the decay reaches the given level in decibels. For example
  tauToDecayTime(tau, -60.0) returns the time required to fall to -60 dB for the given value of
  tau. */
  RS_INLINE double rsTauToDecayTime(double tau, double decayLevel);

  /** Converts a time-stamp given in whole notes into seconds according to a tempo measured
  in beats per minute (bpm). */
  RS_INLINE double rsWholeNotesToSeconds(double noteValue, double bpm);

  //-----------------------------------------------------------------------------------------------
  // implementation:

  RS_INLINE double rsAmp2dB(double amp)
  {
    return 8.6858896380650365530225783783321 * log(amp);
  }

  RS_INLINE double rsAmp2dBWithCheck(double amp, double lowAmplitude)
  {
    if( amp >= lowAmplitude )
      return rsAmp2dB(amp);
    else
      return rsAmp2dB(lowAmplitude);
  }

  RS_INLINE double rsBeatsToSeconds(double beat, double bpm)
  {
    return (60.0 / bpm)*beat;
  }

  RS_INLINE double rsDB2amp(double dB)
  {
    return exp(dB * 0.11512925464970228420089957273422);
  }

  RS_INLINE double rsCubicFadeIn(double x)
  {
    //return x*(1+x*(1-x)); // slope = 1 at x = 0
    return x*(x*((PI/2-2)*x+(3-PI))+(PI/2));
  }

  RS_INLINE double rsCubicFadeOut(double x)
  {
    //return x*x*(x-2)+1;   // slope = -1 at x = 1
    return x*x*((2-PI/2)*x+(PI/2-3))+1;
  }

  RS_INLINE double rsFreqToPitch(double freq, double masterTuneA4)
  {
    return 12.0 * rsLog2(freq / masterTuneA4) + 69.0;
  }

  RS_INLINE double rsPitchOffsetToFreqFactor(double pitchOffset)
  {
    return exp(0.057762265046662109118102676788181 * pitchOffset);
    //return pow(2.0, pitchOffset/12.0); // naive, slower but numerically more precise
  }

  RS_INLINE double rsPitchToFreq(double pitch)
  {
    return 8.1757989156437073336828122976033 * exp(0.057762265046662109118102676788181 * pitch);
    //return 440.0*( pow(2.0, (pitch-69.0)/12.0) ); // naive, slower but numerically more precise
  }

  RS_INLINE double rsPitchToFreq(double pitch, double masterTuneA4)
  {
    return masterTuneA4 * 0.018581361171917516667460937040007
      * exp(0.057762265046662109118102676788181 * pitch);
  }

  RS_INLINE double rsSecondsToBeats(double timeInSeconds, double bpm)
  {
    return timeInSeconds * (bpm/60.0);
  }

  RS_INLINE double rsTauToDecayTime(double tau, double decayLevel)
  {
    static const double k = LN10 / 20.0;
    return -k*tau*decayLevel;
  }

  RS_INLINE double rsWholeNotesToSeconds(double noteValue, double bpm)
  {
    return (240.0/bpm) * noteValue;
  }

}

#endif
