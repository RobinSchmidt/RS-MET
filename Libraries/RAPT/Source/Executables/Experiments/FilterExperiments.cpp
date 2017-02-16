#include "FilterExperiments.h"

//template <class T1, class T2>
//void rsScale(T1 *buffer, int length, T2 scaleFactor)
//{
//  for(int n = 0; n < length; n++)
//    buffer[n] *= scaleFactor;
//}

//-------------------------------------------------------------------------------------------------

void ladderResonanceManipulation()
{
  // We create two ladder filters with the same cutoff frequency, one with and one without 
  // resonance and plot the difference of the states of the two ladder filters when we pass
  // a sawtooth wave into them.

  // parameters:
  int   N   =  2000;      // number of samples
  float fs  = 44100;      // sample rate
  float fc  =  500;       // cutoff frequency
  float res =    0.99f;   // resonance
  float fIn =   80;       // input frequency

  // create and set up the filters:
  typedef RAPT::LadderFilter<float, float> LDR;  // for convenience

  LDR withReso;
  withReso.setSampleRate(fs);
  withReso.setCutoff(fc);
  withReso.setResonance(res);
  withReso.setMode(LDR::LP_24);

  LDR noReso;
  noReso.setSampleRate(fs);
  noReso.setCutoff(fc);
  noReso.setResonance(0);
  noReso.setMode(LDR::LP_24);

  // create time axis and input signal;
  vector<float> t(N), x(N);
  createTimeAxis(N, &t[0], fs);
  createWaveform(&x[0], N, 1, fIn, fs, 0.f, false);

  // filter input signal and record the difference between the state variables of teh two filters:
  vector<float> y0(N), y1(N), y2(N), y3(N), y4(N);
  vector<float> z0(N), z1(N), z2(N), z3(N), z4(N);
  vector<float> r0(N), r1(N), r2(N), r3(N), r4(N);
  vector<float> offset(N);
  float y[5], z[5];
  float dummy;
  for(int n = 0; n < N; n++)
  {
    // let the filter compute outputs (we don't actuualy need them):
    dummy = withReso.getSample(x[n]);
    dummy = noReso.getSample(  x[n]);

    // obtain the filter states:
    withReso.getState(y);
    noReso.getState(z);

    // record the states and state-differences:
    y0[n] = y[0];
    y1[n] = y[1];
    y2[n] = y[2];
    y3[n] = y[3];
    y4[n] = y[4];

    z0[n] = z[0];
    z1[n] = z[1];
    z2[n] = z[2];
    z3[n] = z[3];
    z4[n] = z[4];

    r0[n] = y[0] - z[0];
    r1[n] = y[1] - z[1];
    r2[n] = y[2] - z[2];
    r3[n] = y[3] - z[3];
    r4[n] = y[4] - z[4];

    // scale difference-states with appropriate factor to make the sinusoids have the same 
    // amplitude:
    r1[n] *= (float)SQRT2;      // s2^1, s2: sqrt(2)
    r2[n] *= 2;                 // s2^2
    r3[n] *= 2*(float)SQRT2;    // s2^3
    r4[n] *= 4;                 // s2^4

    // According to the model, this value would be zero because r0 and r4 are 180° out of phase. 
    // Any deviation from 0 is a shortcoming of the model and considered an offset to the idealized
    // value:
    offset[n] = r0[n] + r4[n];
  }

  // plot:
  GNUPlotter plt;
  plt.addDataArrays(N, &t[0], &x[0]);
  //plt.addDataArrays(N, &t[0], &y0[0], &y1[0], &y2[0], &y3[0], &y4[0]);
  //plt.addDataArrays(N, &t[0], &z0[0], &z1[0], &z2[0], &z3[0], &z4[0]);
  //plt.addDataArrays(N, &t[0], &r0[0], &r1[0], &r2[0], &r3[0], &r4[0]);
  plt.addDataArrays(N, &t[0], &r0[0], &r4[0]); // should be 180° out of phase
  plt.addDataArrays(N, &t[0], &offset[0]);
  //plt.addDataArrays(N, &t[0], &y0[0], &z0[0], &r0[0]);
  plt.plot();

  // Observations:
  // The states of the imaginary "pure-resonance" difference filters (resonant - nonResonant)
  // represent (decaying) sinusoids, where each of the successive states has the same sinusoid
  // multiplied by 1/sqrt(2) and phase-shifted by -45 degrees with respect to the stage before. 
  // We assume that at any time instant the resonance waveform has instantaneous amplitude a and 
  // phase p. We must then have:
  //
  // r4 = a             * sin(p)
  // r3 = a * sqrt(2)^1 * sin(p + 1*pi/4)
  // r2 = a * sqrt(2)^2 * sin(p + 2*pi/4) = a * 2 * sin(p + pi/2)
  // r1 = a * sqrt(2)^3 * sin(p + 3*pi/4)
  // r0 = a * sqrt(2)^4 * sin(p + 4*pi/4) = a * 4 * sin(p + pi) 
  //
  // That means, given the states r4,r2 we should be able to figure out the instantaneous amplitude
  // a and phase p. We just divide the state r2 by two and together with r4 we can take it as the 
  // sine and cosine part of our sine from which amplitude and phase can be computed. Having a and
  // p, we can do whatever we want to them (for example sync the phase to an input signal) and then 
  // compute our new states r0,...,r4 using the formulas above. Having done that, we may update the 
  // states of the resonant filter according to y0 = z0 + r0, ..., y4 = z4 + r4.  This should give 
  // us a ladder filter in which we can freely mess with the instantaneous phase of the resonance 
  // without needing to post-process the resonance signal.
  // hmm...well, all of this holds only approximately. Maybe we should estimate the sine parameters
  // not from r2,r4 but from r3,r4 (so we use the 2 most filtered outputs) and maybe we should 
  // compute the difference between the actual states and the idealized states according to the 
  // model and add this difference back after modification of the states, such that when we don't
  // apply any modification, we really don't apply any modification (when recomuting the states
  // according to the formulas)

  // Note: it works only when the compensation gain is applied at the input and when the filter
  // is linear.
}


/** N: num samples, x: abscissa-values (mostly time), y: ordinate values, avg: average - the 
output, width: length/width/range of the support of the filter, weightFunc: normalized weighting
function - should have a support in the range -1..+1 */
template<class T>
void movingAverage(int N, T* x, T* y, T* avg, T width, T (*weightFunc)(T))
{
  T w2  = 0.5 * width;                         // half width
  T w2r = 1 / w2;                              // reciprocal of half width
  T dist;                                      // distance
  int k;                                       // inner loop index
  for(int n = 0; n < N; n++){                  // outer loop over all points
    T wgt = weightFunc(0);                     // weight
    T sw  = wgt;                               // sum of weights
    T swv = wgt * y[n];                        // sum of weighted values
    k = n-1;                                   // immediate left neighbour
    while(k >= 0 && dist = (x[n]-x[k]) <= w2){ // left side loop
      wgt  = weightFunc(dist * w2r);           // compute weight for distance
      sw  += wgt;                              // accumulate weight sum
      swv += wgt * y[k];                       // accumulate weighted values
      k--; }                                   // jump to next neighbour
    k = n+1;
    while(k < N && dist = (x[k]-x[n]) <= w2){  // right side loop
      wgt  = weightFunc(dist * w2r);
      sw  += wgt;
      swv += wgt * y[k];
      k++; }
    avg[n] = swv / sw; }
}

void nonUniformMovingAverage()
{



  int dummy = 0;
}
