#ifndef RAPT_BANDWIDTHCONVERTER_H
#define RAPT_BANDWIDTHCONVERTER_H

/** This is a class that provides a bunch of static functions to convert back and forth between
various parametrizations of the center-frequency and bandwidth of a bandpass (or peaking) filter.
There are a whole lot of possible ways to parametrize a bandpass filter such as via its lower and
upper bandedge frequencies, center-frequency and bandwidth (where the latter van be given either
in Hz or in octaves), center-frequency and Q value, and many more. Each of these parametrizations
may be more or less convenient in different contexts and it's sometimes necessary to convert
between the different parametrizations - that's what this class is for.

We define:
fl: lower bandedge frequency  fl = fc/k
fu: upper bandedge frequency  fu = fc*k
fc: center frequency          fc = sqrt(fl*fu)
ba: absolute bandwidth        ba = fu-fl
br: relative bandwidth        br = ba/fc
Q:  quality factor            Q  = 1/br
bo: bandwidth in octaves      bo = log2(fu/fl)
k:  bandedge factor           k  = 0.5*br+sqrt(0.25*br*br+1)


\todo maybe move this file into the filters section  */

class rsBandwidthConverter
{

public:

  /** Converts a lower and upper bandedge frequency (given by fl and fu respectively) to a
  corresponding center frequency. */
  template<class T>
  static T bandedgesToCenterFrequency(T fl, T fu);

  /** Converts a lower and upper bandedge frequency (given by fl and fu respectively) to a
  corresponding absolute bandwidth in Hz. */
  template<class T>
  static T bandedgesToAbsoluteBandwidth(T fl, T fu);

  /** Given a relative bandwidth br, this function computes the factor k by which the center
  frequency must be multiplied to get the upper bandedge frequency (or divided to get the lower
  bandedge frequency) such that fu = fc * k, fl = fc / k. */
  template<class T>
  static T relativeBandwidthToBandedgeFactor(T br);

  /** Given an absolute bandwidth in Hz and a center frequency fc, this function computes the
  lower and upper bandedge frequencies fl and fu. */
  template<class T>
  static void absoluteBandwidthToBandedges(T bw, T fc, T *fl, T *fu);

  /** Converts an absolute bandwidth bw given in Hz to the corresponding relative bandwidth,
  assuming a center frequency of fc. */
  template<class T>
  static T bandwidthAbsoluteToRelative(T bw, T fc);

  /** Converts an absolute bandwidth bw given in Hz to the corresponding quality factor Q,
  assuming a center frequency of fc. */
  template<class T>
  static T absoluteBandwidthToQ(T bw, T fc);

  /** Converts a lower and upper bandedge frequency (given by fl and fu respectively) to a
  corresponding bandwidth in octaves. */
  template<class T>
  static T bandedgesToBandwidthInOctaves(T fl, T fu);

  /** Converts an absolute bandwidth bw given in Hz to the corresponding value in octaves,
  assuming a center frequency of fc. */
  template<class T>
  static T absoluteBandwidthToOctaves(T bw, T fc);

  /** Assuming that we have an N-th order lowpass or (2*N)-th order bandpass Butterworth
  prototype filter with cutoff/bandedge frequency/ies defined at the half-power point (i.e.
  magnitude-squared is 0.5 at the cutoff/bandedges), this function computes a scaling factor for
  the cutoff (or bandwidth), such that the filter will have a magnitude of g at the cutoff point,
  when the filter is applied "numPasses" times (the same filter may be run over the signal
  several times). */
  template<class T>
  static T multipassScalerButterworth(int numPasses, int order = 1, T g = SQRT2_INV);
    // maybe implement a scaling function that, instead of matching at a specified gain point
    // match the area under the curve - requires to evaluate the integral of the Butterworth-power
    // magnitude (squared) response - wolfram alpha can solve this integral for specific values
    // for numPasses (M) and order (N), but not for general N,M - maybe we need to tell it to
    // assume N and M to be positive integers.....
    // the integral of 1/(1+x^N) from 0 to inf is pi / (N*sin(pi/N)) (try with sympy gamma for 
    // N=20 or N=21)
    // to find the inetgral, we must first find the partial frcation expansion of 
    // 1/(x^N+1)^M - i think, the roots are given by the N-th roots of -1, each with 
    // multiplicity M


  // \todo write a function: absoluteBandwidthToRingingTime


};

#endif
