#ifndef romos_FilterDesign_h
#define romos_FilterDesign_h

//namespace romos
//{

// todo: re-implement this using RAPT::rsLadderFilter, RAPT::rsBiquad - consolidate all formulas 
// there

// maybe use struct rsLadderCoeffs, struct rsBiquadCoeffs, to conveniently pass around pointers
// maybe have also struc biquadParams, ladderParams, etc - then just do
// RAPT::rsLadder::computeCoeffs(params, coeffs)

//-------------------------------------------------------------------------------------------------

// The following functions are for calculating biquad coefficients for filters realizing the 
// difference equation: y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]
// the coeffs array is assumed to be of the form: 
// coeffs[0] = b0, coeffs[1] = b1, coeffs[2] = b2, coeffs[3] = a1, coeffs[4] = a2

/** Assigns the biquad coefficients to produce a filter that passes the signal through 
unchanged. */
void biquadBypassCoeffs(double *coeffs);

/** First order lowpass via bilinear transform. */
void biquadLowpassCoeffsBilinear1(double *coeffs, double f);

/** First order highpass via bilinear transform. */
void biquadHighpassCoeffsBilinear1(double *coeffs, double f);

/** First order low shelving filter via bilinear transform. g is the linear amplitude gain. */
void biquadLowShelfCoeffsBilinear1(double *coeffs, double f, double g);

/** First order high shelving filter via bilinear transform. g is the linear amplitude gain. */
void biquadHighShelfCoeffsBilinear1(double *coeffs, double f, double g);

/** First order allpass filter. f is the -90 degree frequency in Hz. */
void biquadAllpassCoeffsBilinear1(double *coeffs, double f);

/** Second order lowpass via bilinear transform (RBJ cookbook design). */
void biquadLowpassCoeffsBilinear2(double *coeffs, double f, double q);

/** Second order highpass via bilinear transform (RBJ cookbook design). */
void biquadHighpassCoeffsBilinear2(double *coeffs, double f, double q);

/** Second order bandpass with constant skirt gain, peak gain = q (RBJ cookbook design). */
void biquadBandpassConstSkirtCoeffs(double *coeffs, double f, double q);

/** Second order bandpass with constant peak gain of unity (RBJ cookbook design). */
void biquadBandpassConstPeakCoeffs(double *coeffs, double f, double q);

/** Second order bandreject (RBJ cookbook design). */
void biquadBandrejectCoeffs(double *coeffs, double f, double q);

/** Second order peak/dip filter (RBJ cookbook design). g is the linear amplitude gain. */
void biquadPeakCoeffs(double *coeffs, double f, double q, double g);

/** Second order low shelving filter (RBJ cookbook design). g is the linear amplitude gain. */
void biquadLowShelfCoeffsBilinear2(double *coeffs, double f, double q, double g);

/** Second order high shelving filter (RBJ cookbook design). g is the linear amplitude gain. */
void biquadHighShelfCoeffsBilinear2(double *coeffs, double f, double q, double g);

/** Second order allpass filter (RBJ cookbook design). q determines slope of the phase response at 
the -180 degree frequency. */
void biquadAllpassCoeffsBilinear2(double *coeffs, double f, double q);

/** Computes ladder filter coefficients for a filter that implements the difference equations:
y0[n] = x[n]  - k * y4[n-1];
y1[n] = y0[n] + a1*(y0[n]-y1[n]);
y2[n] = y1[n] + a1*(y1[n]-y2[n]);
y3[n] = y2[n] + a1*(y2[n]-y3[n]);
y4[n] = y3[n] + a1*(y3[n]-y4[n]);
y[n]  = c0*y0[n] + c1*y1[n] + c2*y2[n] + c3*y3[n] + c4*y4[n];
where the "coeffs"-array is organized as follows: 
coeffs[0]=a1, [1]=k, [2]=c0, [3]=c1, [4]=c2, [5]=c3, [6]=c4
"f" is the cutoff frequency, mode is the filter-mode, "r" is the amount of resonance (1.0 means 
self-oscillation) and "ag" is the amount of automatic gain compensation that counteracts the gain 
loss at high resonances (1.0 means full compensation) */
void ladderCoeffs(double *coeffs, int mode, double f, double r, double ag);

//}

#endif
