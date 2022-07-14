// todo: merge with the other FilterExperiments.cpp file

using namespace RAPT;

// Move these into the prototypes and later maybe into the library:

/** Computes coefficients for 2nd order (biquad) allpass for the difference equation:
  
  y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] - a1*y[n-1] - a2*y[n-2]

The formula is from DAFX (1st Ed), pg 41 and implements a bilinear transfrom design (verify!).

ToDo: document what wb is. It determines the steepness of the phase response and therefore the
bandwidth of a derived bandpass or bandreject filter. But how exactly? Is this an absolute 
bandwidth in Hz? Or is it somehow in terms of Q? What would be the formula to convert from a 
desired bandwidth in octaves to wb? Compare to RBJ allpass formula. They use also the BLT, so 
maybe they are equivalent? Or maybe not? RBJ writes something about prewarping center-freq and
bandwidth but I guess, one could also choose to match lower and upper band-edges or maybe do
something else. */
template<class T> 
void coeffsAllpassDAFX(T wc, T wb, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  T d = -cos(wc);
  T t =  tan(T(0.5)*wb);
  T c =  (t-T(1)) / (t+T(1));
  *b0 = -c;
  *b1 = d * (T(1)-c);
  *b2 = T(1);
  *a1 = *b1;
  *a2 = *b0;

  // The a-coeffs are obtained by reversing the array of b-coeffs (generally for allpass filters)
}

/** Computes coeffs for a biquad bandpass derived from an allpass via: 

  y[n] = (x[n] + a[n]) / 2

where a[n] is the output of the allpass filter. See DAFX (1st Ed), pg 43. */
template<class T> 
void coeffsBandpassDAFX(T wc, T wb, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  coeffsAllpassDAFX(wc, wb, b0, b1, b2, a1, a2);
  *b0  = T(1) + *b0;
  *b0 *= T(0.5);
  *b1 *= T(0.5);
  *b2 *= T(0.5);
}

/** Computes coeffs for a biquad bandstop derived from an allpass via: 

  y[n] = (x[n] - a[n]) / 2

where a[n] is the output of the allpass filter. See DAFX (1st Ed), pg 43. */
template<class T> 
void coeffsBandstopDAFX(T wc, T wb, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  coeffsAllpassDAFX(wc, wb, b0, b1, b2, a1, a2);
  *b0  = T(1) - *b0;
  *b0 *= T(0.5);
  *b1 *= T(0.5);
  *b2 *= T(0.5);
}

/** RBJ cookbook allpass for analog prototype transfer function:
  H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)  */
template<class T> 
void coeffsAllpassRBJ(T w0, T Q, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  T s, c; rsSinCos(w0, &s, &c);   // s = sin(w0), c = cos(w0)
  T a = s / (T( 2) * Q);          // alpha
  T k = 1 / (T( 1) + a);          // 1/a0, scaler for all coeffs
  *b0 = k * (T( 1) - a);
  *b1 = k * (T(-2) * c);
  *b2 = T(1);
  *a1 = *b1;
  *a2 = *b0;
}

/** RBJ cookbook bandpass with constant peak gain for analog prototype transfer function:
  H(s) = (s/Q) / (s^2 + s/Q + 1) */
template<class T> 
void coeffsBandpassPeakRBJ(T w0, T Q, T* b0, T* b1, T* b2, T* a1, T* a2)
{
  T s, c; rsSinCos(w0, &s, &c);   // s = sin(w0), c = cos(w0)
  T a = s / (T( 2) * Q);          // alpha
  T k = 1 / (T( 1) + a);          // 1/a0, scaler for all coeffs
  *b0 =  k * a;
  *b1 =  T(0);
  *b2 = -k * a;
  *a1 =  k * T(-2) * c;
  *a2 =  k * (T(1) - a);
}

//coeffsBandpassPeakRBJ(wc, Q, &b0, &b1, &b2, &a1, &a2);

/** RBJ cookbook formula for converting from a bandwidth bw given in octaves to the quality factor
"Q" of a digital filter using the bilinear transform (BLT). */
template<class T>  
T bandwidthToDigitalQualityRBJ(T w0, T bw)
{
  return T(1) / (T(2)*sinh(log(2)/2 * bw * w0/sin(w0)));
}



void bandpassAndNotch()
{
  // Under construction...so far, perfect reconstruction seems to work but I did not yet check the
  // frequency responses of the bands and they are almost certain to be still a mess because I have
  // not yet figured out how to compute wb from bw (see code comments below).

  // We create a bandpass filter and a complementary notch filter by subtracting the bandpass 
  // signal from the original. Eventually, the goal is to use this as a building block to design a
  // bandpass filter bank with perfect reconstruction capabilities. In the DAFX book page 44, we
  // see that a biquad bandpass has a phase response of zero at the center frequency (it goes from
  // +90° at DC to -90° at fs/2), so subtraction of the bandpass signal from the original should
  // indeed yield a notch...although in the book itself, they create both bandpass and notch from
  // adding or subtracting an allpass signal from the original (and dividing by 2). The allpass
  // has a phase response from 0° at DC to -360° at fs/2, going through -180° at the "cutoff". I
  // find that a bit confusing - is it possible to create the notch in either of the two ways:
  // adding the allpass (and dividing by 2) or subtracting the bandpass? ...Yes! Let: 
  //
  //   x                      input 
  //   y_a                    allpass output
  //   y_b = (x - y_a) / 2    bandpass output
  //   y_n = (x + y_a) / 2    notch output as in DAFX book
  //       =  x/2 + y_a/2     ...some algebra
  //       =  x - (x-y_a)/2   ...some more
  //       =  x - y_b         notch output as I want to see it
  //
  // So the signal flow diagram may look like this:
  //
  //                  -                  -                  -    
  //  x -----> BFF1 ---+--------> BPF2 ---+--------> BPF3 ---+------   ...and so on
  //      |            ^     |            ^     |            ^
  //      |            |     |            |     |            |
  //      --------------     --------------     --------------
  //
  // However, this is only the way I like to think about it conceptually. The actual implementation
  // below uses an allpass for each stage and creates the BP and BR signals via the add/subtract 
  // formulas above. Also we have added an inital highpass stage there.
  //
  // Notes:
  // -The band ouptuts are the outputs of the bandpasses
  // -At the adders, we produce the notch outputs which are only intermediate signals except maybe
  //  for the last one which we may then split into lowpass and highpass parts and provide as 
  //  additional output bands
  // -The filters should perhaps be ordered by descending center frequencies such that the total 
  //  delay is largest in low freqs and smallest in the high freqs because delay matters more in 
  //  the highs.
  //
  // Questions:
  // -Maybe it does indeed make sense to implement the BPF/BRF pair as an allpass and then do the
  //  add/subtract? Maybe that could be an optimization in terms of memory usage because allpasses
  //  only need 2 coeffs rather than 5
  // -What if we use higher order bandpasses (or allpasses) instead of 2nd order? Will we get 
  //  wiggles in the responses of the BPFs and or BRFs? 
  // -What about sign inverting the notch signals? Could that helpwith the phase/delay response?
  //  But no - that will mess up the reconstruction...unless we take the inversion into account in
  //  the reconstruction which would be easy to do - just invert again. Experiment with this, Plot
  //  phase delay and group delay responses for various configurations of inversions. Maybe invert
  //  only every other notch. Do whatever gives the least amount of delay.

  using Real   = double;
  using Vec    = std::vector<Real>;
  using Mat    = rsMatrix<Real>;
  using Biquad = rsBiquadDF1<Real, Real>;
  using FOB    =  rsFirstOrderFilterBase<Real, Real>;

  // Setup:
  // -The 1st BP passes 4k..16k, the 2nd 1k..4k, the 3rd 250..1000, the 4th 62.5..250
  // -The lowpass passes freqs below 62.5 and the highpass freqs above 16k
  Real fs = 44100;                      // Sample rate
  Vec  fc = { 8000, 2000, 500, 125 };   // Bandpass center frequencies in Hz
  Real B = 1.5;
  Vec  bw = {    B,    B,   B,   B };   // Bandwidths in octaves
  //Real fl = 62.5;                       // Lowpass cutoff in Hz
  Real fh = 16000;                      // Highpass cutoff in Hz
  //int noiseLength   = 512;
  //int impulseLength = 512;
  int N = 512;                          // Signal length in samples
  int fftOversample = 8;                // Oversampling factor for the FFT for spectral plots


  // Test:
  //fc = { 4000, 250 }; B  = 3.0; bw = { B, B}; 
  fc ={ 1000 }; bw = { 0.2 };



  // Create and set up the filters:
  int numBPFs  = (int)fc.size();
  int numBands = numBPFs + 2;           // = num bandpass outputs plus low and high
  std::vector<Biquad> apfs(numBPFs);    // Allpass filters
  //Biquad lpf;                           // Lowpass
  Biquad hpf;                           // Highpass
  Real b0, b1, b2, a1, a2;              // Temporary biquad coeffs
  for(int i = 0; i < numBPFs; i++)
  {
    Real wc = 2*PI*fc[i]/fs;

    // Old:
    //Real wb = 2*wc;               // preliminary - should later be computed from bw[i]...somehow
    //coeffsAllpassDAFX(wc, wb, &b0, &b1, &b2, &a1, &a2);

    // New:
    Real Q = bandwidthToDigitalQualityRBJ(wc, bw[i]);
    coeffsAllpassRBJ(wc, Q, &b0, &b1, &b2, &a1, &a2);

    apfs[i].setCoefficients(b0, b1, b2, -a1, -a2);   // uses the other sign convention
  }
  FOB::coeffsHighpassBLT(2*PI*fh/fs, &b0, &b1, &a1);
  hpf.setCoefficients(b0, b1, 0.0, a1, 0.0);
  //FOB::coeffsLowpassBLT( 2*PI*fl/fs, &b0, &b1, &a1);
  //lpf.setCoefficients(b0, b1, 0.0, a1, 0.0);


  // Helper function to produce the M output signals from the given input signal:
  auto splitIntoBands = [&](const Vec& x)
  {
    int M = numBands;
    int N = (int) x.size();
    Mat Y(M, N);
    int i, n;

    // Maybe factor out into resetFilters:
    //lpf.reset();
    hpf.reset();
    for(i = 0; i < (int) apfs.size(); i++)
      apfs[i].reset();

    // Apply bandsplitting filters:
    for(n = 0; n < N; n++) 
    {
      Real tol = 1.e-14;
      Real err;

      // Apply highpass to obtain topmost band:
      Real res    = x[n];                      // Residual is initially equal to input
      Real sum    = 0.0;                       // Total sum of outputs is initially zero
      Real hpOut  = hpf.getSample(res);        // Apply highpass - we slice off highs first
      Y(M-1, n)   = hpOut;                     // Record highpass output
      sum        += hpOut;                     // Accumulate what we have sliced off already
      res         = x[n] - hpOut;              // Residual of highpass filter

      // Apply the bandpasses. They should fill the rows with indices M-2 down to 1. The first stage
      // stage gets the residual from the highpass as input and subsequent stages get the bandstop 
      // output of of the previous stage as input. The final residual is considered the lowpass
      // band:
      for(i = 0; i < numBPFs; i++) 
      {
        Real apOut   = apfs[i].getSample(res); // Allpass output of stage i
        Real bpOut   = 0.5 * (res + apOut);    // Bandpass output of stage i
        Real bsOut   = 0.5 * (res - apOut);    // Bandstop output of stage i
        Real test    = bpOut + bsOut;          // Should equal res. That means the BP and BR sum to
        err          = res - test;             // their input. 
        rsAssert(rsAbs(err) <= tol);           // ...that seems to work indeed.
        Y(M-2-i, n)  = bpOut;                  // Record bandpass output into final output
        sum         += bpOut;
        test         = res - bpOut;
        res          = bsOut;                  // Should be the same - yes, looks good
        err          = res - test;
        rsAssert(rsAbs(err) <= tol);
      }
      Y(0, n) = res;                           // The last residual is the low band
      sum += res; 

      // After slicing off frequency content band by band and accumulating it into our sum 
      // variable, at the very end, the sum should be equal to x[n] up to roundoff. That is the
      // prefect reconstruction condition which we check here:
      err = x[n] - sum;
      rsAssert(rsAbs(err) <= tol);
    } 

    return Y;
    // Actually, we don't really use the lpf here. Instead, the lowpass output is taken to be the 
    // final residual. Maybe get rid of the object...but maybe in some other implementation (one 
    // that slices off frequency content starting at low frequencies), we'll use it, so maybe don't
    // delete it yet.
  };


  // Create and set up plotter object fro plotting the frequency responses:
  SpectrumPlotter<Real> plt;
  plt.setFftSize(fftOversample*N);
  plt.setLogFreqAxis(true);


  // Create a couple of example input signals,split them into bands and plot results:
  Mat Y;
  Vec noise   = createNoise(  N, -1.0, +1.0);
  Vec impulse = createImpulse(N);
  //Y = splitIntoBands(noise);   plotMatrixRows(Y);

  Y = splitIntoBands(impulse); 
  //plotMatrixRows(Y);
  plt.plotDecibelSpectraOfRows(Y);

  // Observations:
  // -The purple impulse resonse looks strange - it has only positive values but is very erratic.
  //  Is this the highpass response? Or the highest bandpass? But it's the only one that's only 
  //  positive which could mean, it's actually the lowpass response? Figure out - it's strange!
  //  ...with the new RBJ formulas, it looks different and actually goes slightly below zero. I 
  //  think, it's indeed the highest bandpass response
  // -The magnitude responses look totally useless. Looks like I have created a splitter that has
  //  the nice perfect reconstruction property but kinda useless per band frequency responses.
  // -Most of the frequency content seems to end up in the final residual (= "lowpass") band.

  // Conclusion:
  // -So far, the idea seems not to lead to a useful bandpass-based frequency splitter. We should 
  //  probably stick to the LP/HP tree in case of IIR bandpasses and for a bandpass bank, we'll
  //  probably need linear-phase FIR filters.

  // ToDo:
  // -Maybe try a different algorithm for a 3-band splitter:
  //  -Extract the middle band first by a bandpass.
  //  -The residual is splitted further into low/high by a pair of complementary LPF/HPF whose
  //   split frequency equals the bandpass center frequency.
  // -Maybe to create more bands, the same idea can now be applied recursively to the so extracted
  //  low/mid/high bands...it could be especially interesting to split the low and high bands 
  //  further - either by another such 3-way splitter or by a two way splitter. We could combine 
  //  such 2-way and 3-way splitters in various ways to build multi-way splitters
  // -Write a helper function to reconstruct the full signal from the individual bands. It should 
  //  add the matrix rows into a single row, i.e. each column is summed to give a single signal 
  //  value. Maybe it should allow partial reconstructions
  // -In setting up the filters, we need to find a correct formula to compute the wb variable. 
  //  Currently, we use a preliminary (wrong) formula. this doesn't affect the perfect 
  //  reconstruction, though - it only messes up the frequency responses of our band filters.
  //  Check out the RBJ allpass design formulas for that - maybe also the bandpass (const peak 
  //  gain) formulas
  // -Verify that subtracting a bandpass signal from the original (i.e. previous resiudal) gives 
  //  the same results as we get with the allpass add/sub technique. Maybe using a bandpass may be
  //  even more efficient? But maybe using an allpass is better in terms of modulation? Try it!
  // -Plot the magnitude- and phase-responses of all the band outputs. I think, we should bandpass 
  //  responses with additional notches above the main lobe. The first BP has no notch, the 2nd has
  //  one 1, the 3rd has 2 and so on -> check, if that is actually the case
  // -Plot also responses of partial reconstructions, i.e. reconstructions that use only a subset
  //  of all bands.
  // -Plot the impulse responses, step-responses, delay, etc.
  // -Try using more filters with narrower bandwidths (maybe 1, 1/2, 1/3, 1/6, 1/12 oct) and 
  //  investigate how this changes the responses - especially with respect to delay.
  // -Try reversing the order of the filters, starting at low freqs and going up. I think, this 
  //  should delay the high frequencies more - which is undesirable which is why by default, the
  //  high bands should come first.
  // -Maybe try to somehow split off the highpass part first and investigate, what difference that
  //  makes

  int dummmy = 0;
}

void bandSplittingTwoWay()
{
  // user parameters:
  float sampleRate = 44100;
  float splitFreq  = sampleRate/4; // 11025

  // set up the splitter:
  float omega = 2*float(PI)*splitFreq/sampleRate;
  rsTwoBandSplitter<float, float> splitter;
  splitter.setOmega(omega);

  // plot poles/zeros and magnitude responses:
  FilterPlotter<float> plt;
  // ...
}

void bandSplittingThreeWay()
{
  // We first split a mid signal off using a 2nd order bandpass filter, then obtain the difference
  // between original and bandpass signal which should give a bandreject signal. Then we split this
  // bandreject signal further using a complementary lowpass/highpass pair with split frequency 
  // equal to the bandpass center frequency. What we expect to see is a clean bandpass magnitude
  // response and the resulting lowpass/highpass responses should feature an additional notch at
  // that frequency. Their "cutoffs" will be determined implicitly by the bandwidth given to the 
  // bandpass filter.

  using Real   = double;
  using Vec    = std::vector<Real>;
  using Mat    = rsMatrix<Real>;

  // Setup:
  Real fs = 44100;     // Sample rate
  Real fc = 1000;      // Bandpass center frequency in Hz
  Real bw = 2.5;       // Bandwidth in octaves
  int  N  = 1024;      // Signal length in samples
  int  K  = 4*N;       // Number of FFT bins for plot

  // Compute intermediate variables, create and set up filters:
  Real wc = 2*PI*fc/fs;
  Real Q = bandwidthToDigitalQualityRBJ(wc, bw);
  Real b0, b1, b2, a1, a2;                             // Temporary biquad coeffs
  coeffsBandpassPeakRBJ(wc, Q, &b0, &b1, &b2, &a1, &a2);
  rsBiquadDF1<Real, Real> bpf;
  bpf.setCoefficients(b0, b1, b2, -a1, -a2);
  rsTwoBandSplitter<Real, Real> splitter;
  splitter.setOmega(wc);

  // Create input signal, split into 3 bands and plot results:
  Mat Y(3, N);
  Vec x = createImpulse(N);
  for(int n = 0; n < N; n++)
  {
    // Apply the filters:
    Real yB = bpf.getSample(x[n]);            // Bandpass output
    Real yN = x[n] - yB;                      // Notch output
    Real yL = splitter.getLowpassSample(yN);  // Lowpass output
    Real yH = yN - yL;                        // Highpass output

    // Record lowpass, bandpass and highpass outputs:
    Y(0, n) = yL;
    Y(1, n) = yB;
    Y(2, n) = yH;

    // Check perfect reconstruction
    Real xR  = yL + yB + yH;                  // Reconstructed input
    Real err = x[n] - xR;
    rsAssert(rsAbs(err) < 1.e-14);
  }

  // Plot results:
  //plotMatrixRows(Y);  // Impulse responses
  SpectrumPlotter<Real> plt;
  plt.setFftSize(K);
  plt.setLogFreqAxis(true);
  plt.plotDecibelSpectraOfRows(Y);

  // Observations:
  // -We do indeed see a nice bandpass response and punctured lowpass and highpass responses
  // -Plot is not normalized - perhaps the plotting code is buggy. Or maybe it's the division by 
  //  the FFT size which is not appropriate here? The freq-axis also looks wrong
  // -When the bandwidth gets small (e.g. 0.1), the freq responses show ripples. I think, it's 
  //  because of the truncation of the impulse response after a finite length
  // -When bw = 0, the bandpass becomes an oscillator (a2 = 1) and the input gain (b0,b1,b2) 
  //  becomes zero. In an implementation for production use, it may make more sense to have the 
  //  zeros after the poles (i.e. use DF2 or TDF1) for this very reason - when the user dials in a
  //  zero bandwidth coming from some nonzero value, we don't want the oscillation to persist. 
  //  This is a very nice feature because the 3-way splitter then reduces to our old 2-way 
  //  splitter, using only low and high bands. The mid-band just becomes zero for such a setting.

  // ToDo:
  // -Maybe add a parameter midSeparation going from 0 to 1. It should scale yB in:
  //    yB = midSeparation * bpf.getSample(x[n])
  //  can be used to fade between 2-way and 3-way modes.
  //  -what if we dial it beyond the 0..1 range?
  // -Can we also have some other seperation parameter that somehow puts the whole signal into the
  //  mid channel? Maybe it would be nice if that would happen for our existing midSeparation when
  //  it's at -1?...but no - that makes no sense. Maybe a parameter that when at 0 gives full 
  //  seperation, at -1 reduces to a 2-way, i.e. puts nothing in the mid-band and at +1 puts 
  //  everything into the mid band?
  // -Try to avoid the bilinear frequency cramping by using prescribed Nyquist gain designs for
  //  bandpass and lowpass.
  // -Try to figure out hwo to do splitters with steeper slopes like we already did for the 2-way
  //  splitter

  int dummmy = 0;
}

void bandSplittingMultiWay()
{
  // Uses rsMultiBandSplitter to split a test signal (a unit impulse or noise) into various bands.

  // user parameters:
  int numSamples = 128;
  float sampleRate = 44100;
  //vector<float> splitFreqs = { 100, 300, 1000, 3000, 10000 };
  //vector<float> splitFreqs = { 80, 150, 250, 500, 1000, 2000, 4000, 8000 };
  vector<float> splitFreqs = { 125, 250, 500, 1000, 2000, 4000, 8000 };  // 8 bands
  //vector<float> splitFreqs = { 60, 90, 135, 200, 300, 450, 675, 1000, 1500, 2200, 3000, 4500, 6000, 9000, 13000 };  // 16 bands

  // set up dsp objects and signal arrays:
  int n, k;
  rsMultiBandSplitter<float, float> splitter;
  splitter.setSplitFrequencies(splitFreqs);
  splitter.setSplitMode(rsMultiBandSplitter<float, float>::ACCUMULATE_INTO_LOWPASS);
  //splitter.setSplitMode(rsMultiBandSplitter<float, float>::ACCUMULATE_INTO_HIGHPASS);
  //splitter.setSplitMode(rsMultiBandSplitter<float, float>::BINARY_TREE);
  int numBands = splitter.getNumActiveBands();
  std::vector<float> x(numSamples);
  std::vector<std::vector<float>> y(numBands);
  for(k = 0; k < numBands; k++)
    y[k].resize(numSamples);
  RAPT::rsNoiseGenerator<float> ng;
  std::vector<float> tmp(numBands);  // holds the array of bandpass output samples

  // produce input and output signals:
  x[0] = 1;
  for(n = 0; n < numSamples; n++)
  {
    //x[n] = ng.getSample();  // uncomment to use noise as test signal
    splitter.processSampleFrame(x[n], &tmp[0]);
    for(k = 0; k < numBands; k++)
      y[k][n] = tmp[k];

    // test, if the band outputs recombine correctly:
    float xr = 0.f;
    for(k = 0; k < numBands; k++)
      xr += y[k][n];
    float error = x[n] - xr;
    rsAssert(fabs(error) < 1.e-7);
  }

  // write outputs to files:
  rosic::writeToMonoWaveFile("Input.wav", &x[0], numSamples, 44100, 16);
  for(k = 0; k < numBands; k++)
    writeToMonoWaveFile("Output"+to_string(k)+".wav", &y[k][0], numSamples, (int) sampleRate, 16);

  //FilterPlotter<float> plt;

  // Plot impulse responses:
  GNUPlotter plt;
  for(k = 0; k < numBands; k++)
    plt.addDataArrays(numSamples, &y[k][0]);
  plt.plot();

  // Plot magnitude responses:
  SpectrumPlotter<float> splt;
  splt.setFftSize(8*numSamples);
  splt.setLogFreqAxis(true);
  //splt.plotSpectra(y);  // doesn't compile because y has the wrong format
  // ToDo: ass a function to SpectrumPlotter that can deal with a vector of vectors or use rsMatrix
  // for the band outputs

}

//void append(std::string& str, std::string& appendix


int updateBandSplitStrings(std::vector<std::string>& str, int start, int length)
{
  static int numCalls = 0;
  numCalls++;

  if(length == 1)
    return numCalls;

  int L = length/2;
  str[start+L] = str[start] + "H";
  str[start]   = str[start] + "L";  

  // recursion:
  updateBandSplitStrings(str, start+L, L); // upper branch (highpass)
  updateBandSplitStrings(str, start,   L); // lower branch (lowpass)

  return numCalls;
}
void bandSplittingTreeAlgo()
{
  // We test the strategy to build a binary tree of splits symbolically by using strings that
  // represent the path of partial signal along the splitter tree, for example, "LHL" means
  // the signal went through a lowpass then a highpass then another lowpass

  int numBands = 8; // algo currently assumes a power of 2
  std::vector<std::string> str1(numBands);
  int numCalls = updateBandSplitStrings(str1, 0, numBands); 
    // numCalls should be 2*numBands-1 (for numBands a power of two)

  // ok, the recursive function seems to work - now let's try to convert it into an iterative algo
  std::vector<std::string> str2(numBands);
  int inc = numBands; // nextPowerOfTwo(numBands) for non-powers-of-2?
  int numIts = 0;
  while(inc > 1)
  {
    int pos = 0;
    while(pos < numBands)
    {
      str2[pos+inc/2] = str2[pos] + "H"; // maybe call inc/2 offset
      str2[pos]       = str2[pos] + "L";
      pos += inc;
      numIts++;
    }
    inc /= 2;
  }
  // ok - this seems to work, too. In the actual splitter, we need to initialize y[0] with the 
  // input and then run exactly that type of iteration using y[pos] as input to the two-way 
  // splitters, y[pos] (also) as lowpass output and y[pos+inc/2] as highpass output

  bool ok = (str2 == str1);
  // return ok;  // todo: make function return a bool
}

void bandSplitFreqResponses()
{
  // Computes and plots frequency response curves of the multiband bandsplitter outputs. The 
  // response is calculated in two ways: (1) FFT of the impulse-response and (2) using the built-in
  // getMagnitude... functions. The prupose is mainly to compare the two to verify that 
  // implementation of the latter is correct.

  // maybe turn into unit test

  typedef std::vector<float> Vec;

  // user parameters:
  int N = 2048;     // should be power of two for easy FFT
  float sampleRate = 44100;
  //Vec splitFreqs = { };
  //Vec splitFreqs = { 10000 };
  //Vec splitFreqs = { 7000, 14000 };
  Vec splitFreqs = { 5000, 10000, 15000 };
  //Vec splitFreqs = { 4000, 8000, 12000, 16000 };
  //Vec splitFreqs = { 3000, 6000, 9000, 12000, 15000, 18000 };

  // create and set up the splitter:
  rsMultiBandSplitter<float, float> splitter;
  splitter.setSampleRate(sampleRate);
  splitter.setSplitFrequencies(splitFreqs);
  splitter.setSplitMode(rsMultiBandSplitter<float, float>::ACCUMULATE_INTO_LOWPASS);
  //splitter.setSplitMode(rsMultiBandSplitter<float, float>::ACCUMULATE_INTO_HIGHPASS);
  //splitter.setSplitMode(rsMultiBandSplitter<float, float>::BINARY_TREE);

  // create and set up the fourier transformer:
  typedef rsFourierTransformerRadix2<float> FT;
  rsFourierTransformerRadix2<float> ft;
  ft.setBlockSize(N);
  //ft.setNormalizationMode(FT::NORMALIZE_ON_FORWARD_TRAFO);
  ft.setNormalizationMode(FT::NORMALIZE_ON_INVERSE_TRAFO);

  //
  GNUPlotter plt;
  int numBands = splitter.getNumActiveBands();
  Vec tmp(numBands);  // array of bandpass output samples
  Vec fftFreqs(N);    // fft bin frequencies
  Vec impResp(N);     // impulse response of current band
  Vec magFFT(N);      // magnitude response computed by FFT
  Vec magTF(N);       // magnitude response computed by transfer function
  int i, k, n;        // indices for band, fft-bin and sample
  for(k = 0; k < N; k++)
    fftFreqs[k] = ft.binIndexToFrequency(k, N, sampleRate);
  for(i = 0; i < numBands; i++)
  {
    splitter.reset();

    // compute impulse response of i-th band:
    splitter.processSampleFrame(1, &tmp[0]);
    impResp[0] = tmp[i];
    for(n = 1; n < N; n++) {
      splitter.processSampleFrame(0, &tmp[0]);
      impResp[n] = tmp[i];
    }

    // get FFT magnitudes from impulse response:
    ft.getRealSignalMagnitudes(&impResp[0], &magFFT[0]);

    // get band magnitude reponse from transfer function:
    for(k = 0; k < N; k++)
      magTF[k] = splitter.getBandMagnitudeAt(i, fftFreqs[k]);

    // now, the contents of magFFT and magTF should be the same up to errors due to 
    // truncation of impulse-reponse before FFT and roundoff - plot them:
    //if(i == 3) // just to pick out one for test
    {
      plt.initialize();
      plt.plotFunctionTables(N, &fftFreqs[0], &magFFT[0], &magTF[0]);
    }
    // -The last 2 are correct, the others only similar
    // -It plots the freq-range up to fs rather than fs/2. Is that intentional? I don't think so.
  }

  // ToDo:
  // -Instead of making a tree, maybe do a salami-sclicing strategy, shaving off more and more
  //  high (or low) frequency content from the signal
}

void complementaryFiltersIIR()
{
  // Experiment to figure out pole/zero placements in the s-domain to obtain a high/low IIR 
  // splitter with perfect reconstruction...

  // Maybe useful:
  analyzeComplementaryFilter( complementaryLowpass1p1z()   );
  analyzeComplementaryFilter( complementaryLowpass2p2z()   );

  // Probably useless:
  //analyzeComplementaryFilter( complementaryLowpass2p3z()   );  // highly resonant
  //analyzeComplementaryFilter( complementaryLowpass4p4z1t() );  // unstable
  //analyzeComplementaryFilter( complementaryLowpass4p4z()   );  // weird
  //analyzeComplementaryFilter( complementaryLowpass4p5z()   );  // unstable

  // -Freq-axis scaling seems wrong - see comments in analyzeComplementaryFilter 
}


template<class T>
rsFilterSpecificationZPK<T> sos2zpk(
  const T* b0, const T* b1, const T* b2, const T* a1, const T* a2, int numBiquads)
{
  rsFilterSpecificationZPK<T> zpk;
  zpk.z.reserve(numBiquads);
  zpk.p.reserve(numBiquads);

  using Poly = rsPolynomial<T>;
  using Cmp  = std::complex<T>;
  Cmp r1, r2;    // for the 2 complex roots
  zpk.k = T(1);

  for(int i = 0; i < numBiquads; i++)
  {
    // add zero(s) of i-th stage:
    if(b2[i] == T(0)) {
      if(b1[i] != T(0))
        zpk.z.push_back(-b1[i]/b0[i]); }
    else {
      Poly::rootsQuadraticComplex(Cmp(b2[i]), Cmp(b1[i]), Cmp(b0[i]), &r1, &r2);
      zpk.z.push_back(r1);
      zpk.z.push_back(r2);  }

    // add pole(s) of i-th stage:
    if(a2[i] == T(0)) {
      if(a1[i] != T(0))
        zpk.p.push_back(-a1[i]); }
    else {
      Poly::rootsQuadraticComplex(Cmp(a2[i]), Cmp(a1[i]), Cmp(T(1)), &r1, &r2);
      zpk.p.push_back(r1);
      zpk.p.push_back(r2);  }

    // accumulate gain of i-th stage:
    zpk.k *= b0[i]; // is that correct?
  }

  return zpk;
}
// maybe move to rsFilterSpecificationZPK as static member: fromBiquads

void engineersFilterRingResp1()
{
  // We plot the magnitude response of a filter together with a new kind of response that attempts
  // to show the amount of ringing of the filter as function of frequency. Ringing is in this 
  // context seen as a sine-like oscillation with a (more or less) stable frequency in the filter's
  // impulse- and step-response. It is desired that the function shows a strong peak at the 
  // resonance frequency of a (2-pole) resonator, at the ringing frequency of an allpass, at the 
  // cutoff frequency of a (steep) lowpass, etc. - it shoould peak at whatever frequency we expect 
  // strong ringing for the given filter at hand. Clearly, the fact that allpass filters can 
  // feature ringing shows that looking at the magnitude response alone will be not enough. 
  // Likewise, the fact that a steep linear- or zero-phase lowpass filter features ringing shows 
  // that looking at the phase response alone won't be enough either. We instead look at a 
  // combination of the two. In such cases, it's often useful to take a step back and look and 
  // real and imaginary parts instead because they are on a more equal footing. Observing that 
  // "ringing" occurs whenever "something happens" in the frequency response, we look at the 
  // absolute value of the derivative of the frequency response (as a complex function of a real 
  // variable).

  // See this thread: https://www.kvraudio.com/forum/viewtopic.php?f=33&t=569114 

  // filter parameters:
  double fs  = 44100;  // samplerate in Hz
  double fc  = 100;    // center or cutoff frequency in Hz
  double bw  = 2.0;    // bandwidth in octaves
  int    ord = 5;      // prototype order

  using EF   = rsEngineersFilter<double, double>;
  using PTD  = rsPrototypeDesigner<double>;
  using IIRD = rsInfiniteImpulseResponseDesigner<double>;

  EF flt;
  //flt.setApproximationMethod(PTD::BUTTERWORTH);
  flt.setApproximationMethod(PTD::ELLIPTIC);
  //flt.setApproximationMethod(PTD::CHEBYCHEV);
  //flt.setApproximationMethod(PTD::INVERSE_CHEBYCHEV);
  //flt.setApproximationMethod(PTD::BESSEL);
  //flt.setApproximationMethod(PTD::PAPOULIS);
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setBandwidth(bw);
  //flt.setMode(IIRD::LOWPASS);
  //flt.setMode(IIRD::HIGHPASS);
  flt.setMode(IIRD::BANDPASS);
  //flt.setMode(IIRD::BANDREJECT);
  //flt.setMode(IIRD::PEAK);
  //flt.setMode(IIRD::LOW_SHELV);
  //flt.setMode(IIRD::HIGH_SHELV);
  flt.setGain(10.0);
  flt.setRipple(1.0);
  flt.setStopbandRejection(20.0); 
  flt.setPrototypeOrder(ord);

  // Plot different filter responses:
  //plotImpulseResponse(      flt, 10000, 1.0);
  //plotFrequencyResponse(    flt,  5000, 10.0, 1000.0, fs, true);
  //plotFrequencyResponseReIm(flt, 5000, 10.0, 1000.0, fs, true);
  plotMagAndRingResponse(   flt,  5000, 10.0, 1000.0, fs, true);  // experimental


  /*
  // temporary throw-away code, for generating a pic for the kvr forum:
  flt.setApproximationMethod(PTD::INVERSE_CHEBYCHEV);
  flt.setSampleRate(44100);
  flt.setFrequency(1000);
  flt.setMode(IIRD::LOWPASS);
  flt.setStopbandRejection(80.0); 
  flt.setPrototypeOrder(10);
  plotImpulseResponse(   flt, 1200, 1.0, true, -200.0);
  plotMagAndRingResponse(flt, 5000, 100.0, 10000.0, 44100.0, true); 
  // This suggests a way to estimate the ringing time numerically: just fit a straight line to
  // the peaks of this function. This will give a decay time constant in dB/sample from which the 
  // ringing time can be computed
  // ...maybe move this code into another function
  */



  // Observations:
  // -Notation: H(f): complex frequency response, H'(f): derivative of H(f) with respect to f, 
  //  R_i(f): i-th attempt of defining a "ringing response"
  // -The function R_1(f) := abs(H'(f)) does indeed show peaks at the bandedge frequencies.
  // -The function R_2(f) := abs(H'(f)) * f seems to make the plot symmetrical on a log-freq axis
  // -I think, R_1 is proportional to the absolute ringing time (in seconds) and R_2 is 
  //  proportional to the relative ringing time (in cycles) -> verify!
  // -An elliptic bandpass or bandreject filter does indeed show strong spikes in the so defined 
  //  "ringing response" around the passband edge frequencies

  // -A narrow 1st order bandpass with bw = 0.2 has a broader ringing response peak at the center
  //  frequency than the corresponding notch with otherwise equal settings, the height of the 
  //  peaks is the same though. With a wide bandpass (bw >= 2), the single ringing maximum splits 
  //  into two. The effect is very pronounced at bw = 3. This splitting happens at around bw = 5.2
  //  for the corresponding notch.

  // ToDo:
  // -Try to find a function that, in the case of a resonator, represents the decay time of a given 
  //  freq. Maybe measure it by driving the filter at the frequency into steady state, turn the 
  //  input signal off (smoothly, to avoid a strong switch-off transient!) and measure the tail 
  //  length. That technique may apply to other filter types, too.
  // -try leaving out the sqrt
  // -what happens, when we chain a bunch of (equal) filters? I guess, using N filters of same type
  //  in series just raises the response to the power of N?, just like with magnitude? Or will it 
  //  multiply by N just like with phase?
  // -test with allpasses, pure delay, FIR filters (windowed sinc, boxcar)
  // -test a narrow bandpass vs a notch - with the same Q, their ring-responses should be similar
  // -write the findings up into a paper, maybe "Filter Ringing as Function of Frequency"

  // Questions:
  // -Is H'(w) just the imaginary part of the complex derivative of H(s) at s = j*w? I think so. If
  //  this is the case, it will be easier to evaluate the desired derivative analytically for a 
  //  given rational H(s)...but wait - no: we actually want a complex derivative into the imaginary
  //  direction, the above would give the imaginary direction of a complex derivative. Can we 
  //  define such a derivative more generally: given a function w = f(z), consider the 1D 
  //  derivative at a given fixed real part (here zero) into the direction of the imaginary axis. 
  //  It's like a vector valued partial derivative...maybe the 1st row or column of the Jacobi 
  //  matrix when we consider f(z) as a 2D vector field: 
  //    f(z) = u(x,y) + j*v(x,y)  ->  f(x,y) = (u(x,y), v(x,y))
  //  could other combinations of partial derivatives reveal further information? (magnitude of)
  //  derivative along the real axis? divergence? curl? maybe the magnitude of the gradient can 
  //  reveal pole (or zero) proximity? ...but wait: it's a vector field, the gradient is a matrix.
  //  maybe look at the gradient of the magnitude (squared). We need an analytical framework that
  //  can deal with bivariate rational functions H(x,y) = N(x,y) / D(x,y). Create a class
  //  rsBivariateRationalFunction using two rsBivariatePolynomial members. It should be able to 
  //  compute partial derivatives, gradients, etc. ...maybe also path integrals
  // -hmm..actually if H(z) is analytic at a given z0, it does not matter from which direction we
  //  approach z0, so the derivative of H(z) with respect to w (approaching z0 along the imaginary
  //  axis) does not need to be distinguished from the derivative approaching from any other 
  //  direction and all the vector-analytic combinations of derivatives may simplify.
  // -What about lookign at the 2nd derivative of the phase response? Could this be some sort of 
  //  measure of dispersion? ...but i think, that's what the group delay is for
  // -Can evaluation of the transfer function along the real axis reveal something interesting? 
  //  Maybe how strong the filter responds to exponentially decaying functions as function of the 
  //  decay? can this be related to some concept of transient response, maybe? ...like: responds 
  //  strongly to quick transients and weakly to slow transients or vice versa?
  // -What about trying to find a scalar potential for the vector field defined by the transfer 
  //  function? Maybe we need to take the complex conjugate because of the conjugation involved in 
  //  forming "Polya vcetor fields"?

  // Ideas:
  // -what other kinds of responses can we derive from the transfer-function, frequency-response or
  //  the time responses?
  //  -Maybe define a drop-response as function of two frequencies w0,w1 as 
  //   abs(H(j*w1)) - abs(H(j*w2)) which measures, how much the magnitude response drops off between 
  //   the two given frequencies
  //  -Can we find a potential function of the transfer-function and interpret it somehow? maybe 
  //   also by considering differences of values at various points (along the w-axis). These 
  //   differences represent line integrals of some kind
}

void engineersFilterRingResp2()
{
  // We plot the ringing response of a series connection of 2 high-shelving filters in order to 
  // see, if it depends on the actual magnitude...tbc...

  // User parameters:
  double fs    = 44100;  // samplerate in Hz

  // Filter 1:
  double frq1  = 50.0;   // frequency
  double gain1 = 2.0;    // gain as factor
  int    ord1  = 15;     // prototype order

  // Filter 2:
  double frq2  = 200.0;
  double gain2 = 2.0;
  int    ord2  = 15;

  using EF   = rsEngineersFilter<double, double>;
  using PTD  = rsPrototypeDesigner<double>;
  using IIRD = rsInfiniteImpulseResponseDesigner<double>;

  // Create and setup filter 1:
  EF flt1;
  flt1.setApproximationMethod(PTD::BUTTERWORTH);
  flt1.setMode(IIRD::HIGH_SHELV);
  flt1.setSampleRate(fs);
  flt1.setFrequency(frq1);
  flt1.setGain(rsAmp2dB(gain1));
  flt1.setPrototypeOrder(ord1);

  // Create and setup filter 2:
  EF flt2;
  flt2.setApproximationMethod(PTD::BUTTERWORTH);
  flt2.setMode(IIRD::HIGH_SHELV);
  flt2.setSampleRate(fs);
  flt2.setFrequency(frq2);
  flt2.setGain(rsAmp2dB(gain2));
  flt2.setPrototypeOrder(ord2);

  // Create and setup combined filter:
  int numBiquads = (ord1 + ord2 + 2) / 2;  // todo: take into account order doubling of BP,PK,etc.
  rsBiquadCascade<double, double> flt(numBiquads);
  flt.chainWith(flt1);
  flt.chainWith(flt2);

  // Plot ringing responses:
  //plotMagAndRingResponse(flt1, 5000, 10.0, 1000.0, fs, true);
  //plotMagAndRingResponse(flt2, 5000, 10.0, 1000.0, fs, true);
  plotMagAndRingResponse(flt, 5000, 10.0, 1000.0, fs, true);
  plotImpulseResponse(   flt, 5000, 1.0);

  // Plot responses of an allpass version of the filter:
  flt.turnIntoAllpass();
  plotMagAndRingResponse(flt, 5000, 10.0, 1000.0, fs, true);
  plotImpulseResponse(   flt, 5000, 1.0);


  // Observations:
  // -With gain1 = 2.0, gain2 = 2.0:
  //  -When we don't divide by the magnitude, the 2nd peak is indeed roughly twice as high as the
  //   first. But that's really only very rough - the 1st peak is actually quite a bit higher than 
  //   half of the 2nd.
  //  -When we do divide by the magnitude, both peaks seem to be the same height. However, dividing 
  //   by the magnitude seems to be a bad idea in general due to ptotential division-by-zero issues
  //   which are apparent in elliptic filters (and probably notches, too). Dividing by 1+mag 
  //   instead of just mag itself fixes the division by zero, but then the heights are still 
  //   unequal here
  // -With gain1 = 2.0, gain2 = 0.5:
  //  -When we don't divide by the magnitude, both peak have equal height of around 2.2 and 
  //   equal shape.
  // -Turning the filter into an allpass also yields a symmetric response. Whether we divide by the
  //  magnitude or not does not matter in this case (obviously, since the magnitude is 1 
  //  everywhere). But what would that do to FIR filters?
  //  -the ringing amount of the allpass is around 1 (in 0.99..1.015)
  //  -the impulse response looks initially more like chirping than ringing (i.e. a frequency 
  //   sweepdown) and the later section looks rather irregular
  //
  // Conclusions:
  // -Due to the potential issues with division by zero, it does not seem to be a good idea to 
  //  divide by the magnitude response - at least not, when the filters features zeros. For allpole
  //  filters, it could make sense, but we want a measure that is universally applicable.
  //
  // Questions:
  // -What about the 2nd derivative? Maybe instead of dividing by the magnitude, we could divide by 
  //  a sum of magnitude and 2nd derivative? Maybe it happens that the 2nd derivative is nonzero 
  //  whenever the magnitude is zero and that could avoid the division-by-zero problems? 
  //  -> plot 2nd derivative
  // -What about the integral? could that reveal something useful, too? -> plot it!
  //  -> maybe plot (magnitudes of) derivatives and integrals of various orders and try to figure
  //  out if these plots can reveal anything useful
}

void engineersFilterRingResp()
{
  engineersFilterRingResp1();
  engineersFilterRingResp2();
}

void engineersFilterFreqResps()
{
  // We plot frequency responses of engineer's filter for various orders. The main aim is actually
  // to demonstrate how the plotting works..

  int minOrder  = 1;      // minimum filter order
  int numOrders = 5;      // number of orders
  int orderInc  = 1;      // spacing between the orders
  int numFreqs  = 501;    // number of frequency samples
  double fc     = 1000;   // cutoff frequency
  double fs     = 44100;  // sample rate

  using EF   = rsEngineersFilter<double, double>;
  using PTD  = rsPrototypeDesigner<double>;
  using IIRD = rsInfiniteImpulseResponseDesigner<double>;
  using ZPK  = rsFilterSpecificationZPK<double>;

  EF flt;
  flt.setApproximationMethod(PTD::BUTTERWORTH);
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setMode(IIRD::LOWPASS);

  FilterPlotter<double> plt;
  double *a1, *a2, *b0, *b1, *b2;
  a1 = flt.getAddressA1();
  a2 = flt.getAddressA2();
  b0 = flt.getAddressB0();
  b1 = flt.getAddressB1();
  b2 = flt.getAddressB2();
  for(int i = 0; i < numOrders; i++)
  {
    int order = minOrder + i*orderInc;
    flt.setPrototypeOrder(order);
    ZPK zpk = sos2zpk(b0, b1, b2, a1, a2, flt.getNumStages());
    zpk.sampleRate = fs;
    plt.addFilterSpecificationZPK(zpk);
    int dummy = 0;
  }

  plt.plotFrequencyResponses(numFreqs, 20.0, fs/2, true);
  //plt.plotFrequencyResponses(numFreqs, 20.0, fs/2, true, true, true, false);


  // Observations:
  // -the colors need to be adjusted

  int dummy = 0;
}

void firstOrderFilters()
{
  typedef rsOnePoleFilter<double, double> FLT;

  double sr = 44100;
  //double fc = sr/4;  // halfband
  double fc = sr/8;    // quarterband

  FLT flt;
  flt.setSampleRate(sr);
  flt.setCutoff(fc);
  //flt.setShelvingGain(1.5);
  flt.setShelvingGain(0.75);
  //flt.setMode(FLT::LOWPASS_IIT);
  //flt.setMode(FLT::HIGHPASS_MZT);
  //flt.setMode(FLT::ALLPASS_BLT);
  //flt.setMode(FLT::LOWPASS_BLT);
  //flt.setMode(FLT::HIGHPASS_BLT);
  flt.setMode(FLT::LOWSHELV_BLT);
  flt.setMode(FLT::HIGHSHELV_BLT);

  // translate object variables into filter spec as needed by the plotter:
  RAPT::rsFilterSpecificationBA<double> spec;
  spec.sampleRate = sr;
  spec.a.resize(2);
  spec.a[0] = 1;
  spec.a[1] = -flt.a1;  // the filter class uses the other sign convention
  spec.b.resize(2);
  spec.b[0] = flt.b0;
  spec.b[1] = flt.b1;

  // make some plots:
  spec.sampleRate = 1;  // without it, the plots are messed up -> debug this
  showFilterPlots(spec);
}

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
  typedef RAPT::rsLadderFilter<float, float> LDR;  // for convenience

  LDR withReso;
  withReso.setSampleRate(fs);
  withReso.setCutoff(fc);
  withReso.setResonance(res);
  withReso.setMode(LDR::Mode::LP_24);

  LDR noReso;
  noReso.setSampleRate(fs);
  noReso.setCutoff(fc);
  noReso.setResonance(0);
  noReso.setMode(LDR::Mode::LP_24);

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


/** Moving average filter for non-uniformly spaced samples. N: num samples, x: abscissa-values
(mostly time), y: ordinate values, avg: average - the output (must be distinct from y), width: 
length/width/range of the support of the filter, weightFunc: normalized weighting function, 
should have a support in the range -1..+1, this is the filter kernel as continuous function. 
\todo: maybe let the user pass a functor for the weighting function instead of a function 
pointer (maybe it's possible to write a functor-wrapper for ordinary functions with automatic
type casting such that ordinary functions can still be passed as usual?) */
template<class Tx, class Ty>
void movingAverage(int N, Tx* x, Ty* y, Ty* avg, Tx width, Ty (*weightFunc)(Tx))
{
  Tx w2  = 0.5f * width;                       // half width
  Tx w2r = 1 / w2;                             // reciprocal of half width
  Tx dist;                                     // distance
  int k;                                       // inner loop index
  for(int n = 0; n < N; n++){                  // outer loop over all points
    Ty wgt = weightFunc(0);                    // weight
    Ty sw  = wgt;                              // sum of weights
    Ty swv = wgt * y[n];                       // sum of weighted values
    k = n-1;                                   // immediate left neighbour
    while(k >= 0 && (dist = x[n]-x[k]) <= w2){ // left side loop
      wgt  = weightFunc(dist * w2r);           // compute weight for distance
      sw  += wgt;                              // accumulate weight sum
      swv += wgt * y[k];                       // accumulate weighted values
      k--; }                                   // jump to next neighbour
    k = n+1;
    while(k < N  && (dist = x[k]-x[n]) <= w2){ // right side loop
      wgt  = weightFunc(dist * w2r);
      sw  += wgt;
      swv += wgt * y[k];
      k++; }
    avg[n] = swv / sw; }
}
// Weighting functions for nonuniform MA. Maybe implement more, see here:
// https://en.wikipedia.org/wiki/Kernel_(statistics)#Kernel_functions_in_common_use
// actually, for use in movingAverage, the checks against > 1 are superfluous bcs the loops there
// ensure already that x <= 1. maybe we can wrap all this stuff into a class 
// NonUniformMovingAverage. We could have functions like setNominalWidth (the width used for 
// uniform weighting), setWeightingFunction(T (*func)(T), T area) and for functions with unknown
// area just setWeightingFunction(T (*func)(T)) which finds the area by numerical integration
// and then calls setWeightingFunction(T (*func)(T), T area) with the numerically computed area
// other functions to try: (1-x^a)/(1+x^a), a = 1, 1/2, 1/3, 2, 
float uniform(float x)       // rename to uniform, maybe scale by 0.5
{
  if(fabs(x) > 1)
    return 0;
  return 1;
}
float triangular(float x)      // half-area: 0.5 (integral from 0 to 1)
{
  x = fabs(x);
  if(x > 1)
    return 0;
  return 1 - x;
}
float rationalTent(float x)   // half-area: log(4)-1
{
  x = fabs(x);
  if(x > 1)
    return 0;
  return (1-x) / (1+x);
}
float parabolic(float x)      // half-area: 2/3  (Epanechikov)
{
  if(fabs(x) > 1)
    return 0;
  return 1 - x*x;
}
float cubicBell(float x)
{
  return RAPT::rsPositiveBellFunctions<float>::cubic(fabs(x));
}
float quinticBell(float x)
{
  return RAPT::rsPositiveBellFunctions<float>::quintic(fabs(x));
}
float hepticBell(float x)
{
  return RAPT::rsPositiveBellFunctions<float>::heptic(fabs(x));
}
void nonUniformMovingAverage()
{
  // user parameters:
  static const int N = 200;  // number of samples
  float minDist = 1.f;       // minimum distance between successive samples
  float maxDist = 10.f;      // maximum ...
  float minY    = -10.f;     // minimum y value
  float maxY    = +10.f;     // maximum
  float width   =  30.f;     // filter support width for box, others are scaled by their area
  int   seed    = 0;

  // create input data:
  float x[N], y[N];
  float t = 0.f; // time
  float dt;      // time delta between samples
  x[0] = t;
  y[0] = (float)round(RAPT::rsRandomUniform(minY, maxY, seed)); 
  for(int n = 1; n < N; n++){
    //dt = (float)round(rsRandomUniform(minDist, maxDist));
    dt = (float)RAPT::rsRandomUniform(minDist, maxDist);
    t += dt;
    x[n] = t;
    y[n] = (float)RAPT::rsRandomUniform(minY, maxY); 
    //y[n] = (float)round(rsRandomUniform(minY, maxY)); 
  }

  // create filtered versions:
  float a = 1.f / float(log(4.f)-1); // reciprocal of area under rationalTent weighting function
                                     // ...the definite integral from 0 to 1
  float y0[N], y1[N], y2[N], y3[N], y4[N],
    yR[N],  // rational
    yP[N];  // parabolic

  // 0,1,2,3,4 is the smoothness of the weight function

  movingAverage(N, x, y, y0,  width,      uniform);
  movingAverage(N, x, y, y1,  width*2,    triangular);
  movingAverage(N, x, y, yR,  width*a,    rationalTent); 
  movingAverage(N, x, y, yP,  width*1.5f, parabolic);
  movingAverage(N, x, y, y2,  width*2,    cubicBell);
  movingAverage(N, x, y, y3,  width*2,    quinticBell);
  movingAverage(N, x, y, y4,  width*2,    hepticBell);

  // For all weighting functions other than the box, we scale the support width by factor to
  // equal to the reciprocal of the half of the area under the function to make the plots more 
  // easily comparable. 

  // plot:
  GNUPlotter plt;
  //plt.addDataArrays(N, x, y, y0, y1, y2, y3, y4);
  //plt.addDataArrays(N, x, y0, y1, yR);
  //plt.addDataArrays(N, x, y0, y1, y2, y3, y4, yR);
  //plt.addDataArrays(N, x, y0, y2, y4);
  plt.addDataArrays(N, x, y, y0, y1, yR, yP);
  plt.plot();

  // Observations: The rationalTent weighting seems best in terms of smoothing. It shows the 
  // smallest deviations from the average, i.e. the curve is the most "inside" of all (stays
  // closer to the mean than the others). 
  // ToDo: maybe try to find even better weighting functions. Maybe the smoothness is related to 
  // the area? smaller area -> more width widening -> better smoothing? maybe because the jitter
  // gets averaged out better?
}

void nonUniformOnePole1()
{
  // This test compares the output of the non-uniform one-pole filter to that of a regular uniform
  // one-pole. The sample spacing for the non-uniform filter *is* actually chosen uniformly for 
  // this comparison to make sense, but of course the algorithm in the non-unform filter is a 
  // different one - one that generalizes to the actual non-uniform case. Both filters get a chunk 
  // of noise as input signal and they should produce the same output.

  int N = 20;
  double fs = 1.0;   // sample rate
  double fc = 0.1;   // cutoff freq

  // Tests:
  // 1: compare to the output of a regular uniform one-pole - choose the time-stamps for the 
  //    non-uniform to be equal to the sample-numbers - the non-uniform filter should match the 
  //    output of the uniform one ...done: works
  // 2: look at the non-uniform impulse response - compare to analytic result
  // 3: maybe compute what the paper calls "ground" truth by oversampling - don't use random
  //    intervals for the nonuniform case but put the sample instants at times where we have
  //    actual oversampled data

  std::vector<double> x(N), yu(N), yn(N);                       // input- and output signals
  RAPT::rsArrayTools::fillWithRandomValues(&x[0], N, -1.0, 1.0, 1);  // fill input signal array

  // apply uniform filter:
  rsOnePoleFilter<double, double> fltUni;
  fltUni.setSampleRate(fs);
  fltUni.setCutoff(fc);
  fltUni.setMode(fltUni.LOWPASS_IIT);
  int n;
  for(n = 0; n < N; n++)
    yu[n] = fltUni.getSample(x[n]);

  // apply non-uniform filter:
  rsNonUniformOnePole<double> fltNonUni;
  typedef rsNonUniformOnePole<double>::NormalizeMode NM;
  //fltNonUni.setNormalizationMode(NM::noNormalization);
  //fltNonUni.setNormalizationMode(NM::spatiallyVariantScaling);
  fltNonUni.setNormalizationMode(NM::piecewiseResampling);
  fltNonUni.setOmega(2*PI*fc/fs);
  fltNonUni.reset();   // todo: figure out, if it makes a difference, which formula is used there
  for(n = 0; n < N; n++)
    yn[n] = fltNonUni.getSample(x[n], 1.0);

  GNUPlotter plt;
  plt.addDataArrays(N, &x[0]);
  plt.addDataArrays(N, &yu[0]);
  plt.addDataArrays(N, &yn[0]);
  plt.plot();

  // Observations:
  // -Comparison with uniform filter, noise input:
  //  -when calling getSample2, it seems to make no difference, whether or not we add Phi for the
  //   comparison with the uniform filter
  // -Comparison of non-uniform impulse response with analytic one:
  //  ...make a separate function for this test - then we may get rid of the time-axis here
}

void nonUniformOnePole2()
{
  // This test compares the impulse response of a non-uniform one pole to the ground truth which
  // is computed analytically. The samples of the filter's impulse reponse should fall on the 
  // continuous function defined by the analytic, continuous impulse-response

  int Nf = 50;   // number of samples taken from the filter
  int Nc = 500;  // number of dense (pseudo-continuous) samples for analytic target solution
  double dtMin = 0.2;   // minimum time-difference between non-uniform samples
  double dtMax = 1.8;   // maximum ..
  double fc    = 0.015; // cutoff freq
  double fs    = 1.0;   // sample rate
  double x = 1;         // 0: impulse response, 1: step response



  double wc    = 2*PI*fc/fs;   // normalized radian cutoff frequency


  std::vector<double> tf(Nf), yf(Nf);  // time-axis and output of filter
  std::vector<double> tc(Nc), yc(Nc);  // times axis and values for pseudo-continuous response

  // create time-stamp array for non-uniform sampling:
  typedef RAPT::rsArrayTools AR;
  tf = randomSampleInstants(Nf, dtMin, dtMax, 0);

    // compute non-uniform filter output
  rsNonUniformOnePole<double> flt;
  typedef rsNonUniformOnePole<double>::NormalizeMode NM;
  //flt.setNormalizationMode(NM::noNormalization);         // matches impulse response exactly
  //flt.setNormalizationMode(NM::spatiallyVariantScaling); // erratic around desired values
  flt.setNormalizationMode(NM::piecewiseResampling);       // matches step response exactly
  flt.setOmega(wc);
  flt.reset();
  yf[0] = flt.getSample(1.0, 1.0);
  for(int n = 1; n < Nf; n++)
    yf[n] = flt.getSample(x, tf[n]-tf[n-1]);

  // compute uniform filter output:
  std::vector<double> yu(Nf);
  rsOnePoleFilter<double, double> fltUni;
  fltUni.setSampleRate(fs);
  fltUni.setCutoff(fc);
  fltUni.setMode(fltUni.LOWPASS_IIT);
  yu[0] = fltUni.getSample(1.0);
  for(int n = 1; n < Nf; n++)
    yu[n] = fltUni.getSample(x);

  // create densely sampled (pseudo-continuous) impulse- or step response:
  double tMax = rsLast(tf);
  AR::fillWithRangeLinear(&tc[0], Nc, 0.0, tMax);
  if(x == 0) {  // plot continuous impulse response
    //double scaler = wc;  // wrong - i think, this normalizes the integral, not the sum?
    double scaler = fltUni.getB0();
    for(int n = 0; n < Nc; n++)
      yc[n] = scaler * exp(-tc[n]*wc);  // tau = 1/wc 
  }
  else { // if x == 1, plot the step-response:
    for(int n = 0; n < Nc; n++)
      yc[n] = 1 - exp((-tc[n]-1)*wc);   // verify formula
  }
  // https://en.wikipedia.org/wiki/RC_time_constant
  // http://www.dspguide.com/ch19/2.htm


  GNUPlotter plt;
  plt.addDataArrays(Nc, &tc[0], &yc[0]);
  plt.addGraph("index 0 using 1:2 with lines lw 2 lc rgb \"#808080\" notitle");
  plt.addDataArrays(Nf, &tf[0], &yf[0]);
  plt.addGraph("index 1 using 1:2 with points pt 7 ps 0.8 lc rgb \"#000000\" notitle");
  //plt.addDataArrays(Nf, &yu[0]);
  //plt.addGraph("index 2 using 1 with points pt 7 ps 0.8 lc rgb \"#008000\" notitle");
  plt.setPixelSize(1000, 250);
  plt.setGrid(false, false);
  plt.plot();

  // Observations:
  // -without normalization, the impulse response looks perfect but the step response looks 
  //  disturbed
  // -with piecewise resampling normalization, the step response looks perfect but the impulse
  //  response looks disturbed
  // -with spatially variant scaling, both look somewhat disturbed - it seems to be a compromise
  //  between good impulse-response and good step response (todo: verify, if it's actually
  //  correctly implemented)
  // -it seems like piecewise resampling is better than time-variant scaling (unless my 
  //  implementation of the latter is faulty - verify it!)
  // -actually, the impulse response with piecewise resampling looks "almost good" - the signal 
  //  just seems to be overall a bit too small - except for the first sample, which is exact
  //  -could this be an issue about the initial state?

  // ToDo:
  // -implement highpass..but how would the continuous highpass look like? a delta function minus 
  //  the lowpass response....but how would we represent that?
  // -maybe we can do a highpass by applying a bi-directional lowpass and then subtract the result
  //  from the original
  // -figure out, how the deviation of the step-response in the non-normalized case depends on the
  //  dt values - i think, for small dt (when samples are dense) it shoots above the target and when 
  //  samples are sparse, it is below the target ...maybe a scaling by 1/dt helps? -> i tried to
  //  scale by dt and 1/dt - and nope, that doesn't make any sense - but maybe adding something
  //  could help - but actually that's what Phi in the piecewise resampling method does
  // -try, if the complex implementation gives the same results, when the resonance(?) frequency is
  //  set to 0Hz
}

void nonUniformComplexOnePole()
{
  // We create a non-uniform decaying sine filter and compare its output to a uniform version
  // and a pseudo-continuous version (which is computed as an oversampled uniform filter). The
  // output samples of the non-unifrom filter should match the underlying pseudo-continuous output.


  // decaying sine parameters:
  double amplitude = 1.0;   // overall amplitude
  double phase     = 45;    // start-phase in degrees
  double decay     = 50;    // number of samples to decay to A/e
  double freq      = 0.05;  // normalized resonance frequency (freq/sampleRate)
  double in        = 0;     // 0: impulse response, 1: step response

  // sampling parameters:
  int Nf = 100;           // number of samples taken from the filter
  int oversampling = 10;  // oversampling factor for pseudo-continuous signal
  double dtMin = 0.2;     // minimum time-difference between non-uniform samples
  double dtMax = 1.8;     // maximum ..


  // create uniformly sampled impulse-response:
  typedef RAPT::rsArrayTools AR;
  double a[3], b[2];
  a[0] = 1.0;
  rsDampedSineFilterCoeffs(2*PI*freq, amplitude, decay, rsDegreeToRadiant(phase),
    &b[0], &b[1], &a[1], &a[2]);
  std::vector<double> x(Nf), yu(Nf);
  AR::fillWithValue(&x[0], Nf, in); x[0] = 1;  // create impulse or step as input signal
  AR::filter(&x[0], Nf, &yu[0], Nf, b, 1, a, 2);

  // create non-uniformly sampled impulse-response:
  std::vector<double> t(Nf), yn(Nf);
  t = randomSampleInstants(Nf, dtMin, dtMax, 0);
  std::complex<double> r, p;  // residue and pole
  rsDampedSineFilterResidueAndPole(b[0], b[1], a[1], a[2], &r, &p);
  rsNonUniformComplexOnePole<double> flt;
  flt.setCoeffs(r, p);
  yn[0] = 2*flt.getSampleReal(1.0, 1.0);
  for(int n = 1; n < Nf; n++)
    yn[n] = 2*flt.getSampleReal(in, t[n]-t[n-1]);

  // create pseudo-continuous impulse response (via oversampling):
  int Nc = Nf * oversampling;
  rsDampedSineFilterCoeffs(2*PI*freq/oversampling, amplitude, decay*oversampling,
    rsDegreeToRadiant(phase),  &b[0], &b[1], &a[1], &a[2]);
  std::vector<double> tc(Nc), yc(Nc);
  AR::fillWithValue(&yc[0], Nc, in); yc[0] = 1;  // create impulse or step as input signal
  //if(in == 1) AR::scale(&yc[0], Nc, 1./oversampling);  // ad hoc - not sure, if correct
  AR::fillWithRangeLinear(&tc[0], Nc, 0.0, (Nc-1.0)/oversampling);
  AR::filter(&yc[0], Nc, &yc[0], Nc, b, 1, a, 2);
  //AR::scale(&yc[0], Nc, yu[0]/yc[0]);  // ad hoc - nope!
  // pseudo-continuous step response is still wrong - why? ...looks scaled and shifted



  // plot:
  GNUPlotter plt;
  plt.addDataArrays(Nc, &tc[0], &yc[0]);  // pseudo-continuous data
  plt.addGraph("index 0 using 1:2 with lines lw 2 lc rgb \"#808080\" notitle");
  plt.addDataArrays(Nf, &t[0], &yn[0]);   // non-uniformly sampled data
  plt.addGraph("index 1 with points pt 7 ps 0.8 lc rgb \"#000000\" notitle");
  plt.addDataArrays(Nf, &yu[0]);          // uniformly sampled data
  plt.addGraph("index 2 with points pt 7 ps 0.8 lc rgb \"#0000ff\" notitle");
  plt.setGrid(false, false);
  plt.setPixelSize(1000, 250);
  plt.plot();
}

// next: 

// -figure out, how to implement a general non-uniform biquad 
//  -a uniform biquad can be represented as a unit-delayed two-pole-one-zero plus a (weighted)
//   input
//  -the gastal paper does not really address, how to implement the FIR part - but maybe the julia
//   code does? ...figure out, how the plots in figure 5 were created (impulse responses of
//   butterworth, cauer, etc.) - that should be applicable to our case
//  -maybe we should look at the continuous time transfer biquad function and obtain an expression
//   for the continuous time impulse response? ..and maybe work out the impulse-invariant transform
//   for a general analog biquad and then take a z-transform of a sampled version thereof
//   -> use sage for this
//  -maybe, for the time being, restrict ourselves to all-pole filters (this will include 
//   butterworth, bessel, gaussian, papoulis, halpern, cheby-1) to avoid this problem
//   -> we may also do a proper impulse-invariant transform with such filters

void nonUniformAllpole()  // rename - it now also includes filters with zeros (elliptic)
{
  // Test for a high-order non-uniform allpole filter (Butterworth, Bessel, etc.)

  int N = 1000;           // number of samples
  double d     = 0.9;     // amount of randomness in the sample spacing
  double dtMin = 1-d;     // minimum time-difference between non-uniform samples
  double dtMax = 1+d;     // maximum ..
  double fc    = 0.01;    // cutoff freq
  double x     = 1;       // 0: impulse response, 1: step response

  std::vector<int> orders = { 1,2,3,4,5,6,7,8 };

  typedef rsNonUniformFilterIIR<double>::ApproximationMethod AM;
  rsNonUniformFilterIIR<double> flt;
  //flt.setApproximationMethod(AM::butterworth);
  //flt.setApproximationMethod(AM::bessel);
  flt.setApproximationMethod(AM::gaussian);
  //flt.setApproximationMethod(AM::papoulis);
  //flt.setApproximationMethod(AM::halpern);
  //flt.setApproximationMethod(AM::elliptic);
  flt.setFrequency(fc);
  //flt.setOrder(order);
  // flt.setType(FT::lowpass);

  GNUPlotter plt;

  typedef std::vector<double> Vec;
  Vec h(N);  // impulse response
  Vec t = randomSampleInstants(N, dtMin, dtMax, 0);

  for(size_t i = 0; i < orders.size(); i++)
  {
    flt.setOrder(orders[i]);
    flt.reset();

    h[0] = flt.getSample(1.0, 1.0);
    for(int n = 1; n < N; n++)
      h[n] = flt.getSample(x, t[n]-t[n-1]);

    plt.addDataArrays(N, &t[0], &h[0]);
  }
  
  plt.setPixelSize(1000, 300);
  plt.plot();

  // Observations:
  // -without normalization: 
  //  -impulse-responses look good, but the step responses are a total mess
  // -piecewise resampling: 
  //  -step responses look good - but the first order ones approach a value != 1 (slightly above) 
  //   higher order filters don't have that problem
  //  -impulse responses look also ok, but the first sample of the first order responses seems too
  //   large
  // -spatially variant scaling:
  //  -both responses are a total mess - but that may not be surprising since we do the 
  //   normalization per stage and the paper says it must be done once for the whole filter

  // todo: check, why the first order filters have wrong DC gain 
  // -check initial conditions
  // -check bahavior at uniform sample rate -> the problem persists! -> look at the coeffs, compare
  //  to simple first order filter
  // -it seems when order = 1, it approaches the value of s (1.038 for gaussian) - maybe just 
  //  divide everything by s? ...or: divide the final result by the sum of s values? 
  //  ..seems to work!


  // -check gaussian cutoff normalization (maybe make it available in EngineersFilter/ToolChain
  //  and check mag-resp there
  // -allow filters with zeros (elliptic, etc.)
  //  -has been implemented but doesn't work - todo: try to decompose a uniform elliptic filter
  //   into partial fractions - or: try some more general filter that has a somewhat longer
  //   FIR part ...maybe 6th order numerator and 4th order denominator
  //  -it seems like for odd orders the 0th polynomial coeff is always 0 and for even order it's
  //   always 1 - but that seems to be wrong - we have an offset in the step response for even 
  //   order filters
  //  -but wait - what's the actual number of *finite* zeros in elliptic filters? could there be 
  //   one zero at infinity in certain cases? yes - there is, but that seems to be taken into
  //   account already
  //  -is it actually correct to directly add the polynomial part of the laplace trafo? what does
  //   a polynomial part in the s-domain actually mean? the inverse laplace trafo of a constant is
  //   a dirac impulse
  //  -> multiplying the FIR part by k fixes the offset - but the step responses still look strange
  //  todo: compare with uniform EngineersFilter
  //  -i think, it my due to the way i'm scaling - compare to paarmann, page 195 and 206 - the 
  //   even ones have a non-unity DC gain - i think, if i scale my even filters appropriately,
  //   the step response will look better - maybe allow for an optional setting that normalizes the
  //   filters like in paarmann - no overshoot over unit gain

}

void nonUniformBiquad()
{
  // We design a biquad filter (via RBJ cookbook formulas), separate it into a parallel connection 
  // of two complex one-pole filters (via partial fraction expansion) and apply it to non-uniformly
  // sampled data. We compare the resulting samples to a pseudo-continuous version of the same data
  // which we obtain via an oversampled, uniform filter.


  int Nf = 50;             // number of samples taken from the filter
  int oversampling = 10;   // oversampling factor for pseudo-continuous signal


  int Nc = Nf * oversampling; // number of samples for oversampled signal
}

// todo: for testing the complex bandpass filters later, maybe try to separate two sinusoids of
// different frequencies and/or an sinusoid buried in white noise

// next step: use a pair of complex one-poles to implement a biquad. can we get a formula for the 
// analytic biquad implese response? if not, obtain the pseudo-continuous impulse-response by means
// of taking an oversampled discrete time response 
// ...the partial fraction expansion of the biquad can be computed analytically ...somewhere, i have 
// already done this...
// ...maybe then try an N-th order butterworth 
// filter

void nonUniformBiDirectional()
{
  // We test our non-uniform bidirectional filter on a mix of 3 sinusoids with the goal to 
  // separate them.

  // user parameters:
  int N = 1000;       // number of samples
  double d  = 0.9;    // amount of randomness in the sample spacing
  double fs = 500;    // average sample rate
  double f1 =   5;    // first frequency in Hz
  double f2 =  10;    // second frequency in Hz
  double f3 =  20;    // third frequency in Hz
  double fc =   7;    // filter cutoff/center freq
  int order     = 20; // filter order 
  int numPasses = 1;  // number of filter passes

  // test: normalize the sample-rate to 1 and scale all frequencies accordingly:
  double s = 1/fs;
  fs *= s; f1 *= s; f2 *= s; f3 *= s; fc *= s;
  // not needed anymore since the dtScaler has been implemented in rsNonUniformFilterIIR
  // ...now s can be anything from very small to very large, the filter doesn't care

  // create input signals:
  typedef std::vector<double> Vec;
  double dtMin = (1-d)/fs;           // minimum time-difference between non-uniform samples
  double dtMax = (1+d)/fs;           // maximum ..
  Vec xn1(N), xn2(N), xn3(N), xn(N); // the 3 non-uniform sines and their sum
  Vec xu1(N), xu2(N), xu3(N), xu(N); // same signals uniformly sampled
  Vec tn = randomSampleInstants(N, dtMin, dtMax, 0);
  Vec tu = rsLinearRangeVector( N, 0.0, (N-1)/fs);
  // tn should also end at N-1 to avoid artifact at right boundary:
  RAPT::rsArrayTools::scale(&tn[0], N, tu[N-1]/tn[N-1]);
  for(int n = 0; n < N; n++) {
    xn1[n] = sin(2*PI*f1*tn[n]);
    xn2[n] = sin(2*PI*f2*tn[n]);
    xn3[n] = sin(2*PI*f3*tn[n]);
    xn[n]  = xn1[n]+xn2[n]+xn3[n];
    xu1[n] = sin(2*PI*f1*n/fs);
    xu2[n] = sin(2*PI*f2*n/fs);
    xu3[n] = sin(2*PI*f3*n/fs);
    xu[n]  = xu1[n]+xu2[n]+xu3[n];
  }

  // create uniformly and non-uniformly sampled filtered signals:
  Vec yu(N), yn(N);
  typedef RAPT::rsBiDirectionalFilter BDF;
  BDF::applyButterworthLowpass(&xu[0],         &yu[0], N, fc, fs, order, numPasses);
  BDF::applyButterworthLowpass(&xn[0], &tn[0], &yn[0], N, fc,     order, numPasses);

  // plot results:
  GNUPlotter plt;
  //plt.addDataArrays(N, &tn[0], &xn[0]);
  //plt.addDataArrays(N, &tu[0], &xu[0]);

  //plt.addDataArrays(N, &tu[0], &xu1[0]);  // lowest sine only
  plt.addDataArrays(N, &tu[0], &yu[0]);   // output (should be mostly the lowest sine)

  //plt.addDataArrays(N, &tn[0], &xn1[0]);  // lowest sine
  plt.addDataArrays(N, &tn[0], &yn[0]);   // output (should be mostly the lowest sine)

  plt.plot();

  //rsPlotVector(rsDifference(tn)); // just to inspect to dt values




  // Observations:
  // -isolating the lowest sine via lowpass works
  // -todo: implement highpass, bandpass and bandreject and use them to isolate the upper or middle
  //  sine ...maybe implement peaking and shelving types, too - although, i currently don't see any
  //  use for them...maybe to boost certain modulating frequencies in the sine model?
  // -todo: let the uniform filter use impulse-invariant transform too - then, the outputs should
  //  be exactly the same when d=0



  // -todo: figure out, what the "best" sample-rate is in terms of numerical precision - tweak the
  //  s-factor above -> it seems like when going down like s=4/fs, s=2/fs, s=1/fs, the last is 
  //  the first, for which the non-uniformly sampled signal is not above the uniformly sampled one
  //  s=0.5/fs also still works, but 0.25/fs triggers a singluar-matrix assert (all tested with 4th
  //  order), s=1.5 looks really good - could this be a hint to pi/2? is this plausible for an 
  //  optimal spacing? why would pi show up? ...but maybe the comparison should be done again when 
  //  we have impulse invariant design in place for the uniform filter...
  //  -there is now a dtScale parameter in rsNonUniformFiterIIR that normalizes everything with
  //   respect to the sample-rate - but may try to figure out, if the selected operating
  //   point is really chosen optimally 




  // Conclusion:
  // We should normalize the average time-delta between samples to be of the order of unity and 
  // scale the cutoff frequency of the non-uniform filter accordingly.
}

void smoothingFilterOrders()
{
  // We plot the step responses of the rsSmoothingFilter for various orders.

  static const int numOrders = 4;   // number of filters with different orders
  bool expSpacing = false;          // if true, orders are 1,2,4,8,.. else 1,2,3,4,..
  int orderIncrement = 2;           // or 1,3,5,.. or 1,4,7,...
  static const int N = 300;         // number of samples
  float fs  = 1.0f;                  // sample rate
  float tau = 101.0f;                  // time constant

  // create and set up the smoother:
  rsSmoothingFilter<float, float> smoother;
  //rsSmoothingFilterDD smoother;
  smoother.setTimeConstantAndSampleRate(tau, fs);
  //smoother.setNumSamplesToReachHalf(tau*fs);
  //smoother.setShape(rsSmoothingFilterFF::FAST_ATTACK);
  smoother.setShapeParameter(0.5f);

  // compute step responses:
  int order = 1;
  float y[numOrders][N];
  for(int i = 0; i < numOrders; i++)
  {
    smoother.setOrder(order);
    if(expSpacing) // update order for next iteration
      order *= 2;
    else
      order += orderIncrement;

    smoother.reset();
    y[i][0] = smoother.getSample(1.f);
    for(int n = 1; n < N; n++)
    {
      //y[i][n] = smoother.getSample(0.f);   // impulse response
      y[i][n] = smoother.getSample(1.f); // step response
    }
  }

  // plot:
  float t[N];
  createTimeAxis(N, t, fs);
  GNUPlotter plt;
  for(int i = 0; i < numOrders; i++)
    plt.addDataArrays(N, t, y[i]); 
  plt.plot();

  // Observations:
  // The step responses of the different orders are comparable in terms of overall transition time.
  // However, they do not meet in a common point. From a user's perspective, it would perhaps be
  // most intuitive, if he could just set up the time-instant, where the step-response goes through
  // 0.5. No matter what the order is, the time instant where it passes through 0.5 should remain 
  // fixed. To achieve that, i think, we need a closed form expression for the step-response and 
  // then set that expression equal to 0.5 and solve for the time-constant. If it's too hard to 
  // find such a formula, we could also create a table by just reading off the time-instants where
  // the various step responses go through 0.5 and then use that value as divider.

  // ToDo:
  // Measure the sample instants where the the step responses pass through 0.5 for some normalized
  // filter setting and create tables from that. These tables will eventually be used for scaling. 
  // It should be a 2D table, 1st dimension is the order (say 1..32) and 2nd dimension is the value
  // of the shape parameter (maybe from 0 to 4 in 32 steps?) the 2nd dimenstion can use linear
  // interpolation, the 1st dimension doesn't need that because inputs are integers anyway.
  // The values of the table should be overall scalers for all the time constants.


  // try higher orders - like 100 - see what kind of shaped is approached for order -> inf
}

void smoothingFilterTransitionTimes()
{
  // We plot the step response for smoothing filters of a given order (and shape/asymmetry) for
  // different transition time constants. What we should see is time-stretched/compressed versions
  // of the same response. But the nonworking tuning tables suggest that somthing might be wrong
  // with that assumption, so we check.

  // user parameters:
  static const int N = 400;           // number of samples
  int order = 5;                      // order of the filters
  float asymmetry = 1.f;              // asymmetry parameter
  static const int numFilters = 5;    // number of transition times
  float transitionTimes[numFilters] = { 51.f, 101.f, 151.f, 200.1f, 251.f }; // transition times in samples


  // create and set up the smoother:                                                             
  rsSmoothingFilter<float, float> smoother;
  smoother.setShapeParameter(asymmetry);
  smoother.setOrder(order);

  // create the step responses:
  float y[numFilters][N];
  for(int i = 0; i < numFilters; i++)
  {
    smoother.setNumSamplesToReachHalf(transitionTimes[i]);
    smoother.reset();
    for(int n = 0; n < N; n++)
      y[i][n] = smoother.getSample(1.f);
  }

  // plot:
  GNUPlotter plt;
  for(int i = 0; i < numFilters; i++)
    plt.addDataArrays(N, y[i]); 
  plt.plot();

  // Observations:
  // The intuitive assumption that scaling all the time constants by the same factor has the effect
  // of stretching/compressing the step-response by that amount turns out to be false. It seems to
  // be true for 1st order filters but for higher orders, especially when there's asymmetry, the 
  // actually observed step repsonses deviate from that simple rule more and more. 
  // Why is that? Is it because the design formula uses impulse-invariant transform and should use
  // step-invariant? ...or maybe even bilinear?

  // see here: http://web.cecs.pdx.edu/~tymerski/ece452/Chapter4.pdf
  // http://homes.esat.kuleuven.be/~maapc/Sofia/slides_chapter8.pdf
  // http://dspcan.homestead.com/files/IIRFilt/zfiltsii.htm
}




template<class T>
void reflectRoots(complex<T>* roots, int N) // rename to mirrorFirstHalf, make generic, add to library
{
  for(int i = 0; i <= N/2; i++)
    roots[N-1-i] = conj(roots[i]);
}
template<class T>
void removeInfiniteValues(vector<complex<T>>& z)
{
  for(size_t i = 0; i < z.size(); i++) {
    if(RAPT::rsIsInfinite(z[i])) {   // as rs prefix
      z.erase(z.begin() + i);
      i--; }}
}
template<class T>
rsFilterSpecificationZPK<T> getFilterSpecificationZPK(RAPT::rsPrototypeDesigner<T>& pd)
{
  int nz = pd.getNumFiniteZeros();
  int np = pd.getNumFinitePoles();
  vector<complex<T>> z(np);          // np is correct, we may get infinite zeros which will be..
  vector<complex<T>> p(np);          // ..removed later
  pd.getPolesAndZeros(&p[0], &z[0]); // returns only the non-redundant upper halfplane poles/zeros
  reflectRoots(&p[0], np);           // create full pole/zero set by reflection
  reflectRoots(&z[0], np);
  removeInfiniteValues(z);
  rsFilterSpecificationZPK<T> spec;
  spec.p = p;
  spec.z = z;
  spec.sampleRate = RS_INF(T);
  //spec.gain  = pd.getGain();
  return spec;
}
void prototypeDesign()
{
  //rsEngineersFilterFF ef;

  /*
  // something is wrong with odd order bessel filters (exept for order = 1), looks like it never 
  // worked in rapt - ive gone back to the commit of 21.11.2017, 00:06:28 with message:
  // implemented FilterPlotter<T>::plotMagnitude
  // code for debugging :
  // compare float/double versions:
  rsPrototypeDesignerD pdD;
  complex<double> polesD[5], zerosD[5];
  pdD.setApproximationMethod(pdD.BESSEL);
  pdD.setOrder(3);
  pdD.getPolesAndZeros(polesD, zerosD);

  rsPrototypeDesignerF pdF;
  complex<float> polesF[5], zerosF[5];
  pdF.setApproximationMethod(pdF.BESSEL);
  pdF.setOrder(3);
  pdF.getPolesAndZeros(polesF, zerosF);

  int dummy = 0;
  // ok, yes - that's the problem: the double precision poles are correct (they match with what the 
  // rosic version computes) but the single precision float poles are wrong. It's probably 
  // something about the root finder
  // in rsPrototypeDesigner<T>::makeBesselLowShelv the pickNonRedundantPolesAndZeros fills the p 
  // member array differently for double and float versions, the pTmp array is the same in both
  // cases

  // in rsZeroNegligibleImaginaryParts the "real" pole has some small negative imaginary part which 
  // doesn't get zeroed out and in the subsequent call to rsOnlyUpperHalfPlane, the real pole
  // gets thrown away - we need a different threshold for float than we use for double
  // ...ok - seems to be fixed.
  */

  typedef float Real;
  typedef RAPT::rsPrototypeDesigner<Real> PD;


  // min and max filter order to plot:
  int minOrder = 1;
  int maxOrder = 10;

  // create and set up prototype designer:
  PD pd;
  pd.setApproximationMethod(PD::GAUSSIAN);
  //pd.setApproximationMethod(PD::BESSEL);
  //pd.setApproximationMethod(PD::BUTTERWORTH);
  //pd.setApproximationMethod(PD::PAPOULIS);
  //pd.setApproximationMethod(PD::HALPERN);
  //pd.setApproximationMethod(PD::CHEBYCHEV);
  //pd.setApproximationMethod(PD::INVERSE_CHEBYCHEV);
  //pd.setApproximationMethod(PD::ELLIPTIC);

  //pd.setPrototypeMode(PD::LOWSHELV_PROTOTYPE); // comment fo lowpass
  pd.setGain(+6.02f); // needs to be nonzero for plots

  pd.setPassbandRipple(1); 
  pd.setStopbandRejection(20);

  // create plotter, add filter specs for the desired orders to it and plot:
  FilterPlotter<Real> plt;
  for(int i = minOrder; i <= maxOrder; i++)
  {
    pd.setOrder(i);
    rsFilterSpecificationZPK<Real> spec = getFilterSpecificationZPK(pd);
    plt.addFilterSpecificationZPK(spec);
  }
  //plt.plotPolesAndZeros(600);
  plt.plotMagnitude(1000, 0, 3, false, false);

  // issues:

  // it seems, even order elliptic prototypes have an overall gain equal to the reciprocal of the 
  // linear stopband rejection (passband ripple seems to have no effect), odd order elliptics have
  // a different gain and it seems to depend on the ripple (might be computed by evaluating DC 
  // gain). Papoulis design has also a wrong DC gain -> add overall gain to the prototype designer
  // in EngineersFilter, it's not a problem, because i have a gain normalization step at the very
  // end of the design pipeline - but it would be desirable to have correct gains at all stages
  // of the design process rather than renormalizing at the end (which is sort of dirty)

  // todo: test the gaussian filter design
}


// this is for testing the rewrite attempt - which i then decided to give up:

// helper function to convert form raw arrays of poles and zeros to the FilterSpecificationZPK 
// structure used by the plotter
template<class T>
rsFilterSpecificationZPK<T> analogPrototypeSpecZPK(
  int N, std::complex<T>* za, std::complex<T>* pa, T k)
{
  vector<complex<T>> z = RAPT::toVector(za, N);
  vector<complex<T>> p = RAPT::toVector(pa, N);
  reflectRoots(&p[0], N);           // create full pole/zero set by reflection
  reflectRoots(&z[0], N);
  removeInfiniteValues(z);
  rsFilterSpecificationZPK<T> spec;
  spec.p = p;
  spec.z = z;
  spec.sampleRate = RS_INF(T); // indicates analog (s-plane) poles and zeros
  spec.k = k;
  return spec;
}

void poleZeroPrototype()
{
  // tests the new implementation of analog filter prototype design

  typedef float Real;
  typedef PoleZeroPrototype<Real> PZP;
  static const int minOrder = 4;   // minimum filter order to plot
  static const int maxOrder = 4;  // maximum filter order to plot
  Real G0 = 0.f;    // use G0=0, G=1 for lowpass and G0=2, G=8 for low-shelf
  Real G  = 1.f;     

  PZP pzp;
  pzp.setApproximationMethod(PZP::BUTTERWORTH);

  // maybe the new version should use optimized versions of halpernT2, papoulisL2 from 
  // Prototypes.h/cpp (avoids dynamic memory allocation)

  Real k; 
  std::complex<Real> p[maxOrder], z[maxOrder]; // (maxOrder+1)/2 should be enough
  FilterPlotter<Real> plt;
  for(int n = minOrder; n <= maxOrder; n++)
  {
    pzp.setOrder(n);
    pzp.getPolesZerosAndGain(p, z, &k);
    rsFilterSpecificationZPK<Real> spec = analogPrototypeSpecZPK(n, z, p, k);
    plt.addFilterSpecificationZPK(spec);
  }
  plt.plotPolesAndZeros(600);
  //plt.plotMagnitude(500, 0, 3, false, false);

  // something is wrong...

  //int dummy = 0;
}

void seriesConnectionDecay()  
{
  // We compare the impulse responses of a series connection of two first order lowpass filters 
  // with the impulse response of both filters alone.
  // ..well, we actually also look at a parallel connection - rename the function - maybe also look
  // at the difference

  int    N  = 501;
  double fs = 44100;
  double f1 = 100;
  double f2 = 200;
  double f3 = 200;



  RAPT::rsOnePoleFilter<double, double> flt;
  flt.setMode(flt.LOWPASS_IIT);
  //flt.setMode(flt.ALLPASS_BLT);
  flt.setSampleRate(fs);

  using Vec = std::vector<double>;

  flt.setCutoff(f1);
  Vec y1 = impulseResponse(flt, N, 1.0);
  flt.setCutoff(f2);
  Vec y2    = impulseResponse(flt, N, 1.0);
  Vec ys12  = filterResponse( flt, N, y1);   // serial 1->2
  Vec yp12  = y1 + y2;                       // parallel 1+2
  Vec yd12  = y1 - y2;                       // difference - doesn't start at zero
  flt.setCutoff(f3);
  Vec y3    = impulseResponse(flt, N, 1.0);
  Vec ys123 = filterResponse( flt, N, ys12); // serial 1->2->3


  // maybe compute the number of samples, after which they fall below a given threshold 
  // (like -60dB) ...just look at the output signals

  double minDb = -100;
  Vec db1    = ampToDb(y1,    minDb);
  Vec db2    = ampToDb(y2,    minDb);
  Vec db3    = ampToDb(y3,    minDb);
  Vec dbs12  = ampToDb(ys12,  minDb);
  Vec dbp12  = ampToDb(yp12,  minDb);
  Vec dbs123 = ampToDb(ys123, minDb);


  //rsPlotVectors(y1,  y2,  ys12,  yp12);
  //rsPlotVectors(db1, db2, dbs12, dbp12);
  rsPlotVectors(db1, db2, db3, dbs12, dbs123);
  //rsPlotVectors(y1, y2, ys, yp, yd);

  // todo: 
  // -find the T60 experimentally
  // -what happens, if we put a 3rd filter in series
  // -maybe, as a first approximation, use the sum of the T60s of the individual stages?
  // -the x-offset of the serial connection with respect to the longer filter seems to be equal to 
  //  the x-coordinate of the point, where 1st, 2nd and the serial filter meet - so, to find it,
  //  we'd have to equate a1*exp(-t/tau1) = a2*exp(-t/tau2)
  //  -> log(a1) - t/tau1 = log(a2) - t/tau2 -> solve for t
  // -but how does this generalize to more than two filters? ..it seems to be also the location of
  //  the peak of the serial filter - perhaps that should be used as the offset in general
  // -so in general, use the decay-time of the longest tail and add to that the location of the 
  //  peak-amplitude?
  //
  // https://www.kvraudio.com/forum/viewtopic.php?f=33&t=533696
}

void quantileFilterElongation()
{
  // We try to produce a sample that a length L+1 *would have* produced with a length L filter.

  double q = 0.2;    // quantile
  int    L = 2;      // length of non-elongated filter

  using Vec = std::vector<double>;
  Vec x;
  // x = Vec({ 2,4,6,8 });
  // x = Vec({ 0,2,4,6,8,10,12,14 });
  //x = Vec({2,4,6,8,6,4,2,0,-2,-4,-6,-8,-6,-4,-2,0,2,4,6,8,6,4,2,0,-2,-4,-6});
  x = rsRandomIntVector(200, -9, +9, 0);
  int N = (int) x.size();   // number of samples

  rsQuantileFilterCore<double> flt;
  flt.setMaxLength(32);

  // produce a target signal with a filter that actually *is* one sample longer than L:
  flt.setLengthAndQuantile(L+1, q);
  Vec t(N); // t for target
  for(int n = 0; n < N; n++)
    t[n] = flt.getSample(x[n]);


  rsDelayBuffer<double> dly;
  dly.setCapacity(L+1);
  dly.setLength(L+1);

  flt.reset();
  flt.setLengthAndQuantile(L, q);
  Vec y(N);  // normal output of filter
  Vec z(N);  // output of elongated filter - should match t
  for(int n = 0; n < N; n++)
  {
    dly.getSample(x[n]);       // feed the delayline
    y[n] = flt.getSample(x[n]);
    z[n] = flt.getElongatedOutput(dly[L]);
  }
  Vec err = t-z;

  double xOld = dly[L];
  //double branch;
  double tmp  = flt.getElongatedOutput(xOld);

  //rsPlotVectors(t, z, err);
  rsPlotVectors(x, t, z, err);

  // Observations: it works

  // todo: move explanation into comment in readOutputWithOneMoreInput
  // -if p1==p and xOld falls into the large heap, we just need to use xS = S0; xL = L0; as usual
  //  because xOld the large heap can actually admit for xOld as additional sample, because it's
  //  one sample longer - no data would have to be moved from large to small. If, on the other 
  //  hand, xold falls into the small heap, we have to take into account that a datapoint from
  //  the front of the small heap would have to be moved over to the large heap.
  // -if p1==p+1, the situation is reversed: when xOld falls into the small heap, we use 
  //  xS = S0; xL = L0; as usual and if xOld falls inot the large heap, some additional logic is
  //  required. try L=5, q= 0.6
}

void quantileFilterSweep()
{
  // We test the non-integer length quantile filter by making a continuous sweep of the length 
  // (and/or the quantile, but that's not really relevant for this test) and write the result to 
  // a wavefile. For comparison, we also generate a signal with a filter that is restricted to
  // integer lengths.

  // user parameters:
  double fs          = 44100;    // sample rate
  int    N           = 200000;   // number of samples
  double minLength   = 2.0;      // minimum length in sweep
  double maxLength   = 20.0;     // maximum length in sweep
  double minQuantile = 0.5;
  double maxQuantile = 0.5;

  // Create lengs and quantile sweep and in- and output signals:
  using Vec = std::vector<double>;
  using AT  = rsArrayTools;
  Vec lengths(N), quantiles(N);
  AT::fillWithRangeExponential(&lengths[0],   N, minLength,   maxLength);
  AT::fillWithRangeExponential(&quantiles[0], N, minQuantile, maxQuantile);
  Vec x = rsRandomVector(N, -1.0, 1.0, 0);  // input to the filter
  Vec y1(N), y2(N);                         // rough and smooth sweep

  // Create filter and produce the non-smooth sweep - it's not smooth because we do not yet 
  // assign the delay buffer, which is a requirement to make a smooth sweep work:
  rsQuantileFilterCore2<double> flt;
  flt.setMaxLength((int)ceil(maxLength));
  for(int n = 0; n < N; n++) {
    flt.setLengthAndQuantile(lengths[n], quantiles[n]);
    y1[n] = flt.getSample(x[n]); }

  // Create and assign the delay-buffer and produce a smooth sweep:
  rsDelayBuffer<double> delayLine;
  delayLine.setCapacity((int)ceil(maxLength) + 1); // +1 for the 1 extra elongation sample 
  flt.setDelayBuffer(&delayLine);
  flt.reset();
  for(int n = 0; n < N; n++) {
    delayLine.getSample(x[n]);
    flt.setLengthAndQuantile(lengths[n], quantiles[n]);
    y2[n] = flt.getSample(x[n]); }

  // plot both sweeps:
  //rsPlotVectors(y1, y2);

  // write both sweeps to files:
  rosic::writeToMonoWaveFile("QuantileFilterSweepRough.wav",  &y1[0], N, 44100);
  rosic::writeToMonoWaveFile("QuantileFilterSweepSmooth.wav", &y2[0], N, 44100);
}

void quantileFilterDelay()
{
  // Creates a sinewave and applies an rsQuantilefilter to it with a given length and various 
  // settings of the quantile and plots the results together with a suitable delayed version of the
  // input. 

  int    N  = 500;    // number of samples for plot
  double L  = 251;    // length of filter
  double lo = 1.0;    // lowpass gain
  double hi = 0.0;    // highpass gain
  double P  = 100.0;  // period of input wave

  // quantiles:
  using  Vec = std::vector<double>;
  //Vec q({0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0});
  Vec q({0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0});
  //Vec q({0.0, 0.25, 0.5, 0.75, 1.0});
  //Vec q({0.5});

  // Create and set up rsQuantileFilter object:
  rsQuantileFilter<double> flt;
  int maxLength = (int) ceil(L);
  flt.setSampleRateAndMaxLength(1.0, maxLength);
  flt.setFrequency(1.0/L);
  flt.setLowpassGain(lo);
  flt.setHighpassGain(hi);
  flt.updateInternals();        // because we want to retrieve the delay below

  // Create input waveform and delayed version of it:
  Vec x(N); createWaveform(&x[0], N, 0, 1./P, 1.0);
  rsDelayBuffer<double> dly(maxLength);
  double d = flt.getDelayInSamples();
  Vec xd(N);
  for(int n = 0; n < N; n++) {
    dly.getSample(x[n]);
    xd[n] = dly[d]; }

  // create filtered outputs for various quantile settings and plot them along with the 
  // delayed input:
  GNUPlotter plt;
  plt.addDataArrays(N, &xd[0]);  
  Vec y(N);
  for(size_t i = 0; i < q.size(); i++)
  {
    flt.setQuantile(q[i]);
    flt.reset();
    for(int n = 0; n < N; n++)
      y[n] = flt.getSample(x[n]);
    plt.addDataArrays(N, &y[0]);
  }
  plt.plot();

  // Observations:
  // -when lo = hi = 1, the output matches the delayed input, as it should
  // -L=50,lo=1,hi=0: higher and lower quantiles bulge out the top or bottom halfwaves respectively
  //  and narrow the other halfwaves - the corresponding "highpass" waves (lo=0,hi=1) look quite 
  //  interesting
  // -L=P=100,lo=1,hi=0: the outputs become contants - increasing L to 120, we see weird 
  //  small-scale wiggles in the output (why? that seems strange! is that a bug?)
  //  -the effect is even more pronounce for odd L, like 121 or 151 - the output stays contant over
  //   3 samples, leading to stairsteps - with L = 251, the steps are 5 samples wide
}

// nove to test uitilities:
std::vector<double> applyDelay(const std::vector<double>& x, double delay)
{
  size_t N = x.size(); // maybe optionally use N + ceil(delay)
  rsDelayBuffer<double> dly((int)ceil(delay));
  std::vector<double> y(N);
  for(size_t n = 0; n < N; n++) {
    dly.getSample(x[n]);  
    y[n] = dly[delay];   }
  return y;
}

void quantileFilterDual()
{
  double fs = 44100;  // sample rate
  int    N  = 1000; // number of samples
  int    L  = 100;    // filter length in samples (can we make this a double, too?)
  double q  = 1.0;   // filter quantile, 0.0: minimum, 0.5: median, 1.0: maximum
  double f1 = fs/2;
  double f2 = fs/500;  // fs/256 is a nice end value for a sweep

  f1 = f2 = fs/51;

  // Create and set up filter:
  using QF  = rsDualQuantileFilter<double>;
  using Vec = std::vector<double>;
  using AT  = rsArrayTools;
  double maxLength = ceil(rsMax(1/f1, 1/f2)); // required maximum length
  QF flt;
  flt.setSampleRateAndMaxLength(fs, maxLength);
  flt.setQuantile(q);
  flt.setLowpassGain(0.5);
  flt.setHighpassGain(0.0);
  //flt.setFeedback(0.0);    // later - maybe
  flt.setCore2Complementary();
  flt.updateInternals();  // so we have a non-dirty state to look at

  // Create input signal:
  //Vec x = rsRandomVector(N, -1.0, +1.0, 0);  // try sinusoids, too
  //Vec x = rsRandomIntVector(N, 1, 9, 0);  // looks not random at all - very periodic!
  //Vec x = rsRandomIntVector(N, 0, 16, 0);   // periodic with period 16
  //Vec x = rsRandomIntVector(N, 1, 99, 0);
  // todo: use a irwin-hall generator - we can see it better when the distribution is more skewed
  // toward the middle
  //createWaveform(&x[0], N, 0, 1.0/L, 1.0);  // sine wave
  //createSineSweep(&x[0], N, 0.0/L, 2.0/L);
  Vec x(N); AT::fillWithImpulse(&x[0], N, 1.0, N/2);  // for testing the delay
  //Vec x = createCrackle(N, 0.02);

  // Create frequency sweep, delayed input and output signal:
  Vec f(N), y(N), xd(N);
  AT::fillWithRangeExponential(&f[0], N, f1, f2);
  double delay;   // maybe use later to compare output with delayed input
  rsDelayBuffer<double> dly((size_t)ceil(maxLength*fs));
  for(int n = 0; n < N; n++) {
    flt.setFrequency(f[n]);
    flt.setCore2Complementary();
    flt.updateInternals();
    delay = flt.getDelayInSamples();
    dly.getSample(x[n]);
    xd[n] = dly[delay];
    y[n] = flt.getSample(x[n]); }


  /*
  // create a median-filtered version of x (obsolete):
  int nS = L/2;    // floor division
  int nL = L-nS;
  Vec t(N);
  //for(int n = nS; n < N-nL; n++)
  //  t[n+nS] = rsArrayTools::median(&x[n-nS], nS+nL); // why t[n+nS]?

  // 7: 36,67,76,41,82,45,74 -> 36,41,45,67,74,76,82 -> 67
  // 6: 36,67,76,41,82,45    -> 36,41,45,67,76,82    -> (45+67)/2 = 56 ..but occurs at sample 101 - why?
  // so, for even lengths, t is lagging one sample

  // let's try it with the naive implementation of rsQuantileFilter
  rsQuantileFilterNaive<double> fltN(nS, nL);
  Vec z(N);
  //for(int n = 0; n < N; n++)
  //  z[n] = fltN.getSampleMedian(x[n]);
  */


  //rosic::writeToMonoWaveFile("QuantileFilterInput.wav",  &x[0], N, 44100);
  //rosic::writeToMonoWaveFile("QuantileFilterOutput.wav", &y[0], N, 44100);
  //std::cout << "Files written.";

  //rsPlotVectors(x, y);
  rsPlotVectors(x, xd, y);
  //rsPlotVectors(y);
  //rsPlotVectors(x, y, t);
  //rsPlotVectors(y, t);
  //rsPlotVectors(t, z);
  //rsPlotVectors(t, z, y);
  //rsPlotVectors(z, y);  // for even lengths, z is the better reference - t has a delay there
  //rsPlotVectors(x, y);



  // ToDo: 
  // -check, if the computed delay corresponds to what we see in the plots - ok - looks good

  // maybe try it with a square-wave with period 100, set the length also to 100 - this is an even
  // number, so we should get an interpolation coeff of 0.5 - i think, the output should be a 
  // triangle wave

  // how about using 4 of them in series and building the feedback loop around the whole thing?

  // Observations:
  // -the outputs generally have a sort of "sample-and-hold" character
  // -when sweeping the freq, there are audible artifacts that sound like switching to a new
  //  settings - could the be because of the rounding? it's especially obvious for very small 
  //  lengths - yes, it happens most obviously when the (rounded) length L switches from 2 to 3 
  //  -> maybe we need indeed implement some sort of support for non-integer length
  // -the highpass signal has strong components at the Nyquist freq when q=0 or q=1...this seems
  //  to have disappeared - may have been due to the false delay computation (delay was off by
  //  a factor of 1/2
  // 

  // ToDo:
  // -setMaxLength call takes really long - figure out why and fix it! ah - it was because it
  //  wants its input in seconds and i assumed samples and passed an unreasonably large value
  // -test, if the delay computation is correct - use an impulse at sample 100 - a non-delayed
  //  maximum should have a 1 at samples 98,99,100,101,102 - that's our reference. the filter
  //  output will have ones at 100..104, so the computed delay should be 2 
  //  -it's most easy to test when using loGain = hiGain = 1, in which case the filter produces
  //   a pure delay
  // -benchmark against naive implementation
  // -try putting a regular bandpass in front
  // -test modulation 
  //  -switch done - looks good
  //  -try sweeping the freq exponentially
  // -try a sine-sweep as input - looks weird, but that's expected
  // -figure out if we should use a scale factor when converting from frequency to length (it 
  //  seems mathematically natural to just use length = 1/freq, but maybe that's not a good choice
  //  in practice) - maybe look at frequency responses of moving-average filters to find the
  //  cutoff freq as function of filter length
  // -how about introducing feedback?
  // -can the filter length somehow be made a non-integer, too? what would it mean to take the 
  //  median over 10.7 samples? maybe something like 0.3 * median_over_10 + 0.7 * median_over_11?
  //  -maybe the double-heap should use the ceil and the calulation would not only use the front
  //   element but also one of the children? ...or should a non-integer length somehow go into 
  //   the calculation of the interpolation weight w? or should we use a modified quantile? we 
  //   compute the non-integer read position q as:
  //     T q = quantile * sampleRate * (*L - 1);
  //   maybe we can just use this formula as is but with the non-integer L?
  //  ...but no: the artifacts also appear when q = 0 or q = 1, so adjusting q seems not to be the 
  //  solution - we really need soem sort of formula/algo that computes the quantile over L and 
  //  over L+1 samples and crossfade between them - how can we do this? maybe we need to look not 
  //  only at the front samples of the heaps but also to the first child-nodes *and* at at a sample
  //  at L+1 to figure out, where fits in? between the middle samples
  // -support of non-integer length would allow for a smoother modulation of the cutoff freq, so
  //  it may be a worthwhile feature
  // -i think, this filter could be especially interesting to create filtered noises, it tends
  //  to some sort of "sample-and-hold" behavior - maybe try in conjunction with the TriModalNoise
  // -maybe re-introduce the second core and give it similar parameters - it may make sense to 
  //  form a linear combination of min and max, for example - but maybe do this in a subclass
  //  rsDualQuantileFilter

  //  ...maybe it makes sense to give presets in ToolChain and its submodules an Info field which
  //  can be edited in the plugins - there, we could store a text that explains what this preset
  //  does and hwo it should be used - a short info could be displayed in the info field when
  //  the mouse is over the preset field and a longer info coul be displayed in a popup text
  //  editor on right-click

}

void quantileFilterResonant()
{
  // Test the pseudo-resonance for the quantile filter that is introduced by using an additional 
  // min-max filter

  double fs = 44100;  // sample rate
  double f1 = 20000;  // filter frequency 1
  double f2 = 20;     // filter frequency 2
  double r  = 1.0;    // resonance mix
  int    N  = 200000; // number of samples

  rsQuantileFilterResonant<double> flt;
  flt.setSampleRate(fs);
  flt.setMaxLength(1.0);   // allows for 1 Hz as lowest cutoff (reciprocal of minFreq)
  //flt.setFrequency(f);
  flt.setResonanceMix(r);

  // Create input signal:
  using Vec = std::vector<double>;
  Vec x = rsRandomVector(N, -0.5, +0.5, 0);  // try sinusoids, too

  rsOnePoleFilter<double, double> lpf;
  lpf.setSampleRate(fs);
  lpf.setMode(rsOnePoleFilter<double, double>::LOWPASS_IIT);
  lpf.setCutoff(f2); // maybe use something else
  for(int n = 0; n < N; n++)
    x[n] = lpf.getSample(x[n]);


  //rsArrayTools::cumulativeSum(&x[0], &x[0], N);

  Vec f(N);
  rsArrayTools::fillWithRangeExponential(&f[0], N, f1, f2);

  // produce output signal:
  Vec y(N);
  for(int n = 0; n < N; n++)
  {
    flt.setFrequency(f[n]);
    y[n] = flt.getSample(x[n]);
  }

  //rsPlotVectors(x, y);
  rosic::writeToMonoWaveFile("ResoQuantileIn.wav",  &x[0], N, 44100);
  rosic::writeToMonoWaveFile("ResoQuantileOut.wav", &y[0], N, 44100);

  // Observations:
  // -when creating a sweep, the resonance gets louder for lower frequencies - i think, it's 
  //  because the min/max values are taken over a longer interval, so the amplitude tends to 
  //  increase
  // -maybe the length over which min and max are taken should be a fixed (adjustable) parameter
  //  that is independent from the filter frequency setting - at least, when the resonance 
  //  frequency is determined by the bandpass - but that tends to create weird aliasing (alike) 
  //  artifacts at high resonance frequencies
  // -for high frequencies, the resonance turns itself into noise
  //  -maybe it can be counteracted by applying a highpass and amplification that somehow track the 
  //   filter frequency (and become neutral at lower frequencies)...yes - that seems like a good 
  //   idea - try f_hp = rsMax(0, f - 2000)....but the appropriate amplification factor may depend
  //   on the input signal - for white noise, 1 may be appropriate, for brown noise, something
  //   proportional to frequency is more appropriate

  // Ideas:
  // -maybe use a bandpass and use max when the output is >= 0 and min otherwise
  //  -done - that seems promising for a tuned pseudo-resonance
  //  -the bandwidth can be adjusted by the user

  return;
}


void quantileFilter()
{
  //quantileFilterElongation();  // tests producing the length L+1 output by length L filter
  //quantileFilterSweep();  // tests non-integer length quatile filter
  //quantileFilterDelay();
  //quantileFilterDual();  // tests the dual-quantile filter (with highpass mode, etc.)
  quantileFilterResonant();
}

template<class T, int M>  // M: simd vector size
void simdFilter()
{
  // We try to simd-vectorize a filter by packing M samples into the input...that doesn't work yet

  int N = 5*M;               // number of samples, must be divisible by M
  float b0 = 2.f, b1 = 3.f;   // filter coeffs

  using simdVec = rsSimdVector<T, M>;
  //using simdVec = rsFloat32x4;

  using VecN    = std::vector<float>;
  VecN x = rsConvert(rsRandomIntVector(N, -9, 9, 0), float());
  VecN y(N);

  // filter the regular, scalar way:
  float x1 = 0.f;
  for(int n = 0; n < N; n++)
  {
    float x0 = x[n];
    y[n] = b0*x0 + b1*x1;
    x1 = x0;
  }

  // Vectorize x:
  std::vector<simdVec> vx(N/M), vy(N/M);
  for(int n = 0; n < N/M; n++) 
  {
    vx[n][0] = x[n+0*M];
    vx[n][1] = x[n+1*M];
    vx[n][2] = x[n+2*M];
    vx[n][3] = x[n+3*M];

    // Better:
    //vx[n].set(x[n+0*M], x[n+1*M], x[n+2*M], x[n+3*M]);
    // But this does not compile when RS_USE_SSE is not defined because ony in the explicit 
    // specialization we hae such a set function that takes 4 parameters. In the fallback version,
    // it will perhaps need to be implemented by means of variadic templates
  }

  // Filter vectorized x:
  simdVec vx1 = 0.f;
  for(int n = 0; n < N/M; n++)
  {
    simdVec vx0 = vx[n];
    vy[n] = b0*vx0 + b1*vx1;
    vx1 = vx0;
  }

  // Devectorize y:
  VecN z(N);
  for(int n = 0; n < N/M; n++)
  {
    z[M*n+0] = vy[n][0];
    z[M*n+1] = vy[n][1];
    z[M*n+2] = vy[n][2];
    z[M*n+3] = vy[n][3];
  }

  // ääähhhmmm...nope - this is not how it works!

  // -Maybe try it with an exponential decay filter
  // -The M decimated filters should have their time constants shortened by a factor of M. Maybe we
  //  should take the M-th power of the original recursion coeff.
  // -Try to generate the correct impulse-response first. We may need to compute M initial states 
  //  in a direct (non-vectorized) way
  // -To make the decimation impulse-invariant with respect to the original filter, we may need to
  //  average the input signal in some way...maybe it whould work like:
  //  MA -> exp-decay-filter -> extract
  // -Maybe we need special coeffs to feed the inputs into the M filters, like if an impulse 
  //  occured at tap k (k=0,..,M-1), the filter at tap i (i=0,..,M-1) should weight it the the 
  //  recurson coeff to the power of k-i or i-k or something. With "tap" i mean here the different
  //  phases of the filter (like in a polyphase filter)
  // -Consider the difference equation: y[n] = x[n] + c*y[n-1]. Assume that a new vector of M=4 
  //  input samples X[0],X[1],X[2],X[3] arrives and we want to find the coeffs by which each of 
  //  these values should go into our 4 partial filters. Let the coefficient a_ij denote the coeff
  //  by which X[i] goes into partial filter j. I think we need:
  //    a_00 = a_11 = a_22 = a33 = c^0 = 1    (no delay or advance)
  //    a_01 = a_12 = a_23       = c^1        (simulates delayed arrival by 1 sample)
  //    a_10 = a_21 = a_32       = c^-1       (simulates arrival 1 sample earlier)
  //    a_02 = a_13              = c^2
  //    a_20 = a_31              = c^-2
  //    a_03                     = c^3
  //    a_30                     = c^-3
  //  multiplication by these factors simulates the effect of the sample arriving with an 
  //  appropriate time delay (or time advance). Maybe we also need an additional factor of 1/M to 
  //  compensate for the summation. Maybe that matrix should be called something like decimation 
  //  matrix: it applies a matrix to M subsequent input samples to establish the input vector to 
  //  the filter that operates at a decimated sample rate. This matrix is specific to the filter. 




  rsPlotVectors(x, y);



  int dummy = 0;
}

template void simdFilter<float, 4>();
