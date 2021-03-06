using namespace RSLib;


//=================================================================================================

// x: input, y: output, N: length, processor: the filter/processor to be applied (needs to 
// implement double getSample(double), P: padding length, numPasses: number of forward/backward 
// passes maybe move to another file
template<class T>
void rsApplyBiDirectionally(double *x, double *y, int N, T &processor, int P, int numPasses)
{
  // create a buffer containing the signal with (pre- and post) zero padding to allow the 
  // processor/filter to ring out at the ends:
  int M = N+2*P;
  double *tmp = new double[M];
  rsFillWithZeros(tmp, P);
  rsCopyBuffer(x, &tmp[P], N);
  rsFillWithZeros(&tmp[P+N], P);

  // apply processor (multipass, bidirectionally):
  int p, n;
  for(p = 1; p <= numPasses; p++)
  {
    // forward pass:
    for(n = 0; n < M; n++)
      tmp[n] = processor.getSample(tmp[n]);

    // backward pass:
    for(n = M-1; n >= 0; n--)
      tmp[n] = processor.getSample(tmp[n]);
  }

  // copy result to output and clean up:
  rsCopyBuffer(&tmp[P], y, N);
  delete[] tmp;

  // todo: Actually, we don't really need pre- and post padding. Using just post padding and 
  // applying first all forward passes and then all backward passes should give the same result.
}

int rsBiDirectionalFilter::getPaddingLength(double bw, double fs)
{
  return rsCeilInt(10 * fs / bw); 
  // factor 10 is ad hoc - experiment to find optimal factor, maybe the formula should include 
  // the number of passes as well?
}

void rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(double *x, double *y, int N, double fc, 
  double bw, double fs, int numPasses, double gc)
{
  // compute desired bandwidth for single-pass filter in octaves:
  bw *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, 1, gc);  
  double bo = rsBandwidthConverter::absoluteBandwidthToOctaves(bw, fc);

  // create and set up the filter:
  rsStateVariableFilter flt;  // maybe use a biquad later
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setBandwidth(bo);
  flt.setMode(rsStateVariableFilter::BANDPASS_PEAK);

  // apply filter:
  int P = getPaddingLength(bw, fs);
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

void rsBiDirectionalFilter::applyButterworthBandpassBwInHz(double *x, double *y, int N, double fc, 
  double bw, double fs, int order, int numPasses, double gc)
{
  // compute desired bandwidth for single-pass filter in octaves:
  bw *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, order, gc);
  double bo = rsBandwidthConverter::absoluteBandwidthToOctaves(bw, fc);

  // create and set up the filter:
  rsEngineersFilter flt;
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setBandwidth(bo);
  flt.setMode(rsInfiniteImpulseResponseDesigner::BANDPASS);
  flt.setApproximationMethod(rsPrototypeDesigner::BUTTERWORTH);
  flt.setPrototypeOrder(order);

  // apply filter:
  int P = rsCeilInt(numPasses * flt.getRingingTimeEstimate(rsDB2amp(-100.0)));
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

void rsBiDirectionalFilter::applyButterworthLowpass(double *x, double *y, int N, double fc, 
  double fs, int order, int numPasses, double gc)
{
  // compute desired lowpass cutoff:
  fc *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, order, gc);

  // create and set up the filter:
  rsEngineersFilter flt;
  flt.setSampleRate(fs);
  flt.setFrequency(fc);
  flt.setMode(rsInfiniteImpulseResponseDesigner::LOWPASS);
  flt.setApproximationMethod(rsPrototypeDesigner::BUTTERWORTH);
  flt.setPrototypeOrder(order);

  // apply filter:
  int P = rsCeilInt(numPasses * flt.getRingingTimeEstimate(rsDB2amp(-100.0)));
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

void rsBiDirectionalFilter::applyButterworthHighpass(double *x, double *y, int N, double fc, 
  double fs, int order, int numPasses, double gc)
{
  applyButterworthLowpass(x, y, N, fc, fs, order, numPasses, gc);
  rsSubtract(x, y, y, N);


  //// quick and dirty code duplication from applyButterworthLowpass - refactor

  //// compute desired lowpass cutoff:
  //fc *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, order, gc);

  //// create and set up the filter:
  //rsEngineersFilter flt;
  //flt.setSampleRate(fs);
  //flt.setFrequency(fc);
  //flt.setMode(rsInfiniteImpulseResponseDesigner::HIGHPASS); // this line is the only difference to lowpass
  //flt.setApproximationMethod(rsPrototypeDesigner::BUTTERWORTH);
  //flt.setPrototypeOrder(order);

  //// apply filter:
  //int P = rsCeilInt(numPasses * flt.getRingingTimeEstimate(rsDB2amp(-100.0)));
  //rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}



void rsBiDirectionalFilter::applyLowpass(double *x, double *y, int N, double fc, double fs, 
  int numPasses, double gc)
{
  // compute desired cutoff for single-pass filter:
  fc *= rsBandwidthConverter::multipassScalerButterworth(2*numPasses, 1, gc);

  // create and set up the filter:
  rsOnePoleFilter flt;
  flt.setMode(rsOnePoleFilter::LOWPASS); // gives impulse invariant lowpass design - maybe switch
                                         // to bilinear later
  flt.setSampleRate(fs);
  flt.setCutoff(fc);

  // apply filter:
  int P = getPaddingLength(fc, fs);
  rsApplyBiDirectionally(x, y, N, flt, P, numPasses);
}

//=================================================================================================

bool rsZeroCrossingFinder::isUpwardCrossing(double *x, int n)
{
  if(x[n] < 0.0 && x[n+1] >= 0.0)
    return true;
  return false;
}

int rsZeroCrossingFinder::numUpwardCrossings(double *x, int N)
{
  int nz = 0;
  for(int n = 0; n < N-1; n++)
  {
    if( isUpwardCrossing(x, n) )
      nz++;
  }
  return nz;
}

std::vector<int> rsZeroCrossingFinder::upwardCrossingsInt(double *x, int N)
{
  int Nz = numUpwardCrossings(x, N);  // number of zero crossings in x
  int nz = 0;                         // index of the current zero crossing
  std::vector<int> z;
  z.resize(Nz);
  for(int n = 0; n < N-1; n++)
  {
    if(isUpwardCrossing(x, n))
    {
      z[nz] = n;
      nz++;
    }
  }
  return z;
}

std::vector<rsFractionalIndex> rsZeroCrossingFinder::upwardCrossingsIntFrac(double *x, 
  int N, int p)
{
  std::vector<int> zi = upwardCrossingsInt(x, N);
  int Nz = (int)zi.size();
  std::vector<rsFractionalIndex> z;   // the positions of the zero crossings
  z.resize(Nz);
  double *a = new double[2*p+2];      // polynomial coefficients for interpolant 
  double nf;                          // fractional part of zero-crossing sample-index
  int q;                              // p, restricted to avoid access violation
  int n;                              // integer part of zero-crossing sample-index 
  for(int nz = 0; nz < Nz; nz++)
  {
    n  = zi[nz];
    nf = x[n]/(x[n]-x[n+1]); // zero of linear interpolant
    q  = rsMin(p, n, N-n-2);
    if(q > 0)
    {
      // refine linear zero estimate by Newton iteration on a higher order interpolating 
      // polynomial using the zero of the linear interpolant as initial guess:
      rsInterpolatingPolynomial(a, -q, 1, &x[n-q], 2*q+2);
      nf = getRootNear(nf, a, 2*q+1, 0.0, 1.0);
    }
    z[nz].intPart  = n;    // store integer part
    z[nz].fracPart = nf;   // store fractional part
  }
  delete[] a;
  return z;
}

std::vector<double> rsZeroCrossingFinder::upwardCrossings(double *x, int N, int p)
{
  std::vector<rsFractionalIndex> z = upwardCrossingsIntFrac(x, N, p);
  std::vector<double> zd;
  zd.resize(z.size());
  for(unsigned int i = 0; i < z.size(); i++)
    zd[i] = z[i].intPart + z[i].fracPart; // combine integer and fractional part into a "double" value
  return zd;
}

std::vector<double> rsZeroCrossingFinder::bandpassedUpwardCrossings(double *x, int N, double fc, 
  double bw, double fs, int np, int p)
{
  double *y = new double[N];
  rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y, N, fc, bw, fs, np);
  std::vector<double> z = upwardCrossings(y, N, p);
  delete[] y;
  return z;
}

//=================================================================================================

// move to ArrayFunctions:
double rsMaxPosition(double *buffer, int N)
{
  int i = rsMaxIndex(buffer, N);

  if(i < 1 || i >= N-1)
    return i;

  double x[3] = { -1, 0, 1 };
  double a[3];
  fitQuadratic(a, x, &buffer[i-1]);

  double offset = 0;
  if( abs(2*a[2]) >= abs(a[1]) ) // abs(offset) shall be <= 1
    offset = -a[1] / (2*a[2]);   // maximum of the parabola (zero-crossing of its slope)

  return i; // for test
  //return i + offset;  // also test
  //return i - offset;
  // figure out first, if it has to be + or - (i think +), then check, if the added constant in 
  // the delta computation in line 302 should be 0.5, 1, 1.5 (or something else)
}

rsCycleMarkFinder::rsCycleMarkFinder(double sampleRate, double minFundamental, double maxFundamental)
{
  fs   = sampleRate;
  fMin = minFundamental;
  fMax = maxFundamental;
}

void rsCycleMarkFinder::refineCycleMarksByCorrelation(double *x, int N, std::vector<double>& cm, 
  double f0)
{
  double* y = new double[N];
  if(correlationHighpass > 0)
    rsBiDirectionalFilter::applyButterworthHighpass(x, y, N, f0*correlationHighpass, fs, 4, 1);
  else
    rsCopyBuffer(x, y, N);

  int maxLength = 5000; // preliminary - use something based on the maximum time-delta between the 
                        // cycle-marks in cm array

  int left = (int) cm[0];
  std::vector<double> cl(maxLength), cr(maxLength), corr(2*maxLength-1);
  for(unsigned int i = 1; i < cm.size(); i++)
  {
    int right      = (int) cm[i];
    int halfLength = rsFloorInt(correlationLength * 0.5 * (right-left));
    int length     = 2*halfLength;

    // prepare buffers for correlation computation:
    cl.resize(length); // maybe use raw arrays instead of vectors, avoid memory re-allocations
    cr.resize(length);
    corr.resize(2*length-1);
    rsCopySection(y, N, &cl[0],  left-halfLength, length);
    rsCopySection(y, N, &cr[0], right-halfLength, length);

    // use a matched filter to do the correlation, 
    // see here: https://en.wikipedia.org/wiki/Matched_filter
    rsReverse(&cl[0], length);                            // reversed left cycle is used as "impulse response"
    rsConvolve(&cr[0], length, &cl[0], length, &corr[0]); // OPTIMIZE: use FFT convolution
    double delta = rsMaxPosition(&corr[0], 2*length-1) - length + 1.5; // is +1.5 correct? was found by trial and error
    //delta = rsMaxPosition(&corr[0], 2*length-1) - length + 1.0; // test

    // overwrite array entry and update for next iteration:
    cm[i] = right + delta;
    left  = (int) cm[i];
  }

  delete[] y;
}

std::vector<double> rsCycleMarkFinder::findCycleMarks(double *x, int N)
{
  // Get initial estimate of fundamental by using an autocorrelation based algorithm at the center
  // of the input signal:
  double f0 = rsInstantaneousFundamentalEstimator::estimateFundamentalAt(x, N, N/2, fs, fMin, fMax);
    // here, we have a mutual dependency between rsInstantaneousFundamentalEstimator and 
    // rsCycleMarkFinder - maybe break up by dragging estimateFundamentalAt out of the class

  double bw = bandPassWidth*f0; // absolute bandwidth
  double *y = new double[N];
  rsBiDirectionalFilter::applyConstPeakBandpassBwInHz(x, y, N, f0, bw, fs, bandpassSteepness);

  // Find cycle marks:
  std::vector<double> z; 
  if(algo == F0_ZERO_CROSSINGS)
    z = rsZeroCrossingFinder::upwardCrossings(y, N, precision);
  else
  { // use CYCLE_CORRELATION as cycleMarkAlgo
    z = rsZeroCrossingFinder::upwardCrossings(y, N, 0); // initial estimates
    refineCycleMarksByCorrelation(x, N, z, f0);
  }
  // todo: maybe work with the more precise version that splits integer and fractional parts
  // of the zero-crossings

  delete[] y;
  return z;
}

//=================================================================================================

// converts raw value from the cosine-generator into the window value - todo: write a class
// rsWindowFunctionIterator (subclass of rsSineIterator) that wraps that:
RS_INLINE double cosineToWindow(double c)
{
  // exact Blackman coefficients, modified such that we don't need to precompute c2=2*c*c-1:
  const double a0 = 6508.0/18608.0;
  const double a1 = 9240.0/18608.0;
  const double a2 = 2860.0/18608.0;
  return a0 + a1*c + a2*c*c;
}
// rename to cosineToBlackman
// todo: maybe write a general function that recursively computes cosines of successive multiples 
// of an angle from the cosine of that angle. This recursion is based on the trigonometric 
// identity:
// cos(a)*cos(b) = (cos(a-b) + cos(a+b))/2 which implies
// cos(w)*cos(w) = (cos(0)   + cos(2*w))/2, so cos(2*w) = 2*cos^2(w) - cos(0)
// the general recursion formula is:
// c[0] = cos(0*w) = 1
// c[1] = cos(1*w) = c
// c[n] = cos(n*w) = a*c[n-1]-c[n-2] with a = 2*cos(w) = 2*c
// from linear combinations these successive cosine values, various windows can be created, 
// see https://en.wikipedia.org/wiki/Window_function#Higher-order_generalized_cosine_windows


RS_INLINE void sincInterpolatorLoop(int mMin, int mMax, double &tf, rsSineIterator &sinIt, 
  rsSineIterator &wndIt, double &y, double *&x, int &ti, double &ws)
{
  double w;
  for(int m = mMin; m <= mMax; m++)
  {
    w   = sinIt.getValue() * cosineToWindow(wndIt.getValue()) / (m-tf);
    ws += w;
    y  += w * x[ti+m];
  }
}
RS_INLINE void sincInterpolatorLoopNoStretch(int mMin, int mMax, double &tf, double &s,
  rsSineIterator &wndIt, double &y, double *&x, int &ti, double &ws)
{
  double w;
  for(int m = mMin; m <= mMax; m++)
  {
    w  = s * cosineToWindow(wndIt.getValue())  / (m-tf);
    ws += w;
    y  += w * x[ti+m];      
    s  *= -1.0;
  }
}
double rsResampler::signalValueViaSincAt(double *x, int N, double t, double sincLength, 
  double stretch)
{
  int L  = (int) floor(sincLength);
  int ti = (int) floor(t);         // integer part of t
  if( ti < 0 || ti >= N )
    return 0.0;
  double tf = t - ti;              // fractional part of t    
  double y  = 0.0;                 // output value
  double s;                        // sine value
  double ws = 0.0;                 // sum of tap weights
  int mMin = -rsMin(L/2-1, ti);    // minimum shift
  int mMax = +rsMin(L/2, N-ti-1);  // maximum shift

  // optimized loop for stretch == 1.0 (used for downward transpositions):
  rsSineIterator wndIt(2*PI/sincLength, 2*PI*(mMin-tf)/sincLength+PI/2);
  if( stretch == 1.0 )
  {
    s = sin(PI*(mMin-tf)/stretch) / PI;
    if( tf > EPS )
      sincInterpolatorLoopNoStretch(mMin, mMax, tf, s, wndIt, y, x, ti, ws);
    else
    {
      // split loop into a part for m < 0 and for m > 0, the case for m == 0 is handled separately 
      // because it would lead to division-by-zero error:
      sincInterpolatorLoopNoStretch(mMin, -1, tf, s, wndIt, y, x, ti, ws);
      wndIt.getValue(); // for internal state-update
      ws += 1.0;
      y  += x[ti];
      sincInterpolatorLoopNoStretch(1, mMax, tf, s, wndIt, y, x, ti, ws);
    }
    return y/ws;
  }

  // general case with stretch (could handle special case above also, but less efficiently):
  rsSineIterator sinIt(PI/stretch, PI*(mMin-tf)/stretch, 1.0/PI);
  if( tf > EPS )
    sincInterpolatorLoop(mMin, mMax, tf, sinIt, wndIt, y, x, ti, ws);
  else
  {
    sincInterpolatorLoop(mMin, -1, tf, sinIt, wndIt, y, x, ti, ws);
    sinIt.getValue();
    wndIt.getValue();
    ws += 1.0/stretch;
    y  += x[ti]/stretch;
    sincInterpolatorLoop(1, mMax, tf, sinIt, wndIt, y, x, ti, ws);
  }
  return y/ws;
}
// \todo: make a version of this function for periodic signals, maybe the loop should look like:
//  for(int m = sincLength/2; m <= sincLength/2; m++)
//    y += x[(ti+m)%N] * rsWindowedSinc(m-tf, sincLength, stretch); 
// ...but we need to check, if this is correct. This version could be used for looped 
// sample-playback (for example, in a sampler) or for wavetable-oscillators

void rsResampler::transposeLinear(double *x, int xN, double *y, int yN, double factor)
{
  int    nw;         // write position
  double nr = 0.0;   // read position
  int    nri;        // integer part of nr
  double nrf;        // fractional part of nr
  for(nw = 0; nw < yN; nw++)
  {
    nri = (int) floor(nr);

    if( nri+1 >= xN )
      break;  // end of input reached - todo: precompute the maximum nw and use it in the loop
              // condition to get rid of that conditional here

    nrf = nr - nri;
    y[nw] = (1.0-nrf)*x[nri] + nrf*x[nri+1];
    nr += factor;
  }
  rsFillWithZeros(&y[nw], yN-nw); // fill tail with zeros
}

void rsResampler::transposeSinc(double *x, int xN, double *y, int yN, double factor,                            
  double sincLength, bool antiAlias)
{
  double stretch = 1.0; 
  if( antiAlias == true )
    stretch = rsMax(1.0, factor);

  int    nw;         // write position
  double nr = 0.0;   // read position
  for(nw = 0; nw < yN; nw++)
  {
    y[nw] = signalValueViaSincAt(x, xN, nr, sincLength, stretch);
    nr += factor;
  }
  rsFillWithZeros(&y[nw], yN-nw);
}

void rsResampler::shiftSinc(double *x, double *y, int N, double amount, double sincLength)
{
  if( x == y )
  {
    double *tmp = new double[N];
    shiftSinc(x, tmp, N, amount, sincLength);
    rsCopyBuffer(tmp, y, N);
    delete[] tmp;
  }
  else
  {
    for(int n = 0; n < N; n++)
      y[n] = rsResampler::signalValueViaSincAt(x, N, n-amount, sincLength, 1.0);
  }
}

//=================================================================================================

void rsTimeWarper::timeWarpSinc(double *x, int xN, double *y, double *w, int yN,                             
  double minSincLength, double maxLengthScaler, bool antiAlias)
{
  if( antiAlias == false )
  {
    for(int n = 0; n < yN; n++)
      y[n] = rsResampler::signalValueViaSincAt(x, xN, w[n], minSincLength, 1.0);
  }
  else
  {
    double tOld = w[0];
    for(int n = 0; n < yN; n++)
    {
      double speed   = rsAbs(w[n]-tOld);
      double stretch = rsMax(1.0, speed);
      //double length  = minSincLength * rsMin(stretch, maxLengthScaler); // old, probably buggy
      double length  = rsMax(minSincLength, minSincLength*rsMin(stretch, maxLengthScaler));

      // maybe better, do just: 
      // length = rsLimitToRange(minSincLength*stretch, minSincLength, maxSincLength);
      // where maxSincLength is precomputed as minSincLength*maxLengthScaler

      y[n] = rsResampler::signalValueViaSincAt(x, xN, w[n], length, stretch);
      tOld = w[n];
    }
  }
}

void rsTimeWarper::invertMonotonousWarpMap(double *w, int N, double *wi)
{
  //bool cubic = true; // make this a parameter for the function

  int M = (int) ceil(w[N-1]); // length of inverse map 
  int n, i, iOld = 0;
  for(n = 0; n < M; n++)
  {
    // find index i, such that ty[i] <= n and ty[i+1] >= n:
    for(i = iOld; i < N-1; i++)
    {
      if( w[i] <= n && w[i+1] >= n )
        break;
    }
    wi[n] = rsInterpolateLinear(w[i], w[i+1], (double)i, (double)(i+1), (double)n);


    /*
    if( cubic == true && i >= 1 && i < N-2 )
    {
      wi[n] = rsInterpolateCubicHermite(w[i-1], w[i], w[i+1], w[i+2], (double)(i-1), (double)i, 
        (double)(i+1), (double)(i+2), (double)n);
    }
    else
      wi[n] = rsInterpolateLinear(w[i], w[i+1], (double)i, (double)(i+1), (double)n);
    */
    //wi[n] = rsInterpolateLinear(w[i], (double)i, w[i+1], (double)(i+1), (double)n);
      // Can be optimized: (i+1)-i == 1 - always, the function computes the value, but maybe it's 
      // not worth it. On the other hand, maybe the quality could be improved by using cubic 
      // hermite instead of linear interpolation.

    iOld = i;
  }
  // Remark: I have tried to use cubic hermite interpolation instead of linear, but it doesn't make
  // a difference quality wise (the spectral artifacts do not decrease measurably), so linear seems
  // to be appropriate. Moreover, with linear interpolation a forward/backward mapping roundtrip of
  // any time-value should be an indentity operation (up to roundoff error) whcih seems a desirable
  // property.
}

 
int rsTimeWarper::getPitchModulatedLength(double *r, int N)
{
  double s = 0.0;
  for(int n = 1; n < N; n++)
    s += 1.0 / r[n-1];
  return (int) ceil(s);
  // Remark: the last readout speed value r[N-1] is irrelevant because it would formally determine
  // the time-instant, where the (nonexistent) next-to-last sample x[N] would have to be written
}

void rsTimeWarper::applyPitchModulation(double *x, double *r, int N, double *y,
  double minSincLength, double maxLengthScaler, bool antiAlias)
{
  rsVariableSpeedPlayer vsp;
  vsp.setInputAndSpeed(x, r, N);
  vsp.getOutput(y, minSincLength, maxLengthScaler, antiAlias);
}

//=================================================================================================

rsVariableSpeedPlayer::rsVariableSpeedPlayer()
{
  init();
}

rsVariableSpeedPlayer::~rsVariableSpeedPlayer()
{
  clear();
}

void rsVariableSpeedPlayer::setInputAndSpeed(double *input, double *r, int length)
{
  clear();

  x  = input;
  Nx = length;

  // Create inverse warping map (i.e., the map, that assigns for each input sample-index n in x a 
  // (noninteger) time-instant t in y where the sample value x[n] should be written, such that
  // y(t) = x[n]:
  wi = new double[Nx];
  wi[0] = 0.0; // y(0.0) = x[0]
  for(int n = 1; n < Nx; n++)
    wi[n] = wi[n-1] + 1.0 / r[n-1];

  // Obtain desired warping map by inverting the wi map:
  Ny = (int) ceil(wi[Nx-1]);
  w  = new double[Ny]; 
  rsTimeWarper::invertMonotonousWarpMap(wi, Nx, w);
}

double rsVariableSpeedPlayer::warpTime(double tx)
{
  return rsInterpolateClamped(wi, Nx, tx);
}

double rsVariableSpeedPlayer::unwarpTime(double ty)
{
  return rsInterpolateClamped(w, Ny, ty);
}

void rsVariableSpeedPlayer::getOutput(double *y, double minSincLength, double maxLengthScaler,
  bool antiAlias)
{
  rsTimeWarper::timeWarpSinc(x, Nx, y, w, Ny, minSincLength, maxLengthScaler, antiAlias);
}

vector<double> rsVariableSpeedPlayer::getOutput(double minSincLength, double maxLengthScaler, 
  bool antiAlias)
{
  vector<double> y(Ny);
  getOutput(&y[0], minSincLength, maxLengthScaler, antiAlias);
  return y;
}

vector<double> rsVariableSpeedPlayer::getTimeWarpMapXY()
{
  vector<double> map(Nx);
  rsCopyBuffer(wi, &map[0], Nx);
  return map;
}

vector<double> rsVariableSpeedPlayer::getTimeWarpMapYX()
{
  vector<double> map(Ny);
  rsCopyBuffer(w, &map[0], Ny);
  return map;
}

vector<double> rsVariableSpeedPlayer::invertSpeeds(vector<double>& speeds)
{
  rsVariableSpeedPlayer vsp;
  vsp.setInputAndSpeed(nullptr, &speeds[0], (int)speeds.size());
  vector<double> map = vsp.getTimeWarpMapYX();
  int N = (int)map.size();
  for(int n = 0; n < N-1; n++)
    map[n] = 1.0 / (map[n+1]-map[n]);
  map[N-1] = 1;   // value irrelevant (not used in map computation)
  return map;
}

vector<double> rsVariableSpeedPlayer::applyPlaybackSpeed(vector<double>& x, vector<double>& s)
{
  rsAssert(x.size() == s.size());
  rsVariableSpeedPlayer vsp;
  vsp.setInputAndSpeed(&x[0], &s[0], (int)x.size());
  return vsp.getOutput();
}

void rsVariableSpeedPlayer::init()
{
  x = w = wi = nullptr;
  Nx = Ny = 0;
}

void rsVariableSpeedPlayer::clear()
{
  delete[] wi;
  delete[] w;
  w = wi = nullptr;
}

//=================================================================================================

void rsPitchFlattener::setInput(double *x, double *f, int N, double ft)
{
  // x: input, f: instantaneous frequencies, N: length of x and f, ft: target frequency

  // set default target frequency, if 0 is passed:
  if(ft == 0.0)
    ft = rsMean(f, N);

  // create temporary read-speed array:
  double *r = new double[N];
  for(int n = 0; n < N; n++)
    r[n] = ft / f[n];

  // setup inherited variable speed player and clean up temporary array:
  rsVariableSpeedPlayer::setInputAndSpeed(x, r, N);
  delete[] r;
}

//=================================================================================================

void rsPhaseLockedCrossfader::setInputs(double *in1, double *f1, int len1, double *in2,
  double *f2, int len2, double ft)
{
  x1 = in1;
  x2 = in2;
  N1 = len1;
  N2 = len2;
  if(ft == 0.0)
    ft = 0.5 * (rsMean(f1, N1) + rsMean(f2, N2));
  pf1.setInput(x1, f1, N1, ft);
  pf2.setInput(x2, f2, N2, ft);
}

void rsPhaseLockedCrossfader::setFlattenedCrossfade(double start, double end, double shift)
{
  // adjust crossfade-start such that it is an integer with respect to x1 and crossfade-end to
  // be an integer with respect to x2 for seamless splicing:
  cs1   = (int)round(pf1.unwarpTime(start));
  start = pf1.warpTime(cs1);
  cs2   = pf2.unwarpTime(start-shift);
  ce2   = (int)round(pf2.unwarpTime(end-shift));
  end   = pf2.warpTime(ce2) + shift;
  ce1   = pf1.unwarpTime(end);
  this->shift = shift;
}

int rsPhaseLockedCrossfader::getCrossfadeOutputLength()
{
  double cl1 = ce1 - cs1;              // crossfade length in x1
  double cl2 = ce2 - cs2;              // crossfade length in x2
  double cly;                          // crossfade length in y 
  //cly = 0.5*(cl1+cl2);               // arithmetic mean
  cly = 1 / (0.5*(1/cl1 + 1/cl2));     // harmonic mean
  return (int)round(cly);

  // Using arithmetic mean, the vibrato speed during crossfade will be between the vibrato speeds
  // of x1, x2. Using the harmonic mean, the output frequency will be the arithmetic mean between 
  // the input frequencies of x1, x2. The geometric mean is in between. Maybe we can use a 
  // generalized mean and expose the exponent as user parameter.
  // https://en.wikipedia.org/wiki/Generalized_mean
}

vector<double> rsPhaseLockedCrossfader::getFlattenedSignal1()
{
  return pf1.getOutput();
}
vector<double> rsPhaseLockedCrossfader::getFlattenedSignal2()
{
  return pf2.getOutput();
}

vector<double> rsPhaseLockedCrossfader::getOutput()
{
  // splice together heading section of x1, crossfade section and trailing section of x2:
  vector<double> yc = getCrossfadeSection();
  int L  = (int)yc.size();
  int s2 = ce2+1;             // start in x2
  int N  = cs1 + L + N2-s2;   // output length
  vector<double> y(N);
  int n;
  for(n = 0; n < cs1; n++)
    y[n] = x1[n];
  for(n = cs1; n < cs1+L; n++)
    y[n] = yc[n-cs1];
  for(n = cs1+L; n < N; n++)
    y[n] = x2[n-cs1-L+s2];
  return y;
}

vector<double> rsPhaseLockedCrossfader::getTimeWarpMapXY1()
{
  return pf1.getTimeWarpMapXY();
}
vector<double> rsPhaseLockedCrossfader::getTimeWarpMapYX1()
{
  return pf1.getTimeWarpMapYX();
}
vector<double> rsPhaseLockedCrossfader::getTimeWarpMapXY2()
{
  return pf2.getTimeWarpMapXY();
}
vector<double> rsPhaseLockedCrossfader::getTimeWarpMapYX2()
{
  return pf2.getTimeWarpMapYX();
}

vector<double> rsPhaseLockedCrossfader::getCrossfadeSection()
{
  computeReadoutTimes();
  int    L     = (int)t1.size();
  double scale = 1.0 / (L-1);
  double y1, y2, c;
  vector<double> yc(L);
  for(int n = 0; n < L; n++)
  {
    c  = scale * n;
    y1 = rsResampler::signalValueViaSincAt(x1, N1, t1[n], 64, 1.0);
    y2 = rsResampler::signalValueViaSincAt(x2, N2, t2[n], 64, 1.0);
    yc[n] = (1-c)*y1 + c*y2;
  }
  return yc;
  // maybe later include optional anti-aliasing
}

void rsPhaseLockedCrossfader::computeReadoutTimes()
{
  double cl1 = ce1 - cs1;                        // crossfade length in x1
  double cl2 = ce2 - cs2;                        // crossfade length in x2
  int    L   = getCrossfadeOutputLength();       // crossfade length in output
  double scl = 1.0 / (L-1);                      // scaler to map 0..L-1 to 0..1
  double t, c;                                   // t: normalized time 0..1, c(t)
  double a = 1, b = 0;                           // parameters of c(t) = a*t + b*t^2
  bool sweep = true;                             // maybe make this a member variable
  if(sweep == true)
  {
    double s0 = L/cl1;                           // slope of c(t) at t=0
    double s1 = L/cl2;                           // slope of c(t) at t=1
    a = s0;
    b = 0.5*(s1-a);
  }
  t1.resize(L);
  t2.resize(L);
  for(int n = 0; n < L; n++)
  {
    t = scl * n;
    c = a*t + b*t*t;
    double tx1 = pf1.warpTime(cs1+c*cl1);
    double tx2 = pf2.warpTime(cs2+c*cl2)+shift;
    double txw = (1-t)*tx1 + t*tx2;              // aggree on warped readout time
    t1[n] = pf1.unwarpTime(txw);
    t2[n] = pf2.unwarpTime(txw-shift);
  }

  // Notes: 
  // The boolean "sweep" variable could at some point be exposed as user option. If sweep mode
  // is off and the two input signals have different fixed frequencies, we will see a fixed output
  // frequency during the crossfade which is intermediate between the two input frequencies. Which 
  // frequency exactly that is, is determined by the return value of getCrossfadeOutputLength - if 
  // it returns the harmonic mean between both input crossfade-lengths, the frequency during 
  // crossfade will be the arithmetic mean between both input frequencies (and vice versa).
  // The computation of txw as weighted mean between tx1 and tx2: txw = (1-t)*tx1 + t*tx2 could 
  // use also "c" instead of "t" as weight - or some combination of t and c and we could also try 
  // weighted means different from the arithmetic mean. All of these modifications will lead to 
  // different behavior when both input signals have vibratos with different speeds.
}

//=================================================================================================

void rsInstantaneousFundamentalEstimator::estimateReliability(double *x, int N, 
  const std::vector<double>& z, double *r)
{
  int Nz = (int) z.size(); // number of zero crossings in z
  int nz;                  // index of zero-crossing
  int sl, sr, el, er;      // start and end, left and right
  double c;                // cross-correlation value

  sr = rsCeilInt( z[0]);
  er = rsFloorInt(z[1]);
  for(nz = 2; nz < Nz; nz++)
  {
    sl = sr;
    el = er;
    sr = el + 1;
    er = (int) z[nz];

    double test = z[nz]; // for debug
      // seems like the last zero crossing is too close to the end of the signal such that the
      // zero crossing detection accesses invalid array indices - or something

    //c  = rsCrossCorrelation(&x[sl], el-sl+1, &x[sr], er-sr+1);
    c  = rsStretchedCrossCorrelation(&x[sl], el-sl+1, &x[sr], er-sr+1);
    rsFillWithValue(&r[sl], el-sl+1, c);
  }

  // extend first and last reliability value to the start and end of the r-array:
  rsFillWithValue(&r[el], N-el, c);
  sl = rsCeilInt(z[0]);
  rsFillWithValue(r, sl, r[sl]);
}

void rsInstantaneousFundamentalEstimator::measureInstantaneousFundamental(double *x, double *f, 
  int N, double fs, double fMin, double fMax, double *r, int cycleMarkAlgo)
{
  rsCycleMarkFinder cmf(fs, fMin, fMax); // todo: maybe set it up - or maybe have it a member and allow client code to set it up
  std::vector<double> z = cmf.findCycleMarks(x, N);

  // old:
  //std::vector<double> z = rsCycleMarkFinder::findCycleMarks(x, N, fs, fMin, fMax, cycleMarkAlgo, 3);
    // 3 is the precision - may be tweaked

  // Measure periods. They are given by the distance between successive zero-crossings. The time 
  // instant at which we consider this period measurement to be effective is midway between the 
  // two zero-crossings:
  int Nz = (int) z.size();
  double *p  = new double[Nz-1];
  double *tp = new double[Nz-1];
  int n;
  for(n = 0; n < Nz-1; n++)
  {
    p[n]  = z[n+1]-z[n];
    tp[n] = 0.5*(z[n+1]+z[n]);
  }

  // Interpolate period measurements up to samplerate and convert to frequencies
  double *tn = new double[N];
  rsFillWithIndex(tn, N);
  rsInterpolateSpline(tp, p, Nz-1, tn, f, N, 1);
  for(n = 0; n < N; n++)
    f[n] = fs/f[n];

  // maybe a better way of interpolation would be to compute an instantaneous phase for each
  // zero-crossing/pitch-mark (each would have a value of phase[n] = 2*pi*n), interpolate these
  // phases, and obtain the instantaneous frequencies as difference - this would ensure a "phase
  // coherent interpolation" of the instantaneous frequencies (i.e. the instantaneous phase for a
  // resynthesized signal will not tend to drift depending on the specifics of the interpolation 
  // scheme)
  // we should not really extrapolate the head and tail values with the first and last polynomial
  // repeating the last value will probably be better

  // If more actual data values are needed before the interpolation, we could find the zeros of
  // the derivatives which would give values halfway in between our current values. If the envelope
  // is exponential, the distance between maxima also equals the period (see Physics of Musical 
  // Instruments, 2nd Ed, p.12)

  // assessment of reliability, if desired:
  if( r != nullptr )
    estimateReliability(x, N, z, r);

  // Cleanup:
  delete[] p;
  delete[] tp;
  delete[] tn;
}

double rsInstantaneousFundamentalEstimator::estimateFundamentalAt(double *x, int N, int n, 
  double fs, double fMin, double fMax)
{
  double pMax = fs/fMin; // maximum detectable period

  // Compute desired Length of autocorrelation sequence. It should be significantly longer than the
  // lag that corresponds to our maximum detectable period because the accuracy of the measured 
  // values degrades for lags that are close to the end of the sequence due to less samples used in 
  // the averaging. We need it to be odd so we can use a chunk that is centered around n.
  double k = 2.0;
  int L = (int) ceil(k*pMax);  
  if( rsIsEven(L) )
    L += 1;

  // Get a chunk from the signal, typically centered around n, but the center is shifted at the 
  // start and end of the input signal, such that we don't need to zero-extend the input 
  // (conceptually):
  double *y = new double[L];            // chunk from input signal
  int ns = rsMin(rsMax(0, n-L/2), N-L); // start of the chunk
  rsCopySection(x, N, y, ns, L); 

  // measure frequency:
  double r;  // reliability
  double f = rsAutoCorrelationPitchDetector::estimateFundamental(y, L, fs, fMin, fMax, &r);

  // cleanup and return:
  delete[] y;
  return f;
}

//=================================================================================================

void RSLib::rsSineQuadraturePart(double *x, double *y, int N, double f, double fs, bool backward)
{
  rsOnePoleFilter flt;
  flt.setSampleRate(fs);
  flt.setCutoff(f);
  flt.setMode(rsOnePoleFilter::ALLPASS);
  if( backward == true )
  {
    for(int n = N-1; n >= 0; n--)
      y[n] = flt.getSample(x[n]);
  }
  else
  {
    for(int n = 0; n < N; n++)
      y[n] = flt.getSample(x[n]);
  }
}

// used internally by rsSineEnvelopeViaQuadrature, rsSineEnvelopeViaAmpFormula
void rsSmoothSineEnvelope(double *y, int N, double f, double fs, double s)
{
  double a[3], b[3]; // filter coeffs
  if( s > 0.0 && f/s < 0.5*fs )
  {
    rsBiquadDesigner::calculateFirstOrderLowpassCoeffs(b[0], b[1], b[2], a[1], a[2], 1.0/fs, f/s);
    rsNegate(a, a, 3);
    a[0] = 1.0;
    rsFilter(y, N, y, N, b, 1, a, 1);
    rsReverse(y, N);
    rsFilter(y, N, y, N, b, 1, a, 1);
    rsReverse(y, N);
     // do this with a 1st order filter as well, maybe use rsBiDirectionalFilter (lowpass version)
  }
}

void RSLib::rsSineEnvelopeViaQuadrature(double *x, double *y, int N, double f, double fs, 
  double s)
{
  // get 90� phase shifted version:
  rsSineQuadraturePart(x, y, N, f, fs, true);

  // obtain magnitude of complex phasor:
  for(int n = 0; n < N; n++)
    y[n] = sqrt(x[n]*x[n] + y[n]*y[n]);

  // optional smoothing:
  rsSmoothSineEnvelope(y, N, f, fs, s);
}

void RSLib::rsSineEnvelopeViaAmpFormula(double *x, double *y, int N, double f, double fs, 
  double s)
{
  double w = 2*PI*f/fs;
  double a, p;
  for(int n = 0; n < N-1; n++)
  {
    rsSineAmplitudeAndPhase(x[n], x[n+1], w, &a, &p);
    y[n] = a;
  }

  // for the last value, we look backwards:
  rsSineAmplitudeAndPhase(x[N-1], x[N-2], w, &a, &p);
  y[N-1] = a;

  // optional smoothing:
  rsSmoothSineEnvelope(y, N, f, fs, s);
}

void RSLib::rsEnvelopedSine(double *y, int N, double f, double fs, double p, double *a)
{
  rsSineIterator sinIt(2*PI*f/fs, p, 1.0);
  for(int n = 0; n < N; n++)
    y[n] = a[n] * sinIt.getValue();
}

double rsUnwrappedPhaseForCatchSweep(double p0, double pk, double wk, int k, 
  int sweepDirection = 0)
{
  // compute phase that would be reached at k, if we were running a fixed frequency sine at w from 
  // 0 to k:
  double pkw = p0 + k*wk;

  // find the unwrapped target phase pk at k - this is given by pk plus/minus some integer multiple
  // of 2*PI, such that we get most closely to pkw:
  double pku = pk;
  int n = 0;
  while( pku < pkw )
  {
    pku = pk + n*2*PI;
    n++;
  }
  if( sweepDirection > 0 || (sweepDirection == 0 && pku-pkw > PI) )
    pku -= 2*PI;

  return pku;
}

void RSLib::rsEnvelopedPhaseCatchSweep(double *y, int k, double p0, double pk, double wk, 
  double *a, int sweepDirection)
{
  // Notation: p0, pk, w0, wk: instantaneous phases and normalized radian frequencies at sample 0
  // and sample k, respectively, wn: inst. freq. at sample n = 0...k-1. We have
  // pk = p0 + sum_{n=0}^{k-1} wn and we assume a linear frequency sweep such that 
  // wn = w0 + n * dw, where dw = (wk-w0)/k. Inserting this into the pk formula, we get:
  // pk = p0 + sum_{n=0}^{k-1} (w0 + n*dw) = p0 + k*w0 + sn*dw, where sn = sum_{n=0}^{k-1} n
  // This can be solved for w0 = (pk - p0 - sn*wk/k) / (k - sn/k)

  pk = rsUnwrappedPhaseForCatchSweep(p0, pk, wk, k, sweepDirection);
  double snk = (double) rsSum(0, k-1) / (double) k;
  double w0  = (pk-p0-snk*wk)/(k-snk); // normalized radian frequency at start
  double dw  = (wk-w0)/k;              // increment for normalized radian frequency
  double wn  = w0;                     // instantaneous normalized radian frequency
  double pn  = p0;                     // instantaneous phase
  for(int n = 0; n < k; n++)
  {
    y[n] = a[n] * sin(pn);
    wn   = w0 + n*dw;
    pn  += wn;
  }

  // todo: factor out a function to compute an array of w[n] - mainly for plotting the frequency 
  // sweep and comparison between linear and parabolic sweep (to be implemented)
}

void RSLib::rsPhaseCatchSweep(double *y, int k, double p0, double pk, double wk, 
  int sweepDirection)
{
  rsFillWithValue(y, k, 1.0);
  rsEnvelopedPhaseCatchSweep(y, k, p0, pk, wk, y, sweepDirection);
}

void RSLib::rsEnvelopedPhaseCatchSine(double *y, int N, double f, double fs, double p0, 
  double pk, int k, double *a, int sweepDirection)
{
  // 1st part (sweeping frequency):
  double w = 2*PI*f/fs;
  rsEnvelopedPhaseCatchSweep(y, k, p0, pk, w, a, sweepDirection);

  // 2nd part (fixed frequency):
  rsEnvelopedSine(&y[k], N-k, f, fs, pk, &a[k]);
}

void RSLib::rsRecreateSine(double *x, double *y, int N, double fx, double fy, double fs, double p0, 
  double smooth)
{
  rsSineEnvelopeViaAmpFormula(x, y, N, fx, fs, smooth);
  rsEnvelopedSine(y, N, fy, fs, p0, y);
}

void RSLib::rsRecreateSineWithPhaseCatch(double *x, double *y, int N, double fx, double fy, 
  double fs, double p0, double pk, int k, double smooth, int sd)
{
  rsSineEnvelopeViaAmpFormula(x, y, N, fx, fs, smooth);
  rsEnvelopedPhaseCatchSine(y, N, fy, fs, p0, pk, k, y, sd);
}

/*
void RSLib::rsRecreateSine(double *x, double *y, int N, double fx, double fy, double fs, double p0, 
  double smooth)
{
  // get the input sinusoid's envelope (we may use y as temporary buffer for this)
  rsSineEnvelopeViaAmpFormula(x, y, N, fx, fs, smooth);

  // generate the sine at the desired frequency and multiply it by the envelope:
  double w = 2*PI*fy/fs;
  double p = p0;
  for(int n = 0; n < N; n++)
  {
    y[n] *= sin(p);
    p    += w;
  }

  // optimize using rsSineOscillator

  // factor out the sine creation part - should take an amplitude array as input
}
*/