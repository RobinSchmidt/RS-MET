
// a little helper function that is supposed to reorder a buffer (of delayline-lengths) that is
// initially in ascending order in a way that maximizes the average crossfeedback between
// delaylines with most different lengths (when the diffusion parameter is neither 0 nor 1)
// we implement it as template function outside the class because we need to apply it to the vector
// of output-gains as well
// maybe try the same with a buffer of initially descending lengths - this promotes stroger
// crossfeedback between the longer delaylines (wheareas we currently have stronger crossfeedback
// between shorter delaylines - hmm but that could actually be desirable - higher modes denser)
template<class T>
void shuffleBuffer(T *buffer, int length)
{
  int k  = length;
  int s0 = k/2;
  int s;
  while( s0 > 0 )
  {
    s = s0;
    while( s < length )
    {
      RAPT::rsArrayTools::reverse(&buffer[s], k/2);
      s += k;
    }
    s0 /= 2;
    k  /= 2;
  }
}
// ToDo: 
// -document, how this algo works..it seems to reverse the upper half of the buffer and then 
//  recursively applies the same algo to both half-buffers?
// -make it a member function
// -maybe try other shufflings:
//  -swap first with second half, then do the same recursively on both half-buffers
//   -maybe combine this with reversal
//  -swap highest quarter with 2nd-to-lowest, recurse to half-buffers


// Construction/Destruction:

FeedbackDelayNetwork::FeedbackDelayNetwork()
{
  delayLines         = nullptr;
  readIndices        = nullptr;
  writeIndices       = nullptr;
  delayLineLengths   = nullptr;
  relativeDelayTimes = nullptr;

  //numDelayLines = 4;   // 4 is the minimum, sounds kinda shattery, gritty - not so nice
  //numDelayLines = 8;   // 8 has already a nice reverberant character
  //numDelayLines   = 16;  // 16 sounds smoother, denser, noisier
  //numDelayLines = 32;  // 32 sounds still denser
  //numDelayLines = 64;  // still denser, initial section acquires a snare'ish character
  numDelayLines = 128; // initial "electronic snare" character becomes very apparent
  //numDelayLines = 256; // snare sound becomes higher pitched
  // ...conclusion: 
  // -16 or 32 seems to be the sweet spot for reverb, 8 and 64 are also usable
  // -larger sizes acquire this weird pitch-sweepdown electronic drum sound
  //  -maybe that can be avoided by using only a subset of the delaylines for output
  //  -> experiment with output vectors containing zeros, maybe an output-density parameter
  //  can be defined
  // -maybe try to use not one N=64 network but a parallel connection of 4 N=16 networks with
  //  different delay-sttings - maybe that creates more density and also avoids the pitch-sweepdown
  // -try different formulas for the relative delay times, currently, apower rule is used

  sampleRate         = 44100.0;
  referenceDelayTime = 20.0/1000.0;  // 20 ms
  allocateMemory();
  setDiffusion(100.0);
  //setDiffusion(0.0);
  reset();
}

FeedbackDelayNetwork::~FeedbackDelayNetwork()
{
  freeMemory();
}

// Setup:

void FeedbackDelayNetwork::setDiffusion(double newDiffusion)
{
  diffusion = newDiffusion;

  // hmm...seems like we should have a nonlinear mapping - we need more precision in the lower
  // range

  double d = 0.01*diffusion;
  //d = d*d*d;  // nah - that's not yet a good mapping

  //d = mapLinearToRational(d, +0.7);

  //d = d*d;


  /*
  // found empirically to provide a perceptually uniform mapping - still between 75% and 100%
  // there's too little difference - maybe use a cubic spline and enforce steeper slope at
  // (x,y)=(1,1) - in any event the spline should go through (x,y)=(0.5, 0.25) and provide higher
  // resolution around that point (slope should be shallow there)
  if( d < 0.5 )
    d =         0.875*d - 0.75 * d*d;
  else
    d = 0.875 - 2.625*d + 2.75 * d*d;
  */

  // i think, this required mapping comes about because the percived diffusion depends on the ratio
  // a/b which has singularities quickly towards the ends of of the interval ...maybe some tan or
  // atan based formula could help?

  RAPT::rsSinCos(d*PI/4, &b, &a);

  // \todo - when introducing modulation for d (or the angle phi), check if it is better to
  // modulate the angle before or after the non linear mapping. if the modulation is applied before
  // the mapping, we should wrap the modulated value into 0...1 before the sin/cos
}


// Others:

void FeedbackDelayNetwork::reset()
{
  memset(writeIndices, 0, numDelayLines*sizeof(int));
  for(int i = 0; i < numDelayLines; i++)
    memset(delayLines[i], 0, delayLineLengths[i]*sizeof(double));
  setupReadIndices();
}

// code is obsolete - we have better code in rapt Transforms now - move to prototypes for reference
// ...but before doing so, do performance tests - it may be more efficient, even though it doesn't
// work in place
void FeedbackDelayNetwork::fastGeneralizedHadamardTransform(
  double *x, int N, int log2N, double *y, double a, double b, double c, double d)
{
  int  N2 = N/2;
  bool xContainsResult = true;
  for(int i = 1; i <= log2N; i++)
  {
    for(int j = 0; j < N2; j++)
    {
      y[j]    = a*x[2*j] + b*x[2*j+1];
      y[j+N2] = c*x[2*j] + d*x[2*j+1];
      // maybe precompute 2*j, 2*j+1 and measure performance - but the compiler should be able
      // to do that optimization itself
      // maybe this computation lends itself well to SSE2-optimization (do the 4 multiplications at
      // once)
    }
    RAPT::rsSwap(x, y);
    xContainsResult = !xContainsResult;
  }
  if( !xContainsResult )
    memcpy(y, x, N*sizeof(double));  // try rsArrayTools::copy instead, measure performance
}
// After the algo, x contains garbage - that's not nice
// todo: try to write the algorithm in a way that works in place - maybe it's sufficient to pull 
// out the x[2*j], x[2*j+1] into temp variables? ...but i actually don't think so -> figure out!
// ..i think it would work, if we also write y into slots 2*j, 2*j+1 - but then we would end up
// with a differently ordered output - but that may actually not be an issue in this context, but 
// it could also be counteracted by a subsequent reordering algo, like in the FFT with the 
// bit-reversed ordering...maybe it's even the same re-ordering here?
// maybe the algo should be moved to some appropriate place in rapt, too

// Here is a python implemenation of the FWHT in place:
// https://en.wikipedia.org/wiki/Fast_Walsh%E2%80%93Hadamard_transform
// ToDo: 
// -implement the algo from here and from wikipedia in rapt
// -make a unit test for both by comparing the results with explicit matrix multiplication of 
//  matrices created by the Sylvester construction (use the Kronecker-product in rsMatrix for
//  this)
// -there is alredy a unit test in testFastGeneralizedHadamardTransform -> extend this


void FeedbackDelayNetwork::fastInverseGeneralizedHadamardTransform(
  double *x, int N, int log2N, double *work, double a, double b, double c, double d)
{
  double s = 1.0 / (a*d - b*c);
  fastGeneralizedHadamardTransform(x, N, log2N, work, s*d, -s*b, -s*c, s*a);
}
// also obsolete now

void FeedbackDelayNetwork::setupRelativeDelayTimes()
{
  double dMax = 2.4;

  //dMax = 5;        // maybe vary this between 1.5...2.5 via a user parameter "Shape"

  //dMax = sqrt((double) numDelayLines + 1);
  //dMax = pow((double) numDelayLines, 2.0/3.0);

  //dMax = (double) numDelayLines / 5.0;

  dMax = (double) numDelayLines / 5.0;
  //dMax = 5;
  //dMax = 1.5;
  //dMax = 2;
  //dMax = 3;
  //dMax = 4;
  //dMax = 5;
  dMax = 8;
  // dMax = 1.5: grouping of primary reflections impulses
  // dMax = 2.0: most compact, i guess
  // dMax = 3.0: seems a goo all around value
  // dMax = 5.0: some later primary reflections are already bathed in noise
  // dMax = 8.0: primary reflectiosn becoem distinguishable (with N=16)

  // Power rule:
  for(int i = 0; i < numDelayLines; i++)
  {
    double D = dMax;
    double p = (double)i / (double)(numDelayLines-1);

    //relativeDelayTimes[i] = 1.0 + p * (D - 1.0);   // linear rule
    relativeDelayTimes[i] = pow(D, p);               // power rule
    //relativeDelayTimes[i] = 1.0 + sqrt(p) * (D - 1.0);   // square-root rule
    //relativeDelayTimes[i] = 1.0 + cbrt(p) * (D - 1.0);   // cube-root rule
    //relativeDelayTimes[i] = 1.0 + p*p * (D - 1.0);        // square rule
    //relativeDelayTimes[i] = 1.0 + p*p*p * (D - 1.0);        // cube rule
  }
  // with N = 128, dMax = 3
  // -Power:  has a clear pitch slide-down
  // -Linear: has a tonal component at the start
  // -rules of the form: 1 + p^k * (D-1) have a slide-down effekt for k > 1 and slide-up for
  //  k < 1 and are tonal for k = 1. The p^k shape is the key here - maybe try a shape in a form
  //  of a saturation curve or a sinh-like curve..should be neither up-down or down-up
  // -this can actually be used to synthesize (snare) drum sounds, maybe the initial section
  //  should be lowpassed a bit
  //  -to avoid the spikey character at the start, maybe try other input signals like rectangular 
  //   impulses, triangular, bipolar, etc - maybe make an impulse-generator class that contains the 
  //   unit-impulse as special case
  // -with larger dMax, the pitch slides further down and later impulses are more bathed in noise,
  //  but when it gets too large (8), it begins to lose the pitch-slide effekt and becomes more 
  //  gritty/noisy


  shuffleBuffer(relativeDelayTimes, numDelayLines);
  // put an if(shuffleLengths) around, or make a case-statement switching on an "ordering" member
  // with shuffling, the initial impulses are not so cramped towards the beginning of the response
  // so it seems beneficial..but should be optional. maybe let the user select between various
  // shuffling algos

  allocateDelayLines();

  // or maybe use a formula that that is a weighted sum between this exponential rule and a linear
  // rule where the shape parameter determines the weights: shape = 0 - linear, very regular
  // echo-pattern, shape = 1 - exponential, more irregular. ....maybe call this parameter
  // "Regularity" or something
  // or maybe use a formula based on mode frequencies of a rectangular room - maybe generalize
  // the formula (replace the squares and square-root with an adjustable exponent or something)
}

void FeedbackDelayNetwork::setupReadIndices()
{
  for(int i = 0; i < numDelayLines; i++)
  {
    readIndices[i] = writeIndices[i] - (delayLineLengths[i] - 1);
    if( readIndices[i] < 0 )
      readIndices[i] += delayLineLengths[i];
  }
}

void FeedbackDelayNetwork::allocateMemory()
{
  freeMemory();
  readIndices        = new int[numDelayLines];
  writeIndices       = new int[numDelayLines];
  delayLineLengths   = new int[numDelayLines];
  relativeDelayTimes = new double[numDelayLines];
  setupRelativeDelayTimes(); // calls allocateDelayLines
}

void FeedbackDelayNetwork::allocateDelayLines()
{
  freeDelayLines();
  delayLines = new double*[numDelayLines];
  for(int i = 0; i < numDelayLines; i++)
  {
    delayLineLengths[i] = ceilInt(referenceDelayTime*relativeDelayTimes[i]*sampleRate);
    delayLines[i]       = new double[delayLineLengths[i]];
      // \todo: maybe use powers of 2 for the lengths to enable wraparound via masking
  }
}

void FeedbackDelayNetwork::freeMemory()
{
  delete[] readIndices;        readIndices        = nullptr;
  delete[] writeIndices;       writeIndices       = nullptr;
  delete[] delayLineLengths;   delayLineLengths   = nullptr;
  delete[] relativeDelayTimes; relativeDelayTimes = nullptr;
  freeDelayLines();
}

void FeedbackDelayNetwork::freeDelayLines()
{
  if( delayLines != nullptr )
  {
    for(int i = 0; i < numDelayLines; i++)
      delete[] delayLines[i];
    delete[] delayLines;
    delayLines = nullptr;
  }
}

// \todo inline this function later - during development, we want to have it non-inlined and in the
// cpp file to decrease build times
void FeedbackDelayNetwork::processFrame(double *inOutL, double *inOutR)
{
  //using AT = RAPT::rsArrayTools;

  static const int bufSize = 256; // should be >= numDelayLines
  double delayLineOuts[bufSize];

  int i;

  for(i = 0; i < numDelayLines; i++)
  {
    // increment read/write indices with wraparound:
    writeIndices[i] = (writeIndices[i] + 1) % delayLineLengths[i];
    readIndices[i]  = (readIndices[i]  + 1) % delayLineLengths[i];
    // optimize: let the buffer-lengths (not actual delayline lengths) be powers of 2 and do the 
    // wrap-around via bit-mask

    // pull out the outputs of the delaylines:
    delayLineOuts[i] = delayLines[i][readIndices[i]];


    // apply damping-filters to the delayline-outputs:
    //delayLineOuts[i] = dampingFilters[i].getSample(delayLineOuts[i]); // reactivate later


    // stuff the new input into the delaylines via the injection-vectors:
    //delayLines[i][readIndices[i]] = injectionVectorL[i] * *inOutL + injectionVectorR[i] * *inOutL;
    delayLines[i][writeIndices[i]] = 0.5 * (*inOutL + *inOutR);
  }


  // compute wet signal for left and right channel by using output vectors:
  // 1 0 -1 0 1 0 -1 0.... and 0 1 0 -1 0 1 0 -1
  double wetL = 0.0;
  double wetR = 0.0;


  /*
  double sign = 1.0;
  for(i = 0; i < numDelayLines; i+=2)
  {
    wetL += sign * delayLineOuts[i];
    wetR += sign * delayLineOuts[i+1];
    sign *= -1.0;
  }
  */

  // output vectors: left: [1 -0.5 -1 0.5], right: [0.5 1 -0.5 -1] ...and then repeated. these
  // output vectors are orthogonal (but maybe not orthonormal - fix this):
  for(i = 0; i < numDelayLines; i += 4)
  {
    wetL +=         delayLineOuts[i]
            - 0.5 * delayLineOuts[i+1]
            -       delayLineOuts[i+2]
            + 0.5 * delayLineOuts[i+3];
    wetR +=   0.5 * delayLineOuts[i]
            +       delayLineOuts[i+1]
            - 0.5 * delayLineOuts[i+2]
            -       delayLineOuts[i+3];
  }
  // ToDo: 
  // -check, if left and right impulse responses are uncorrelated
  // -this output vector implies that 4 delaylines and a 4x4 feedback matrix is the smallest 
  //  supported case - which makes sense
  // -use a factor that compensates for the number of delaylines N, maybe 1/sqrt(N), because
  //  currently, the impulse response gets louder, the more delaylines we use

  // Apply feedback matrix via a generalized fast Hadamard transform:
  double feedback = 0.9;    // 0.9 currently hardcoded decay
  RAPT::rsFGHT(delayLineOuts, numDelayLines, a, b, -b, a); 
  for(i = 0; i < numDelayLines; i++)
    delayLines[i][writeIndices[i]] += feedback * delayLineOuts[i];
  // todo: 
  // -compute feedback factor from a desired decay time RT60
  // -maybe we should take the average delay-time as basis for this calculation

  // preliminary - later use members:
  double dryGain = 0.0;
  double wetGain = 1.0; // * 1.0/sqrt((double)numDelayLines);
    // compensates for the energy-addition of the delayline-outputs

  // mix dry and wet:
  *inOutL = dryGain*(*inOutL) + wetGain*wetL;
  *inOutR = dryGain*(*inOutR) + wetGain*wetR;
}
