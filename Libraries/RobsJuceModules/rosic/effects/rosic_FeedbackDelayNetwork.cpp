
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
  numDelayLines   = 16;  // 16 sounds smoother, denser, noisier
  //numDelayLines = 32;  // 32 sounds still denser
  //numDelayLines = 64;  // still denser, initial section acquires a snare'ish character
  //numDelayLines = 128; // initial "electronic snare" character becomes very apparent
  //numDelayLines = 256; // snare sound becomes higher pitched
  // ...conclusion: 
  // -16 or 32 seems to be the sweet spot for reverb, 8 and 64 are also usable
  // -larger sizes acquire this weird pitch-sweepdown electronic drum sound
  //  -maybe that can be avoided by using only a subset of the delaylines for output
  //  -> experiment with output vectors containing zeros, maybe an output-density parameter
  //  can be defined
  // -maybe try to use not one N=64 network but a parallel connection of 4 N=16 networks with
  //  different delay-sttings - maybe that creates more density and also avoids the pitch-sweepdown
  // -try different formulas for the relative delay times, currently, a power rule is used

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

  double d = 0.01*diffusion;  // 0..100 -> 0..1
  //d = d*d*d;  // nah - that's not yet a good mapping
  //d = mapLinearToRational(d, +0.7);
  //d = d*d;

  // I think, we need an odd function that maps 0 to 0, 1 to 1, 2 to 2, and has adjustable slope s 
  // at x = 1 (should be s < 1). rationale: the points d = -2 and d = 2 both lead to 3 diffusion 
  // and the parameter should be periodic: we may wrap to -2 when d > 2 and vice versa - it has a 
  // period of 4
  //   f( x) = a1*x +   a3*x^3 +   a5*x^5
  //   f'(x) = a1   + 3*a3*x^2 + 5*a5*x^4
  // so:
  //   f( 1) = 1 = a1   +   a3   +   a5
  //   f( 2) = 2 = a1*2 +   a3*8 +   a5*32
  //   f'(1) = s = a1   + 3*a3   + 5*a5 
  // then adjust s by ear, test with wraparound

  // var("x s a1 a3 a5")
  // eq1 = a1 +     a3   + a5    == 1
  // eq2 = a1*2 +   a3*8 + a5*32 == 2
  // eq3 = a1   + 3*a3 + 5*a5    == s
  // solve([eq1,eq2,eq3],[a1,a3,a5])
  //
  // a1 == -2/3*s + 5/3, a3 == 5/6*s - 5/6, a5 == -1/6*s + 1/6
  // a1 == -2*s/3 + 5/3, a3 == 5*s/6 - 5/6, a5 == -1*s/6 + 1/6
  //
  // with that function, the range between 0..100 is mapped differently than between 200..100 which
  // may or may not be desirable - the 100..200 range has a different character, so a different 
  // curve (that's not some sort of mirror image of the 0..100 curve) may actually be desirable
  // -> listening tests are needed


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

  // i think, this required mapping comes about because the perceived diffusion depends on the 
  // ratio a/b which has singularities quickly towards the ends of of the interval ...maybe some 
  // tan or atan based formula could help?

  RAPT::rsSinCos(d*PI/4, &b, &a);

  // \todo - when introducing modulation for d (or the angle phi), check if it is better to
  // modulate the angle before or after the non-linear mapping. if the modulation is applied before
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
// But before doing so, do performance tests - it may be more efficient, even though it doesn't
// work in place. Maybe let the number of delaylines be a compile-time parameter, i.e. a template
// parameter. Then, the compiler might be able to unroll the loop.
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
  // actually, we can simplify this, based on whether log2N is even of odd and get rid of that
  // flip-flopping flag
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

/*
void FeedbackDelayNetwork::fastInverseGeneralizedHadamardTransform(
  double *x, int N, int log2N, double *work, double a, double b, double c, double d)
{
  double s = 1.0 / (a*d - b*c);
  fastGeneralizedHadamardTransform(x, N, log2N, work, s*d, -s*b, -s*c, s*a);
}
// also obsolete now
*/

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
  dMax = 3;
  //dMax = 4;
  //dMax = 5;
  //dMax = 6;
  //dMax = 8;
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

  using LT = RAPT::rsLinearTransforms;
  if(use3x3) LT::kronecker3x3(delayLineOuts, numDelayLines, seedMat3x3);
  else       LT::kronecker2x2(delayLineOuts, numDelayLines, a, b, -b, a); // use seedMat2x2

  for(i = 0; i < numDelayLines; i++)
    delayLines[i][writeIndices[i]] += feedback * delayLineOuts[i];
  // todo: 
  // -compute feedback factor from a desired decay time RT60
  // -maybe we should take the average delay-time as basis for this calculation
  // -try the symmetric variant of the FGHT using a,b,b,-a - this here with a,b,-b,a is the 
  //  antisymmetric variant. only the symmetric variant is actually a special case of the Hadamard
  //  trafo...but are both actually equivalent for some other setting of diffusion and maybe some 
  //  permutaion of the delaylines? no - i think, it can't because a,-a are on the diagonal, 
  //  meaning they are responsible for self feedback - no way can we turn an a,-a into an a,a, so
  //  the self-feedback is ways same-signed or always different-signed - it's really a different
  //  family of feedback matrices - so: provide an option to switch between both variants, maybe 
  //  the option could be named: symmetric...or maybe the name could be based on whether it's
  //  a pure rotation or a rotation + reflection
  // -idea: can we also start with a 3x3 seed matrix and would that give us more freedom? ...maybe
  //  3x3 matrix kroneckered with itself, giving a 27x27 matrix could be interesting as well


  // preliminary - later use members:
  double dryGain = 0.0;
  double wetGain = 1.0; // * 1.0/sqrt((double)numDelayLines);
    // compensates for the energy-addition of the delayline-outputs

  // mix dry and wet:
  *inOutL = dryGain*(*inOutL) + wetGain*wetL;
  *inOutR = dryGain*(*inOutR) + wetGain*wetR;
}


/*

ToDo:

-the memory management is totally silly - way to many re-allocations
 -allocate memory for all delaylines once in the constructor (via an allocate function) and then
  use pointers into that memory

Ideas:

-allow a different set of a,b values (diffusion) for each stage, maybe arrange it such that the
 feedback between longer lines spreads faster and between shorter lines more slowly, per roundtrip,
 such that the total time for the build-up is roughly the same in all lines

-maybe allow a fully general product of matrices that are multiplied together via the Kronecker 
 product..maybe call it KronVerb, KronReverb, KronyVerb

-maybe we could introduce a permutation as part of the feedback, maybe that could be time-varying
 to introduce further randomness

-try a Kronecker trafo with complex coeffs. oh - and figure out for which choices of a,b,c,d the
 resulting matrix is unitary - will it be when using a matrix with c=-b, d=a or c=b, d=-a. 
 oh - i checked with the 2x2 seed matrix it is unitary indeed does this imply the higher order 
 matrices are unitary too? yes, i think preserving unitariness this is a general feature of the 
 Kronecker product. 

-Try using dual numbers or hyperbolic numbers...dual numbers could indeed be interesting in an FDN
 because they model derivatives...maybe that leads to some sort of highpass feature? ..can we also 
 have a number system tha simulates integration?

-Try arbitrary a,b,c,d and rescale the output to have the same length as the input vector - this 
 makes the system nonlinear (really?) - maybe in an interesting way?

-Maybe it could be useful to do incomplete Kronecker trafos, i.e. let the outer while(h < N) loop
 only run up to N/2 or N/4, etc? ..we could pass a "divider" parameter or somehow pass a 
 loopLimit parameter

-When using the 3x3 version, we have 3 density parameters - we could modulate them by 3 different
 modulation frequencies (mod signals are obtained by bandpassing the rectified (and possibly 
 leveled) output - or maybe the leveling should come after the bandpass)
 -maybe the 3x3 rotation should be specified not in terms of Euler angles but by an axis and an 
  angle..or maybe in some other way -> figure out how 3D rotations are most conveniently 
  parametrized in graphics
 -maybe we could make 3 groups of delaylines (like short, medium, long or somehow sorted by 
  harmonic relationships) and have the 3 angles control the intermixing between any pair of groups

-to avoid the initial predelay (the delay before the 1st reflection arrives), maybe use multiple
 output taps (that are not part of the feeback loop)...maybe it's sufficient to let only one 
 delayline have multiple output taps (perhaps the shortest...or the longest, if we want these 
 output spikes to intermix with the build-up of the reflections)...or maybe use a dedicated 
 delayline for the early reflections, maybe the early reflections should be spatialized via
 HRTF processing

-try a different ordering of the delaylines using the shuffling function in the cpp file.
 this is supposed to be better for the diffusion parameter: short delaylines would mostly
 crossfeed long ones and vice versa. if this is used, shuffle the output gains in the same way
 (-> introduce output vectors as members and throw the same function at them) ...hmm but maybe
 not - it seems that shuffled delayline-length together with non-shuffled output vector work well
 - especiall with low diffusion

-introduce modulation of the angle "phi" by rectifying the FDNs output signal (maybe use
 the mono-sum), pass it into a pair of filters tuned to the modulation-frequency and are 90
 degrees out of phase (use the ModalFilter class for this). Two out-of-phase signals are needed in
 order to obatin the instantaneous envelope of the modulation signal in order to divide by it to
 obtain an amplitude-normalized modulation signal from either of the two filter outputs (perhaps
 use the sine-component). 

-For the output vector: it doesn't seem desirable to just repeat the first 4 values multiple times, 
 because that seems to lead to constructive interference in the left channel and destructive 
 interference in the right channel between the 1st and 2nd order reflections. Maybe permute the 
 next 4 values (in the same way for left and right output vector). It also seems a bit unnatural 
 that they have only 2 values for the amplitude. Maybe start with two random vectors and 
 orthogonalize them via Gram-Schmidt. Maybe do the same for the injection vectors. Maybe after 
 orthogonalization, do a (joint) sorting such that the values where (aL^2 + aR^2) are greatest come 
 first where aL,aR are the output weights for left and right. We want the first feq reflections be 
 loudest.

-Provide a function compute the poles. The eigenvalues of the feedback matrix can be easily 
 computed, see here:
 https://en.wikipedia.org/wiki/Kronecker_product
 https://en.wikipedia.org/wiki/Kronecker_product#Abstract_properties
 ...we can easily compute the eigenvalues! the eigenvalues of a Kronecker product of two matrices
 A,B are given by the products of their eigenvalues...that may be helpful in figuring out 
 expressions for the poles of the FDN and maybe even to place the resonance frequencies of the
 FDN deliberately. I think, the eigenvalues of a 2D rotation matrix is just the complex number on
 the unit circle corresponding to the rotation angle (and the eigenvectors are not real). But what 
 about the delays? They have impact on the poles, too. I think, the number of poles may be equal to
 the product of the delaytimes in samples? Is that correct?

-Maybe apply a nonlinear waveshaping function in the feedback loop. The idea is to modify the 
 effective feedback strength depending on the signal level. Maybe something that has a slope < 1 
 around zero and approaches unit slope at higher values - may lead to a sort of cutting off the 
 tail early as in gated reverb. But maybe it's better to do this based on an envelope detector.

-Implement a pre-diffuser based on a series of allpass delays. Try short lengths with prime numbers
 like 2,3,5,7. Maybe their feedback coeff should be related to the length (lower for longer ones). 
 Or maybe use post-diffuser with different delays for left and right channel, for example:
 L: 5,13,17,31, R: 7,11,19,29 - the corresponding L/R pairs are chosen to be twin primes and the 
 total delay is the same for left and right: 5+13+17+31 = 7+11+19+29 = 66. Or maybe include these 
 allpasses also in the feedback loop - maybe at the input of the delaylines but after the feedback 
 point. Maybe the feedback coeff of each allpass should be inversely proportional to the allpass 
 length/delay, or some power of it (with exponent <= 1). Maybe try to also use a post-diffuser. It 
 shouldn't make a difference for an LTI-FDN, but if we go non-LTI, it may matter

-Compute delayline lengths from geometric considerations. For example, use the modal distribution 
 of rectangular rooms. The formula for the modal frequencies is:
   f ~ sqrt( nx^2/Lx^2 + ny^2/Ly^2 + nz^2/Lz^2 ) where nx,ny,nz are the 3 independent modal indices
 (including zero!) and Lx,Ly,Lz are the dimensions of the room. The delays should then be 
 reciprocals of these (i think).

-Maybe try powers of the golden ration (reduced to the octave [1,2)

-Instead of using the 2 produced outputs for left/right, try using them for mid/side. Maybe 
 introduce a stereo width parameter as well

-Introduce delay-based panning of the middle delaylines, as explained here:
 https://www.kvraudio.com/forum/viewtopic.php?f=33&t=123095

-For modulation: Schlecht suggests that modulation should be confined to mid- and high frequencies 
 (above 400 Hz). This can't be done by modulating he matrix but maybe we could modulate the 
 effective delayline lengths by putting a (modulated, 1st or 2nd order) allpass in series which 
 only affects the delay of high-frequencies. He suggests to split the signal and use a modulated 
 feedback matrix only for the high frequencies. We could implement using multiplpe feedback 
 matrices in parallel efficiently by using simd - with rsFloat32x4, we can compute 4 different 
 feedback matrices at the same time, so we could use it for both: early vs late decay and 
 modulated vs unmodulated

-At the moment, we use a common damping factor for all delaylines. That may be justified when 
 there's a lot of crosstalk between the delaylines (i.e. high feedback diffusion), but when the 
 matrix is more like parallel bank of combs, really each delayline should have its own attenuation 
 coeff

-Would it make sense to replace the delaylines with bidirectional delayline as used in waveguide 
 modeling? ..i.e. a feedback waveguide network (FWN)



Resources:

Julius Smith:
https://www.dsprelated.com/freebooks/pasp/Feedback_Delay_Networks_FDN.html
https://www.dsprelated.com/freebooks/pasp/FDN_Reverberation.html

Digital delay networks for designing artificial reverberators (Jean-Marc Jot, Antoine Chaigne):
https://www.researchgate.net/publication/243779004_Digital_delay_networks_for_designing_artificial_reverberators

Maximally diffusive yet efficient feedback delay networks for artificial reverberation (Davide Rocchesso)
https://www.researchgate.net/publication/3342318_Maximally_diffusive_yet_efficient_feedback_delay_networks_for_artificial_reverberation

Sebastian Schlecht:
Feedback Delay Networks in Artificial Reverberation and Reverberation Enhancement (PhD Thesis):
https://www.researchgate.net/publication/322951473_Feedback_Delay_Networks_in_Artificial_Reverberation_and_Reverberation_Enhancement
https://dafx2020.mdw.ac.at/proceedings/papers/DAFx2020_paper_53.pdf (with matlab toolbox)

Choosing Optimal Delays for Feedback Delay Networks:
http://pub.dega-akustik.de/DAGA_2014/data/articles/000025.pdf

there was also some paper by Miller Puckette..


About reverb in general:
https://valhalladsp.com/2011/01/21/reverbs-diffusion-allpass-delays-and-metallic-artifacts/
https://www.kvraudio.com/forum/viewtopic.php?t=564078


Sound material for quality assessment:
https://tech.ebu.ch/publications/sqamcd
https://tech.ebu.ch/docs/tech/tech3253.pdf

*/