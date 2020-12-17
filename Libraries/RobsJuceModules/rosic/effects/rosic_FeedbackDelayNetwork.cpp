//#include "rosic_FeedbackDelayNetwork.h"
//using namespace rosic;

// a little helper function that is supposed to reorder a buffer (of delayline-lengths) that is
// initially in ascending order in a way that maximizes the average crossfeedback between
// delaylines with most different lengths (when the diffusion parameter is neither 0 nor 1)
// we implement it as template function outside the class because we need to apply it to the vector
// of output-gains as well
// maybe try the same with a buffer of initially descending lengths - this promotes stroger
// crossfeedback between the longer delaylines (whearas we currently have stronger crossfeedback
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


// Construction/Destruction:

FeedbackDelayNetwork::FeedbackDelayNetwork()
{
  // some test code for development:
  int buf8[8];
  RAPT::rsArrayTools::fillWithIndex(buf8, 8);
  shuffleBuffer(buf8, 8);

  int buf16[16];
  RAPT::rsArrayTools::fillWithIndex(buf16, 16);
  shuffleBuffer(buf16, 16);

  int buf32[32];
  RAPT::rsArrayTools::fillWithIndex(buf32, 32);
  shuffleBuffer(buf32, 32);







  delayLines         = NULL;
  readIndices        = NULL;
  writeIndices       = NULL;
  delayLineLengths   = NULL;
  relativeDelayTimes = NULL;

  //numDelayLines      = 2;
  //log2NumDelayLines  = 1;
  //numDelayLines      = 4;
  //log2NumDelayLines  = 2;
  //numDelayLines      = 8;
  //log2NumDelayLines  = 3;
  numDelayLines      = 16;
  log2NumDelayLines  = 4;
  //numDelayLines      = 32;
  //log2NumDelayLines  = 5;
  //numDelayLines      = 64;
  //log2NumDelayLines  = 6;
  //numDelayLines      = 128;
  //log2NumDelayLines  = 7;

  //numDelayLines      = 1024;
  //log2NumDelayLines  = 10;


  sampleRate         = 44100.0;
  referenceDelayTime = 20.0/1000.0;  // 20 ms
  //referenceDelayTime = 0.5/1000.0;  // 2 ms

  allocateMemory();
  //setDiffusion(100.0);
  setDiffusion(0.0);
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

void FeedbackDelayNetwork::fastGeneralizedHadamardTransform(
  double *x, int N, int log2N, double *y, double a, double b, double c, double d)
{
  //double *y               = (double*) alloca(N*sizeof(double));
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
    memcpy(y, x, N*sizeof(double)); // try rosic::copy instead, measure performance
     // huh - shouldn't we copy from y to x? ...but maybe it works because the pointers are also
     // swapped ->that's confusing, clean this up - maybe just rename the flag into something
     // neutral
}
// todo: try to write the algorithm in a way that works in place - maybe it's sufficient to pull 
// out the x[2*j], x[2*j+1] into temp variables? ...but i actually don't think so -> figure out!
// ..i think it would work, if we also write y into slots 2*j, 2*j+1 - but then we would end up
// with a differently ordered output - but that may actually not be an issue in this context, but 
// it could also be counteracted by a subsequent reordering algo, like in the FFT with the 
// bit-reversed ordering...maybe it's even the same re-ordering here?
// maybe the algo should be moved to some appropriate place in rapt, too
// 

void FeedbackDelayNetwork::fastInverseGeneralizedHadamardTransform(
  double *x, int N, int log2N, double *work, double a, double b, double c, double d)
{
  double s = 1.0 / (a*d - b*c);
  fastGeneralizedHadamardTransform(x, N, log2N, work, s*d, -s*b, -s*c, s*a);
}

void FeedbackDelayNetwork::setupRelativeDelayTimes()
{
  double dMax = 2.4; // hmmm..actually, this is not really the length of the longest delaylin
                     // ...find a better name

  //dMax = 5;        // maybe vary this between 1.5...2.5 via a user parameter "Shape"

  //dMax = sqrt((double) numDelayLines + 1);
  //dMax = pow((double) numDelayLines, 2.0/3.0);

  //dMax = (double) numDelayLines / 5.0;

  dMax = (double) numDelayLines / 5.0;
  dMax = 5;
  dMax = 3;


  for(int i = 0; i < numDelayLines; i++)
    relativeDelayTimes[i] = pow(dMax, (double) i / (double) numDelayLines);


  shuffleBuffer(relativeDelayTimes, numDelayLines);
    // put an if(shuffleLengths) around, or make a case-statement switching on an "ordering" member

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
  delete[] readIndices;        readIndices        = NULL;
  delete[] writeIndices;       writeIndices       = NULL;
  delete[] delayLineLengths;   delayLineLengths   = NULL;
  delete[] relativeDelayTimes; relativeDelayTimes = NULL;
  freeDelayLines();
}

void FeedbackDelayNetwork::freeDelayLines()
{
  if( delayLines != NULL )
  {
    for(int i = 0; i < numDelayLines; i++)
      delete[] delayLines[i];
    delete[] delayLines;
    delayLines = NULL;
  }
}

// \todo inle this function later - during development, we want to have it non-inlined and in the
// cpp file to decrease build times
void FeedbackDelayNetwork::processFrame(double *inOutL, double *inOutR)
{
  double delayLineOuts[1024];
  double fhtWorkspace[1024];  
  // later use a members with dynamic allocation, also 1024 is excessive






  int i;

  for(i = 0; i < numDelayLines; i++)
  {
    // increment read/write indices with wraparound:
    writeIndices[i] = (writeIndices[i] + 1) % delayLineLengths[i];
    readIndices[i]  = (readIndices[i]  + 1) % delayLineLengths[i];

    // pull out the outputs of the delaylines:
    delayLineOuts[i] = delayLines[i][readIndices[i]];


    // apply damping-filters to the delayline-outputs:
    //delayLineOuts[i] = dampingFilters[i].getSample(delayLineOuts[i]); // reactivate later


    // stuff the new input into the delaylines via the injection-vectors:
    //delayLines[i][readIndices[i]] = injectionVectorL[i] * *inOutL + injectionVectorR[i] * *inOutL;
    delayLines[i][writeIndices[i]] = 0.5 * (*inOutL + *inOutR);



    //int dummy = 0;
    // huh? when "i" gets larger than 131, it wraps around to zero, so we end up in an endless loop for
    // FDNs larger than 128. WTF?

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
    //double c = 1.0 / (0.25*i+1);
    double c = 1.0;
    wetL += c * (        delayLineOuts[i]
                 - 0.5 * delayLineOuts[i+1]
                 -       delayLineOuts[i+2]
                 + 0.5 * delayLineOuts[i+3]);
    wetR += c * (  0.5 * delayLineOuts[i]
                 +       delayLineOuts[i+1]
                 - 0.5 * delayLineOuts[i+2]
                 -       delayLineOuts[i+3]);
  }



  // apply feedback via matrix:
  fastGeneralizedHadamardTransform(delayLineOuts, numDelayLines, log2NumDelayLines, fhtWorkspace, 
    a, b, -b, a);
  for(i = 0; i < numDelayLines; i++)
    delayLines[i][writeIndices[i]] += 0.9 * delayLineOuts[i];

  // preliminary - later use members:
  double dryGain = 1.0;
  double wetGain = 1.0; // * 1.0/sqrt((double)numDelayLines);
    // compensates for the energy-addition of the delayline-outputs

  // mix dry and wet:
  *inOutL = dryGain*(*inOutL) + wetGain*wetL;
  *inOutR = dryGain*(*inOutR) + wetGain*wetR;
}
