using namespace RSLib;

// Construction/Destruction:

rsPhonoFilter::rsPhonoFilter()
{
  mode = PRE_EMPHASIS;
  setSampleRate(44100.0);
  reset();                 
}

// Setup:

void rsPhonoFilter::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  calcCoeffs();
}

void rsPhonoFilter::setMode(int newMode)
{
  mode = newMode;
  calcCoeffs();
}

// Audio Processing:
    
void rsPhonoFilter::processBlock(float *in, float *out, int blockSize)
{
  for(int n = 0; n < blockSize; n++)
    out[n] = (float) getSample(in[n]);
}

// Misc:

void rsPhonoFilter::reset()
{
  y2 = y1 = x2 = x1 = 0.0;
}

double rsPhonoFilter::getMagnitudeAt(double frequency)
{
  return biquadMagnitudeAt(b0, b1, b2, a1, a2, 2*PI*frequency/sampleRate);
}

void rsPhonoFilter::calcCoeffs()
{  
  double bb0[2], bb1[2], aa1[2]; // first order filter coeffs for the two stages

  // model 1st stage:
  double f[3], m[3], w[3];
  f[0] = 0.0;    
  f[1] = 1000.0; 
  f[2] = sampleRate/3;  
  m[0] = prototypeMagnitude1(f[0]);
  m[1] = prototypeMagnitude1(f[1]);
  m[2] = prototypeMagnitude1(f[2]);
  for(int i = 0; i < 3; i++)
    w[i] = 2*PI*f[i]/sampleRate;
  magnitudeMatchedOnePoleCoeffs(bb0[0], bb1[0], aa1[0], w, m);

  // model 2nd stage:
  m[0] = prototypeMagnitude2(f[0]);
  m[1] = prototypeMagnitude2(f[1]);
  m[2] = prototypeMagnitude2(f[2]);
  magnitudeMatchedOnePoleCoeffs(bb0[1], bb1[1], aa1[1], w, m);

  // consolidate the two 1st order filters into a single biquad:
  twoOnePolesToBiquad(bb0, bb1, aa1, b0, b1, b2, a1, a2);

  // invert, if desired:
  if( mode == DE_EMPHASIS )
    invertBiquad(b0, b1, b2, a1, a2);
}

double rsPhonoFilter::prototypeMagnitudeAt(double frequency)
{
  return prototypeMagnitude1(frequency) * prototypeMagnitude2(frequency);
}

double rsPhonoFilter::prototypeMagnitude1(double frequency)
{
  return unnormalizedMagnitude1(frequency) / unnormalizedMagnitude1(1000);
}

double rsPhonoFilter::prototypeMagnitude2(double frequency)
{
  return unnormalizedMagnitude2(frequency) / unnormalizedMagnitude2(1000);
}
  
double rsPhonoFilter::unnormalizedMagnitude1(double frequency)
{
  // component values:
  double C1 = 1.0e-6;   // 1 uF
  double R1 = 3.18e+3;  // 3.18 kOhm
  double R2 = 353.33;   // 353.33 Ohm

  // evaluate and return transfer-function magnitude at s = j*w, w=2*PI*frequency:
  rsComplexDbl s  = rsComplexDbl(0.0, 2*PI*frequency);
  rsComplexDbl H1 = (s + 1.0/(R1*C1)) / (s + (R1+R2)/(R1*R2*C1));
  return H1.getRadius();
}
 
double rsPhonoFilter::unnormalizedMagnitude2(double frequency)
{
  double C2 = 10.0e-9;  // 10 nF
  double R3 = 7.5e+3;   // 7.5 kOhm
  double R4 = 353.33;   // 353.33 Ohm
  rsComplexDbl s  = rsComplexDbl(0.0, 2*PI*frequency);
  rsComplexDbl H2 = (s + 1.0/(R3*C2)) / (s + (R3+R4)/(R3*R4*C2));
  return H2.getRadius();
}
