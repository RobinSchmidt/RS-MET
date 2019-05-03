// Construction/Destruction:

template<class TSig, class TPar>
rsPhonoFilter<TSig, TPar>::rsPhonoFilter()
{
  mode = PRE_EMPHASIS;
  setSampleRate(44100.0);
  reset();                 
}

// Setup:

template<class TSig, class TPar>
void rsPhonoFilter<TSig, TPar>::setSampleRate(TPar newSampleRate)
{
  sampleRate = newSampleRate;
  calcCoeffs();
}

template<class TSig, class TPar>
void rsPhonoFilter<TSig, TPar>::setMode(int newMode)
{
  mode = newMode;
  calcCoeffs();
}

// Audio Processing:
    
template<class TSig, class TPar>
void rsPhonoFilter<TSig, TPar>::processBlock(TSig *in, TSig *out, int blockSize)
{
  for(int n = 0; n < blockSize; n++)
    out[n] = (float) getSample(in[n]);
}

// Misc:

template<class TSig, class TPar>
void rsPhonoFilter<TSig, TPar>::reset()
{
  y2 = y1 = x2 = x1 = 0.0;
}

template<class TSig, class TPar>
TPar rsPhonoFilter<TSig, TPar>::getMagnitudeAt(TPar frequency)
{
  return biquadMagnitudeAt(b0, b1, b2, a1, a2, 2*PI*frequency/sampleRate);
}

template<class TSig, class TPar>
void rsPhonoFilter<TSig, TPar>::calcCoeffs()
{  
  TPar bb0[2], bb1[2], aa1[2]; // first order filter coeffs for the two stages

  // model 1st stage:
  TPar f[3], m[3], w[3];
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

template<class TSig, class TPar>
TPar rsPhonoFilter<TSig, TPar>::prototypeMagnitudeAt(TPar frequency)
{
  return prototypeMagnitude1(frequency) * prototypeMagnitude2(frequency);
}

template<class TSig, class TPar>
TPar rsPhonoFilter<TSig, TPar>::prototypeMagnitude1(TPar frequency)
{
  return unnormalizedMagnitude1(frequency) / unnormalizedMagnitude1(1000);
}

template<class TSig, class TPar>
TPar rsPhonoFilter<TSig, TPar>::prototypeMagnitude2(TPar frequency)
{
  return unnormalizedMagnitude2(frequency) / unnormalizedMagnitude2(1000);
}
  
template<class TSig, class TPar>
TPar rsPhonoFilter<TSig, TPar>::unnormalizedMagnitude1(TPar frequency)
{
  // component values:
  TPar C1 = 1.0e-6;   // 1 uF
  TPar R1 = 3.18e+3;  // 3.18 kOhm
  TPar R2 = 353.33;   // 353.33 Ohm

  // evaluate and return transfer-function magnitude at s = j*w, w=2*PI*frequency:
  std::complex<TPar> s  = std::complex<TPar>(0.0, 2*PI*frequency);
  std::complex<TPar> H1 = (s + 1.0/(R1*C1)) / (s + (R1+R2)/(R1*R2*C1));
  return abs(H1);
}
 
template<class TSig, class TPar>
TPar rsPhonoFilter<TSig, TPar>::unnormalizedMagnitude2(TPar frequency)
{
  TPar C2 = 10.0e-9;  // 10 nF
  TPar R3 = 7.5e+3;   // 7.5 kOhm
  TPar R4 = 353.33;   // 353.33 Ohm
  std::complex<TPar> s  = std::complex<TPar>(0.0, 2*PI*frequency);
  std::complex<TPar> H2 = (s + 1.0/(R3*C2)) / (s + (R3+R4)/(R3*R4*C2));
  return abs(H2);
}
