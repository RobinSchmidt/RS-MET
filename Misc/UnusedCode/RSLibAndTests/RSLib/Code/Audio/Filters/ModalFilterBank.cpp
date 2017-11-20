using namespace RSLib;

//=================================================================================================
// class rsTwoPoleFilter:

rsTwoPoleFilter::rsTwoPoleFilter()
{
  a1 = a2 = 0;
  g  = 1;
  reset();
}

void rsTwoPoleFilter::setFrequencyAndDecay(double w, double d)
{
  double P = exp(-1/d);  // pole radius
  a1 = -2*P*cos(w);
  a2 = P*P;    
}

void rsTwoPoleFilter::setOutputGain(double newGain)
{
  g = newGain;
}

double rsTwoPoleFilter::getMagnitudeAt(double w)
{
  return biquadMagnitudeAt(1, 0, 0, a1, a2, w);
    // optimize: a simpler (less general) formula may be used
}

void rsTwoPoleFilter::reset()
{
  y1 = y2 = 0;
}

//=================================================================================================
// helper functions for damped-sine filter design:

//void RSLib::modalParamsToFilterCoeffs(double w, double A, double d, double p, 
//                                      double *g, double *b1, double *a1, double *a2)
//{
//  // implementation is obsolete - replace by new one in rsDampedSineFilter - still in the
//  // Experiments project. ah, and we should implement the filter in terms of b0, b1, a1, a2 
//  // instead g, b1, a1, a2. the new design is more flexible and the coefficient computations are
//  // simpler
//
//  // calculate intermediate variables:
//  double P, pp, ri, R;
//  P  = exp(-1.0/d);
//  p  = rsWrapToInterval(rsDegreeToRadiant(p), 0, 2*PI);
//  pp = p-PI/2;
//  ri = 0.5*tan(pp);
//  R  = rsSqrt(0.25+ri*ri);
//
//  // todo: use rsSinCos, double-angle formula:
//  double c1, c2, s1, s2;
//  c1 = cos(w);
//  c2 = cos(2*w);
//  s1 = sin(w);
//  s2 = sin(2*w);
//
//  // calculate coefficients:
//  *a2 = P*P;
//  *a1 = -2*P*c1;
//  *b1 = -(P/2)*(2*(1-c2)*ri+s2)/s1; 
//  *g  = A/(2*R);
//  if( p > PI )
//    *g = -*g;
//}

void RSLib::rsDampedSineFilter(double w, double A, double d, double p, double *b0, double *b1, 
  double *a1, double *a2)
{
  double cw, sw, cp, sp, P;
  rsSinCos(w, &sw, &cw);
  rsSinCos(p, &sp, &cp);
  P   = exp(-1.0/d);        // = exp(-alpha), pole radius
  *a1 = -2*P*cw;            // = -2*P*cos(w)
  *a2 = P*P;                // = P^2
  *b0 = A*sp;               // = A*sin(p)
  *b1 = A*P*(sw*cp-cw*sp);  // = A*P*sin(w-p) via addition theorem
}

double RSLib::findDecayScalerLess1(double c)
{
  if( c <= 0.0 || c >= 1.0 )
  {
    rsError("Function assumes 0 < c < 1");
    return 1.0;
  }

  // precomputations:
  double kp  = 1/c;                                    // location of the the peak of g(k)
  double k   = 1 + 2*(kp-1);                           // initial guess for the zero of g(k)
  double eps = std::numeric_limits<double>::epsilon(); // relative tolerance
  int    i   = 0;                                      // iteration counter

  // Newton iteration:
  double g, gp;      // g(k), g'(k)
  double kOld = 2*k; // ensure to enter the loop
  while( fabs(k-kOld) > k*eps && i < 1000 )
  {
    kOld = k;
    g    = log(k) + c*(1-k); // g(k)
    gp   = 1/k - c;          // g'(k)
    k    = k - g/gp;         // Newton step
    i++;                     // count iteration
  }
 
  return k;

  // \todo: check this function in the range 0 <= c < 1, if all works well and the iteration count 
  // is always low, get rid of the iteration counter - it serves a purpose only during development

  // \todo find a refined formula for the initial guess by plotting the output against the input in
  // the range 0...1 and fit a polynomial (or other suitable function) to the data
}

void RSLib::expDiffScalerAndTau2(double tau1, double tp, double *tau2, double *scaler)
{
  if( tp >= tau1 )
  {
    rsError("assumes tp < tau1");
    *tau2   = tau1;
    *scaler = 1.0;
    return;
  }
  double a1 = 1/tau1;
  double c  = a1 * tp;
  double k  = findDecayScalerLess1(c);
  double a2 = k*a1;
  double hp = exp(-a1*tp) - exp(-a2*tp); // peak height

  *tau2   = 1/a2;
  *scaler = 1/hp; 
}

//=================================================================================================
// class rsModalFilter:

rsModalFilter::rsModalFilter()
{
  b0 = 1.0;
  b1 = a1 = a2 = 0.0;
  reset();
}

void rsModalFilter::setModalParameters(double frequency, double amplitude, double decayTime, 
                                       double startPhase, double sampleRate)
{
  rsDampedSineFilter(2*PI*frequency/sampleRate, amplitude, decayTime*sampleRate, startPhase,
    &b0, &b1, &a1, &a2);  
}

void rsModalFilter::copyCoefficientsFrom(const rsModalFilter &other)
{
  b0 = other.b0;
  b1 = other.b1;
  a1 = other.a1;
  a2 = other.a2;
}

double rsModalFilter::getDecayTime(double sampleRate)
{
  return -1.0 / (log(rsSqrt(a2))*sampleRate);
}

double rsModalFilter::getMagnitudeAt(double w)
{
  return biquadMagnitudeAt(b0, b1, 0, a1, a2, w);
}

void rsModalFilter::reset()
{
  x1 = y1 = y2 = 0.0;
}

void rsModalFilter::processBlock(double in[], double out[], int blockSize)
{
  // under construction - still it is less efficient to call the block-based version compared to 
  // the sample-based version, also it is not yet entirely correct (the scale factor should be 
  // rendered into the state-variables - we need to swicth from g,b1 to b0,b1 coeffs

  out[0] = in[0] + b1*x1 - a1*y1 - a2*y2;
  x1 = in[0];  
  y2 = y1;
  y1 = out[0];

  if( blockSize < 2 )
    return;

  out[1] = in[1] + b1*x1 - a1*y1 - a2*y2;
  x1 = in[1];  
  y2 = y1;
  y1 = out[1];

  if( blockSize < 3 )
    return;

  for(int n = 2; n < blockSize; n++)
  {
    out[n] = in[n] + b1*in[n-1] - a1*out[n-1] - a2*out[n-2];
    //out[n] = g*y;  // we need to use b0, b1 instead of g, b1
  }

  x1 = in[blockSize-1];
  y2 = out[blockSize-2];
  y1 = out[blockSize-1];

  //scale(out, blockSize, g);

  /*
  for(int n = 0; n < blockSize; n++)
  {
    //out[n] = getSample(in[n]);
    double y = in[n] + b1*x1 - a1*y1 - a2*y2;
    x1 = in[n];
    y2 = y1;
    y1 = y;
    out[n] = g*y;
  }
  */
}

//=================================================================================================
// class rsNonlinearModalFilter:

rsNonlinearModalFilter::rsNonlinearModalFilter()
{  
  a  = 0.0; 
  cr = 1.0;
  ci = 0.0;

  startPhase    = 0.0;
  amplitude     = 1.0;
  phaseModByIm  = 0.0;
  phaseModByAbs = 0.0;

  reset();
}

void rsNonlinearModalFilter::setModalParameters(double frequency, double amplitude, double decayTime,                                  
                                                double startPhase, double sampleRate)
{
  // compute recursion coeff:
  double w     = 2*PI*frequency/sampleRate;
  double alpha = 1.0 / (decayTime*sampleRate); 
  double r     = exp(-alpha);
  a.setRadiusAndAngle(r, w);

  this->amplitude  = amplitude;
  this->startPhase = startPhase;

  // compute output weights:
  rsSinCos((PI/180.0)*startPhase, &cr, &ci);
  //cr *= amplitude;
  //ci *= amplitude;
}

void rsNonlinearModalFilter::setAmplitude(double newAmplitude)
{
  amplitude = newAmplitude;
}

void rsNonlinearModalFilter::setPhaseModulation(double newPhaseModulation)
{
  // todo: rename into phasemodualtionByIm, try to find a mor intuitive parametrization in terms
  // of the overtone spectrum (maybe a formula may be derived that relates a spectral slope to the 
  // modulation index or something) ...or maybe a slope (in dB/oct) can be related to the overall
  // level in dB ...so we may dial in a value in (dB/oct)/dB = dB^2/oct
  phaseModByIm = newPhaseModulation;
}

void rsNonlinearModalFilter::copyCoefficientsFrom(const rsNonlinearModalFilter &other)
{
  a  = other.a;
  cr = other.cr;
  ci = other.ci;
}

void rsNonlinearModalFilter::reset()
{
  z = 0.0;
}

//=================================================================================================
// class ModalFilterWithAttack:

void rsModalFilterWithAttack::setModalParameters(double frequency, double amplitude, 
                                                 double attackTime, double decayTime, 
                                                 double startPhase, double sampleRate, 
                                                 double detuneFactor)
{
  rsAssert(attackTime < decayTime);  // attackTime >= decayTime will not work (because of math)

  double tau1, tau2, scaler;
  tau1 = decayTime;
  expDiffScalerAndTau2(tau1, attackTime, &tau2, &scaler);
  amplitude *= scaler;
  modalFilter1.setModalParameters(frequency,              amplitude, tau1, startPhase, sampleRate);
  modalFilter2.setModalParameters(frequency*detuneFactor, amplitude, tau2, startPhase, sampleRate);

  // \todo "detuneFactor" does not really work well because when attack and decay are very similar, 
  // the amplitude explodes (due to a high value of "scaler") - either remove this parameter or 
  // find a way to alleviate this (maybe the amplitude excess can be computed - math has to be 
  // worked out) the detune is supposed to introduce some roughness into the transient by detuning 
  // the quickly decaying sinusoid
}

double rsModalFilterWithAttack::getLength(double decayLevel, double sampleRate)
{
  double td = modalFilter1.getDecayTime(sampleRate);       // decay time
  double a1 = 1.0 / td;
  double a2 = 1.0 / modalFilter2.getDecayTime(sampleRate);
  double tp = (log(a1)-log(a2))/(a1-a2);                   // peak time
  return tp + rsTauToDecayTime(td, decayLevel);
}

void rsModalFilterWithAttack::reset()
{
  modalFilter1.reset();
  modalFilter2.reset();
}

//=================================================================================================
// class ModalFilterWithAttack2:

rsModalFilterWithAttack2::rsModalFilterWithAttack2()
{
  //b0 = 1.0;
  b1 = b2 = b3 = a1 = a2 = a3 = a4 = 0.0;
  reset();
}

void rsModalFilterWithAttack2::setModalParameters(double frequency, double amplitude, 
                                                  double attackTime, double decayTime, 
                                                  double startPhase, double sampleRate)
{
  rsAssert(attackTime < decayTime);  // attackTime >= decayTime will not work (because of math)

  double tau1, tau2, scaler;
  double w = 2*PI*frequency/sampleRate;
  double p = startPhase;
  tau1 = decayTime;
  expDiffScalerAndTau2(tau1, attackTime, &tau2, &scaler);

  // compute coefficients for a parallel connection of filters:
  double a = amplitude * scaler;
  double b10, b11, a11, a12;
  double b20, b21, a21, a22;
  rsDampedSineFilter(w, a, tau1*sampleRate, p, &b10, &b11, &a11, &a12);
  rsDampedSineFilter(w, a, tau2*sampleRate, p, &b20, &b21, &a21, &a22);

  // convert parallel connection into a single 4th order filter:
  //b0 = b10-b20; // b0 comes out as zero
  b1 = b11-b21+a21*b10-a11*b20;
  b2 = a21*b11-a11*b21+a22*b10-a12*b20;
  b3 = a22*b11-a12*b21;
  a1 = a21+a11;
  a2 = a22+a11*a21+a12;
  a3 = a11*a22+a12*a21;
  a4 = a12*a22;
}

void rsModalFilterWithAttack2::reset()
{
  w1 = w2 = w3 = w4 = 0.0;
}

//=================================================================================================
// class ModalFilterBank:

// construction/destruction:

rsModalFilterBank::rsModalFilterBank()
{
  sampleRate         = 44100.0;
  referenceFrequency = 440.0;
  referenceAttack    = 0.1;
  referenceDecay     = 1.0;
  numModes           = maxNumModes;
  modalFilters.reserve(maxNumModes);
  for(int m = 0; m < maxNumModes; m++)
    modalFilters.push_back(rsModalFilterWithAttack());
}

rsModalFilterBank::~rsModalFilterBank()
{

}

// parameter settings:

void rsModalFilterBank::setSampleRate(double newSampleRate)
{
  sampleRate = newSampleRate;
  calculateModalFilterCoefficients();
}

void rsModalFilterBank::setReferenceFrequency(double newFrequency)
{
  referenceFrequency = newFrequency;
  calculateModalFilterCoefficients();
}

void rsModalFilterBank::setReferenceAttack(double newAttack)
{
  referenceAttack = newAttack;
  calculateModalFilterCoefficients();
}

void rsModalFilterBank::setReferenceDecay(double newDecay)
{
  referenceDecay = newDecay;
  calculateModalFilterCoefficients();
}

void rsModalFilterBank::setModalParameters(rsVectorDbl newFrequencies, rsVectorDbl newAmplitudes, 
                                           rsVectorDbl newAttackTimes, rsVectorDbl newDecayTimes, 
                                           rsVectorDbl newStartPhases)
{
  frequencies = newFrequencies;
  amplitudes  = newAmplitudes;
  attackTimes = newAttackTimes;
  decayTimes  = newDecayTimes;
  startPhases = newStartPhases;
  calculateModalFilterCoefficients();
}

// inquiry:

double rsModalFilterBank::getLength(double decayLevel)
{
  double max = 0.0;
  double tmp;
  for(int m = 0; m < decayTimes.dim; m++)
  {
    tmp = modalFilters[m].getLength(decayLevel, sampleRate);
    if( tmp > max )
      max = tmp;
  }
  return max;
  // actually, we could just look it up in our attckTime/decayTime members instead of letting the
  // filters compute their lengths
}

// audio processing:

void rsModalFilterBank::processBlock(double in[], double out[], int blockSize)
{
  for(int n = 0; n < blockSize; n++)
    out[n] = getSample(in[n]);

  // \todo optimize this
}

// others:

void rsModalFilterBank::resetModalFilters()
{
  for(unsigned int m = 0; m < modalFilters.size(); m++)
    modalFilters[m].reset();
}

void rsModalFilterBank::calculateModalFilterCoefficients()
{
  int nm = rsMin(numModes, frequencies.dim, amplitudes.dim, decayTimes.dim);
  nm = rsMin(nm, startPhases.dim);
  for(int m = 0; m < nm; m++)
  {
    modalFilters[m].setModalParameters(
      referenceFrequency * frequencies.v[m], 
      amplitudes.v[m],
      referenceAttack * attackTimes.v[m],
      referenceDecay * decayTimes.v[m], 
      startPhases.v[m], 
      sampleRate); 
  }
}

// static member functions:

double rsModalFilterBank::modeDecayTime(double f, double fc, double p)
{
  double k = pow(fc, -p);
  return (1-k) / ((1-k*2) + pow(f/fc, p));
}

rsVectorDbl rsModalFilterBank::randomModePhases(const rsVectorDbl &a, 
                                                double randomness, 
                                                int seed)
{
  // The formulas used here were obtained by requiring a1*sin(p1) + a2*sin(p2) = k and solving 
  // for p2. For the last value in case of an odd number of modes, it is: a*sin(p) = k.
  // In order to let one mode cancel the other at n = 0, we use k = 0. We need to ensure that the 
  // mode for which we adjust the phase has at least the same amplitude as the mode for which we
  // assume a fixed phase in order to be able to cancel in all cases for the phase - this is why 
  // there is this if(a[n] < a[n+1]) thing.

  rsRandomUniform(0.0, 1.0, seed);
  rsVectorDbl p(a.dim);
  double k = 0.0;    // target value for the sample at n = 0
  int    N = p.dim;  // upper limit for loop below
  if( rsIsOdd(p.dim) )
  {
    if( a[p.dim-1] == 0.0 )
      p[p.dim-1] = 0.0;
    else
      p[p.dim-1] = asin(k/a[p.dim-1]);
    N--;
  }
  for(int n = 0; n < N; n += 2)
  {
    if( fabs(a[n]) == 0.0 && fabs(a[n+1]) == 0.0 )
      continue;
    if( fabs(a[n]) < fabs(a[n+1]) )
    {
      p[n]   = rsRandomUniform(0.0, 2*PI*randomness);
      p[n+1] = asin((k-a[n]*sin(p[n]))/a[n+1]);
    }
    else // roles reversed
    {
      p[n+1] = rsRandomUniform(0.0, 2*PI*randomness);
      p[n]   = asin((k-a[n+1]*sin(p[n+1]))/a[n]);      
    }
  }
  rsScale(p.v, p.dim, 360.0/(2*PI));
  return p;
}

rsVectorDbl rsModalFilterBank::modeDecayTimes(rsVectorDbl f, double fc, double p)
{
  rsVectorDbl d(f.dim);
  for(int n = 0; n < d.dim; n++)
    d[n] = modeDecayTime(f[n], fc, p);
  return d;
}

rsVectorDbl rsModalFilterBank::scaleAtIntervals(rsVectorDbl v, 
                                                int startIndex, 
                                                int interval, 
                                                double scaler)
{
  rsVectorDbl r = v;
  for(int n = startIndex; n < r.dim; n += interval)
    r[n] *= scaler;
  return r;
}

