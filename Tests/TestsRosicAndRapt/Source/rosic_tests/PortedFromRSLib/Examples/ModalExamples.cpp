#include "ModalExamples.h"

// try to get rid of these prototype declarattions:
std::vector<double> rsLinearRangeVector(     int N, double min, double max);
//std::vector<double> rsExponentialRangeVector(int N, double min, double max);
std::vector<double> rsRandomVector(          int N, double min, double max, int seed = 0);
std::vector<double> rsApplyFunction(const std::vector<double>& v, double p, 
  double (*f) (double, double));


// convenience functions (should go to and included from TestUtilities or something):
/*
std::vector<double> rsApplyFunction(const std::vector<double>& x, double p, 
  double (*f) (double, double))
{
  std::vector<double> y(x.size());
  for(int i = 0; i < x.size(); i++)
    y[i] = f(x[i], p);
  return y;
}
*/

// \todo: let this function accept a relative path:
template<class T>
void writeImpulseResponseToFile(const char *path, T &module, int N, int fs, int numBits)
{  
  double *h = new double[N];  
  getImpulseResponse(module, h, N);
  rosic::writeToMonoWaveFile(path, h, N, fs, numBits);
  delete[] h;
}




void createModalFilterExamples()
{
  std::cout << "Rendering rsModalFilter example 1\n";

  // user parameters:
  int    N   = 88200;  // number of samples 
  double fs  = 44100;  // samplerate in Hz
  double td  = 0.2;    // decay time constant in seconds
  double f   = 55;     // frequency in Hz
  double phs = 0;      // phase in degrees
  double A   = 1.0;    // amplitude as raw factor
  double pm  = 8.0;    // phase-modulation

  // synthesize the signal:
  rsNonlinearModalFilterDD nmf;
  nmf.setModalParameters(f, A, td, phs, fs);
  nmf.setPhaseModulation(pm);
  //writeImpulseResponseToFile("d:\\TmpData\\ModalFilterTest1.wav", nmf, N, (int) fs, 16);
  writeImpulseResponseToFile("ModalFilterTest1.wav", nmf, N, (int) fs, 16); 


  // create a "noise-response":
  N = 3*N;
  nmf.reset();
  double *y = new double[N];  
  RAPT::rsRandomUniform(-1.0, 1.0, 0);
  for(int n = 0; n < N; n++)
    y[n] = nmf.getSample(RAPT::rsRandomUniform(-0.01, 0.01));
  rsNormalize(y, N, 1.0);
  //rosic::writeToMonoWaveFile("d:\\TmpData\\ModalFilterNoiseTest1.wav", y, N, (int) fs, 16);
  rosic::writeToMonoWaveFile("ModalFilterNoiseTest1.wav", y, N, (int) fs, 16);
}







// gives the relative mode decay time for mode with relative frequency f given a (relative)
// cutoff frequency fc and an ultimate slope of the decay-time function (with respect to f)
// determined by the exponent p
// d(f) = a / (b + (f/fc)^p) where a and b are adjusted such that d(f=1)=1 and d(f=fc)=1/2
// fc must be > 1
/*
double modeDecayTime(double f, double fc, double p)
{
  double k = pow(fc, -p);
  return (1-k) / ((1-k*2) + pow(f/fc, p));
}
*/

// Creates a vector of decay-times for the modes that have relative frequecies given in the 
// vector f. The rule for the decay time for a mode with relative frequency f is:
// d(f) = a / (b + (f/fc)^p) where the constants a and b are adjusted such that 
// d(f=1)=1 and d(f=fc)=1/2, -> fc > 1 must be satisfied (hmm - maybe fc != 1 suffices)
/*
rsVectorDbl modeDecayTimes(rsVectorDbl f, double fc, double p)
{
  rsVectorDbl d(f.dim);
  for(int n = 0; n < d.dim; n++)
    d[n] = modeDecayTime(f[n], fc, p);
  return d;
}
*/

/** Scales the values at a given interval by the given scaler starting at the given startIndex. 
This is useful for scaling the decay times of even or odd harmonics, for example. */
/*
rsVectorDbl scaleAtIntervals(rsVectorDbl v, int startIndex, int interval, double scaler)
{
  rsVectorDbl r = v;
  for(int n = startIndex; n < r.dim; n += interval)
    r[n] *= scaler;
  return r;
}
*/

// move into class rsModalParameterGenerator
double combAmplitude(double frequency, double notchDistance, double notchOffset = 0, 
                     double ampFloor = 0, double shape = 1)
{
  double c = PI/notchDistance;
  double a = fabs(sin(c*(frequency-notchOffset)));
  a  = pow(a, shape);
  a  = ampFloor + (1-ampFloor)*a;
  return a;
}

std::vector<double> applyCombWeighting(std::vector<double> v, std::vector<double> f, 
  double notchDistance, double notchOffset = 0, double ampFloor = 0, double shape = 1)
{
  for(size_t i = 0; i < f.size(); i++)
    v[i] = v[i] * combAmplitude(f[i], notchDistance, notchOffset, ampFloor, shape);
  return v;
}


std::vector<double> pseudoHarmonicRatios12TET(int numPartials) // only 0..20
{
  double tmp[21];
  //long double s = pow(2.0, 1.0/12.0); // basis
  double s = pow(2.0, 1.0/12.0); // basis

                         // #    ratio
  tmp[0]  = pow(s,  0);  //  1    1.0
  tmp[1]  = pow(s, 12);  //  2    2.0
  tmp[2]  = pow(s, 19);  //  3    2.9966141537533639
  tmp[3]  = 4.0;         //  4    4.0
  tmp[4]  = pow(s, 28);  //  5    5.0396841995794937
  tmp[5]  = pow(s, 31);  //  6    5.9932283075067279
  tmp[6]  = pow(s, 34);  //  7    7.1271897451227169
  tmp[7]  = 8.0;         //  8    8.0
  tmp[8]  = pow(s, 38);  //  9    8.9796963864749877
  tmp[9]  = pow(s, 40);  // 10   10.079368399158989
  tmp[10] = pow(s, 42);  // 11   11.313708498984765
  tmp[11] = pow(s, 43);  // 12   11.986456615013459
  tmp[12] = pow(s, 44);  // 13   12.699208415745600
  tmp[13] = pow(s, 46);  // 14   14.254379490245437
  tmp[14] = pow(s, 47);  // 15   15.101989002907104
  tmp[15] = 16.0;        // 16   16.0  
  tmp[16] = pow(s, 49);  // 17   16.951409509748732
  tmp[17] = pow(s, 50);  // 18   17.959392772949979
  tmp[18] = pow(s, 51);  // 19   19.027313840043551
  tmp[19] = pow(s, 52);  // 20   20.158736798317982
  tmp[20] = pow(s, 53);  // 21   21.357437666720561

  std::vector<double> r = rsLinearRangeVector(numPartials, 1.0, numPartials); // relative frequencies
  for(int n = 0; n < RAPT::rsMin(numPartials, 21); n++)
    r[n] = tmp[n];
  return r;
}

std::vector<double> ratios12TET(int numPartials)
{  
  double tmp[21];
  std::vector<double> r(numPartials);
  //long double s = pow(2.0, 1.0/12.0); // basis
  double s = pow(2.0, 1.0/12.0); // basis
                         //  #    ratio
  tmp[0]  = pow(s,  0);  //  1    1.0
  tmp[1]  = pow(s, 12);  //  2    2.0
  tmp[2]  = pow(s, 19);  //  3    2.9966141537533639
  tmp[3]  = 4.0;         //  4    4.0
  tmp[4]  = pow(s, 28);  //  5    5.0396841995794937
  tmp[5]  = pow(s, 31);  //  6    5.9932283075067279
  tmp[6]  = pow(s, 34);  //  7    7.1271897451227169
  tmp[7]  = 8.0;         //  8    8.0
  tmp[8]  = pow(s, 38);  //  9    8.9796963864749877
  tmp[9]  = pow(s, 40);  // 10   10.079368399158989
  tmp[10] = pow(s, 42);  // 11   11.313708498984765
  tmp[11] = pow(s, 43);  // 12   11.986456615013459
  tmp[12] = pow(s, 44);  // 13   12.699208415745600
  tmp[13] = pow(s, 46);  // 14   14.254379490245437
  tmp[14] = pow(s, 47);  // 15   15.101989002907104
  tmp[15] = 16.0;        // 16   16.0  
  tmp[16] = pow(s, 49);  // 17   16.951409509748732
  tmp[17] = pow(s, 50);  // 18   17.959392772949979
  tmp[18] = pow(s, 51);  // 19   19.027313840043551
  tmp[19] = pow(s, 52);  // 20   20.158736798317982
  tmp[20] = pow(s, 53);  // 21   21.357437666720561
  // get rid of this code-duplication

  for(int n = 0; n < RAPT::rsMin(numPartials, 21); n++)
    r[n] = tmp[n];

  for(int n = 21; n < numPartials; n++)
    r[n] = pow(s, 54+n-21);

  // maybe use a formula:
  // b  = pow(2.0, 1.0/12.0); // basis - can be generalized
  // fn = b^k where k = round(logB(n*f0, b)...or maybe use a special rounding:

  // roundedExponent(target, basis)
  // a  = logB(target, basis);
  // af = floor(a);
  // ac = ceil(a);
  // yf = pow(basis, yf);
  // yc = pow(basis, yc);
  // if(fabs(target-yc) < fabs(target-yf))
  //   return ac;
  // else
  //   return af;


  //double dbg[200];
  //rsCopyBuffer(r.v, dbg, 200);

  return r;
}

std::vector<double> stiffStringRatios(double frequency, double sampleRate, double inharmonicity)
{
  double B = inharmonicity;

  static const int maxNumPartials = 1000;
  double tmp[maxNumPartials];

  int i, n;
  int numPartials = 0;
  double s = 1.0 / sqrt(1+B); // to scale 1st entry to unity
  for(i = 0; i < maxNumPartials; i++)
  {
    n      = i+1;
    tmp[i] = s * n * sqrt(1+B*n*n);
    if( frequency*tmp[i] >= sampleRate/2 )
      break;
    numPartials++;
  }

  //rsVectorDbl v(numPartials, tmp);
  std::vector<double> v(numPartials);
  RAPT::rsArrayTools::copy(tmp, &v[0], numPartials);

  return v;
}
// function superseded by rsModalParameterGenerator::getFrequencies

/*
struct rsModalBankParameters
{
  double frequency;    // reference frequency in Hz
  double attack;       // attack time in seconds
  double decay;        // decay time in seconds

  rsVectorDbl f;       // relative frequencies
  rsVectorDbl g;       // gains
  rsVectorDbl d;       // relative decay times
  rsVectorDbl a;       // relative attack times
  rsVectorDbl p;       // start-phases
};
*/

/*
rsModalBankParameters modalParametersGuitar55Hz()
{
  int numPartials = 400;  
  rsModalBankParameters p;

  p.frequency = 55.0;
  p.attack    = 0.1;
  p.decay     = 1.0;

  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.d = rsModalFilterBank::modeDecayTimes(p.f, 4.0, 0.95);
  p.a = p.d;
  p.p = rsRandomVector(numPartials, 0.0, 2*PI, 0);

  return p;
}
rsModalBankParameters modalParametersGuitar110Hz()
{
  int numPartials = 200;  
  rsModalBankParameters p;

  p.frequency = 110.0;
  p.attack    = 0.02;
  p.decay     = 0.9;

  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.d = rsModalFilterBank::modeDecayTimes(p.f, 2.5, 0.93);
  p.a = p.d;
  p.p = rsRandomVector(numPartials, 0.0, 2*PI, 0);

  return p;
}
*/

rsModalBankParametersD modalParametersOrganBass55Hz()
{
  int numPartials = 400;  
  rsModalBankParametersD p;

  p.frequency = 55.0;
  p.attack    = 0.15;
  p.decay     = 0.5;

  p.f = pseudoHarmonicRatios12TET(numPartials);
  p.g = rsApplyFunction(p.f, -0.5, &pow);
  p.g[1] = p.g[3] = p.g[7] = 1.0;
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 3.0, 1.8); 
  p.a = p.d;
  p.d = rsModalFilterBankDD::scaleAtIntervals(p.d, 3, 3, 0.5);
  p.p = rsRandomVector(numPartials, 0.0, 2*PI, 0);  

  return p;
}


rsModalBankParametersD modalParametersPiano110Hz()
{
  rsModalBankParametersD p;

  p.frequency = 110.0;
  p.attack    = 0.02;
  p.decay     = 1.2;

  p.f = stiffStringRatios(110.0, 44100.0, 0.00075);



  // preliminary (from guitar):

  p.g = rsModalFilterBankDD::modeDecayTimes(p.f, 35.0, 3.0);


  p.g = applyCombWeighting(p.g, p.f, 7);
  //p.g = applyCombWeighting(p.g, p.f, 10.0/3.0);
  //p.g = applyCombWeighting(p.g, p.f, 10.0/1.0);


  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 35.0, 1.0);
  //p.d = linearRangeVector(p.f.dim, 1.0, 1.0);

  p.a = p.d;
  //p.a = linearRangeVector(p.f.dim, 0.1, 0.1);
  p.p = rsRandomVector((int)p.f.size(), 0.0, 60.0, 0);



  // we need comb-filter effects in the amplitude array - try a comb with a zero at DC and the
  // next zero on the 7th partial


  // it sounds somewhat piano-like, but still a bit "buzzy", especially in the attack

  return p;
}


/*
Sounddesign Notes:
-with stretched partial ratios and sligned intial phases, zap-like transients can be created
-also, with attacks that decrease with frequency (i.e. high partial reach their peak amplitude 
 earlier)

*/

void createModalFilterBankExamples()
{
  std::cout << "Rendering rsModalFilterBank examples...\n";

  // global parameters:
  double sampleRate      = 44100.0;
  double phaseRandomness =     1.0;
  double truncationLevel =   -60.0;
  double fadeCycles      =     5.0;


  //truncationLevel = -10.0;  // for preview
  //truncationLevel = -40.0;  // for preview


  /*
  // gong'ish - but not quite:
  f = pseudoHarmonicRatios12TET(numPartials);
  g = applyFunction(f, -0.7,  &pow);
  d = rsModalFilterBank::modeDecayTimes(f, 20.0, 0.95);
  a = d;                               
  a[1] *= 0.25;
  a[3] *= 0.25;
  */

  // it seems like the modal synth doesn't admit zero mode amplitude? check this out...

  rsModalBankParametersD p;

  //p = modalParametersGuitar55Hz();
  //p = modalParametersGuitar110Hz();
  p = modalParametersOrganBass55Hz();
  //p = modalParametersPiano110Hz();


  rsModalFilterBankDD mfb;
  mfb.setSampleRate(sampleRate);
  mfb.setModalParameters(p.f, p.g, p.a, p.d, p.p);
  mfb.setReferenceFrequency(p.frequency);
  mfb.setReferenceAttack(p.attack);
  mfb.setReferenceDecay(p.decay);

  int    numSamples     = (int) (mfb.getLength(truncationLevel)*sampleRate);
  double fadeOutTime    = fadeCycles/p.frequency;
  int    numFadeSamples = (int) (sampleRate * fadeOutTime);
  numSamples += numFadeSamples;

  mfb.reset();
  double *y = new double[numSamples];  
  y[0] = mfb.getSample(1.0);
  for(int n = 1; n < numSamples; n++)
    y[n] = mfb.getSample(0.0);


  //scale(y, numSamples, 0.25);
  rsNormalize(y, numSamples, 1.0);
  RAPT::rsFadeOut(y, numSamples-numFadeSamples-1, numSamples-1);
  //rosic::writeToMonoWaveFile("d:\\TmpData\\ModalSynthTest.wav", y, numSamples, (int) sampleRate, 16);
  rosic::writeToMonoWaveFile("ModalSynthTest.wav", y, numSamples, (int) sampleRate, 16);
}

template<class T>
rsModalParameterGenerator<T>::rsModalParameterGenerator()
{
  tmp.resize(maxNumPartials);
}

template<class T>
rsModalBankParameters<T> rsModalParameterGenerator<T>::getModalParameters()
{
  rsModalBankParameters<T> mp;

  mp.frequency = frequency;
  mp.gain      = amplitude;
  mp.attack    = attackTime;
  mp.decay     = decayTime;

  getFrequencies(mp.f);
  getAmplitudes( mp.g, mp.f);
  getPhases(     mp.p, mp.f);
  getDecayTimes( mp.d, mp.f);
  getAttackTimes(mp.a, mp.f, mp.d);

  return mp;
}


template<class T>
void rsModalParameterGenerator<T>::getFrequencies(std::vector<T>& f)
{
  int numPartials = 0;

  // for harmonic partials or stiff-string inharmonicity (todo: make a switch to allow for other
  // partial frequency ratios):
  T B = inharmonicity;
  T s = 1.0 / sqrt(1+B); // to scale 1st entry to unity
  for(int i = 0; i < maxNumPartials; i++) {
    int n = i+1;
    tmp[i] = s * n * sqrt(1+B*n*n);
    if(frequency*tmp[i] >= sampleRate/2)
      break;
    numPartials++;
  }

  f.resize(numPartials);
  RAPT::rsArrayTools::copy(&tmp[0], &f[0], numPartials);
}

template<class T>
void rsModalParameterGenerator<T>::getPhases(std::vector<T>& p, const std::vector<T>& f)
{
  p.resize(f.size());
  prng.setSeed(phaseRandomSeed);
  //prng.setOrder(phaseRandomShape);
  prng.setRange(-PI*phaseRandomness, +PI*phaseRandomness);
  for(size_t i = 0; i < p.size(); i++) {
    p[i] = prng.getSample();
  }
}

template<class T>
T rsModalParameterGenerator<T>::modeDecayTime(T f, T fc, T p)
{
  if(p == 0.0) 
    return 1.0;
  T k = pow(fc, -p);
  if(k == 1.0)
    return 1.0;
  return (1-k) / ((1-k*2) + pow(f/fc, p));
}

template<class T>
T rsModalParameterGenerator<T>::combAmplitude(T frequency, T notchDistance, T amount,
  T notchOffset, T shape)
{
  T c = PI/notchDistance;
  T a = fabs(sin(c*(frequency-notchOffset)));
  a = pow(a, shape);
  a = (1-amount) + amount*a;
  return a;
}

template<class T>
void rsModalParameterGenerator<T>::getAmplitudes(std::vector<T>& a, const std::vector<T>& freqs)
{
  a.resize(freqs.size());
  for(size_t i = 0; i < a.size(); i++) {
    T f   = freqs[i];
    a[i]  = amplitude * pow(f, -ampSlope1);
    a[i] *= modeDecayTime(f, ampCutoff, ampSlope2);  // rename this function
    a[i] *= combAmplitude(f, ampCombHarmonic, ampCombAmount);
    if(RAPT::rsIsEven(i+1))
      a[i] *= evenAmpScale;
  }
}

template<class T>
void rsModalParameterGenerator<T>::getDecayTimes(std::vector<T>& d, const std::vector<T>& freqs)
{
  d.resize(freqs.size());
  for(size_t i = 0; i < d.size(); i++) 
  {
    T f = freqs[i];
    d[i]  = pow(f, -decaySlope1);
    d[i] *= modeDecayTime(f, decayCutoff, decaySlope2);
    d[i] *= combAmplitude(f, decayCombHarmonic, decayCombAmount);
    if(RAPT::rsIsEven(i+1))
      d[i] *= evenDecayScale;
  }
}

template<class T>
void rsModalParameterGenerator<T>::getAttackTimes(std::vector<T>& a, 
  const std::vector<T>& f, const std::vector<T>& d)
{
  a.resize(f.size());

  // attack should perhaps be a function of frequency plus some function of decay
  // attack = a*func(freq) + b*func(decay) maybe a*freq^p + b*decay^q or maybe
  // a*(1/freq)^p + b*decay^q

  // preliminary - make them adjustable members:
  T attackFreqCoeff  = 0.0;
  T attackFreqPower  = 1.0;
  T attackDecayCoeff = 0.1;
  T attackDecayPower = 1.0;

  for(size_t i = 0; i < a.size(); i++) {
    a[i]  = attackFreqCoeff  * pow(f[i], -attackFreqPower);
    a[i] += attackDecayCoeff * pow(d[i],  attackDecayPower);
    a[i]  = RAPT::rsMin(a[i], attackDecayRatioLimit*d[i]);
  }
}

void createPiano1()
{
  // Trying to create a piano like sound, features:
  // -partials are slightly inharmonic
  // -even partials decay at different rate than odd nd may have different amplitude
  // -harmonics expose comb-filter like profile

  rsModalParameterGenerator<double> mpg;

  // frequency:
  mpg.setSampleRate(44100);
  mpg.setFrequency(100);
  mpg.setInharmonicity(0.0);

  // amplitude:
  mpg.setAmplitude(1.0);
  mpg.setAmpSlope1(1.0);
  mpg.setAmpCutoff(2.0);
  mpg.setAmpSlope2(0.0);
  mpg.setEvenAmpScale(1.0);
  //mpg.setAmpCombHarmonic(7.0);
  //mpg.setAmpCombAmount(1.0);

  // phase:
  mpg.setPhaseRandomness(1.0);
  mpg.setPhaseRandomSeed(0);

  // decay;
  mpg.setDecay(0.5);
  mpg.setDecaySlope1(0.5);
  mpg.setDecayCutoff(2.0);
  mpg.setDecaySlope2(0.0);
  mpg.setEvenDecayScale(0.5);
  //mpg.setDecayCombHarmonic(7.0);
  //mpg.setDecayCombAmount(0.5);

  // attack:
  mpg.setAttack(0.1);

  rsModalBankParametersD mp = mpg.getModalParameters();

  // plot:
  //GNUPlotter plt;

  // ...

  int dummy = 0;
}




void createBass1()
{
  SampleMapGeneratorModal g;
  g.setName("Bass1");

  int numPartials = 800; 
  rsModalBankParametersD p;
  
  // A1, 55 Hz (key 33)
  p.frequency = 55.0;
  p.attack    =  0.15;
  p.decay     =  0.5;
  p.f = pseudoHarmonicRatios12TET(numPartials);
  p.g = rsApplyFunction(p.f, -0.5,  &pow);
  p.g[1] = p.g[3] = p.g[7] = 1.0;
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 3.0, 1.8); 
  p.a = p.d;
  p.d = rsModalFilterBankDD::scaleAtIntervals(p.d, 3, 3, 0.5);
  p.p = rsRandomVector(numPartials, 0.0, 2*PI, 0);    
  g.setModalParametersForKey(33, p);

  // A2, 110 Hz (key 45)
  p.frequency = 110.0;
  p.attack    =   0.07;
  p.decay     =   0.45;
  p.f = pseudoHarmonicRatios12TET(numPartials);
  p.g = rsApplyFunction(p.f, -0.5,  &pow);
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 2.3, 2.2);
  p.a = p.d;
  p.d = rsModalFilterBankDD::scaleAtIntervals(p.d, 3, 3, 0.5);
  p.p = rsRandomVector(numPartials, 0.0, 2*PI, 0);    
  g.setModalParametersForKey(45, p);

  // for preview and sounddesign:
  //g.setKeyRangeToRender(21, 21);
  //g.setKeyRangeToRender(33, 33);
  //g.setKeyRangeToRender(45, 45);
  //g.setKeyRangeToRender(57, 57);
  //g.setKeyRangeToRender(69, 69);
  //g.setKeyRangeToRender(81, 81);
  //g.setKeyRangeToRender(93, 93);

  //g.setKeyRangeToRender(45, 46);
  g.setKeyRangeToRender(33, 45);

  //g.setTruncationLevel(-10);
  //g.setTruncationLevel(-1);
  //g.setTruncationLevel(-3);
  //g.setTruncationLevel(-20);

  g.generateSampleMap(true);
}


void createGong1()
{
  SampleMapGeneratorModal g;
  g.setName("Gong1");

  int numPartials = 71;    // 71: number of modes for 55 Hz
  rsModalBankParametersD p;
  

  // A1, 55 Hz (key 33)
  numPartials = 71;
  p.frequency = 55.0;
  p.attack    =  0.5;
  p.decay     =  0.8;
  p.f = ratios12TET(numPartials);
  p.g = rsModalFilterBankDD::modeDecayTimes(p.f, 16.0, 3.0);
  p.g[1] *= 0.5;
  p.a = std::vector<double>(numPartials);
  rsFillWithRangeExponential(&p.a[0],  16,              0.05, 0.5);
  rsFillWithRangeExponential(&p.a[16], numPartials-16,  0.5,  0.6);
  p.d = std::vector<double>(numPartials);
  rsFillWithRangeExponential(&p.d[0],  16,              1.5, 1.0);
  rsFillWithRangeExponential(&p.d[16], numPartials-16,  1.0, 0.8);
  p.p = rsRandomVector(numPartials, 0.0, 2*PI, 0);    
  g.setModalParametersForKey(33, p);


  // for preview and sounddesign:
  //g.setKeyRangeToRender(21, 21);
  g.setKeyRangeToRender(33, 33);
  //g.setKeyRangeToRender(45, 45);
  //g.setKeyRangeToRender(57, 57);
  //g.setKeyRangeToRender(69, 69);
  //g.setKeyRangeToRender(81, 81);
  //g.setKeyRangeToRender(93, 93);

  //g.setKeyRangeToRender(45, 46);
  //g.setKeyRangeToRender(33, 45);

  //g.setTruncationLevel(-10);
  //g.setTruncationLevel(-1);
  //g.setTruncationLevel(-3);
  //g.setTruncationLevel(-20);

  g.generateSampleMap(true);
}



void createPluck1()
{
  SampleMapGeneratorModal g;
  g.setName("Pluck1");

  int numPartials = 800; 
  rsModalBankParametersD p;

  double minPhase   = 0.0;
  double maxPhase   = 360.0;
  int    randomSeed = 0;      // for the randomized phases (standard: 0)

  // A0, 22.5 Hz (key 21)
  p.frequency = 22.5;  // will be overwritten anyway
  p.gain      =  1.0;
  p.attack    =  0.1;
  p.decay     =  1.5;
  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.g[0] *= 0.25;    // why? reduce subsonic fundamental?
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 5.0, 0.95);
  p.a = p.d;
  p.p = rsRandomVector(numPartials, minPhase, maxPhase, randomSeed);
  g.setModalParametersForKey(21, p);

  // A1, 55 Hz (key 33)
  p.frequency = 55.0;  // will be overwritten anyway
  p.gain      =  1.0;
  p.attack    =  0.1;
  p.decay     =  1.0;
  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 4.0, 0.95);
  p.a = p.d;
  p.p = rsRandomVector(numPartials, minPhase, maxPhase, randomSeed);
  g.setModalParametersForKey(33, p);

  // A2, 110 Hz (key 45): 
  p.frequency = 110.0;  // will be overwritten anyway
  p.gain      =   1.0;  // test - later set back to 1.0
  p.attack    =   0.02;
  p.decay     =   0.9;
  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 2.5, 0.93);
  p.a = p.d;
  p.p = rsRandomVector(numPartials, minPhase, maxPhase, randomSeed);
  g.setModalParametersForKey(45, p);

  // A3, 220 Hz (key 57): 
  p.frequency = 220.0;  // will be overwritten anyway
  p.gain      =   1.0;  // test - later set back to 1.0
  p.attack    =   0.02;
  p.decay     =   0.7;
  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 2.7, 0.80);  // before: 2.35, 0.85
  p.a = p.d;
  p.p = rsRandomVector(numPartials, minPhase, maxPhase, randomSeed);
  g.setModalParametersForKey(57, p);

  // A4, 440 Hz (key 69): 
  p.frequency = 440.0;  // will be overwritten anyway
  p.gain      =   1.0;  // test - later set back to 1.0
  p.attack    =   0.01;
  p.decay     =   0.5;
  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 2.9, 0.65);  // before: 2.2, .82
  p.a = p.d;
  p.p = rsRandomVector(numPartials, minPhase, maxPhase, randomSeed);
  g.setModalParametersForKey(69, p);

  // A5, 880 Hz (key 81): 
  p.frequency = 880.0;  // will be overwritten anyway
  p.gain      =   1.0;  // test - later set back to 1.0
  //p.attack    =   0.01;
  p.attack    =   0.02;
  //p.decay     =   0.4;
  p.decay     =   0.2;
  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  //p.g = rsApplyFunction(p.f, -1.0,  &pow);
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 2.4, 0.60); // before 2.0, 0.8
  //p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 2.0, 0.80); 
  p.a = p.d;
  p.p = rsRandomVector(numPartials, minPhase, maxPhase, randomSeed);
  g.setModalParametersForKey(81, p);

  // A6, 1760 Hz (key 93): 
  p.frequency = 1760.0;  // will be overwritten anyway
  p.gain      =    1.0;  // test - later set back to 1.0
  p.attack    =    0.01;
  p.decay     =    0.3;
  p.f = rsLinearRangeVector(numPartials, 1.0, numPartials);
  p.g = rsApplyFunction(p.f, -0.7,  &pow);
  p.d = rsModalFilterBankDD::modeDecayTimes(p.f, 2.0, 0.50); // before: 1.8, 0.75
  p.a = p.d;
  p.p = rsRandomVector(numPartials, minPhase, maxPhase, randomSeed);
  g.setModalParametersForKey(93, p);

  g.setKeyRangeToRender(80, 80);  
  // For production rendering, use 21..93. For preview and sound-design, use smaller range
  // ToDo: maybe use an increment (default 1, 12 means one sample per octave - which is also good 
  // for preview), see SampleMapGenerator::generateAllSamples

  g.setTruncationLevel(-40);
  // For production rendering, use -80 or -60. For preview -40


  g.generateSampleMap(true);

  // -the higher notes sound rather synthetic, maybe we need to add some transient sample that has
  //  frequency content below the fundamental?
  // -and/or maybe mix in a sample an octave or 2 octaves lower to create subharmonics and partials
  //  in between the actual ones - maybe it should have a fast decay to affect only the transient
  //  section
  // -maybe try a little bit of inharmonicity
  // -maybe try mixing in a karplus-strong string sound

  // ToDo:
  // -provide a function g.plotParameters - this should plot the modal parameters as function of 
  //  the key
  // -provide load/save functions for the settings
  // -provide functions to split samples into harmonic, inharmonic, transient, noise (maybe split
  //  harmonic part further into even/odd)
  // -these different parts should go into different groups in the sfz
  // -maybe we should render in 96 kHz and use key-crossfades
  // -Try to render the upper keys at a higher sample-rate with supersonic harmonics. Maybe that 
  //  helps against the synthetic sound. It's plausible that even when downsampling with AA 
  //  lowpass, the transient portion may be different, because the supersonic harmonics decay so 
  //  quickly that they supply a broadband signal at the start of the sample (i think). ..that may
  //  generally an interesting way to render transients: oversampled, supersonic, quickly decaying 
  //  modes - then downsampled with AA filter.
}

void testHighPluck()
{
  // We experiment with the creation of a rather high pitched note for a guitarish plucked string 
  // like sound. The high notes are problematic in the redering above (they sound synthetic), so 
  // here we specifically experiment with ideas trying to fix this...

  double sampleRate  = 44100;
  double fundamental = 1000;   // fundamental frequency of the string
  double minPhase    = 0.0;
  double maxPhase    = 360.0;
  double attack      = 0.02;
  double decay       = 0.3;
  double length      = 3.0;    // in seconds
  double attackFac   = 0.02;   // factor for the attack for the transient
  double decayFac    = 0.10;   // factor for the decay for the transient  ..0.08
  int    randomSeed  = 0;      // for the randomized phases (standard: 0)

  using Vec = std::vector<double>;
  using AT  = RAPT::rsArrayTools;
  using MFB = rsModalFilterBankDD;

  // Create the main signal:
  int numPartials = (int) floor(0.5*sampleRate / fundamental);
  int numSamples  = (int) ceil(sampleRate * length);
  int N = numPartials;
  Vec frq(N), amp(N), att(N), dec(N), phs(N);
  frq = rsLinearRangeVector(N, 1.0, N);
  amp = rsApplyFunction(frq, -0.7,  &pow);          // use variable ampExponent
  dec = MFB::modeDecayTimes(frq, 2.4, 0.30); // use variables
  att = dec;
  phs = rsRandomVector(N, minPhase, maxPhase, randomSeed);
  MFB mfb;
  mfb.setReferenceFrequency(fundamental);
  mfb.setReferenceAttack(attack);
  mfb.setReferenceDecay(decay);
  mfb.setModalParameters(frq, amp, att, dec, phs);
  N = numSamples;
  Vec x1(N);  // the main signal
  x1[0] = mfb.getSample(1.0);
  for(int n = 1; n < N; n++)
    x1[n] = mfb.getSample(0.0);

  // Create a signal with a subharmonic frequency
  int D = 3;  // divider
  double baseFreq = fundamental / D;
  numPartials = (int) floor(0.5*sampleRate / baseFreq);
  N = numPartials;
  frq = rsLinearRangeVector(N, 1.0, N);
  amp = rsApplyFunction(frq, -0.0, &pow);     // use variable ampExponent2
  for(int i = 0; i < D; i++)                  // highpass
    amp[i] = 0.0;   
  for(int i = D-1; i < N; i += D)             // comb
    amp[i] *= 0.0;
  dec = MFB::modeDecayTimes(frq, 20.0, 0.1);  // use variables
  att = dec;
  phs = rsRandomVector(N, minPhase, maxPhase, randomSeed);
  mfb.setReferenceFrequency(baseFreq);
  mfb.setReferenceAttack(attackFac*attack); 
  mfb.setReferenceDecay(decayFac*decay); 
  mfb.setModalParameters(frq, amp, att, dec, phs);
  N = numSamples;
  Vec x2(N); 
  mfb.reset();
  x2[0] = mfb.getSample(1.0);
  for(int n = 1; n < N; n++)
    x2[n] = mfb.getSample(0.0);

  AT::normalize(&x1[0], N, 0.5);
  rosic::writeToMonoWaveFile("HighPluckMain.wav", &x1[0], N, (int)sampleRate);
  AT::normalize(&x2[0], N, 0.5);
  rosic::writeToMonoWaveFile("HighPluckSub.wav", &x2[0], N, (int)sampleRate);
  Vec mix = x1 + x2;
  rosic::writeToMonoWaveFile("HighPluckMix1.wav", &mix[0], N, (int)sampleRate);
  mix = x1 - x2;
  rosic::writeToMonoWaveFile("HighPluckMix2.wav", &mix[0], N, (int)sampleRate);

  // amp-modulation:
  mix = x1 * (1.0*x2 + 1.0); 
  rosic::writeToMonoWaveFile("HighPluckAmpMod.wav", &mix[0], N, (int)sampleRate);
  // todo: avoid aliasing by oversampling by 2: upsample by 2 using sinc interpolation, do the 
  // amp-mod, downsample by 2 again using sinc interpolation...wrap that into a function 
  // rsAmpModulate(const T* carrier, const T* modulator, int N, T* result, 
  //   int modIndex = 1, int sincLength = 512)...or maybe that should just implement ring-mod, 
  // amp-mod can be obtained via original + ringmod, maybe for the ring-mod signal, use a single
  // mode (or a few) all below the fundamental
  // 


  // Test: write to 24 bit wavefile (todo: make a unit test for this and remove the code here):
  //rsWaveFile

  // int24
  RSLib::rsOutputWaveFile wavFile("HighPluckMain24Bit.wav", (int)sampleRate, 24, 1);

  // float32:
  //rsOutputWaveFile wavFile("HighPluckMain32BitFloat.wav", (int)sampleRate, 32, 1, 
  //  rsWaveFile::SampleFormat::IEEE_FLOAT);

  //// float64:
  //rsOutputWaveFile wavFile("HighPluckMain64BitFloat.wav", (int)sampleRate, 64, 1, 
  //  rsWaveFile::SampleFormat::IEEE_FLOAT);

  wavFile.write(&x1[0], N);



  printf("%s", "Rendering HighPluck*.wav done\n");


  // Observations:
  // -when increasing the decayFac, the transient of the mix becomes more dirty/gritty. the 
  //  transient of the amp-modulated becomes just longer
  // - a decayFac of 0.08 ...0.1


  // ToDo:
  // -render transient samples with dividers 2,3,4,5,6,7,8 and use them as mix-in samples that can
  //  be layered in the sfz
  // -they could be generally useful - not only for the Pluck
  // -when rendering sample-sets, maybe don't normalize each sample seperately, instead, normalize
  //  the whole sample set by the same factor - it's inconvenient to have settings like 
  //  -7.353135dB in the sfz file
  // -i think, divider=3 works best
  // -it may actually also be useful for lower keys/frequencies
  // -it may also help to mix in a second transient with a little delay...and a third
  //  ...maybe some sort of feedback-echo or resonant comb-filter effect could be useful
  // -maybe render some sort of broadband bandpass signals that can be combined in various ways
  //  to be added to the signal. for example, a butterworth from 1100-1900, 2100-2900, 3100-3900, 
  //  ...
  // -maybe amplitude-modulate main signal by the transient instead of just adding it
  //  ...done - gives the transient more of a piano-like transient...todo:
  //  subtract the original from the amp-modulated to obtain the amp-mod transient for mixing
  //  in...but: there is aliasing taking place...maybe for amp-mod, we should oversample by 2 to 
  //  avoid it

}
// ToDo: make it possible to write such rendering scripts in python

void insertionSortSound(double *cycle, int cycleLength, double *signal, int signalLength)
{ 
  rsFillWithZeros(signal, signalLength);
  rsCopyBuffer(cycle, signal, cycleLength);
  for(int i = 1; i < cycleLength; i++)
  {
    double tmp = cycle[i];
    int    j   = i-1;
    while( j >= 0 && tmp < cycle[j] )
      cycle[j+1] = cycle[j--];
    cycle[j+1] = tmp;

    if( (i+1)*cycleLength > signalLength )
      break;
    rsCopyBuffer(cycle, &signal[i*cycleLength], cycleLength);
  }
}
void insertionSortBackwardSound(double *cycle, int cycleLength, double *signal, int signalLength)
{ 
  rsFillWithZeros(signal, signalLength);
  rsCopyBuffer(cycle, signal, cycleLength);  
  int k = 0;
  for(int i = cycleLength-2; i >= 0; i--)
  {
    double tmp = cycle[i];
    int    j   = i;
    while( tmp > cycle[j+1] && j+1 < cycleLength )
    {
      cycle[j] = cycle[j+1];
      j++;
    }
    cycle[j] = tmp;
    k++;

    if( (k+1)*cycleLength > signalLength )
      break;
    rsCopyBuffer(cycle, &signal[k*cycleLength], cycleLength);
  }
}
void selectionSortSound(double *cycle, int cycleLength, double *signal, int signalLength)
{ 
  rsFillWithZeros(signal, signalLength);
  rsCopyBuffer(cycle, signal, cycleLength);  
  for(int i = 0; i < cycleLength-1; i++)
  {
    double tmp = cycle[i];
    int    k   = i;
    for(int j = i+1; j < cycleLength; j++)
    {
      if( cycle[j] < tmp )
      {
        k   = j;
        tmp = cycle[j];
      }
    }
    cycle[k] = cycle[i];
    cycle[i] = tmp;

    if( (i+2)*cycleLength > signalLength )
      break;
    rsCopyBuffer(cycle, &signal[(i+1)*cycleLength], cycleLength);
  }
}
void selectionSortBackwardSound(double *cycle, int cycleLength, double *signal, int signalLength)
{ 
  rsFillWithZeros(signal, signalLength);
  rsCopyBuffer(cycle, signal, cycleLength);  
  int ii = 0;
  for(int i = cycleLength-1; i > 0; i--)
  {
    double tmp = cycle[i];
    int    k   = i;
    for(int j = i-1; j >= 0; j--)
    {
      if( cycle[j] > tmp )
      {
        k   = j;
        tmp = cycle[j];
      }
    }
    cycle[k] = cycle[i];
    cycle[i] = tmp;
    ii++;

    if( (ii+1)*cycleLength > signalLength )
      break;
    rsCopyBuffer(cycle, &signal[ii*cycleLength], cycleLength);
  }
}


void createInsertionSortSound()
{
  static const int cycleLength = 512;
  double cycle[cycleLength];

  int numCycles = cycleLength;
  int signalLength = numCycles*cycleLength;
  double *signal = new double[signalLength];
  rsFillWithZeros(signal, signalLength);

  // initialize cycle with random numbers (later, we can use other intializations):
  rsFillWithRandomValues(cycle, cycleLength, -1.0, +1.0, 3); 
    // good seeds (for insertion-sort, at least): 0,3,4,5,6

  //rsCumulativeSum(cycle, cycleLength, 1, true); // old
  rsCumulativeSum(cycle, cycle, cycleLength);

  rsRemoveMean(cycle, cycleLength);
  rsNormalize(cycle, cycleLength, 1.0);


  //insertionSortSound(cycle, cycleLength, signal, signalLength);
  //insertionSortBackwardSound(cycle, cycleLength, signal, signalLength);
  selectionSortSound(cycle, cycleLength, signal, signalLength);
  //selectionSortBackwardSound(cycle, cycleLength, signal, signalLength);

  rosic::writeToMonoWaveFile("SelectionSortSound.wav", signal, signalLength, 44100, 16);

  delete[] signal;
}



/*
void createBubbleSortSound()
{
  static const int cycleLength = 1000;
  double cycle[cycleLength];

  int numCycles = cycleLength;
  int signalLength = numCycles*cycleLength;
  double *signal = new double[signalLength];
  rsFillWithZeros(signal, signalLength);

  // initialize cycle with random numbers:
  rsFillWithRandomValues(cycle, cycleLength, -1.0, +1.0, 0);

  // do the bubble sort - after each completion of the oute loop, a new cycle is written into the
  // output signal:
  int k = 0;
  for(int i = cycleLength-2; i >= 0; i--)
  {
    rsCopyBuffer(cycle, &signal[k*cycleLength], cycleLength);   
    double tmp = cycle[i];
    int j = i;
    while( tmp > cycle[j+1] && j+1 < cycleLength )
    {
      cycle[j] = cycle[j+1];
      j++;
    }
    cycle[j] = tmp;
    k++;
  }
  rsCopyBuffer(cycle, &signal[k*cycleLength], cycleLength);


  writeToMonoWaveFile("BubbleSortSound.wav", signal, signalLength, 44100, 16);
  delete[] signal;
}
*/


