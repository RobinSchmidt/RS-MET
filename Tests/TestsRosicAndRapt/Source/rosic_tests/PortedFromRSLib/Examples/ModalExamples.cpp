
//#include <rs_testing/rs_testing.h> // ToDo: get rid of that! organize the includes properly!

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
-with stretched partial ratios and aligned intial phases, zap-like transients can be created
-also, with attacks that decrease with frequency (i.e. high partials reach their peak amplitude 
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

void createSamplerWaveforms()
{
  // Create mip-mapped multisamples for single-cycle waveforms that can be used in the sampler 
  // engine. Most tables are of length 2048 except the bottom two which are of length 4096 and 
  // 8192 respectively. We want a cycle length that is a power of two and we want the samples
  // to have a rootKey that is exactly a midi note. All midi notes except for the As have 
  // irrational frequencies. So let's choose an A as rootKey. A4 at 440Hz is midi-key 69. When
  // we use key=21 (27.5 Hz) for the table of length 2048, we need to choose a sample-rate of 
  // 56320 to make it all work out. The formula for the frequency is: sampleRate/cycleLength, 
  // so we get 56320/2048 = 27.5

  // ToDo:
  // -Create prototype wave of length 8192
  // -Create various filtered versions moving average filters...the goal is to preserve the 
  //  time domain waveform as closely as possible ...maybe use a 
  //  filter -> decimate -> interpolate procedure for the higher tables
  //  ...maybe to avoid boudary artifacts from the MA filters, we should use 3 cycles for the
  //  filtering and then etract the middle one

  using Vec = std::vector<double>;
  using SWR = StandardWaveformRenderer;
  using namespace RAPT;
  using namespace rosic;

  std::string waveName = "Saw";  // name of the waveform as it appears in the .wav files

  // Render prototype sawtooth wave of length 8192 and obtain the first two decimated versions of
  // it by just suing averages of successive samples:
  Vec w8192(8192);
  SWR::renderSawWaveform(&w8192[0], 8192); //rsPlotVector(w8192);
  rsScale(w8192, 0.8125); //
  // Some nicely representable approximation to the scaling by the Wilbraham-Gibbs constant to 
  // avoid clipping due to Gibb's overshoot 
  // https://en.wikipedia.org/wiki/Gibbs_phenomenon
  // https://mathworld.wolfram.com/Wilbraham-GibbsConstant.html


  Vec w4096 = rsDecimateViaMean(w8192, 2); //rsPlotVector(w4096);
  Vec w2048 = rsDecimateViaMean(w4096, 2); //rsPlotVector(w2048);

  std::string fileName;
  fileName = waveName + "_K9.wav";
  rosic::writeToMonoWaveFile(fileName.c_str(), &w4096[0], 4096, 56320, 16);

  //fileName = waveName + "_K21.wav";
  //rosic::writeToMonoWaveFile(fileName.c_str(), &w2048[0], 2048, 56320, 16);
  // seems like we don't need the 8192 version..


  // From the waveform of length 2048, we create the mip-map using the FFT/iFFT technique:
  MipMappedWaveTableStereo mipMap;
  double* pWave[2];
  pWave[0] = pWave[1] = &w2048[0];
  mipMap.setWaveform(pWave, 2048);
  Vec tmp(mipMap.getTableLength());  // temp buffer for the succesive mip-map levels
  int key = 21;
  for(int i = 0; i < mipMap.getNumLevels(); i++)
  {
    mipMap.copyDataTo(&tmp[0], 0, i);
    fileName = waveName + "_K" + std::to_string(key) +  ".wav";
    rosic::writeToMonoWaveFile(fileName.c_str(), &tmp[0], 2048, 56320, 16);
    key += 12;
    //rsPlotVector(tmp);
  }

  // Observations:
  // -In the 0th mip-map, it seems like the Nyquist freq is missing?
  // -Allow the user to specify a tapering function such that we do not necessarily have to use
  //  hard brickwall filters with all of their strong ringing
  //
  // ToDo:
  // -Figure out (and maybe fix) the missing Nyquist freq in the 0th mip-map level.
  //  -Commenting out the "Truncate the spectrum for the next iteration" part doesn't fix it.
  //  -Maybe it has to do with the weird encoding of DC and Nyquist gain in the 0th spectral bin?
  // -Looks like we could have one more level...but it will consume more memory
  // -Write the different rendered levels to wavefiles with a meaningfully formatted way, such as
  //  Saw_K21 for the w2048 wave.

  int dummy = 0;
}


// get rid - there's a copy of it now in RenderScriptTools.cpp which should be used
std::vector<double> randomizePhases(const std::vector<double>& x, int seed, double amount)
{
  int N = (int)x.size();
  RAPT::rsAssert(RAPT::rsIsPowerOfTwo(N), "Works only for power of 2 lengths");

  // Do forward transform:
  std::vector<double> y(N), a(N), p(N);
  using FT = RAPT::rsFourierTransformerRadix2<double>;   // todo: use Bluestein later
  FT ft;
  ft.setBlockSize(N);
  ft.getRealSignalMagnitudesAndPhases(&x[0], &a[0], &p[0]);

  // Randomize phases:
  RAPT::rsNoiseGenerator<double> ng;
  ng.setSeed(seed);
  ng.setRange(0.0, amount * 2*PI);
  for(int k = 1; k < N; k++) {    // DC is unaffected
    p[k] += ng.getSample();
    p[k] =  RAPT::rsWrapToInterval(p[k], -PI, +PI); }

  // Do inverse transform and return result:
  ft.getRealSignalFromMagnitudesAndPhases(&a[0], &p[0], &y[0]);
  return y;
}


void createBassdrumPsy1Sample(double freqScale = 1.0, bool plot = false)
{
  // Create a bassdrum sample using a weighted sum of exponential envelopes with different decay
  // times for a sine-sweepdown. Create also overtones at twice and thrice the frequency, maybe let
  // them have an attack/decay and an attack envelope respectively

  // Rendering parameters:
  int fs = 44100;
  //int N  = 44100*2;
  int N = rsPowInt(2, 16); // maybe have special function rsPowerOf2
  //bool plot = false;  // if true, we plot stuff, if false, we render a wavefile

  // Frequency envelope parameters:
  double freqDecay1  =   10;    // in ms
  double freqDecay2  =   80;
  double freqDecay3  =  240;
  double freqWeight1 =  6.0;
  double freqWeight2 = -1.0;
  double freqWeight3 =  1.0;
  double freqFloor   =  0.0;
  double freqCeil    =  800;  // an excursion/depth would be better
  //double freqScale   = 1.0;    // 1: fundamental, 2,3,4,etc: overtones

  // Amplitude envelope parameters:
  double ampDecay1  =   20;
  double ampDecay2  =  100;
  double ampDecay3  =  400;
  double ampWeight1 =  0.75;
  double ampWeight2 = -1.0;
  double ampWeight3 =  1.0;


  int plotDecimate = 32;
  if(plot)
  {
    fs /= plotDecimate;
    N  /= plotDecimate;
  }

  // Create frequency envelope:
  using AT = RAPT::rsArrayTools;
  RAPT::rsAttackDecayEnvelope<double> eg1, eg2, eg3; // a simple decay env would suffice

  eg1.setAttackSamples(0);
  eg2.setAttackSamples(0);
  eg3.setAttackSamples(0);
  eg1.setDecaySamples(freqDecay1 * 0.001 * fs);
  eg2.setDecaySamples(freqDecay2 * 0.001 * fs);
  eg3.setDecaySamples(freqDecay3 * 0.001 * fs);
  using Vec = std::vector<double>;
  Vec env1(N), env2(N), env3(N), env(N);
  eg1.noteOn(60, 100);
  eg2.noteOn(60, 100);
  eg3.noteOn(60, 100);
  for(int n = 0; n < N; n++)
  {
    env1[n] = eg1.getSample();
    env2[n] = eg2.getSample();
    env3[n] = eg3.getSample();
    env[n]  =  freqWeight1*env1[n] + freqWeight2*env2[n] + freqWeight3*env3[n];
  }
  if(plot) rsPlotVectors(env1, env2, env3, env);

  // Create the raw sine-sweep:
  AT::normalize(&env[0], N);  // makes it inconvenient to port to realtime
  Vec xL(N), xR(N);
  double phi = 0;  // phase
  for(int n = 0; n < N; n++)
  {
    xL[n] = sin(phi - PI/4);
    xR[n] = sin(phi + PI/4);
    double f = freqFloor + (freqCeil - freqFloor) * env[n];
    double w = freqScale * 2*PI*f/fs;
    phi += w;
  }
  if(plot) rsPlotVectors(xL, xR);

  // Create amplitude envelope:
  eg1.setDecaySamples(ampDecay1 * 0.001 * fs);
  eg2.setDecaySamples(ampDecay2 * 0.001 * fs);
  eg3.setDecaySamples(ampDecay3 * 0.001 * fs);
  eg1.reset(); eg1.noteOn(60, 100);
  eg2.reset(); eg2.noteOn(60, 100);
  eg3.reset(); eg3.noteOn(60, 100);
  for(int n = 0; n < N; n++)
  {
    env1[n] = eg1.getSample();
    env2[n] = eg2.getSample();
    env3[n] = eg3.getSample();
    env[n]  =  ampWeight1*env1[n] + ampWeight2*env2[n] + ampWeight3*env3[n];
  }
  if(plot) rsPlotVectors(env1, env2, env3, env);

  // Apply amplitude envelope and fade-out:
  AT::normalize(&env[0], N);
  for(int n = 0; n < N; n++)
  {
    xL[n] *= env[n];
    xR[n] *= env[n];
  }
  RAPT::rsFadeOut(&xL[0], N-N/8, N-1);
  RAPT::rsFadeOut(&xR[0], N-N/8, N-1);

  // Create ambience sample:
  Vec amb = randomizePhases(xL+xR, 2, 1.0);
  rsNormalize(&amb[0], N, 1.0);

  // Plot final result or write to wvaefile: 
  if(plot)  rsPlotVectors(xL, xR);
  if(!plot)
  {
    std::string fileName = "BassdrumPsy1";
    if(freqScale == 1.0) fileName += "_Prime";
    if(freqScale == 1.5) fileName += "_Fifth1";
    if(freqScale == 2.0) fileName += "_Octave1";
    if(freqScale == 3.0) fileName += "_Fifth2";
    if(freqScale == 4.0) fileName += "_Octave2";
    std::string nameWithExt = fileName + ".wav";
    rosic::writeToStereoWaveFile(nameWithExt.c_str(), &xL[0], &xR[0], N, (int)fs);

    fileName += "_Ambience";
    nameWithExt = fileName + ".wav";
    rosic::writeToMonoWaveFile(nameWithExt.c_str(), &amb[0], N, (int)fs);
  }

  // ToDo:
  // -Write this as a jupyter notebook - we really want quck REPL evaluation for this
  // -Maybe write this as an APE script
  //  -Give the user just a single choice parameter to select the preset, settings are encoded in 
  //   the code...but maybe give some macro-parameters to the user
  //  -but APE does not yet support midi - maybe trigger the drum with an input impulse
  // -Have an envelop for the waveshape. it should control the amount of feedback PM:
  //  y[n] = sin(phi + fb*y[n])  -> nonlinear -> needs implicit solver maybe using y[n-1] as 
  //  initial guess, Newton iteration: f(y) = y - sin(phi + fb*y), f'(y) = 1 - fb*cos(phi + fb*y)
  //  -maybe have also PM by a 2nd independent signal...maybe a sine at twice the 
  //  freq? I want something that turns a sine into a square - the self-feedback leads to a 
  //  sawtooth shape...or maybe make it squarish by waveshaping maybe use something like
  //  tanh(d*x) / d except when d == 0 which is treated as special case returning just x...have a
  //  drive envelope for d
  // -create a difference between a signal with and without the waveshaping and phasemodulation 
  //  applied - the difference gives us an "overtones" signal that we can mix with the fundamental
  //  signal to taste
  // -Maybe remove the normalization of the freq-env. it leads to the effect that increasing 
  //  freqWeight1 increases the overall freq but we want it to affect only the initial freq
  // -Maybe produce stereo output by using +-45° phase shift for left/right
  // -Mix in a very short noise-burst to the attack (might have settings for: decay, HP, LP, color, seed)
  //  -maybe it should have a filter env
  //  -maybe clip or saturate sweep+noise
  // -Apply reverb - maybe using a phase-randomization algorithm: do one big FFT (maybe 
  //  length = 2^16 = 65536), randomize phases, iFFT ...this need its own amp-env
  // -Add an envelope generator to ToolChain that implements a weighted sum of exponential decays
  // -Maybe give ToolChain a rendering functionality

  // -try to raise the amp-env to a time-varying power - with 2 , it should look more like Gaussian
  // -apply distortion, maybe bitcrushing in a time variant manner to add some noisiness to the
  //  transient
  // -for such effects, obtain a difference signal to the original and store the pure effect 
  //  seperately to be mixed in later
  // -use a more sawtooth-like waveform to give it some overtones...however, maybe the 1st and 2nd
  //  should be notched out because 100 or 150 hz are ugly in bassdrums. or: maybe also subtract
  //  from a sample using a sine, obtain difference and pass that through a highpass before mixing
  // -the (gaussian) envelope should be applied after reverb/convolution -> sort of gated reverb


  // ToDo:
  // -render some goasque creaking sounds and make an sfz (use pan-envelopes)...but how? maybe using a 
  //  percussive sound with echo using a delay-length corresponding to a note and a lot of feedback?
  // -render some multiplicative synthesis lead sounds
  // -render some swooshes and reversed sounds ...maybe a noise with filter env
  // -create a bouncy goa-bass that harmonizes well with this bassdrum...maby they need to also to 
  //  the one_shot thing? but the we must take care of interference...maybe only the amp-envs 
  //  should get added on to of each other - the oscs just keep playing
  // -maybe wrap bassdrum and bass into a single instrument using keyranges - the wen can also 
  //  apply some dynamics processing to their mix (if bus_mode is active)
  // -maybe convolve it with an expoentially decaying white noise sound, 
  //  -maybe that noise should have a color envelope from white to brown, say

  
  // -For dubstep sounds, try to use a percussive sample like a snare, set up a loop and modify its 
  //  location via midi-cc, maybe have a second one an octave below going on. maybe change the loop 
  //  length (while adjusting the increment accordingly), maybe also when moving loop_start and 
  //  loop_end simultaneously, adapt the increment to account for the fact that the loop point is
  //  receding or approaching ...but maybe not
  // -it has often very formantish stuff going on, try also comb-filters and sync-stuff
  // -LFOs are important for modulating the position - their frequency should be midi controlled
  // -use MSEGs, attacks should be normal, decays should use "anti-analog" shape
  // -speed should speed up or down in step of simple ratios (like doubling, tripling, halving, ..)
  // -distortion and bitcrushing also helps
  // -Maybe a multiplicative synthesis sample could also be used as basis
  // https://www.youtube.com/watch?v=dknbNrr4EDo&list=RDQM6d8Vubj1zMg&start_radio=1


  int dummy = 0;
}

void createSweepDrummerSamples(int sampleRate, double freqScale)
{
  // We do the same as the function above but this time with less code using the new rsSweepDrummer
  // class. When this is finished, the function above will be obsolete...





}

// move this to RenderScripts:
void createBassdrumPsy1Samples()
{
  //createBassdrumPsy1Sample(1.0, true); return;  // for development


  bool plot = false;
  createBassdrumPsy1Sample(1.0, plot);
  createBassdrumPsy1Sample(1.5, plot);
  createBassdrumPsy1Sample(2.0, plot);
  createBassdrumPsy1Sample(3.0, plot);
  createBassdrumPsy1Sample(4.0, plot);
}

void createMiscSamples()
{
  // Create miscelanneous other samples that are useful as raw material in the sampler engine.

  using Vec = std::vector<double>;
  using AT  = RAPT::rsArrayTools;

  int fs = 44100;
  int N  = 0;       // number of samples
  
  // Generate a unit impulse:
  double one = 1.0;
  rosic::writeToMonoWaveFile("UnitImpulse.wav", &one, 1, fs);
  // -Maybe we should make it a few samples long. Having a wavfile containing just a single value
  //  may be a corner case that some sampler engines won't like? But it actually makes for a nice
  //  unit test to see what an engine does in such an extreme case. OK, my engine handles it 
  //  nicely: it just produces a unit impulse indeed. 
  // -But: we can't create filter blips this way because the region player immediately stops, 
  //  giving the filter no ringtout time. We need a couple of samples of silence after the impulse
  //  and set up a loop over this silence portion.
  //  

  // Generate 5 seconds of white noise with uniform amplitude distribution:
  N = 5*fs;
  Vec x(N);
  RAPT::rsNoiseGenerator<double> ng;
  for(int n = 0; n < N; n++)
    x[n] = ng.getSample();
  rosic::writeToMonoWaveFile("UniformWhiteNoise.wav", &x[0], N, fs);
  // Notes:
  // -Layering this sample with itself with various offsets can be used to obtain Irwin-Hall
  //  distributed white noise.
  // -5 seconds should be long enough to have no noticable repetition pattern and/or 
  //  comb-filtering artifacts when layering several shifted copies.


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

  g.setKeyRangeToRender(21, 93);  
  // For production rendering, use 21..93. For preview and sound-design, use smaller range
  // ToDo: maybe use an increment (default 1, 12 means one sample per octave - which is also good 
  // for preview), see SampleMapGenerator::generateAllSamples

  g.setTruncationLevel(-50);
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
  // -In the sfz file, we could specify a 2nd region with the same same but a very short decay
  //  envelope and add it to the normal sampel. the result is a two-stage-decay or enhanced
  //  transient. making it velocity dependent, we can do things like: higher vel leads to stronger
  //  transients
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


