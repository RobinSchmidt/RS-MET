#include "Misc.h"

#include <rosic/rosic.h> // get rid
using namespace RAPT;
using namespace rosic;

#include "../rosic_tests/PortedFromRSLib/Examples/ModalExamples.cpp"
#include "../rosic_tests/PortedFromRSLib/Examples/SampleMapGenerator.cpp"

// Move to TestInputCreation.h/.cpp in rs_testing/TestTools/Utilities:
std::vector<double> createPluckedString(int numSamples, double frequency, double sampleRate)
{
  std::vector<double> x(numSamples);


  RAPT::rsModalFilterBank<double, double> mfb;

  int numModes = int(0.5*sampleRate/frequency);
  std::vector<double> f(numModes);  // relative frequencies  
  std::vector<double> g(numModes);  // gains
  std::vector<double> a(numModes);  // relative attack times
  std::vector<double> d(numModes);  // relative decay times
  std::vector<double> p(numModes);  // start-phases

  RAPT::rsArrayTools::fillWithRangeLinear(&f[0], numModes, 1.0, double(numModes));
  double amp = 0.1;
  double c = 0.4;
  for(int k = 0; k < numModes; k++) g[k] = amp / pow(k+1.0, c);  // amplitudes follow 1/k^c rule
  p = mfb.randomModePhases(g);
  d = mfb.modeDecayTimes(f, 10, 1.0);
  a = d;

  mfb.setReferenceFrequency(frequency);
  mfb.setReferenceAttack(0.02);
  mfb.setModalParameters(f, g, a, d, p);

  // factor out - maybe let rsModalFilterBank itself be able to produce its impulse response
  // as std::vector for convenience:
  x[0] = mfb.getSample(1.0);
  for(int n = 1; n < numSamples; n++)
    x[n] = mfb.getSample(0.0);

  return x;
}

void sampleTailExtenderTest()
{
  int N = 10*44100;
  double f = 220;
  double fs = 44100;

  // create test signal:
  std::vector<double> x = createPluckedString(N, f, fs);
  rosic::writeToMonoWaveFile("TestPluck.wav", &x[0], N, (int)fs, 16);


  // extend test signal:
  SampleTailExtender ste;
  ste.setSynthesisStartPointInSeconds(1.5);
  ste.setDecayRate(2.2);
  //ste.setCutoffThresholdInDecibels(-70);
  std::vector<double> y = ste.extendSample(x, fs, 9);
  rosic::writeToMonoWaveFile("TestPluckExtended.wav", &y[0], (int)y.size(), (int)fs, 16);
}














template<class T>
void rsPrint(const RAPT::rsMatrixView<T> A)
{
  int M = A.getNumRows();
  int N = A.getNumColumns();

  printf("%s %d %s %d %s", "rsMatrix - rows: ", M, " columns:", N, "\n");
  for(int r = 0; r < M; r++)
  {
    for(int c = 0; c < N; c++)
      printf("%.4f %s", A(r,c), "  ");
    printf("%s", "\n");
  }

  // ToDo:
  //
  // - There isn't any sensible formatting at the moment. We want aligned columns!
}

template<class T>
void rsPrint(const RAPT::rsComplex<T> z)
{
  printf("%s %.5f %s %.5f %s", "z = ", z.re, " + ", z.im, "j\n");
  // Should eventually go into TestUtilities.h. But currently, we get compilation errors when we 
  // put it there. Something needs to be changed about the include structure.

  // ToDo:
  //
  // - Factor out a function that produces a std::string - maybe call it rsToString()
  //
  // - Give the user some parameters to control the format (precision, etc.)
  //
  // - We someday actually want to not directly invoke printf here but intead call another, 
  //   overloaded rsPrint() function on z.re, r.im to make it it universally useful. For example,
  //   re,im could be of type rsFraction. There should be all kinds of different overloads of 
  //   rsPrint (or bette rsToString()) for all sorts of types. The overload resolution should just
  //   keep invoking overloads until it hits one of the base cases, i.e. an implementation for a
  //   primitive data type such as int, float, double, etc.
}

void printTestComplex()
{
  using Complex = rsComplex<double>;

  Complex z(3, 2);

  rsPrint(z);

  // ToDo:
  //
  // - Try it also with float. Make the function a template and call it like 
  //   printTestComplex<float>(), printTestComplex<double>()
}

void printTestMatrix()
{
  using Elem = double;
  using Mat  = rsMatrix<Elem>;

  Mat A(3, 4, {1,2,3, 4,5,6, 7,8,9, 10,11,12});

  rsPrint(A);
}

void printTests()
{
  printTestComplex();
  printTestMatrix();
}





//MemLeakTest memLeakTest;