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











// From: RS-MET\Libraries\RobsJuceModules\rapt\Unfinished\Math\Matrix.cpp
/*
template<class T>
void rsMatrixOld<T>::print()
{
  printf("%s %d %s %d %s", "rsMatrixOld - rows: ", data->numRows, " columns:", data->numColumns, "\n");
  for(int r = 0; r < data->numRows; r++)
  {
    for(int c = 0; c < data->numColumns; c++)
      printf("%.4f %s", data->m[r][c], "  ");
    printf("%s", "\n");
  }
}
// ToDo: comvert code to work with new rsMatrix
*/




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
}

void printTestComplex()
{
  using Complex = rsComplex<double>;

  Complex z(3, 2);

  rsPrint(z);

  // ToDo:
  //
  // - Try it also with float
}

void printTestMatrix()
{

  int dummy = 0;
}

void printTests()
{
  printTestComplex();
  printTestMatrix();
}


//MemLeakTest memLeakTest;