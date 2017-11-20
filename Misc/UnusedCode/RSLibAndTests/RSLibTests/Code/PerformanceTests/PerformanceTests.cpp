#include "RSCore/CoreTests.h"


#include "RSAudio/Synthesis/ModalTests.h"
#include "RSAudio/Analysis/AnalysisTests.h"
#include "RSAudio/MiscAudioTests.h"

#include "RSMath/MathTests.h"

// \todo: get rid of the different subdirectories RSAudio/Analysis, RSAudio/Synthesis - it should
// be sufficient to have a file for each of the submodules



void runTests()
{
  std::string reportString = createPerformanceTestHeader();

  // RSCore:
  //testFlagArray(reportString);

  // RSMath:
  testAbsAndSign(               reportString);
  //testMultinomialCoefficients(reportString);
  //testPrimeSieves(            reportString);
  //testMatrix(                 reportString);
  //testMatrixAddressing(         reportString);

  // RSAudio:
  //testFourierTransformer(reportString);
  //testAutoCorrelationPitchDetector(reportString);
  //testModalFilter(reportString);
  //testModalFilterBank(reportString);
  //testSincInterpolator(reportString);


  // \todo ask, if the report should be saved into a file - if so, create a file that contains the 
  // date/time, platform, compiler, configuration etc. and write the report into it

  std::cout << reportString;
}

int main(int argc, char** argv)
{
  runTests();
  getchar();
  return(EXIT_SUCCESS);
}

