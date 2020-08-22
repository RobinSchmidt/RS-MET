#include "romos_PerformanceTestRunner.h"
using namespace rsTestRomos;

PerformanceTestRunner::PerformanceTestRunner()
{

}

void PerformanceTestRunner::runAllTestsAndPrintResultsToConsole(bool createLogFile)
{
  rosic::rsString report;

  rosic::rsString header  = "RoMoS performance tests\n";
  header               += "CPU cycles required to compute one sample-frame for one voice:\n";
  header               += "Processing function used:         mono/frame  mono/block  poly/frame  poly/block\n";
  header.printToStandardOutput();
  report += header;

  report += runInternalFunctionPerformanceTests();
  report += runFrameworkPerformanceTests();
  report += runAtomicModulePerformanceTests();

  if( createLogFile == true )
    rosic::writeStringToFile("E:/TmpData/RomosPerformanceLog.txt", report.getRawString());
}


rosic::rsString PerformanceTestRunner::runFrameworkPerformanceTests()
{
  rosic::rsString report;
  PerformanceTest *test; 
  
  rosic::rsString header = "Framework:\n";
  header.printToStandardOutput();
  report += header;

  test = new BiquadMacroPerformanceTest();                report += test->runTestsAndGetReport();  delete test;
  test = new IdentityChainPerformanceTest();              report += test->runTestsAndGetReport();  delete test;
  test = new IdentityChainWithFeedbackPerformanceTest();  report += test->runTestsAndGetReport();  delete test;
  test = new AdderChainPerformanceTest();                 report += test->runTestsAndGetReport();  delete test;
  test = new AdderChainWithFeedbackPerformanceTest();     report += test->runTestsAndGetReport();  delete test;

  printf("%s", "\n");

  return report;
}

rosic::rsString PerformanceTestRunner::runAtomicModulePerformanceTests()
{
  rosic::rsString report;
  PerformanceTest *test; 

  rosic::rsString header = "Atomic Modules:\n";
  header.printToStandardOutput();
  report += header;

  test = new FirstOrderFilterPerformanceTest();         report += test->runTestsAndGetReport();  delete test;
  test = new BiquadPerformanceTest();                   report += test->runTestsAndGetReport();  delete test;
  test = new BiquadDesignerPerformanceTest();           report += test->runTestsAndGetReport();  delete test;
  test = new BandlimitedImpulseTrainPerformanceTest();  report += test->runTestsAndGetReport();  delete test;
  test = new SawOscillatorPerformanceTest();            report += test->runTestsAndGetReport();  delete test;

  test = new Formula11PerformanceTest();                report += test->runTestsAndGetReport();  delete test;
  test = new FormulaN1PerformanceTest();                report += test->runTestsAndGetReport();  delete test;
  test = new FormulaNMPerformanceTest();                report += test->runTestsAndGetReport();  delete test;


  printf("%s", "\n");
  return report;
}

double rsTestRomos::dummyFunction(double x)
{
  return x;
}


double rSinh(double x)
{
  long double ex = exp(x);    // e^x
  return (ex*ex - 1.0) / (ex+ex);
}
double rCosh(double x)
{
  long double ex = exp(x);    // e^x
  return (ex*ex + 1.0) / (ex+ex);
}
double rTanh(double x)
{
  long double e2x = exp(2*x); // e^(2x)
  return (e2x - 1.0) / (e2x + 1.0);
}

double rHypot(double x, double y)
{
  long double xL = x;
  long double yL = y;
  return sqrt(xL*xL + yL*yL);
}



rosic::rsString PerformanceTestRunner::runInternalFunctionPerformanceTests()
{
  rosic::rsString report;

  rosic::rsString header = "Internal Functions:\n";
  header.printToStandardOutput();
  report += header;

  ProcessorCycleCounter counter;

  static const int numCalls = 2000;
  static const int numRuns  = 20;

  double x[numCalls];
  double y[numCalls];
  RAPT::rsArrayTools::fillWithIndex(x, numCalls);   // maybe use random values
  RAPT::rsArrayTools::fillWithIndex(y, numCalls);


  double minCyclesPerCall = pow(10.0, 100.0);
  for(int runIndex = 1; runIndex <= numRuns; runIndex++)
  {
    counter.init();
    for(int n = 0; n < numCalls; n++)
    {
      // insert a macro CALL_FUNCTION here ..of a function callFunction

                                                                 //  cycles per call
      y[n] = dummyFunction(x[n]);              //    3.04
      //y[n] = sqrt(x[n]);                       //   24.04

      // trigonometric functions:
      //y[n] = hypot(x[n], y[n]);                  // 1415.58
      //y[n] = rHypot(x[n], y[n]);                 // 26.11
      //y[n] = sin(x[n]);                        //   94.59
      //y[n] = cos(x[n]);                        //   97.54
      //y[n] = tan(x[n]);                        //  127.92
      //sinCos(x[n], &(y[n]), &(y[n]));          //  139.54

      //y[n] = fabs(x[n]);                       //    3.04
      //y[n] = fmod(x[n], y[n]);                 //  125.xx
      //y[n] = sign(x[n]);                       //    8.49
      //y[n] = exp(x[n]);                        //   58.50
      //y[n] = log(x[n]);                        //  123.95

      // hyperbolic functions:
      //y[n] = sinh(x[n]);                       //  368.91
      //y[n] = rSinh(x[n]);                      //   81.07
      //y[n] = cosh(x[n]);                       //  393.xx
      //y[n] = rCosh(x[n]);                      //   81.07
      //y[n] = tanh(x[n]);                       //  215.81
      //y[n] = rTanh(x[n]);                      //   81.07
      //y[n] = tanhApprox(x[n]);                 //   29.63

      // inverse hyperbolic functions:
      //y[n] = asinh(x[n]);                        // 163.13






    }

    double cyclesPerCall = (double) counter.getNumCyclesSinceInit() / (double) numCalls;
    if( cyclesPerCall < minCyclesPerCall && cyclesPerCall > 0.0 )
      minCyclesPerCall = cyclesPerCall;
  }

  printf("%s %.2f %s", "Cycles per call:", minCyclesPerCall ,"\n");


  printf("%s", "\n");
  return report;
}

