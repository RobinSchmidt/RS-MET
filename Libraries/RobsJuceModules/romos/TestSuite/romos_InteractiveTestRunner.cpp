#include "romos_InteractiveTestRunner.h"
using namespace rsTestRomos;

InteractiveTestRunner::InteractiveTestRunner()
{

}

void InteractiveTestRunner::runTests()
{
    
  InteractivePlotTest *plotTest; 


  //plotTest = new BandlimitedImpulseTrainPlotTest();   plotTest->runTestAndPlotResults(); delete plotTest;
  plotTest = new SawOscillatorPlotTest();   plotTest->runTestAndPlotResults(); delete plotTest;
  //plotTest = new DualBlitSawOscillatorPlotTest();   plotTest->runTestAndPlotResults(); delete plotTest;


  //plotTest = new EnvelopeADSRPlotTest();   plotTest->runTestAndPlotResults(); delete plotTest;
  //plotTest = new BiquadDesignerPlotTest(); plotTest->runTestAndPlotResults(); delete plotTest;





  int dummy = 0;
}

