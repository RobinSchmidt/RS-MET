#include "romos_InteractivePlottingTests.h"
using namespace rsTestRomos;



InteractivePlotTest::InteractivePlotTest(const char *testName)
: ProcessingTest(testName)
{
  numFramesToProcess = 801;   // index == 400 -> f = 0 Hz, index == 600 -> f = fs/2
  xAxis              = timeAxis;
}


void InteractivePlotTest::runTestAndPlotResults()
{
  initTest();
  setTestPolyphonic(false);

  clearOutputSignalArrays();
  moduleToTest->resetStateForAllVoices();
  processModuleInFrames();

  plotResult();
  //plotDesiredAndActualOutput(0, 0, numFramesToProcess, 0);
}



//-----------------------------------------------------------------------------------------------------------------------------------------

BandlimitedImpulseTrainPlotTest::BandlimitedImpulseTrainPlotTest()
: InteractivePlotTest("BandlimitedImpulseTrainPlot")
{
  numFramesToProcess = 4000;   
  //moduleToTest       = ModuleFactory::createModule(ModuleTypeRegistry::BANDLIMITED_IMPULSE_TRAIN);
  moduleToTest = moduleFactory.createModule("BandlimitedImpulseTrain");
}
void BandlimitedImpulseTrainPlotTest::fillInputSignalArraysWithTestSignal()  
{
  //initForConstantFreq(500.0);  
  //initForConstantFreq(551.25);
  //initForConstantFreq(100.0);
  //initForConstantFreq(441.0);
  //initForConstantFreq(440.0);

  //initForConstantFreq(4410.0);   // fs/10
  //initForConstantFreq(8820.0);     // fs/5
  //initForConstantFreq(11025.0);    // fs/4
  initForConstantFreq(44100.0/3.0);  

  //initForConstantFreq(4400.0);
  //initForConstantFreq(1000.0);

  //initForLinearFreqSweep(100.0, 200.0);
  //initForLinearFreqSweep(1000.0, 2000.0);
  //initForLinearFreqSweep(100.0, 1000.0);

}
void BandlimitedImpulseTrainPlotTest::plotResult()
{
  RAPT::rsAssert(false, "plotting code needs update");
  //Plotter::plotData(numFramesToProcess, xAxis, outputs[0][0]);
}

void BandlimitedImpulseTrainPlotTest::initForConstantFreq(double freq)
{
  RAPT::rsArrayTools::fillWithValue(inputs[0][0], numFramesToProcess, freq);
  RAPT::rsArrayTools::fillWithValue(inputs[0][1], numFramesToProcess, 0.0);
}
void BandlimitedImpulseTrainPlotTest::initForLinearFreqSweep(double minFreq, double maxFreq)
{
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, minFreq, maxFreq);
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, 0.0);
}






//-----------------------------------------------------------------------------------------------------------------------------------------

SawOscillatorPlotTest::SawOscillatorPlotTest(const char *testName)
: InteractivePlotTest(testName)
{
  numFramesToProcess = 4000;   
  //moduleToTest       = ModuleFactory::createModule(ModuleTypeRegistry::BLIT_SAW_OSCILLATOR);
  moduleToTest = moduleFactory.createModule("BlitOscillator");
}
void SawOscillatorPlotTest::fillInputSignalArraysWithTestSignal()  
{
  double fs = processingStatus.getSystemSampleRate();


  initForConstantFreq(22.04, 0.505); // seems to cause problems becaus it hitsb the cosine branch sometimes

  initForConstantFreq(22.04, 0.6); 

  //initForConstantFreq(224.0, 0.00); // seems to cause problems becaus it hitsb the cosine branch sometimes when margin is 0.01


  //initForConstantFreq(100.0, 0.0);
  //initForConstantFreq(-100.0, 0.0); 

  //initForConstantFreq(100.0, 0.25);
  //initForConstantFreq(100.0, 0.5);
  //initForConstantFreq(100.0, 0.50485);  // prpoblem case - seems to have a linear upward trend
  //initForConstantFreq(100.0, 0.75);
  //initForConstantFreq(100.0, 0.1);
  //initForConstantFreq(100.0, 1.0); 


  //initForConstantFreq(132.0, 0.0);
  //initForConstantFreq(132.0, 0.5);
  //initForConstantFreq(132.0, 0.4); 
  //initForConstantFreq(132.0, 0.3);  // seems to be a weird special case
  //initForConstantFreq(537, 0.57);

 
  //initForConstantFreq(441.0, 0.0);
  //initForConstantFreq(441.0, 0.05);
  //initForConstantFreq(441.0, 0.10);
  //initForConstantFreq(441.0, 0.25);
  //initForConstantFreq(441.0, 0.5);
  //initForConstantFreq(441.0, 0.505);   // problem case  

    
  // some settings cause DC transients, persumably due to Gibbs-ripple which is present in the bandlimited saw but not taken into account
  // by the state initialization (which uses the value of a non-bandlimited saw as desired value for the first sample):

  //initForConstantFreq(100.0, 0.5005);  
  //initForConstantFreq(100.0, 0.500);  
    // maybe we can use a table with desired values of a bandlimited saw


  //initForConstantFreq(440.0);
  //initForConstantFreq(507.2, 0.231);

  //initForConstantFreq(4411.0);
  //initForConstantFreq(4412.0);
  //initForConstantFreq(4413.0);
  //initForConstantFreq(4400.0);
  //initForConstantFreq(1000.0);

  // frequencies which divide the samplerate without remainder seem to cause problems (DC transients):
  //initForConstantFreq(fs/10);   // problem case
  //initForConstantFreq(fs/10, 0.35);

  //initForConstantFreq(fs/5);  
  //initForConstantFreq(fs/5, 0.5);  


                                               // desired first sample
  //initForConstantFreq(fs/2, 0.0);            // +0.5
  //initForConstantFreq(fs/2, 0.5);            // -0.5

  //initForConstantFreq(fs/3, 0.0);            // +0.5
  //initForConstantFreq(fs/3, 0.25);   
  //initForConstantFreq(fs/3, 0.3);   
  //initForConstantFreq(fs/3, 1.0/3.0);        //  0.0
  //initForConstantFreq(fs/3, 2.0/3.0);        // -0.5

  //initForConstantFreq(fs/4, 0.0);            // +0.5
  //initForConstantFreq(fs/4, 0.25);           //  0.0
  //initForConstantFreq(fs/4, 0.5);            //  0.0
  //initForConstantFreq(fs/4, 0.75);           // -0.5


  //initForConstantFreq(fs/4.2, 0.0);          //

  //initForConstantFreq(fs/4.5, 0.0);          // +0.5 (so-so)
  //initForConstantFreq(fs/4.5, 0.25);         //  0.0
  //initForConstantFreq(fs/4.5, 0.50);         //  0.0
  //initForConstantFreq(fs/4.5, 0.75);         // -0.5 (so-so)

  //initForConstantFreq(fs/5, 0.0);              // +0.5
  //initForConstantFreq(fs/5, 0.20);           // +0.25
  //initForConstantFreq(fs/5, 0.25);           // +0.125 (so-so)
  //initForConstantFreq(fs/5, 0.30);           //  0.0
  //initForConstantFreq(fs/5, 0.35);           //  0.0 (so-so)
  //initForConstantFreq(fs/5, 0.40);           //  0.0
  //initForConstantFreq(fs/5, 0.5);            //  0.0
  //initForConstantFreq(fs/5, 0.6);              // -0.25
  //initForConstantFreq(fs/5, 0.7);              // -0.56
  //initForConstantFreq(fs/5, 0.75);      


  // we can obtain these star-values automatically, by using an initial value of 0, performing the integration over one cycle and
  // the taking the negative value of the mean of the integrated cycle








  //initForConstantFreq(fs/3, 0.0);  // heavy DC transient
  //initForConstantFreq(fs/3, 1.0/3.0);
  //initForConstantFreq(fs/3, 0.5);  





  // frequencies near, at and above Nyquist-Freq:
  //initForConstantFreq(20000.0); 
  //initForConstantFreq(fs/2.0); 
  //initForConstantFreq(fs/2, 0.25);    // should disappear - does, but there's a decaying DC...(wrong integrator initialization? seems so)
  //initForConstantFreq(fs/2, 0.5); 
  //initForConstantFreq(fs/2, 0.75);    // should disappear - does, but there's a decaying DC...

  //initForConstantFreq( 22050.1, 0.0);       // above Nyquist - should produce silence
  //initForConstantFreq(-22050.1, 0.0);       // above Nyquist - should produce silence

  //initForConstantFreq(4410.01);  // slightly off, everything is fine - ...mhhh
  //initForConstantFreq(2205.0);     // not OK
  //initForConstantFreq(1102.5);     // not OK
  //initForConstantFreq( 551.25);
  //initForConstantFreq( 551.25, 0.505);
  //initForConstantFreq( 551.25001);

  //initForConstantFreq(8820.0);     // OK


  //initForLinearFreqSweep(100.0, 200.0);
  //initForLinearFreqSweep(1000.0, 2000.0);
  //initForLinearFreqSweep(100.0, 1000.0);

  //initForFreqSwitch(100, 150, 2000);

  //initForFreqSwitch(507.2, 143.7, 1953);  

  //initForFreqSwitch(507.2, 100, 1950);  
    // 507.2 Hz -> 100 Hz @ 2000 samples (44100 Hz samplerate) gives a big DC offset for the samples after 2000 when not re-initializing 
    // the integrator

}
void SawOscillatorPlotTest::plotResult()
{
  double dc = RAPT::rsArrayTools::mean(outputs[0][0], 50);

  RAPT::rsAssert(false, "plotting code needs update");
  //Plotter::plotData(numFramesToProcess, xAxis, outputs[0][0]);
}

void SawOscillatorPlotTest::initForConstantFreq(double freq, double phase)
{
  RAPT::rsArrayTools::fillWithValue(inputs[0][0], numFramesToProcess, freq);
  RAPT::rsArrayTools::fillWithValue(inputs[0][1], numFramesToProcess, phase);

  RAPT::rsArrayTools::fillWithValue(inputs[0][2], numFramesToProcess,  0.5);
  RAPT::rsArrayTools::fillWithValue(inputs[0][3], numFramesToProcess, -1.0);

}
void SawOscillatorPlotTest::initForLinearFreqSweep(double minFreq, double maxFreq)
{
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, minFreq, maxFreq);
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, 0.0);
}
void SawOscillatorPlotTest::initForFreqSwitch(double freq1, double freq2, int switchTimeInSamples)
{
  int length1 = switchTimeInSamples;
  int length2 = numFramesToProcess - length1;

  RAPT::rsArrayTools::fillWithValue( inputs[0][0],          length1, freq1);
  RAPT::rsArrayTools::fillWithValue(&inputs[0][0][length1], length2, freq2);
  RAPT::rsArrayTools::fillWithValue(inputs[0][1], numFramesToProcess, 0.0);
}



DualBlitSawOscillatorPlotTest::DualBlitSawOscillatorPlotTest()
: SawOscillatorPlotTest("DualBlitSawOscillatorPlotTest")
{
  numFramesToProcess = 4000;   

  moduleFactory.deleteModule(moduleToTest);
  //moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::DUAL_BLIT_SAW_OSCILLATOR);
  moduleToTest = moduleFactory.createModule("BlitOscillator");
}

//-----------------------------------------------------------------------------------------------------------------------------------------

EnvelopeADSRPlotTest::EnvelopeADSRPlotTest()
: InteractivePlotTest("EnvelopeADSRPlot")
{
  //moduleToTest = ModuleFactory::createModule(ModuleTypeRegistry::ENVELOPE_ADSR);
  moduleToTest = moduleFactory.createModule("ADSR-Envelope");

  sustainValue = 0.2;
  attackShape  = 0.0;
  decayShape   = 0.0;
  releaseShape = 0.0; 
  timeScale    = 1.0;

  initEnvelopeTimesAllNonZero();

  numFramesToProcess = 1000;
  events             = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, 500);
}
void EnvelopeADSRPlotTest::fillInputSignalArraysWithTestSignal()  
{
  // initialzation of the time-values in the envelope:
  initEnvelopeTimesAllNonZero();
  //initEnvelopeTimesZeroAttack();
  //initEnvelopeTimesZeroDecay();
  //initEnvelopeTimesZeroAttackAndDecay();
  //initEnvelopeTimesZeroRelease();

  // initialization of the events to occur:
  initEventsForNormalCase();
  //initEventsForDoubleAttack1();
  //initEventsForDoubleAttack2();
  //initEventsForDoubleAttack3();
  //initEventsForDoubleAttack4();
  //initEventsForDoubleAttack5();
  //initEventsForDoubleAttack6(); 
  //initEventsForDoubleAttack7();   
  //initEventsForEarlyRelease1();
  //initEventsForEarlyRelease2();
  //initEventsForEarlyRelease3();

  // create the signals for the envelope's input pins:
  RAPT::rsArrayTools::fillWithValue(inputs[0][0], numFramesToProcess, attackSeconds);  
  RAPT::rsArrayTools::fillWithValue(inputs[0][1], numFramesToProcess, attackShape); 
  RAPT::rsArrayTools::fillWithValue(inputs[0][2], numFramesToProcess, decaySeconds);  
  RAPT::rsArrayTools::fillWithValue(inputs[0][3], numFramesToProcess, decayShape); 
  RAPT::rsArrayTools::fillWithValue(inputs[0][4], numFramesToProcess, sustainValue);  
  RAPT::rsArrayTools::fillWithValue(inputs[0][5], numFramesToProcess, releaseSeconds);  
  RAPT::rsArrayTools::fillWithValue(inputs[0][6], numFramesToProcess, releaseShape); 

  //rosic::fillWithValue(inputs[0][4], numFramesToProcess, timeScale); 
  xAxis = timeAxis;
}
void EnvelopeADSRPlotTest::plotResult()
{
  //Plotter::plotData(110, xAxis, outputs[0][0]);  // old

  RAPT::rsAssert(false, "plotting code needs update");
  //Plotter::plotData(numFramesToProcess, xAxis, outputs[0][0]);
}
  
void EnvelopeADSRPlotTest::initEnvelopeTimesAllNonZero()
{
  attackSamples  = 100;
  decaySamples   = 200;
  sustainSamples = 200;
  releaseSamples = 400;
  setTimesInSecondsFromTimesInSamples();
}
void EnvelopeADSRPlotTest::initEnvelopeTimesZeroAttack()
{
  attackSamples  = 0;
  decaySamples   = 200;
  sustainSamples = 200;
  releaseSamples = 400;
  setTimesInSecondsFromTimesInSamples();
}
void EnvelopeADSRPlotTest::initEnvelopeTimesZeroDecay()
{
  attackSamples  = 100;
  decaySamples   = 0;
  sustainSamples = 200;
  releaseSamples = 400;
  setTimesInSecondsFromTimesInSamples();
}
void EnvelopeADSRPlotTest::initEnvelopeTimesZeroAttackAndDecay()
{
  attackSamples  = 0;
  decaySamples   = 0;
  sustainSamples = 200;
  releaseSamples = 400;
  setTimesInSecondsFromTimesInSamples();
}
void EnvelopeADSRPlotTest::initEnvelopeTimesZeroRelease()
{
  attackSamples  = 100;
  decaySamples   = 200;
  sustainSamples = 200;
  releaseSamples = 0;
  setTimesInSecondsFromTimesInSamples();
}
void EnvelopeADSRPlotTest::setTimesInSecondsFromTimesInSamples()
{
  noteLength = attackSamples + decaySamples + sustainSamples;

  double fs = processingStatus.getSystemSampleRate();
  attackSeconds  = attackSamples  / fs;
  decaySeconds   = decaySamples   / fs;
  sustainSeconds = sustainSamples / fs;
  releaseSeconds = releaseSamples / fs;
}

void EnvelopeADSRPlotTest::initEventsForNormalCase()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
}
void EnvelopeADSRPlotTest::initEventsForDoubleAttack1()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
  int secondNoteStart = attackSamples/2;
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(81, 64, secondNoteStart, noteLength));
}
void EnvelopeADSRPlotTest::initEventsForDoubleAttack2()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
  int secondNoteStart = attackSamples;
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(81, 64, attackSamples, noteLength));
}
void EnvelopeADSRPlotTest::initEventsForDoubleAttack3()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
  int secondNoteStart = attackSamples + decaySamples/4;
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(81, 64, secondNoteStart, noteLength));
}
void EnvelopeADSRPlotTest::initEventsForDoubleAttack4()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
  int secondNoteStart = attackSamples + decaySamples;
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(81, 64, secondNoteStart, noteLength));
}
void EnvelopeADSRPlotTest::initEventsForDoubleAttack5()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
  int secondNoteStart = attackSamples + decaySamples + sustainSamples/4;
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(81, 64, secondNoteStart, noteLength));
}
void EnvelopeADSRPlotTest::initEventsForDoubleAttack6()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
  int secondNoteStart = attackSamples + decaySamples + sustainSamples;
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(81, 64, secondNoteStart, noteLength));
}
void EnvelopeADSRPlotTest::initEventsForDoubleAttack7()
{
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteLength);
  int secondNoteStart = attackSamples + decaySamples + sustainSamples + releaseSamples/4;
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(81, 64, secondNoteStart, noteLength));
}
void EnvelopeADSRPlotTest::initEventsForEarlyRelease1()
{
  int noteOffTime = attackSamples/2;
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteOffTime);
}
void EnvelopeADSRPlotTest::initEventsForEarlyRelease2()
{
  int noteOffTime = attackSamples;
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteOffTime);
} 
void EnvelopeADSRPlotTest::initEventsForEarlyRelease3()
{
  int noteOffTime = attackSamples + decaySamples/2;;
  events = TestEventGenerator::generateNoteOnOffPair(81, 64, 0, noteOffTime);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

BiquadDesignerPlotTest::BiquadDesignerPlotTest()
: InteractivePlotTest("BiquadDesignerPlot")
{
  numFramesToProcess = 801;   // index == 400 -> f = 0 Hz, index == 600 -> f = fs/2
  //moduleToTest       = ModuleFactory::createModule(ModuleTypeRegistry::BIQUAD_DESIGNER);
  moduleToTest = moduleFactory.createModule("BiquadDesigner");
  fMin               = -processingStatus.getSystemSampleRate();
  fMax               = +processingStatus.getSystemSampleRate();
}
void BiquadDesignerPlotTest::fillInputSignalArraysWithTestSignal()  // override initTest instead
{
  //initForFreqSweepLowpassBilinear6();  
  //initForFreqSweepHighpassBilinear6();  
  //initForFreqSweepLowpassBilinear12(2.0);  
  //initForFreqSweepLowpassBilinear12(0.000000001);  
  //initForFreqSweepHighpassBilinear12(2.0);  
  //initForFreqSweepHighpassBilinear12(0.000000001);  
  //initForFreqSweepBandpassConstSkirtBilinear(2.0);
  //initForFreqSweepBandpassConstSkirtBilinear(0.000000001);  // investigate limiting case for q -> 0
  //initForFreqSweepBandrejectBilinear(2.0);
  //initForFreqSweepBandrejectBilinear(0.000000001);

  //initForFreqSweepPeakBilinear2(2.0, 2.0);
  //initForFreqSweepPeakBilinear2(2.0, 1.0);                  // g == 1 does not need special treatment   
  //initForFreqSweepPeakBilinear2(0.000000001, 3.0);          // investigate limiting case for q -> 0, g != 0
  //initForFreqSweepPeakBilinear2(2.0, 0.000000001);          // investigate limiting case for g -> 0, q != 0
  //initForFreqSweepPeakBilinear2(0.000000001, 0.000000001);  // investigate limiting case for q -> 0, g -> 0

  //initForFreqSweepLowShelfBilinear2(2.0, 2.0);
  //initForFreqSweepLowShelfBilinear2(2.0, 1.0);
  //initForFreqSweepLowShelfBilinear2(0.000000001, 3.0);          // investigate limiting case for q -> 0, g != 0

  //initForFreqSweepHighShelfBilinear2(2.0, 2.0);

  initForFreqSweepAllpassBilinear2(2.0);
}
void BiquadDesignerPlotTest::plotResult()
{
  RAPT::rsAssert(false, "plotting code needs update");
  //Plotter::plotData(numFramesToProcess, xAxis, outputs[0][0], outputs[0][1], outputs[0][2], outputs[0][3], outputs[0][4]);
}

void BiquadDesignerPlotTest::initForFreqSweepLowpassBilinear6()
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Lowpass, 6 dB/oct, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepHighpassBilinear6()
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Highpass, 6 dB/oct, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepLowpassBilinear12(double q)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Lowpass, 12 dB/oct, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);  
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepHighpassBilinear12(double q)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Highpass, 12 dB/oct, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);  
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepBandpassConstSkirtBilinear(double q)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Bandpass, const. skirt, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepBandrejectBilinear(double q)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Bandreject, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepPeakBilinear(double q, double g)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Peak, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  RAPT::rsArrayTools::fillWithValue(      inputs[0][2], numFramesToProcess, g);  
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepLowShelfBilinear2(double q, double g)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Low Shelf, 2nd order, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  RAPT::rsArrayTools::fillWithValue(      inputs[0][2], numFramesToProcess, g);  
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepHighShelfBilinear2(double q, double g)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "High Shelf, 2nd order, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  RAPT::rsArrayTools::fillWithValue(      inputs[0][2], numFramesToProcess, g);  
  xAxis = inputs[0][0];
}
void BiquadDesignerPlotTest::initForFreqSweepAllpassBilinear2(double q)
{
  ((romos::BiquadDesigner*) moduleToTest)->setParameter("Mode", "Allpass, 2nd order, BLT");
  RAPT::rsArrayTools::fillWithRangeLinear(inputs[0][0], numFramesToProcess, fMin, fMax);  
  RAPT::rsArrayTools::fillWithValue(      inputs[0][1], numFramesToProcess, q);  
  xAxis = inputs[0][0];
}
