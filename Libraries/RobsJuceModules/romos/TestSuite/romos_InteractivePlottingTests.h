#ifndef romos_InteractivePlottingTests_h
#define romos_InteractivePlottingTests_h

//#include "romos_ProcessingTest.h"
//#include "romos_TestModuleBuilder.h"
//#include "romos_TestEventGenerator.h"
//#include "romos_GenerateDesiredOutput.h"

namespace rsTestRomos
{

  /**

  This file contains concrete subclasses of ProcessingTest that are meant to be used to visualize the output signals of modules in a 
  GNUPlot window. 

  */


  class InteractivePlotTest : public ProcessingTest
  {
  public:
    InteractivePlotTest(const char *testName);
    virtual void runTestAndPlotResults();
  protected:

    /** To be overriden by subclasses to plot the relevant results. */
    virtual void plotResult() = 0;

    double *xAxis;
  };




  /** Plots the output waveform of the BandlimitedImpulseTrain module vs a time axis. */
  class BandlimitedImpulseTrainPlotTest : public InteractivePlotTest
  {
  public:
    BandlimitedImpulseTrainPlotTest();
  protected:
    virtual void fillInputSignalArraysWithTestSignal();
    virtual void plotResult();
    virtual void initForConstantFreq(double freq);
    virtual void initForLinearFreqSweep(double minFreq, double maxFreq);
  };



  /** Plots the output waveform of the SawOscillator module vs a time axis. */
  class SawOscillatorPlotTest : public InteractivePlotTest
  {
  public:
    SawOscillatorPlotTest(const char *testName = "SawOscillatorPlot");
  protected:
    virtual void fillInputSignalArraysWithTestSignal();
    virtual void plotResult();
    virtual void initForConstantFreq(double freq, double phase = 0.0);
    virtual void initForLinearFreqSweep(double minFreq, double maxFreq);
    virtual void initForFreqSwitch(double freq1, double freq2, int switchTimeInSamples);
  };


  class DualBlitSawOscillatorPlotTest : public SawOscillatorPlotTest
  {
  public:
    DualBlitSawOscillatorPlotTest();
  };





  /** Plots the ADSR envelope vs a time axis in samples. different settings for the various time values can be checked and different 
  sequences of events can be fed into the envelope generators for checking double-trigger and early-release behavior. */
  class EnvelopeADSRPlotTest : public InteractivePlotTest
  {
  public:
    EnvelopeADSRPlotTest();
  protected:
    virtual void fillInputSignalArraysWithTestSignal();
    virtual void plotResult();

    virtual void initEnvelopeTimesAllNonZero();
    virtual void initEnvelopeTimesZeroAttack();
    virtual void initEnvelopeTimesZeroDecay();
    virtual void initEnvelopeTimesZeroAttackAndDecay();
    virtual void initEnvelopeTimesZeroRelease();
    //virtual void initEnvelopeFractionalAttack();
    virtual void setTimesInSecondsFromTimesInSamples(); // called internally from the setEnvelopeTimes[...] functions

    virtual void initEventsForNormalCase();      // normal case, goes through attack, decay, sustain, release without interruption
    virtual void initEventsForDoubleAttack1();   // second note-on occurs when the envelope is already running (during attack)
    virtual void initEventsForDoubleAttack2();   // second note-on occurs when the envelope is already running (exactly at end of attack)
    virtual void initEventsForDoubleAttack3();   // second note-on occurs when the envelope is already running (during decay)
    virtual void initEventsForDoubleAttack4();   // second note-on occurs when the envelope is already running (exactly at end of decay)
    virtual void initEventsForDoubleAttack5();   // second note-on occurs when the envelope is already running (during sutain)
    virtual void initEventsForDoubleAttack6();   // second note-on occurs when the envelope is already running (exactly at end of sustain)
    virtual void initEventsForDoubleAttack7();   // second note-on occurs when the envelope is already running (during release)
    virtual void initEventsForEarlyRelease1();   // note-off occurs before sustain is reached (during attack)
    virtual void initEventsForEarlyRelease2();   // note-off occurs before sustain is reached (exactly at end of attack)
    virtual void initEventsForEarlyRelease3();   // note-off occurs before sustain is reached (during decay)

    int    attackSamples, decaySamples, sustainSamples, releaseSamples, noteLength;
    double attackSeconds, decaySeconds, sustainSeconds, releaseSeconds;
    double attackShape, decayShape, releaseShape, timeScale;
    double sustainValue;
  };


  /** Plots the biquad coefficients vs. frequency. Inside the test, you can choose the filter mode and fix the other two inputs
  for Q and gain. */
  class BiquadDesignerPlotTest : public InteractivePlotTest
  {
  public:
    BiquadDesignerPlotTest();
  protected:
    virtual void fillInputSignalArraysWithTestSignal();
    virtual void plotResult();

    virtual void initForFreqSweepLowpassBilinear6();
    virtual void initForFreqSweepHighpassBilinear6();
    virtual void initForFreqSweepLowpassBilinear12(         double q);
    virtual void initForFreqSweepHighpassBilinear12(        double q);
    virtual void initForFreqSweepBandpassConstSkirtBilinear(double q);
    virtual void initForFreqSweepBandrejectBilinear(        double q);
    // insert 2nd order allpass here
    virtual void initForFreqSweepPeakBilinear(              double q, double g);
    virtual void initForFreqSweepLowShelfBilinear2(         double q, double g);
    virtual void initForFreqSweepHighShelfBilinear2(        double q, double g);
    virtual void initForFreqSweepAllpassBilinear2(          double q);

    double fMin, fMax;
  };



} // end namespace romos

#endif 
