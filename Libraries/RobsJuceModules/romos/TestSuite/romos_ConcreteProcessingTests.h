#ifndef romos_ConcreteProcessingTests_h
#define romos_ConcreteProcessingTests_h

//#include "romos_ProcessingTest.h"
//#include "romos_TestModuleBuilder.h"
//#include "romos_TestEventGenerator.h"
//#include "romos_GenerateDesiredOutput.h"

namespace rsTestRomos
{

/** This file contains concrete subclasses of ProcessingTest. 

\todo: Maybe split into AtomicProcessingTests and ContainerProcessingTests */


class IdentityTest : public ProcessingTest
{
public:
  IdentityTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

class AdderTest : public ProcessingTest
{
public:
  AdderTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

class Adder3Test : public ProcessingTest
{
public:
  Adder3Test();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

class Adder4Test : public ProcessingTest
{
public:
  Adder4Test();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

class Adder5Test : public ProcessingTest
{
public:
  Adder5Test();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

class SubtractorTest : public ProcessingTest
{
public:
  SubtractorTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

class WrappedAdderTest : public ProcessingTest
{
public:
  WrappedAdderTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

class UnitDelayTest : public ProcessingTest
{
public:
  UnitDelayTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};


class NoiseGeneratorTest : public ProcessingTest
{
public:
  NoiseGeneratorTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
  virtual void connectTestModuleToInputFeederModules() { } // do nothing - avoid to connect the seed input
};



/** A container with 2 inputs and 3 outputs, computing the sum, difference and product of the 
inputs. If it fails, it may hint that something about containers with multiple I/O pins is 
broken. */
class SumDiffProdTest : public ProcessingTest
{
public:
  SumDiffProdTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** Same as SumDiffProdTest but wraps the container that is used there into yet another container. 
If it fails, it may hint that something about container nesting (with multiple I/O pins) is 
broken. */
class WrappedSumDiffProdTest : public ProcessingTest
{
public:
  WrappedSumDiffProdTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** Tests an adder module with dynamic number of inputs (aka summation module). The adder is inside
a container with a single input and a single output. It should have 10 connected input pins and one 
additional disconnected input pin. The test also removes connections to check if the behavior for 
dynamically deleting pins is as expected. */
class WrappedAdderNTest : public ProcessingTest
{
public:
  WrappedAdderNTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
  virtual bool runTest();
  virtual void removeConnection(int index);
  virtual int  getAdderNumInputPins();
  virtual int  getAdderNumConnectedInputPins();
};


/** ... */
class SummedDiffsTest : public ProcessingTest
{
public:
  SummedDiffsTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** A 2-point moving average filter made from basic modules (unit-delay, multiply, etc.). */
class MovingAverageTest : public ProcessingTest
{
public:
  MovingAverageTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** A leaky integrator made from basic modules (unit-delay, multiply, etc.). If this fails, it may 
hint that something about feedback-paths is broken. */
class LeakyIntegratorTest : public ProcessingTest
{
public:
  LeakyIntegratorTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** A container with a single input and single output and two identity modules in between. The 
second identity module is placed such that it is evaluated before the first one such that there is 
an implicit delay in the connection between the two. If this fails, it may hint that something 
about connections with implicit delay is broken. */
class DelayedConnectionTest : public ProcessingTest
{
public:
  DelayedConnectionTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};


/** Similar to LeakyIntegratorTest, but moving one child-module into a position so as to create an 
additional implicit unit-delay inside the feedback path. */
class LeakyIntegratorDoubleDelayTest : public ProcessingTest
{
public:
  LeakyIntegratorDoubleDelayTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** A "filter" that computes the sum, difference and product from a LeakyIntegrator and 
MovingAverage. */
class TestFilter1Test : public ProcessingTest
{
public:
  TestFilter1Test();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};



/** A biquad made from basic modules (unit-delay, multiply, etc.).  */
class BiquadMacroTest : public ProcessingTest
{
public:
  BiquadMacroTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** Test for the atomic biquad module. */
class BiquadAtomicTest : public ProcessingTest
{
public:
  BiquadAtomicTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** Test for a biquad realized by the Formula module. */
class BiquadFormulaTest : public ProcessingTest
{
public:
  BiquadFormulaTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};

/** Test for the 1-in, 1-out formula. */
class Formula1In1OutTest : public ProcessingTest
{
public:
  Formula1In1OutTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};


class Formula_N_1Test : public ProcessingTest
{
public:
  Formula_N_1Test();
  virtual bool runTest() override;
protected:
  //virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};



// factor out a class BiquadTest - the fillDesiredOutputSignalArrays is the same for both classes


/** A module that creates filter impulse response blips at frequencies defined by the note-events 
that are thrown at it. It uses using a NoteFrequency and a NoteGate module as well as some 
others. */
class BlipTest : public ProcessingTest
{
public:
  BlipTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};


/** A monophonic ConstantModule value (of 1) that goes into a polyphonic UnaryMinusModule which 
then goes to the output. Each voice in the polyphonic minus-module should receive a copy of the 
-1 value. And the multiple copies of the -1 should also appear in the output-pins of the container
(unless the container itself monophonic in which case it has only one output voice (index 0) which 
should also contain the -1 (because we don't have a voice combiner in front of the output module).

\todo: make this stuff with the monophonic output actually true */
class MonoToPolyTest : public ProcessingTest
{
public:
  MonoToPolyTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};


class VoiceCombinerTest : public ProcessingTest
{
public:
  VoiceCombinerTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
};


class GateAndKillTest : public ProcessingTest
{
public:
  GateAndKillTest();
protected:
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);
  //virtual void runTest();
};



// next: UnitDelay, MonoToPoly, PolyToMono, 
// In3Out5 -> then with d1 and d4 set to mono, then set to poly, d4 shifted leftward ...try to cover all cases...
// the VoiceAllocatorTests (stealing, etc.)
// then VoiceKiller tests
// then some more sophisticated processing with events (EventHandlingTests maybe?)

} // end namespace romos

#endif 
