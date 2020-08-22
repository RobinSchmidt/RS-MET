#ifndef romos_ProcessingTest_h
#define romos_ProcessingTest_h

//#include "../framework/romos_NoteEvent.h"
//#include "../framework/romos_ModuleFactory.h"
//#include "romos_TestEventGenerator.h"
//#include "romos_UnitTest.h"

namespace rsTestRomos
{

/** This class is the baseclass for all processing test classes (i.e. test, that produce signals) 
for the romos library. These tets produce output signals by invoking the 4 processing functions of
modules (poly/mono, frame/block-wise) and compare the resulting output signals to known correct 
output signals. The baseclass provides all the facilities that are needed for such tests, like 
members which can be used for input- and output- signals, etc.  */

class ProcessingTest : public UnitTest
{

public:

  //-------------------------------------------------------------------------------------------------------------------------------------
  // construction/destruction:

  /** Constructor. */
  ProcessingTest(const char *testName);

  /** Destructor. */
  virtual ~ProcessingTest();

  //-------------------------------------------------------------------------------------------------------------------------------------
  // setup:

  /** Accepts a vector of events that will be triggered during the processing functions. These events should be passed before calling
  runTests. */
  virtual void setEventsToOccurDuringProcessing(const std::vector<romos::NoteEvent> eventsToOccur)
  {
    events = eventsToOccur;
  }

  /** Sets the polyphony of our modueToTest and recursively for all it's child-modules when it is a container. */
  virtual void setTestModulePolyphonyRecursively(bool shouldBePolyphonic);

  /** Sets up all the members for either polyphonic or monophonic processing. */
  virtual void setTestPolyphonic(bool shouldBePolyphonic);

  //-------------------------------------------------------------------------------------------------------------------------------------
  // information output:

  /** Plots the desired output and the actual output fro the given voice and output pin-index. */
  virtual void plotDesiredAndActualOutput(int voiceIndex, int pinIndex, int numFramesToPlot, int startFrame = 0);

  //=====================================================================================================================================

protected:

  //-------------------------------------------------------------------------------------------------------------------------------------
  // running the test:

  /** This function runs the tests and returns true when the tests have passed and false when the 
  tests have failed. */
  virtual bool runTest() override;

  /** Overriden in order to trigger plotting of tareget signal and actual output signal in case of
  failure. */
  virtual void handleTestResult(bool didTestPass) override;

  /** Function to do some setup stuff like connecting the moduleToTest member to some feeder/retriever modules etc. It is called as first
  thing in runTest, and when you override runTest in a subclass, you should most probabaly call it there, too. */
  virtual void initTest();

  /** Coonects the input pins of our moduleToTest with the output pins of a bunch of IdentityModules that we have here in order to feed
  the moduleToTest with the test-signal. This should be called after the moduleToTest has been created and before the tests are actually
  run. */
  virtual void connectTestModuleToInputFeederModules();

  /** Coonects the output pins of our moduleToTest with the input pins of a bunch of IdentityModules that we have here in order to
  retrieve the output signals that the moduleToTest has produced. This should be called after the moduleToTest has been created and
  before the tests are actually run. */
  virtual void connectTestModuleToOutputRetrieverModules();


  /** Checks whether or not the output signals match the desired output signals, either for the 0th voice only (if the argument is
  false) or for all voices up to maxNumVoices (if the argument is true). */
  virtual bool doOutputsMatchDesiredOutputs(bool polyphonic);

  //-------------------------------------------------------------------------------------------------------------------------------------
  // filling of the internal arrays with default data:

  /** Fills the arrays that contain the input signals for all voices and all pins with pseudo-random numbers, using the passed seed for
  pseudo-random number generator. */
  virtual void fillInputSignalArraysRandomly(int seed);

  /** Fills the arrays that contain the input signals for all voices and all pins with zeros.  */
  virtual void clearInputSignalArrays();

  /** Fills the arrays that contain the output signals for all voices and all pins with zeros.  */
  virtual void clearOutputSignalArrays();

  /** Fills the arrays that contain the desired output signals for all voices and all pins with zeros.  */
  virtual void clearDesiredOutputSignalArrays();

  /** Fills the array that contains our time axis for plots with sample indices. */
  virtual void fillTimeAxisWithSampleIndices();

  //-------------------------------------------------------------------------------------------------------------------------------------
  // methods that subclasses can override:

  /** Subclasses can override this to fill the arrays that contain the input signals for all voices and all pins with some special test
  signal to be used for testing the subclass at hand. If you don't override it, the baseclass implementation will fill the input arrays
  with pseudo-random numbers. */
  virtual void fillInputSignalArraysWithTestSignal();

  /** Subclasses can override this function and fill the member "desiredOutputs" there with the correct output signals for all voices
  that the "moduleToTest" is supposed to produce when the module's polyphony is set to "testModuleIsPolyphonic". */
  virtual void fillDesiredOutputSignalArrays(bool testModuleIsPolyphonic);

  //-------------------------------------------------------------------------------------------------------------------------------------
  // processing of the module's output:

  /** Lets our moduleToTest process a number of sample-frames using the module's processFrame function pointer. */
  virtual void processModuleInFrames();

  /** Similar to processModuleInFrames, but uses the processBlock function pointer. The sizes of the individual block will be chosen
  randomly between 1 and the maximum possible blocksize (as determined by the allocated memory for I/O blocks). */
  virtual void processModuleInBlocks();

  /** Processes a number sample-frames starting from startIndex using the frame-based function-pointer. The functions assumes that no
  events occur within this block. Used internally by processModuleInFrames.
  \todo exchange order of the arguments - this needs some care...
  */
  virtual void processModuleInFramesNoEvents(int numFrames, int startIndex);

  /** Similar to processModuleInFramesNoEvents but uses the block-based function pointer. */
  virtual void processModuleInBlocksNoEvents(int numFrames, int startIndex);

  /** Processes a single sample-frame at the given index. Used internally by processModuleInFramesNoEvents. */
  virtual void processFrame(int frameIndex);

  /** Processes a block of sample-frames starting at the given index. Used internally by processModuleInBlocksNoEvents. */
  virtual void processBlock(int blockStart, int blockSize);

  /** Establishes a block of input samples at the moduleToTest's input pins by copying the data from our "inputs" member array into the
  appropriate memory slots of the module's input pins. Used internally by processFrame and processBlock. */
  virtual void establishInputBlock(int blockStart, int blockSize);

  /** Retrieves a block of output samples from the moduleToTest's output pins by copying the data from the memory slots of the module's
  output pins into the appropriate positions in our "outputs" member array. Used internally by processFrame and processBlock. */
  virtual void retrieveOutputBlock(int blockStart, int blockSize);



  //-------------------------------------------------------------------------------------------------------------------------------------
  // validation:



  //-------------------------------------------------------------------------------------------------------------------------------------
  // data:

  static const int maxNumVoices  = 8;
  static const int maxNumInputs  = 10;
  static const int maxNumOutputs = 10;
  static const int maxNumFrames  = 10000;

  romos::Module *moduleToTest;
  romos::Module *inputFeederModules[maxNumInputs];
  romos::Module *outputRetrieverModules[maxNumOutputs];

  //romos::ContainerModule *inputDummyModule;   // used to feed the input into the actually tested module
  //romos::IdentityModule  *outputDummyModules




  //bool monoInFramesPassed, monoInBlocksPassed, polyInFramesPassed, polyInBlocksPassed;

  std::vector<romos::NoteEvent> events;



  int numVoicesToUse;
  //int numInputs;
  //int numOutputs;
  int numFramesToProcess;
  // number of frames to process in this test - we may want to avoid to always fill the arrays up to the end for running the tests faster
  // ->

  double inputs[maxNumVoices][maxNumInputs][maxNumFrames];
  double outputs[maxNumVoices][maxNumOutputs][maxNumFrames];
  double desiredOutputs[maxNumVoices][maxNumOutputs][maxNumFrames];
  double timeAxis[maxNumFrames];

  double tolerance;

};

} // end namespace romos

#endif 
