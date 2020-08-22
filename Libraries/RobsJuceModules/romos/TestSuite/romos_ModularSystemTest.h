#ifndef romos_ModularSystemTest_h
#define romos_ModularSystemTest_h

//#include "../framework/romos_NoteEvent.h"
//#include "../framework/romos_ModuleFactory.h"
//#include "romos_TestEventGenerator.h"
//#include "romos_UnitTest.h"
//#include "../romos.h"

namespace rsTestRomos
{

  /**

  This class is the baseclass for all test classes that testn the modular system as a whole.

  */

  class ModularSystemTest : public UnitTest
  {

  public:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // construction/destruction:

    /** Constructor. */
    ModularSystemTest(const char *testName);

    /** Destructor. */
    virtual ~ModularSystemTest();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // setup:

    /** Accepts a vector of events that will be triggered during the processing functions. These events should be passed before calling
    runTests. */
    virtual void setEventsToOccurDuringProcessing(const std::vector<romos::NoteEvent> eventsToOccur)
    {
      events = eventsToOccur;
    }

    /** In this function, subclasses should create all the child-modules that the TopLevelModule should contain and connect them to the 
    TopLevelModule's input- and output modules, if necessary. */
    virtual void createAndConnectTestChildModules() = 0;

    //-------------------------------------------------------------------------------------------------------------------------------------
    // information output:

    /** Plots the desired output and the actual output for both channels. */
    virtual void plotDesiredAndActualOutput(int numFramesToPlot, int startFrame = 0);

    /** Plots the error (desiredOutput - actualOutput) for both channels. */
    virtual void plotOutputErrors(int numFramesToPlot, int startFrame = 0);



    //=====================================================================================================================================

  protected:

    //-------------------------------------------------------------------------------------------------------------------------------------
    // running the test:

    /** This function runs the tests and returns true when the tests have passed and false when the tests have failed. */
    virtual bool runTest();

    /** Function to do some setup stuff like connecting the moduleToTest member to the inputs/outputs of the TopLevelModule, etc. */
    virtual void initTest();

    /** Checks whether or not the output signals match the desired output signals. */
    virtual bool doOutputsMatchDesiredOutputs();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // filling of the internal arrays with default data:

    /** Fills the arrays that contain the input signals for the TopLevelModule with pseudo-random numbers, using the passed seed for the
    pseudo-random number generator. */
    virtual void fillInputSignalArraysRandomly(int seed);

    /** Fills the arrays that contain the input signals for the TopLevelModule with zeros.  */
    virtual void clearInputSignalArrays();

    /** Fills the arrays that contain the output signals of the TopLevelModule with zeros.  */
    virtual void clearOutputSignalArrays();

    /** Fills the arrays that contain the desired output signals with zeros.  */
    virtual void clearDesiredOutputSignalArrays();

    /** Fills the array that contains our time axis for plots with sample indices. */
    virtual void fillTimeAxisWithSampleIndices();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // methods that subclasses can override:

    /** Subclasses can override this to fill the arrays that contain the input signals with some special test signal to be used for 
    testing the subclass at hand. If you don't override it, the baseclass implementation will fill the input arrays with pseudo-random 
    numbers. */
    virtual void fillInputSignalArraysWithTestSignal();

    /** Subclasses can override this function and fill the member "desiredOutputs" there with the correct output signals that the 
    TopLevelModule is supposed to produce. */
    virtual void fillDesiredOutputSignalArrays();

    //-------------------------------------------------------------------------------------------------------------------------------------
    // processing of the module's output:


    virtual void processOutputSignal(bool useSinglePrecision);

    virtual void handleEvent(romos::NoteEvent eventToHandle);





    //-------------------------------------------------------------------------------------------------------------------------------------
    // data:

    static const int maxSignalLength      = 44100;
    static const int maxNumInputChannels  = 2;
    static const int maxNumOutputChannels = 2;

    int numInputChannels;
    int numOutputChannels;
    int signalLength;
    int blockSize;

    //romos::TopLevelModule *topLevelModule;
    romos::Liberty   modularSystem;

    romos::TopLevelModule    *topLevelModule;           // == modularSystem.topLevelModule - for convenient access
    romos::AudioInputModule  *inModuleL,  *inModuleR;   // pointers to input modules inside the topLevelModule
    romos::AudioOutputModule *outModuleL, *outModuleR;  // pointers to output modules inside the topLevelModule

    std::vector<romos::NoteEvent> events;

    double inputs        [maxNumInputChannels] [maxSignalLength];
    double outputs       [maxNumOutputChannels][maxSignalLength]; 
    double desiredOutputs[maxNumOutputChannels][maxSignalLength];  
    double outputErrors  [maxNumOutputChannels][maxSignalLength]; 
    double timeAxis                            [maxSignalLength];  

    float  inputsFloat   [maxNumInputChannels] [maxSignalLength];
    float  outputsFloat  [maxNumOutputChannels][maxSignalLength]; 


    double tolerance;

  };

} // end namespace romos

#endif 
