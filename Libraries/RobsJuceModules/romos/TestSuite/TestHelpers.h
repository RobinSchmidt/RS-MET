#ifndef TestHelpers_h
#define TestHelpers_h


//#include "../romos.h"
//#include "../../../romos/modules/romos_ModuleCreation.h"

//using namespace rosic;  // get rid!
//using namespace romos;

namespace rsTestRomos
{

// arrays for I/O signals:
static const int maxNumVoices = 8;
static const int maxNumIns    = 10;
static const int maxNumOuts   = 10;
static const int maxNumFrames = 1000;
static const int N            = 50;
extern double x[maxNumVoices][maxNumIns][maxNumFrames];      // inputs
extern double y[maxNumVoices][maxNumOuts][maxNumFrames];     // outputs
extern double d[maxNumVoices][maxNumOuts][maxNumFrames];     // desired outputs
extern double t[maxNumFrames];                               // timeline in samples for plots
extern double *px0[maxNumIns];                               // pointers to the inputs of voice 0
extern double *py0[maxNumOuts];                              // pointers to the outputs of voice 0
extern double *pd0[maxNumOuts];                              // pointers to the desired(?) outputs of voice 0

extern double *px[maxNumVoices][maxNumIns];                  // pointers to the inputs of all voices 
extern double *py[maxNumVoices][maxNumOuts];                 // pointers to the outputs of all voices 
extern double *pd[maxNumVoices][maxNumOuts];                 // pointers to the desired(?) outputs of all voices 

extern double **ppx[maxNumVoices];
extern double **ppy[maxNumVoices];
extern double **ppd[maxNumVoices];
// why extern?


/** Returns a vector containing a note-on and a corresponding note-off (indicated by velocity == 0).
If zero is passed as duration, it will return an empty vector. */
std::vector<NoteEvent> generateNoteOnOffPair(unsigned int key, unsigned int velocity,                                            
  unsigned int deltaFramesForNoteOn, unsigned int durationInFrames);

/** Generates a bunch of simultaneous notes with equal spacing between the individual notes, given 
by noteSpacing. */
std::vector<NoteEvent> generateSimultaneousNotes(unsigned int key, unsigned int velocity, 
  unsigned int deltaFramesForNoteOn, unsigned int durationInFrames, 
  unsigned int numNotes, unsigned int noteSpacing);

/** Merges two vectors of events. The result will be sorted by the time of occurence of the events. Simlutaneous events may occur in any
order wihtin the array (\todo maybe have well defined ordering criterions in these cases, too). */
std::vector<NoteEvent> mergeEvents(const std::vector<NoteEvent> &eventArray1, const std::vector<NoteEvent> &eventArray2);

/** Converts an array of note-on events with corresponding note-off events into an array that contains only the note-ons and another array 
that contains the corresponding durations (in sample-frames). */
void convertNoteEventsToStartsAndDurations(const std::vector<NoteEvent> &events, std::vector<NoteEvent> &noteOns, 
                                           std::vector<int> &durations);

/** Given a vector of note-events and a particular note-on event, this function returns the index of the note-off event that corresponds to
the given note-on event. If none is found, -1 will be returned. */
int findIndexOfMatchingNoteOff(const std::vector<NoteEvent> &events, NoteEvent noteOnEvent);

void sortEventsByTimeOfOccurence(const std::vector<NoteEvent> eventArray);

// get rid - replace with own prng
inline double random(double min, double max)
{
  double tmp = (1.0/RAND_MAX) * rand();  // between 0...1
  return RAPT::rsLinToLin(tmp, 0.0, 1.0, min, max);
}

// generation of test input signals:

void initializeInputSequences(); 
void zeroDesiredOutputs();
void writeInputSequencesToFile();

void generateSilence(int numChannels, int numFrames, double **x);
void generateImpulse(int N, double *x);
void generateRandomSequence(int N, double *x, double xMin, double xMax, int seed);

// computation of correct reference output signals to be matched:
void getDesiredOutputForUnitDelay(int N, double *x, double *d);
void getDesiredOutputForSummedDiffs(    int N, double **inputs, double **desiredOutputs);
void getDesiredOutputForMovingAverage(int N, double *x, double *b0, double *b1, double *desiredOutput);
void getDesiredOutputForLeakyIntegrator(int N, double *x, double *c, double *desiredOutput);
void getDesiredOutputForLeakyIntegratorDoubleDelay(int N, double *x, double *c, double *desiredOutput);
void getDesiredOutputForTestFilter1(int N, double *x, double *b0, double *b1, double *c, double *dSum, double *dDiff, double *dProd);
void getDesiredOutputForBiquad(int N, double *x, double *b0, double *b1, double *b2, double *a1, double *a2, double *y);
void getDesiredOutputForFilterBlip(int N, double frequency, double q, double *desiredOutput);

//void getDesiredOutputForGatedNoteFrequencies(int N, std::vector<NoteEvent> *events, double ***desiredOutputs, 
//                                             bool containerIsPolyphonic, bool noteFreqModuleIsPolyphonic);
  // we need 4 cases: 
  // -NoteFreq poly / output mono -> sum
  // -NoteFreq poly / output poly -> each voice output gets its own frequency output
  // -NoteFreq mono / output mono -> only first voice gets the first frequency
  

// others:

/** Sets the polyphony flag for the passed module and optionally recursively for its child-modules (if any) */
void setModulePolyphony(romos::Module *module, bool shouldBePolyphonic, bool recursivelyForChildren);

/** Checks whether or not y[c][n] == d[c][n] where indices c and n run from 0 to numChannels-1 and numFrames-1 respectively. */
bool checkResult(double **y, double **d, int numChannels, int numFrames, double tolerance);

/** Similar to checkResult but also prints the result of the test to the console. */
bool checkAndPrintResult(double **y, double **d, int numChannels, int numFrames, const char *testName, double tolerance);


/** Block-processes the passed module  ....*/
bool checkBlockProcessingAndPrintResult(romos::Module *module, double ***x, double ***y, double ***d, int maxNumFramesToProcess, 
                                        int numTests, char *testName, double tolerance);
  // is this function redundant now?


/** Checks all 4 processing functions and prints the results of the test to the console. You need to pass the input signals for the voices 
in x, where the indexing is: x[voiceIndex][pinIndex][frameIndex]. Likewise, you need to pass a pointer y where the outputs fo the module 
are to be stored and a pointer d containing the desired outputs.  */
bool checkProcessingFunctionsAndPrintResults(romos::Module *module, int numVoicesToCheck, int numFrames, 
                                             double ***x, double ***y, double ***d, 
                                             double tolerance, char *testName, std::vector<NoteEvent> *events = NULL);

/** Checks the per-frame, monophonic processing function and prints the result to the console. */
bool checkProcessingInFramesMonoAndPrintResult(romos::Module *module, int numFrames, double ***x, double ***y, double ***d, 
                                               double tolerance, char *testName, std::vector<NoteEvent> *events = NULL);

/** Checks the per-block, monophonic processing function and prints the result to the console.. */
bool checkProcessingInBlocksMonoAndPrintResult(romos::Module *module, int numFrames, double ***x, double ***y, double ***d, 
                                               double tolerance, char *testName, std::vector<NoteEvent> *events = NULL);

/** Checks the per-frame, polyphonic processing function and prints the result to the console.. */
bool checkProcessingInFramesPolyAndPrintResult(romos::Module *module, int numVoicesToCheck, int numFrames, 
                                               double ***x, double ***y, double ***d, 
                                               double tolerance, char *testName, std::vector<NoteEvent> *events = NULL);

/** Checks the per-block, polyphonic processing function and prints the result to the console.. */
bool checkProcessingInBlocksPolyAndPrintResult(romos::Module *module, int numVoicesToCheck, int numFrames, 
                                               double ***x, double ***y, double ***d, 
                                               double tolerance, char *testName, std::vector<NoteEvent> *events = NULL);

/** Sets the passed module into poly- or monophonic mode and lets it process a number of sample-frames using the module's processFrame 
function pointer (which, at this point, should resolve to the appropriate processPoly or processMono function). The input-signals and 
memory for the output-signals must be passed as pointers using the indexing convention inputs[voiceIndex][pinIndex][frameIndex], 
likewise for outputs. In the monophonic case, the 0th voice will be used only. */
void processModuleInFrames(romos::Module *module, int numFrames, double ***inputs, double ***outputs, std::vector<NoteEvent> *events, 
                           bool polyphonic);

/** Similar to processModuleInFrames, but uses the processBlock function pointer. The sizes of the individual block will be chosen 
randomly between 1 and the maximum possible blocksize (as determined by the allocated memory for I/O blocks). */
void processModuleInBlocks(romos::Module *module, int numFrames, double ***inputs, double ***outputs, std::vector<NoteEvent> *events, 
                           bool polyphonic);

void exchangeModulePositions(romos::Module *module1, romos::Module *module2);

void randomizeContainment(romos::Module *module);
void printModuleStructure(romos::Module *module, int indent);

}


#endif
