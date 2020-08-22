#include "AutomaticTests.h"

#if defined(_MSC_VER)
#include <io.h>              // works on msc
#elif defined(__APPLE__)
#include <sys/uio.h>         // works on mac
#else
//#include<sys/io.h>           // works on linux
#include<fcntl.h>
#endif
#include<sys/stat.h>

// try to get rid of these global variables - or at least wrap them into a namespace to not clutter
// the global namespace
namespace rsTestRomos
{
double x[maxNumVoices][maxNumIns][maxNumFrames];      // inputs
double y[maxNumVoices][maxNumOuts][maxNumFrames];     // outputs
double d[maxNumVoices][maxNumOuts][maxNumFrames];     // desired outputs
double t[maxNumFrames];                               // timline in samples for plots
double *px0[maxNumIns];                               // pointers to the inputs of voice 0
double *py0[maxNumOuts];                              // pointers to the outputs of voice 0
double *pd0[maxNumOuts];                              // pointers to the outputs of voice 0

double *px[maxNumVoices][maxNumIns];                  // pointers to the inputs of all voices
double *py[maxNumVoices][maxNumOuts];                 // pointers to the outputs of all voices
double *pd[maxNumVoices][maxNumOuts];                 // pointers to the outputs of all voices

double **ppx[maxNumVoices];
double **ppy[maxNumVoices];
double **ppd[maxNumVoices];
}
using namespace rsTestRomos;  // try to get rid of this!

//-----------------------------------------------------------------------------------------------------------------------------------------
// generation of test event-arrays:

std::vector<NoteEvent> generateNoteOnOffPair(unsigned int key, unsigned int velocity,
                                             unsigned int deltaFramesForNoteOn, unsigned int durationInFrames)
{
  std::vector<NoteEvent> result;
  if( durationInFrames == 0 )
    return result; // return an empty vector - note-ons with simultaneous note-offs are discarded
  result.push_back(NoteEvent(deltaFramesForNoteOn,                    key, velocity));
  result.push_back(NoteEvent(deltaFramesForNoteOn + durationInFrames, key, 0));
  return result;
}

std::vector<NoteEvent> generateSimultaneousNotes(unsigned int key, unsigned int velocity,
                                                 unsigned int deltaFramesForNoteOn, unsigned int durationInFrames,
                                                 unsigned int numNotes, unsigned int noteSpacing)
{
  std::vector<NoteEvent> result;
  unsigned int currentKey = key;
  for(unsigned int i = 1; i <= numNotes; i++)
  {
    result = mergeEvents(result, generateNoteOnOffPair(currentKey, velocity, deltaFramesForNoteOn, durationInFrames));
    currentKey += noteSpacing;
  }
  return result;
}

std::vector<NoteEvent> mergeEvents(const std::vector<NoteEvent> &eventArray1, const std::vector<NoteEvent> &eventArray2)
{
  std::vector<NoteEvent> result;
  result.reserve(eventArray1.size() + eventArray2.size());
  unsigned int i;
  for(i = 0; i < eventArray1.size(); i++)
    result.push_back(eventArray1[i]);
  for(i = 0; i < eventArray2.size(); i++)
    result.push_back(eventArray2[i]);
  std::sort(result.begin(), result.end(), noteEventLessByDeltaFrames);
  return result;
}

void convertNoteEventsToStartsAndDurations(const std::vector<NoteEvent> &events, std::vector<NoteEvent> &noteOns, std::vector<int> &durations)
{
  unsigned int i;
  noteOns.clear();
  for(i = 0; i < events.size(); i++)
  {
    if( events[i].getVelocity() != 0 )
      noteOns.push_back(events[i]);
  }

  durations.reserve(noteOns.size());
  for(i = 0; i < noteOns.size(); i++)
  {
    int noteOffIndex = findIndexOfMatchingNoteOff(events, noteOns[i]);
    if( noteOffIndex > -1 )
      durations.push_back(events[noteOffIndex].getDeltaFrames() - noteOns[i].getDeltaFrames());
    else
      durations.push_back(INT_MAX);  // no matching note-off - potentially infinite duration, but INT_MAX is the largest we have
  }
}

int findIndexOfMatchingNoteOff(const std::vector<NoteEvent> &events, NoteEvent noteOnEvent)
{
  unsigned int i;
  unsigned int startIndex = rosic::findElement(events, noteOnEvent) + 1;
  for(i = startIndex; i < events.size(); i++)
  {
    if( events[i].isNoteOff() && events[i].getKey() == noteOnEvent.getKey() )
      return i;
  }
  return -1;
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// generation of test input signals:

void initializeInputSequences()
{
  RAPT::rsArrayTools::fillWithIndex(t, maxNumFrames);

  int v, c, n;

  RAPT::rsRandomUniform(-1.0, 1.0, 1);

  for(v = 0; v < maxNumVoices; v++)
  {
    for(c = 0; c < maxNumIns; c++)
    {
      for(n = 0; n < maxNumFrames; n++)
      {
        x[v][c][n] = RAPT::rsRandomUniform(-1.0, 1.0);
        y[v][c][n] = 0.0;
        d[v][c][n] = 0.0;
      }
      px[v][c] = &(x[v][c][0]);
      py[v][c] = &(y[v][c][0]);
      pd[v][c] = &(d[v][c][0]);
    }
    ppx[v] = &(px[v][0]);
    ppy[v] = &(py[v][0]);
    ppd[v] = &(pd[v][0]);
  }

  for(c = 0; c < maxNumIns; c++)
  {
    px0[c] = &(x[0][c][0]);
    py0[c] = &(y[0][c][0]);
    pd0[c] = &(d[0][c][0]);
  }
}

void zeroDesiredOutputs()
{
  memset(d, 0, maxNumVoices * maxNumOuts * maxNumFrames * sizeof(double));
}

void writeInputSequencesToFile()
{
  // restrict ranges (we don't want to wite the full arrays to the output):
  int numFrames   = RAPT::rsMin(50, maxNumFrames);
  int numVoices   = RAPT::rsMin(3,  maxNumVoices);
  int numChannels = RAPT::rsMin(4,  maxNumIns);

  rosic::rsString s;
  char tmp[8];

  for(int n = 0; n < numFrames; n++)
  {
    for(int v = 0; v < numVoices; v++)
    {
      for(int c = 0; c < numChannels; c++)
      {
        sprintf(tmp, "%+6.3f", x[v][c][n]);
        s += rosic::rsString(tmp) + rosic::rsString(" ");
        //s += rosic::rsString(x[v][c][n]) + rosic::rsString(" ");
      }
      s += rosic::rsString("  ");
    }
    s += rosic::rsString("\n");
  }

  FILE *file = fopen("c:\\tmp\\InputSequences.txt", "w");
  if( file == NULL )
  {
  #ifdef _MSC_VER
    _creat("c:\\tmp\\InputSequences.txt", S_IWRITE);
  #elif defined(__APPLE__)
    //_creat("c:\\tmp\\InputSequences.txt", S_IWRITE);
  #else
    creat("c:\\tmp\\InputSequences.txt", S_IWRITE);
  #endif
    file = fopen("c:\\tmp\\InputSequences.txt", "w");
  }
  if( file != NULL )
  {
    fprintf(file, "%s", s.getRawString());
    fclose(file);
  }
}


void generateSilence(int numChannels, int numFrames, double **x)
{
  for(int c = 0; c < numChannels; c++)
  {
    for(int n = 0; n < numFrames; n++)
      x[c][n] = 0.0;
  }
}

void generateImpulse(int N, double *x)
{
  x[0] = 1.0;
  for(int n = 1; n < N; n++)
    x[n] = 0.0;
}

void generateRandomSequence(int N, double *x, double xMin, double xMax, int seed)
{
  x[0] = RAPT::rsRandomUniform(xMin, xMax, seed);
  for(int n = 1; n < N; n++)
    x[n] = RAPT::rsRandomUniform(xMin, xMax);
}

//-----------------------------------------------------------------------------------------------------------------------------------------
// computation of correct reference output signals to be matched:

void getDesiredOutputForUnitDelay(int N, double *x, double *d)
{
  d[0] = 0.0;
  for(int n = 1; n < N; n++)
    d[n] = x[n-1];
}

void getDesiredOutputForSummedDiffs(int N, double **x, double **d)
{
  double d1, d2, d3, d4;
  for(int n = 0; n < N; n++)
  {
    d1 = x[0][n] - x[1][n];
    d2 = x[1][n] - x[0][n];
    d3 = x[2][n] - x[1][n];
    d4 = x[1][n] - x[2][n];
    d[0][n] = d1 + d2 + d4;
    d[1][n] = d2 + d4;
    d[2][n] = d1 + d3;
    d[3][n] = d1 + d3 + d4;
  }
}

void getDesiredOutputForMovingAverage(int N, double *x, double *b0, double *b1, double *d)
{
  d[0] = b0[0]*x[0] + b1[0]*0.0;  // assumes 0.0 as initial state of x[n-1] buffer
  for(int n = 1; n < N; n++)
    d[n] = b0[n]*x[n] + b1[n]*x[n-1];
}

void getDesiredOutputForLeakyIntegrator(int N, double *x, double *c, double *d)
{
  d[0] = c[0]*x[0] + (1.0-c[0])*0.0;  // assumes 0.0 as initial state of y[n-1] buffer
  for(int n = 1; n < N; n++)
    d[n] = c[n]*x[n] + (1.0-c[n])*d[n-1];
}

void getDesiredOutputForLeakyIntegratorDoubleDelay(int N, double *x, double *c, double *d)
{
  d[0] = c[0] * x[0] + (1.0 - c[0]) * 0.0;  // assumes 0.0 as initial state of y[n-1] buffer
  d[1] = c[1] * x[1] + (1.0 - c[1]) * 0.0;  // assumes 0.0 as initial state of y[n-2] buffer
  for(int n = 2; n < N; n++)
    d[n] = c[n]*x[n] + (1.0-c[n])*d[n-2];
}

void getDesiredOutputForTestFilter1(int N, double *x, double *b0, double *b1, double *c, double *dSum, double *dDiff, double *dProd)
{
  double *dMovAv   = new double[N];
  double *dLeakInt = new double[N];
  getDesiredOutputForMovingAverage(  N, x, b0, b1, dMovAv);
  getDesiredOutputForLeakyIntegrator(N, x, c,      dLeakInt);
  RAPT::rsArrayTools::add(     dMovAv, dLeakInt, dSum,  N);
  RAPT::rsArrayTools::subtract(dMovAv, dLeakInt, dDiff, N);
  RAPT::rsArrayTools::multiply(dMovAv, dLeakInt, dProd, N);
  delete[] dLeakInt;
  delete[] dMovAv;
}

void getDesiredOutputForBiquad(int N, double *x, double *b0, double *b1, double *b2, double *a1, double *a2, double *y)
{
  // assumes 0.0 as initial values of all buffers:
  y[0] = b0[0]*x[0];
  y[1] = b0[1]*x[1] + b1[1]*x[0] - a1[1]*y[0];
  for(int n = 2; n < N; n++)
    y[n] = b0[n]*x[n] + b1[n]*x[n-1] + b2[n]*x[n-2] - a1[n]*y[n-1] - a2[n]*y[n-2];
}

void getDesiredOutputForFilterBlip(int N, double frequency, double q, double *desiredOutput)
{
  //double coeffs[5];  romos::biquadBandpassCoeffs(coeffs, frequency, q);
  double coeffs[5];  romos::biquadBandpassConstSkirtCoeffs(coeffs, frequency, q);

  double *x  = new double[N];  generateImpulse(N, x);
  double *b0 = new double[N];  RAPT::rsArrayTools::fillWithValue(b0, N, coeffs[0]);
  double *b1 = new double[N];  RAPT::rsArrayTools::fillWithValue(b1, N, coeffs[1]);
  double *b2 = new double[N];  RAPT::rsArrayTools::fillWithValue(b2, N, coeffs[2]);
  double *a1 = new double[N];  RAPT::rsArrayTools::fillWithValue(a1, N, coeffs[3]);
  double *a2 = new double[N];  RAPT::rsArrayTools::fillWithValue(a2, N, coeffs[4]);

  getDesiredOutputForBiquad(N, x, b0, b1, b2, a1, a2, desiredOutput);

  delete[] x;
  delete[] b0;
  delete[] b1;
  delete[] b2;
  delete[] a1;
  delete[] a2;
}




//-----------------------------------------------------------------------------------------------------------------------------------------
// others:


void setModulePolyphony(romos::Module *module, bool shouldBePolyphonic, bool recursivelyForChildren)
{
  if( module->isContainerModule() && recursivelyForChildren == true )
    ((ContainerModule*) module)->setPolyphonicRecursively(shouldBePolyphonic);
  else
    module->setPolyphonic(shouldBePolyphonic);
}

bool checkResult(double **y, double **d, int numChannels, int numFrames, double tolerance)
{
  bool result = true;
  for(int c = 0; c < numChannels; c++)
    result &= RAPT::rsArrayTools::almostEqual(y[c], d[c], numFrames, tolerance);
  return result;
}

bool checkAndPrintResult(double **y, double **d, int numChannels, int numFrames, const char *testName, double tolerance)
{
  bool result = checkResult(y, d, numChannels, numFrames, tolerance);
  if( result == false )
    printf("%s %s %s", "!!! ", testName, " failed !!!\n");
  else
    printf("%s %s", testName, "passed\n");
  return result;
}

bool checkBlockProcessingAndPrintResult(romos::Module *module, double ***x, double ***y, double ***d,
                                        int maxNumFramesToProcess, int numTests, char *testName, double tolerance)
{
  // 1st test - per sample processing:
  //processModuleInFramesMono(module, maxNumFramesToProcess, x, y, NULL);
  processModuleInFrames(module, maxNumFramesToProcess, x, y, NULL, false);

  bool resultCorrect = checkResult(y[0], d[0], module->getNumOutputPins(), maxNumFramesToProcess, tolerance);

  // we need to trigger a note in order to mark voice with index 0 as active - todo: maybe leave voice 0 (optionally) always active
  // for effects:
  voiceAllocator.noteOn(69, 64);

  // now the actual block processing test with different (randomized) number of frames per block:
  module->resetStateForAllVoices();
  RAPT::rsRandomUniform(0.0, 1.0, 7); // seed for the total number of samples
  for(int i=1; i<=numTests; i++)
  {
    generateSilence(module->getNumOutputPins(), maxNumFramesToProcess, y[0]);
    int numFramesToProcess = (int) RAPT::rsRandomUniform(0.01*maxNumFramesToProcess, maxNumFramesToProcess, -1);  // total number of samples to process
    module->resetStateForAllVoices();
    processModuleInBlocks(module, numFramesToProcess, x, y, NULL, false);
    resultCorrect = checkResult(y[0], d[0], module->getNumOutputPins(), numFramesToProcess, tolerance);
    if( !resultCorrect )
    {
      voiceAllocator.killVoice(voiceAllocator.noteOff(69));
      printf("%s %s %s", "!!! ", testName, " failed !!!\n");
      return false;
    }
  }

  voiceAllocator.killVoice(voiceAllocator.noteOff(69));
  printf("%s %s", testName, "passed\n");
  return true;
}

bool checkProcessingFunctionsAndPrintResults(romos::Module *module, int numVoicesToCheck, int numFrames,
                                             double ***x, double ***y, double ***d,
                                             double tolerance, char *testName, std::vector<NoteEvent> *events)
{
  bool result = true;

  setModulePolyphony(module, false, true);
  result &= checkProcessingInFramesMonoAndPrintResult(module, numFrames, x, y, d, tolerance, testName, events);
  result &= checkProcessingInBlocksMonoAndPrintResult(module, numFrames, x, y, d, tolerance, testName, events);

  setModulePolyphony(module, true, true);
  result &= checkProcessingInFramesPolyAndPrintResult(module, numVoicesToCheck, numFrames, x, y, d, tolerance, testName, events);
  result &= checkProcessingInBlocksPolyAndPrintResult(module, numVoicesToCheck, numFrames, x, y, d, tolerance, testName, events);

  return result;
}

bool checkProcessingInFramesMonoAndPrintResult(romos::Module *module, int numFrames, double ***x, double ***y, double ***d,
                                               double tolerance, char *testName, std::vector<NoteEvent> *events)
{
  module->resetStateForAllVoices();
  RAPT::rsArrayTools::fillWithZeros(y[0][0], numFrames);

  processModuleInFrames(module, numFrames, x, y, events, false);

  bool result = RAPT::rsArrayTools::equal(d[0][0], y[0][0], numFrames);

  //Plotter::plotData(numFrames, t, d[0][0], y[0][0]);

  if( result == false )
    printf("%s %s %s", "!!! ", testName, " in frames, monophonic failed !!!\n");
  else
    printf("%s %s", testName, " in frames, monophonic passed\n");
  return result;
}

bool checkProcessingInBlocksMonoAndPrintResult(romos::Module *module, int numFrames, double ***x, double ***y, double ***d,
                                               double tolerance, char *testName, std::vector<NoteEvent> *events)
{
  module->resetStateForAllVoices();
  RAPT::rsArrayTools::fillWithZeros(y[0][0], numFrames);

  processModuleInBlocks(module, numFrames, x, y, events, false);

  bool result = RAPT::rsArrayTools::equal(d[0][0], y[0][0], numFrames);

  //Plotter::plotData(100, t, d[0][0], y[0][0]);

  if( result == false )
    printf("%s %s %s", "!!! ", testName, " in blocks, monophonic failed !!!\n");
  else
    printf("%s %s", testName, " in blocks, monophonic passed\n");
  return result;
}

bool checkProcessingInFramesPolyAndPrintResult(romos::Module *module, int numVoicesToCheck, int numFrames,
                                               double ***x, double ***y, double ***d,
                                               double tolerance, char *testName, std::vector<NoteEvent> *events)
{
  module->resetStateForAllVoices();
  int v;

  for(v = 0; v < numVoicesToCheck; v++)
    RAPT::rsArrayTools::fillWithZeros(y[v][0], numFrames);

  processModuleInFrames(module, numFrames, x, y, events, true);

  bool result = true;
  for(v = 0; v < numVoicesToCheck; v++)
    result &= RAPT::rsArrayTools::equal(d[v][0], y[v][0], numFrames);

  //Plotter::plotData(50, t, d[0][0], y[0][0]);

  if( result == false )
    printf("%s %s %s", "!!! ", testName, " in frames, polyphonic failed !!!\n");
  else
    printf("%s %s", testName, " in frames, polyphonic passed\n");
  return result;
}

bool checkProcessingInBlocksPolyAndPrintResult(romos::Module *module, int numVoicesToCheck, int numFrames,
                                               double ***x, double ***y, double ***d,
                                               double tolerance, char *testName, std::vector<NoteEvent> *events)
{
  module->resetStateForAllVoices();
  int v;

  for(v = 0; v < numVoicesToCheck; v++)
    RAPT::rsArrayTools::fillWithZeros(y[v][0], numFrames);

  processModuleInBlocks(module, numFrames, x, y, events, true);

  bool result = true;
  for(v = 0; v < numVoicesToCheck; v++)
    result &= RAPT::rsArrayTools::equal(d[v][0], y[v][0], numFrames);

  //Plotter::plotData(numFrames, t, d[1][0], y[1][0]);
  //Plotter::plotData(50, t, &d[0][0][950], &y[0][0][950]);

  if( result == false )
    printf("%s %s %s", "!!! ", testName, " in blocks, polyphonic failed !!!\n");
  else
    printf("%s %s", testName, " in blocks, polyphonic passed\n");
  return result;
}









void establishInputBlock(romos::Module *module, double ***inputs, int blockStart, int blockSize)
{
  RAPT::rsAssert(false, "Code seems to be out of date" );
  /*
  double *inputAddress = module->getAudioInputAddress();
  if( !module->isPolyphonic() )
  {
    for(unsigned int c = 0; c < module->getNumInputs(); c++)
    {
      for(int n = 0; n < blockSize; n++)
      {
        int offset = module->getInputPinMemoryOffset(n, 0, c);
        inputAddress[offset] = inputs[0][c][blockStart+n];
      }
    }
  }
  else
  {
    int numPlayingVoices           = processingStatus.getNumPlayingVoices();
    const int *playingVoiceIndices = processingStatus.getPlayingVoiceIndices();
    for(int v = 0; v < numPlayingVoices; v++)
    {
      int voiceIndex = playingVoiceIndices[v];
      for(unsigned int c = 0; c < module->getNumInputs(); c++)
      {
        for(int n = 0; n < blockSize; n++)
        {
          int offset = module->getInputPinMemoryOffset(n, voiceIndex, c);
          inputAddress[offset] = inputs[v][c][blockStart+n];
        }
      }
    }
  }
  */
}

void retrieveOutputBlock(romos::Module *module, double ***outputs, int blockStart, int blockSize)
{
  // code may be out of date
  double *outputAddress = module->getOutputPointer(0);
  if( !module->isPolyphonic() )
  {
    for(unsigned int c = 0; c < module->getNumOutputPins(); c++)
    {
      for(int n = 0; n < blockSize; n++)
      {
        int offset = module->getOutputPinMemoryOffset(n, 0, c);
        outputs[0][c][blockStart+n] = outputAddress[offset];
      }
    }
  }
  else
  {
    int numPlayingVoices           = processingStatus.getNumPlayingVoices();
    const int *playingVoiceIndices = processingStatus.getPlayingVoiceIndices();
    for(int v = 0; v < numPlayingVoices; v++)
    {
      int voiceIndex = playingVoiceIndices[v];
      for(unsigned int c = 0; c < module->getNumOutputPins(); c++)
      {
        for(int n = 0; n < blockSize; n++)
        {
          int offset = module->getOutputPinMemoryOffset(n, voiceIndex, c);
          outputs[v][c][blockStart+n] = outputAddress[offset];
        }
      }
    }
  }
}

void processFrame(romos::Module *module, double ***inputs, double ***outputs, int frameIndex)
{
  establishInputBlock(module, inputs, frameIndex, 1);
  module->processSampleFrame();
  retrieveOutputBlock(module, outputs, frameIndex, 1);
}

void processModuleInFramesNoEvents(romos::Module *module, int numFrames, double ***inputs, double ***outputs, int startIndex)
{
  for(int frameIndex = startIndex; frameIndex < startIndex + numFrames; frameIndex++)
    processFrame(module, inputs, outputs, frameIndex);
}

void processModuleInFrames(romos::Module *module, int numFrames, double ***inputs, double ***outputs, std::vector<NoteEvent> *events,
                           bool polyphonic)
{
  voiceAllocator.reset();
  module->setPolyphonic(polyphonic);

  if( events == NULL || events->size() == 0 )
    processModuleInFramesNoEvents(module, numFrames, inputs, outputs, 0);
  else
  {
    int frameIndex       = 0;
    int numEventsHandled = 0;
    while( frameIndex < numFrames )
    {
      NoteEvent e                 = events->at(numEventsHandled);
      int numFramesUntilNextEvent = e.getDeltaFrames() - frameIndex;

      // process chunk until the next event:
      processModuleInFramesNoEvents(module, numFramesUntilNextEvent, inputs, outputs, frameIndex);

      voiceAllocator.noteOn(e.getKey(), e.getVelocity());

      numEventsHandled++;
      frameIndex += numFramesUntilNextEvent;

      // process tail after all events have been handled:
      if( numEventsHandled == events->size() )
      {
        processModuleInFramesNoEvents(module, numFrames - frameIndex, inputs, outputs, frameIndex);
        frameIndex = numFrames;
      }
    }
  }
}

void processBlock(romos::Module *module, double ***inputs, double ***outputs, int blockStart, int blockSize)
{
  establishInputBlock(module, inputs, blockStart, blockSize);
  module->processBlockOfSamples(blockSize);
  retrieveOutputBlock(module, outputs, blockStart, blockSize);
}

void processModuleInBlocksNoEvents(romos::Module *module, int numFrames, double ***inputs, double ***outputs, int startIndex)
{
  int maxBlockSize = processingStatus.getBufferSize();
  int blockStart   = startIndex;
  while( blockStart < startIndex + numFrames )
  {
    int blockSize = (int) ::round(RAPT::rsRandomUniform(1.0, maxBlockSize));
    if( blockStart + blockSize > startIndex + numFrames )
      blockSize = startIndex + numFrames - blockStart;
    processBlock(module, inputs, outputs, blockStart, blockSize);
    blockStart += blockSize;
  }
}

void processModuleInBlocks(romos::Module *module, int numFrames, double ***inputs, double ***outputs, std::vector<NoteEvent> *events,
                           bool polyphonic)
{
  voiceAllocator.reset();
  module->setPolyphonic(polyphonic);

  if( events == NULL || events->size() == 0 )
    processModuleInBlocksNoEvents(module, numFrames, inputs, outputs, 0);
  else
  {
    int maxBlockSize     = processingStatus.getBufferSize();
    int blockStart       = 0;
    int numEventsHandled = 0;
    while( blockStart < numFrames )
    {
      NoteEvent e                 = events->at(numEventsHandled);
      int numFramesUntilNextEvent = e.getDeltaFrames() - blockStart;

      // process chunk until the next event:
      processModuleInBlocksNoEvents(module, numFramesUntilNextEvent, inputs, outputs, blockStart);

      voiceAllocator.noteOn(e.getKey(), e.getVelocity());
      numEventsHandled++;

      blockStart += numFramesUntilNextEvent;

      // process tail after all events have been handled:
      if( numEventsHandled == events->size() )
      {
        processModuleInBlocksNoEvents(module, numFrames - blockStart, inputs, outputs, blockStart);
        blockStart = numFrames;
      }
    }
  }
}







void exchangeModulePositions(romos::Module *module1, romos::Module *module2)
{
  int x1 = module1->getPositionX();
  int y1 = module1->getPositionY();
  module1->setPositionXY(module2->getPositionX(), module2->getPositionY());
  module2->setPositionXY(x1, y1);
}

bool randBool()
{
  return random(0.0, 1.0) > 0.5;
}

void randomizeContainment(romos::Module *module)
{
  romos::ContainerModule *container = dynamic_cast<romos::ContainerModule*> (module);
  if( container != NULL )
  {
    std::vector<romos::Module*> childModules    = container->getChildModules();
    std::vector<romos::Module*> toBeContainerized;
    std::vector<romos::Module*> toBeUnContainerized;
    for(unsigned int i = 0; i < childModules.size(); i++)
    {
      double randomNumber = random(0.0, 1.0);
      if( randomNumber < 1.0/3.0 )
        toBeContainerized.push_back(childModules[i]);
      else if( randomNumber < 2.0/3.0 )
        toBeUnContainerized.push_back(childModules[i]);
      else
      {
        // do nothing
      }
    }
    container->containerizeModules(toBeContainerized);
    container->unContainerizeModules(toBeUnContainerized);
  }
}

void printModuleStructure(romos::Module *module, int indent)
{
  char *spaces = new char[indent+1];
  for(int i=0; i<indent+1; i++)
    spaces[i] = ' ';
  spaces[indent] = '\0';
  printf("%s", spaces);
  delete[] spaces;

  printf("%s", module->getName().c_str() );
  /*
  printf("%s", " - #AI: " );
  printf("%d", module->getNumAudioInputs() );
  printf("%s", ", #AO: " );
  printf("%d", module->getNumAudioOutputs() );
  printf("%s", ", #EI: " );
  printf("%d", module->getNumEventInputs() );
  printf("%s", ", #EO: " );
  printf("%d", module->getNumEventOutputs() );
  */
  printf("%s", "\n" );

  romos::ContainerModule *container = dynamic_cast<romos::ContainerModule*> (module);
  if( container != NULL )
  {
    for(unsigned int i=0; i<container->getNumChildModules(); i++)
      printModuleStructure(container->getChildModule(i), indent+1);
  }
}
