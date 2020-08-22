#include "TestHelpers.h"
using namespace rosic;  // get rid!
using namespace romos;

//wrap into namespace rsTestRomos

bool runUnitTests();

bool testSorting(bool verboseOutput);
//bool testModuleTypeRegistry();
bool testGain(bool verboseOutput);
bool testSumDiff(bool verboseOutput);
bool testWrappedSumDiff(bool verboseOutput);
bool testSummedDiffs(bool verboseOutput);
bool testMovingAverage(bool verboseOutput);
bool testLeakyIntegrator(bool verboseOutput);
bool testLeakyIntegratorDoubleDelay(bool verboseOutput); 
bool testTestFilter1(bool verboseOutput);
bool testBiquadMacro(bool verboseOutput);
bool testBiquadAtomic(bool verboseOutput);
//bool testContainerizationWithConstant(bool verboseOutput);
bool testContainerizationAddedConstants(bool verboseOutput);
bool testPinSorting(bool verboseOutput);

bool testAdderBlock(bool verboseOutput);
bool testBiquadAtomicBlock(bool verboseOutput);

bool testBlip(bool verboseOutput);

//bool testBlipOneNote(bool verboseOutput);
//bool testBlipTwoNotes(bool verboseOutput);

//bool testBlipBlock(bool verboseOutput);
bool testAdderProcessingFunctions(int numVoicesToCheck);
bool testUnitDelayProcessingFunctions(int numVoicesToCheck);
bool testWrappedAdderProcessingFunctions(int numVoicesToCheck);
bool testMonoToPoly(int numVoicesToCheck);
bool testPolyToMono(int numVoicesToCheck);

bool testGatedNoteFrequency(int numVoicesToCheck);

bool testTriggerAndKill(int numVoicesToCheck);


/*
bool testSummedDiffsBlock();
bool testBiquadMacro();
bool testBiquadAtomic();
bool testBiquadAtomicBlock();
bool testProcessAdderBlock();
bool testProcessSubtractorBlock();
//bool testEventHandling();
//bool testContainerizationWithBiquad();
*/

