#include "romos_GlobalFrameworkTests.h"
//using namespace rsTestRomos;

namespace rsTestRomos
{

VoiceAllocatorTest::VoiceAllocatorTest()
  : UnitTest("VoiceAllocatorTest")
{

}
VoiceAllocatorTest::~VoiceAllocatorTest()
{

}
bool VoiceAllocatorTest::runTest()
{
  bool testPassed = true;

  testPassed &= testStealOldestWithoutRetrigger();
  testPassed &= testStealOldestWithRetrigger();

  return testPassed;
}

bool VoiceAllocatorTest::testStealOldestWithoutRetrigger()
{
  romos::VoiceAllocator voiceAlloc;
  voiceAlloc.setNumVoices(3);
  voiceAlloc.setVoiceStealingMode(romos::VoiceAllocator::STEAL_OLDEST_VOICE);
  voiceAlloc.setRetriggerMode(false);
  const int* playingVoiceIndices = voiceAlloc.getPlayingVoiceIndices();

  int  noteOffVoice;
  bool testPassed = true;

  voiceAlloc.noteOn(1, 64);  // should use voice 0
  testPassed &= voiceAlloc.getNumPlayingVoices() == 1;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 0);

  voiceAlloc.resetTriggerFlags();
  testPassed &= areAllNoteOnTriggerFlagsUnchecked(voiceAlloc);

  voiceAlloc.noteOn(2, 64);  // should use voice 1
  testPassed &= voiceAlloc.getNumPlayingVoices() == 2;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= playingVoiceIndices[1] == 1;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 1);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(1, 100);  // should use voice 2
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 2);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(3, 64);  // should use voice 0 (the oldest)
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 0);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(4, 64);  // should use voice 1 (the oldest)
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 1);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(5, 64);  // should use voice 2 (the oldest)
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 2);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOff(4);  // voice 1
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAlloc, 1);
  testPassed &= voiceAlloc.isNoteOn(1)                     == false;
  testPassed &= voiceAlloc.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAlloc.isVoicePlaying(1)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAlloc.getNumPlayingVoices()           == 3;

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(4, 64);  // should use voice 0 (the oldest)
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 0);

  // now, voice 0 and 1 are playing note with key == 4 (but voice 1 has already velocity == 0)

  voiceAlloc.resetTriggerFlags();
  noteOffVoice = voiceAlloc.noteOff(4);
  testPassed &=  noteOffVoice == 0; // voice 0 should have received this because voice 1 is already off
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAlloc, 0);
  testPassed &= voiceAlloc.isNoteOn(0)                     == false;
  testPassed &= voiceAlloc.getNormalizedVelocityOfVoice(0) == 0.0;
  testPassed &= voiceAlloc.isVoicePlaying(0)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAlloc.getNumPlayingVoices()           == 3;

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.killVoice(1);
  testPassed &= voiceAlloc.isNoteOn(1)                     == false;
  testPassed &= voiceAlloc.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAlloc.isVoicePlaying(1)               == false;
  testPassed &= voiceAlloc.getNumPlayingVoices()           == 2;

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(6, 64);  // should use voice 1 (the one which just became available)
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 1);
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= voiceAlloc.getKeyOfVoice(1)      == 6;


  //voiceAlloc.killVoice(1);
  //voiceAlloc.killVoice(2);
  //voiceAlloc.killVoice(0);

  return testPassed;
}
bool VoiceAllocatorTest::testStealOldestWithRetrigger()
{
  romos::VoiceAllocator voiceAlloc;
  voiceAlloc.setNumVoices(3);
  voiceAlloc.setVoiceStealingMode(romos::VoiceAllocator::STEAL_OLDEST_VOICE);
  voiceAlloc.setRetriggerMode(true);
  const int* playingVoiceIndices = voiceAlloc.getPlayingVoiceIndices();

  int  noteOffVoice;
  bool testPassed = true;

  voiceAlloc.noteOn(1, 64);  // should use voice 0
  testPassed &= voiceAlloc.getNumPlayingVoices() == 1;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 0);

  voiceAlloc.resetTriggerFlags();
  testPassed &= areAllNoteOnTriggerFlagsUnchecked(voiceAlloc);

  voiceAlloc.noteOn(2, 64);  // should use voice 1
  testPassed &= voiceAlloc.getNumPlayingVoices() == 2;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= playingVoiceIndices[1] == 1;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 1);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(1, 100);  // should re-use voice 0
  testPassed &= voiceAlloc.getNumPlayingVoices() == 2;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 0);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(3, 64);  // should use voice 2 
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 2);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(4, 64);  // should use voice 1 (the oldest)
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 1);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(5, 64);  // should use voice 0 (the oldest)
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 0);

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOff(4);  // voice 1
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAlloc, 1);
  testPassed &= voiceAlloc.isNoteOn(1)                     == false;
  testPassed &= voiceAlloc.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAlloc.isVoicePlaying(1)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAlloc.getNumPlayingVoices()           == 3;

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(4, 64);  // should re-use voice 1
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 1);

  voiceAlloc.resetTriggerFlags();
  noteOffVoice = voiceAlloc.noteOff(4);
  testPassed &=  noteOffVoice == 1; // voice 0 should have received this because voice 1 is already off
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAlloc, 1);
  testPassed &= voiceAlloc.isNoteOn(1)                     == false;
  testPassed &= voiceAlloc.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAlloc.isVoicePlaying(1)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAlloc.getNumPlayingVoices()           == 3;

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.killVoice(1);
  testPassed &= voiceAlloc.isNoteOn(1)                     == false;
  testPassed &= voiceAlloc.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAlloc.isVoicePlaying(1)               == false;
  testPassed &= voiceAlloc.getNumPlayingVoices()           == 2;

  voiceAlloc.resetTriggerFlags();
  voiceAlloc.noteOn(6, 64);  // should use voice 1 (the one which just became available)
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAlloc, 1);
  testPassed &= voiceAlloc.getNumPlayingVoices() == 3;
  testPassed &= voiceAlloc.getKeyOfVoice(1)      == 6;

  return testPassed;
}
bool VoiceAllocatorTest::areAllNoteOnTriggerFlagsUnchecked(
  const romos::VoiceAllocator& voiceAlloc)
{
  bool result = true;
  for(int i = 0; i < voiceAlloc.getNumVoices(); i++)
    result &= voiceAlloc.getNoteOnTriggerFlag(i) == false;
  return result;
}
bool VoiceAllocatorTest::isNoteOnTriggerFlagCheckedExclusively(
  const romos::VoiceAllocator& voiceAlloc, 
  int voiceIndexThatShouldHaveFlagSet)
{
  bool result = true;
  for(int i = 0; i < voiceAlloc.getNumVoices(); i++)
  {
    if(i == voiceIndexThatShouldHaveFlagSet)
      result &= voiceAlloc.getNoteOnTriggerFlag(i) == true;
    else
      result &= voiceAlloc.getNoteOnTriggerFlag(i) == false;
  }
  return result;
}
bool VoiceAllocatorTest::isNoteOffTriggerFlagCheckedExclusively(
  const romos::VoiceAllocator& voiceAlloc, 
  int voiceIndexThatShouldHaveFlagSet)
{
  bool result = true;
  for(int i = 0; i < voiceAlloc.getNumVoices(); i++)
  {
    if(i == voiceIndexThatShouldHaveFlagSet)
      result &= voiceAlloc.getNoteOffTriggerFlag(i) == true;
    else
      result &= voiceAlloc.getNoteOffTriggerFlag(i) == false;
  }
  return result;
}



TriggerAndKillTest::TriggerAndKillTest()
  : ProcessingTest("TriggerAndKillTest")
{
  triggerAndKillModule = TestModuleBuilder::createTriggerAndKill("TriggerAndKill", 20, 10, true);
  //moduleToTest         = ModuleFactory::createModule(ModuleTypeRegistry::TOP_LEVEL_MODULE);
  moduleToTest = romos::moduleFactory.createModule("TopLevelModule");
  ((romos::ContainerModule*)moduleToTest)->addChildModule(moduleToTest);



  //theSynth             = new ModularSynth();
  //theSynth->getTopLevelModule()->addChildModule(triggerAndKillModule);
}
TriggerAndKillTest::~TriggerAndKillTest()
{
  //delete theSynth;
}
bool TriggerAndKillTest::runTest()
{
  //initTest();

  std::vector<romos::NoteEvent> events = TestEventGenerator::generateNoteOnOffPair(1, 64, 10, 100);
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(2, 64, 20, 100));
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(3, 64, 80, 100));
  events = TestEventGenerator::mergeEvents(events, TestEventGenerator::generateNoteOnOffPair(4, 64, 120, 100));



  /*
  static const int blockSize    = 100;
  static const int signalLength = 1030;

  double signalL[signalLength];
  double signalR[signalLength];

  int blockStart = 0;
  while( blockStart < signalLength - blockSize )
  {
    theSynth->getBlockOfSampleFramesStereo(&signalL[blockStart], &signalR[blockStart], blockSize, events);

    //romos::NoteEvent::updateDeltasAndRemoveObsoleteEvents(events, blockSize);

    blockStart += blockSize;
  }

  int remainingFrames = signalLength - blockStart;
  */



  return false;  // preliminary
}


TopLevelModuleTest::TopLevelModuleTest()
  : UnitTest("TopLevelModuleTest")
{
  //moduleToTest = (TopLevelModule*) ModuleFactory::createModule(ModuleTypeRegistry::TOP_LEVEL_MODULE);
  //moduleToTest = (TopLevelModule*) moduleFactory.createModule("TopLevelModule");
  moduleToTest = romos::moduleFactory.createTopLevelModule();
}
TopLevelModuleTest::~TopLevelModuleTest()
{
  romos::moduleFactory.deleteModule(moduleToTest);
}
bool TopLevelModuleTest::runTest()
{
  return true;  // at the moment, we only check for memory leaks
}

}