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
  romos::VoiceAllocator voiceAllocator;
  voiceAllocator.setNumVoices(3);
  voiceAllocator.setVoiceStealingMode(romos::VoiceAllocator::STEAL_OLDEST_VOICE);
  voiceAllocator.setRetriggerMode(false);
  const int* playingVoiceIndices = voiceAllocator.getPlayingVoiceIndices();

  int  noteOffVoice;
  bool testPassed = true;

  voiceAllocator.noteOn(1, 64);  // should use voice 0
  testPassed &= voiceAllocator.getNumPlayingVoices() == 1;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 0);

  voiceAllocator.resetTriggerFlags();
  testPassed &= areAllNoteOnTriggerFlagsUnchecked(voiceAllocator);

  voiceAllocator.noteOn(2, 64);  // should use voice 1
  testPassed &= voiceAllocator.getNumPlayingVoices() == 2;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= playingVoiceIndices[1] == 1;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 1);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(1, 100);  // should use voice 2
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 2);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(3, 64);  // should use voice 0 (the oldest)
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 0);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(4, 64);  // should use voice 1 (the oldest)
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 1);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(5, 64);  // should use voice 2 (the oldest)
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 2);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOff(4);  // voice 1
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAllocator, 1);
  testPassed &= voiceAllocator.isNoteOn(1)                     == false;
  testPassed &= voiceAllocator.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAllocator.isVoicePlaying(1)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAllocator.getNumPlayingVoices()           == 3;

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(4, 64);  // should use voice 0 (the oldest)
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 0);

  // now, voice 0 and 1 are playing note with key == 4 (but voice 1 has already velocity == 0)

  voiceAllocator.resetTriggerFlags();
  noteOffVoice = voiceAllocator.noteOff(4);
  testPassed &=  noteOffVoice == 0; // voice 0 should have received this because voice 1 is already off
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAllocator, 0);
  testPassed &= voiceAllocator.isNoteOn(0)                     == false;
  testPassed &= voiceAllocator.getNormalizedVelocityOfVoice(0) == 0.0;
  testPassed &= voiceAllocator.isVoicePlaying(0)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAllocator.getNumPlayingVoices()           == 3;

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.killVoice(1);
  testPassed &= voiceAllocator.isNoteOn(1)                     == false;
  testPassed &= voiceAllocator.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAllocator.isVoicePlaying(1)               == false;
  testPassed &= voiceAllocator.getNumPlayingVoices()           == 2;

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(6, 64);  // should use voice 1 (the one which just became available)
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 1);
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= voiceAllocator.getKeyOfVoice(1)      == 6;


  //voiceAllocator.killVoice(1);
  //voiceAllocator.killVoice(2);
  //voiceAllocator.killVoice(0);

  return testPassed;
}
bool VoiceAllocatorTest::testStealOldestWithRetrigger()
{
  romos::VoiceAllocator voiceAllocator;
  voiceAllocator.setNumVoices(3);
  voiceAllocator.setVoiceStealingMode(romos::VoiceAllocator::STEAL_OLDEST_VOICE);
  voiceAllocator.setRetriggerMode(true);
  const int* playingVoiceIndices = voiceAllocator.getPlayingVoiceIndices();

  int  noteOffVoice;
  bool testPassed = true;

  voiceAllocator.noteOn(1, 64);  // should use voice 0
  testPassed &= voiceAllocator.getNumPlayingVoices() == 1;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 0);

  voiceAllocator.resetTriggerFlags();
  testPassed &= areAllNoteOnTriggerFlagsUnchecked(voiceAllocator);

  voiceAllocator.noteOn(2, 64);  // should use voice 1
  testPassed &= voiceAllocator.getNumPlayingVoices() == 2;
  testPassed &= playingVoiceIndices[0] == 0;
  testPassed &= playingVoiceIndices[1] == 1;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 1);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(1, 100);  // should re-use voice 0
  testPassed &= voiceAllocator.getNumPlayingVoices() == 2;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 0);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(3, 64);  // should use voice 2 
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 2);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(4, 64);  // should use voice 1 (the oldest)
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 1);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(5, 64);  // should use voice 0 (the oldest)
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 0);

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOff(4);  // voice 1
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAllocator, 1);
  testPassed &= voiceAllocator.isNoteOn(1)                     == false;
  testPassed &= voiceAllocator.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAllocator.isVoicePlaying(1)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAllocator.getNumPlayingVoices()           == 3;

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(4, 64);  // should re-use voice 1
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 1);

  voiceAllocator.resetTriggerFlags();
  noteOffVoice = voiceAllocator.noteOff(4);
  testPassed &=  noteOffVoice == 1; // voice 0 should have received this because voice 1 is already off
  testPassed &= isNoteOffTriggerFlagCheckedExclusively(voiceAllocator, 1);
  testPassed &= voiceAllocator.isNoteOn(1)                     == false;
  testPassed &= voiceAllocator.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAllocator.isVoicePlaying(1)               == true;   // voice is still playing in release phase (not yet killed)
  testPassed &= voiceAllocator.getNumPlayingVoices()           == 3;

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.killVoice(1);
  testPassed &= voiceAllocator.isNoteOn(1)                     == false;
  testPassed &= voiceAllocator.getNormalizedVelocityOfVoice(1) == 0.0;
  testPassed &= voiceAllocator.isVoicePlaying(1)               == false;
  testPassed &= voiceAllocator.getNumPlayingVoices()           == 2;

  voiceAllocator.resetTriggerFlags();
  voiceAllocator.noteOn(6, 64);  // should use voice 1 (the one which just became available)
  testPassed &= isNoteOnTriggerFlagCheckedExclusively(voiceAllocator, 1);
  testPassed &= voiceAllocator.getNumPlayingVoices() == 3;
  testPassed &= voiceAllocator.getKeyOfVoice(1)      == 6;

  return testPassed;
}
bool VoiceAllocatorTest::areAllNoteOnTriggerFlagsUnchecked(const romos::VoiceAllocator& voiceAllocator)
{
  bool result = true;
  for(int i = 0; i < voiceAllocator.getNumVoices(); i++)
    result &= voiceAllocator.getNoteOnTriggerFlag(i) == false;
  return result;
}
bool VoiceAllocatorTest::isNoteOnTriggerFlagCheckedExclusively(const romos::VoiceAllocator& voiceAllocator, int voiceIndexThatShouldHaveFlagSet)
{
  bool result = true;
  for(int i = 0; i < voiceAllocator.getNumVoices(); i++)
  {
    if(i == voiceIndexThatShouldHaveFlagSet)
      result &= voiceAllocator.getNoteOnTriggerFlag(i) == true;
    else
      result &= voiceAllocator.getNoteOnTriggerFlag(i) == false;
  }
  return result;
}
bool VoiceAllocatorTest::isNoteOffTriggerFlagCheckedExclusively(const romos::VoiceAllocator& voiceAllocator, int voiceIndexThatShouldHaveFlagSet)
{
  bool result = true;
  for(int i = 0; i < voiceAllocator.getNumVoices(); i++)
  {
    if(i == voiceIndexThatShouldHaveFlagSet)
      result &= voiceAllocator.getNoteOffTriggerFlag(i) == true;
    else
      result &= voiceAllocator.getNoteOffTriggerFlag(i) == false;
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